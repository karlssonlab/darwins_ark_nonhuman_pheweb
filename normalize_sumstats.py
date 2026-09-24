#!/usr/bin/env python3
'''
Convert GWAS summary statistics into PheWeb-ready CSV with alleles anchored to the
reference genome.

WHY THIS EXISTS
---------------
Upstream PheWeb requires that the reference allele "must always match the bases on
the reference genome", and warns about exactly the input these files have:

    "If you have an allele1 and allele2, and sometimes one or the other is the
     reference, then you'll need to modify your files."

GCTA's `A1` is the *effect* allele (usually the minor allele), not the alternate
allele. Assigning ref=A2/alt=A1 unconditionally -- as the original
`mlma_to_csv.py` did -- mislabels every variant where the minor allele happens to
be the reference base. Measured against canFam4 on this project's data, that is
**16.2% of variants**, uniformly across chromosomes.

This matters beyond cosmetics: `sites.py` merges variants across phenotypes on
the key `(chrom, pos, ref, alt)`, so allele labelling decides whether two
phenotypes' results for a locus combine into one variant or split into two.

WHAT THIS DOES
--------------
For every variant, look up the actual reference base and assign:

    reference base == A2  ->  ref=A2, alt=A1,  beta unchanged,  af = freq
    reference base == A1  ->  ref=A1, alt=A2,  beta NEGATED,    af = 1 - freq
    neither matches       ->  dropped and counted

The negation is required, not optional: PheWeb defines `beta` as the "effect size
(of alternate allele)". Swapping which allele is `alt` must flip the sign of the
effect, or the direction becomes wrong. Likewise `af` is an alternate-allele
frequency, so it complements when alt changes.

INPUT FORMATS (auto-detected from the header)
---------------------------------------------
  MLMA (tab):    Chr SNP bp A1 A2 Freq b se p
  ITEM (space):  SNP A1 A2 freq BETA se P N CHR POS

In both, `A1` is the effect allele and the frequency column is the A1 frequency.

OUTPUT
------
  chrom,pos,ref,alt,pval,beta,sebeta,af

`af` rather than `maf`: it is directional, and PheWeb derives MAF from it anyway
(`load_utils.get_maf` uses `min(af, 1-af)`). `sebeta` is carried through; the
original mlma_to_csv.py dropped it.

USAGE
-----
    ./normalize_sumstats.py <reference.fa> <input> <output.csv>
    ./normalize_sumstats.py --species cat <reference.fa> <input> <output.csv>

Needs only Python 3 and pysam, so it runs in a cluster conda env as-is.
The FASTA needs a .fai index; one is built automatically if absent.
'''

import argparse
import csv
import gzip
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pysam  # noqa: E402
from pheweb import species as species_module  # noqa: E402

MLMA_HEADER = ['Chr', 'SNP', 'bp', 'A1', 'A2', 'Freq', 'b', 'se', 'p']
ITEM_HEADER = ['SNP', 'A1', 'A2', 'freq', 'BETA', 'se', 'P', 'N', 'CHR', 'POS']

# column index of each value we need, per format
LAYOUTS = {
    'mlma': dict(delim='\t', chrom=0, pos=2, a1=3, a2=4, freq=5, beta=6, se=7, pval=8),
    'item': dict(delim=None, chrom=8, pos=9, a1=1, a2=2, freq=3, beta=4, se=5, pval=6),
}


def detect_layout(header_line):
    '''Identify the input format from its header. Refuses to guess: an unknown
    header means the column positions below would silently read wrong fields.'''
    for name, expected in (('mlma', MLMA_HEADER), ('item', ITEM_HEADER)):
        if header_line.split() == expected:
            return name
    raise SystemExit(
        'unrecognized header:\n  {}\nexpected one of:\n  MLMA: {}\n  ITEM: {}'.format(
            header_line.strip(), ' '.join(MLMA_HEADER), ' '.join(ITEM_HEADER)))


class Reference:
    '''Reference bases, one chromosome held in memory at a time.

    Per-variant `fetch()` calls are far too slow for ~10M rows, and the whole dog
    genome is 2.4GB. Inputs are grouped by chromosome, so caching the current one
    (~120MB at most) gives one load per chromosome.'''

    def __init__(self, fasta_path, profile):
        if not os.path.exists(fasta_path + '.fai'):
            print('building %s.fai ...' % fasta_path)
            pysam.faidx(fasta_path)
        self._fa = pysam.FastaFile(fasta_path)
        self._contigs = set(self._fa.references)
        # pheweb numeric chrom -> assembly name, e.g. dog '39' -> 'X'
        self._name_for_chrom = profile['ucsc_chrom_map']
        self._cur_chrom = None
        self._cur_seq = None
        self.reloads = 0

    def _contig_name(self, chrom):
        base = self._name_for_chrom.get(chrom, chrom)
        for candidate in ('chr' + base, base):
            if candidate in self._contigs:
                return candidate
        return None

    def base_at(self, chrom, pos):
        '''Reference base at 1-based `pos`, uppercased. None if unavailable.'''
        if chrom != self._cur_chrom:
            name = self._contig_name(chrom)
            if name is None:
                self._cur_chrom, self._cur_seq = chrom, None
            else:
                if self._cur_seq is not None:
                    self.reloads += 1
                self._cur_chrom = chrom
                self._cur_seq = self._fa.fetch(name).upper()
        if self._cur_seq is None or pos < 1 or pos > len(self._cur_seq):
            return None
        return self._cur_seq[pos - 1]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('reference', help='reference genome FASTA (.fai built if missing)')
    parser.add_argument('input', help='MLMA or ITEM summary statistics')
    parser.add_argument('output', help='PheWeb-ready CSV to write')
    parser.add_argument('--species', default=species_module.DEFAULT_SPECIES,
                        choices=species_module.known_species())
    args = parser.parse_args()

    profile = species_module.get_profile(args.species)
    valid_chroms = set(profile['chrom_order_list'])
    ref = Reference(args.reference, profile)

    with open(args.input) as fin:
        layout = LAYOUTS[detect_layout(fin.readline())]
        delim = layout['delim']

        n = kept = flipped = unchanged = 0
        dropped_no_match = dropped_chrom = dropped_bad = 0
        # PheWeb's reader raises if an input's chromosomes are out of order, and
        # MLMA files are grouped by chromosome but scrambled (33, 36, 35, 37, ...),
        # so rows have to be regrouped before writing.
        #
        # Spooling each chromosome to its own temp file rather than holding
        # everything in memory: these inputs are ~10M rows, which as Python tuples
        # is ~2.3GB -- close enough to the container's 7.7GB limit to matter, and
        # needlessly heavy on a shared cluster node. This way peak memory is one
        # chromosome (~500k rows), regardless of file size, and it stays correct
        # even if a chromosome appears in several non-contiguous blocks.
        spool_dir = tempfile.mkdtemp(prefix='normsumstats-', dir=os.path.dirname(os.path.abspath(args.output)))
        spools = {}
        try:
            for line in fin:
                fields = line.split(delim) if delim else line.split()
                if len(fields) < 9:
                    continue
                n += 1
                chrom = fields[layout['chrom']].strip()
                if chrom.startswith('chr'):
                    chrom = chrom[3:]
                if chrom == 'X':
                    # map to this species' numeric X before anything else
                    for k, v in profile['ucsc_chrom_map'].items():
                        if v == 'X':
                            chrom = k
                            break
                if chrom not in valid_chroms:
                    dropped_chrom += 1
                    continue
                try:
                    pos = int(fields[layout['pos']])
                    freq = float(fields[layout['freq']])
                    beta = float(fields[layout['beta']])
                    se = float(fields[layout['se']])
                    pval = float(fields[layout['pval']])
                except ValueError:
                    dropped_bad += 1
                    continue
                a1 = fields[layout['a1']].strip().upper()
                a2 = fields[layout['a2']].strip().upper()

                base = ref.base_at(chrom, pos)
                if base == a2:
                    # A2 is the reference: alt is the effect allele, beta as given
                    out_ref, out_alt, out_beta, out_af = a2, a1, beta, freq
                    unchanged += 1
                elif base == a1:
                    # A1 is the reference: alt becomes A2, so the effect flips
                    out_ref, out_alt, out_beta, out_af = a1, a2, -beta, 1.0 - freq
                    flipped += 1
                else:
                    dropped_no_match += 1
                    continue

                spool = spools.get(chrom)
                if spool is None:
                    spool = spools[chrom] = open(os.path.join(spool_dir, chrom), 'w')
                spool.write('%d\t%s\t%s\t%r\t%r\t%r\t%r\n'
                            % (pos, out_ref, out_alt, pval, out_beta, se, out_af))
                kept += 1

            for spool in spools.values():
                spool.close()

            # gzip when asked for it: PheWeb sniffs the magic bytes and reads
            # either transparently, and these files compress ~60%.
            tmp_path = args.output + '.part'
            opener = gzip.open if args.output.endswith('.gz') else open
            with opener(tmp_path, 'wt', newline='') as fout:
                w = csv.writer(fout)
                w.writerow(['chrom', 'pos', 'ref', 'alt', 'pval', 'beta', 'sebeta', 'af'])
                for chrom in profile['chrom_order_list']:
                    path = os.path.join(spool_dir, chrom)
                    if not os.path.exists(path):
                        continue
                    with open(path) as f:
                        rows = [ln.rstrip('\n').split('\t') for ln in f]
                    # PheWeb also requires ascending positions within a chromosome;
                    # ties broken on (ref, alt) so output is deterministic
                    rows.sort(key=lambda r: (int(r[0]), r[1], r[2]))
                    for r in rows:
                        w.writerow([chrom] + r)
                    del rows
            os.rename(tmp_path, args.output)
        finally:
            for spool in spools.values():
                if not spool.closed:
                    spool.close()
            shutil.rmtree(spool_dir, ignore_errors=True)

    def pct(x):
        return 100.0 * x / n if n else 0.0
    print('wrote %s' % args.output)
    print('  rows read            : %d' % n)
    print('  kept                 : %d (%.2f%%)' % (kept, pct(kept)))
    print('    A2 was the reference (beta unchanged) : %d (%.2f%%)' % (unchanged, pct(unchanged)))
    print('    A1 was the reference (beta NEGATED)   : %d (%.2f%%)' % (flipped, pct(flipped)))
    print('  dropped, neither allele matches reference: %d (%.3f%%)' % (dropped_no_match, pct(dropped_no_match)))
    print('  dropped, chromosome not in %s profile    : %d' % (args.species, dropped_chrom))
    print('  dropped, unparseable numbers             : %d' % dropped_bad)
    if ref.reloads > len(valid_chroms):
        print('  NOTE: reloaded chromosomes %d times -- input is not grouped by'
              ' chromosome, so this ran slower than necessary' % ref.reloads)


if __name__ == '__main__':
    main()
