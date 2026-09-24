#!/usr/bin/env python3
'''
Convert a `Gene / Chrom / Start / End / Coding` gene table into the 5-column BED
that PheWeb annotates variants with.

This is the script form of the conversion that produced
`dog_compulsive_disorder_pheweb/cf4_genes_noLOC.bed` from
`Unique.cf4genes2.txt`. Reasoning and measurements are in
docs/items-and-allele-normalization-implementation.md section 1.

Input (tab-separated, header required, CRLF tolerant):

    Gene    Chrom   Start   End     Coding
    ENPP1   1       209085  282880  protein_coding

Output (tab-separated, no header) -- what `get_gene_tuples_with_ensg()` reads:

    chrom   start   end     gene_name       ensg
    1       209084  282880  ENPP1           ENPP1

FOUR TRANSFORMATIONS, each a silent corruption if skipped
---------------------------------------------------------
1. **Coordinates.** The input is 1-based inclusive; BED is 0-based half-open, so
   `start` is decremented and `end` left alone. This was established against
   UCSC's canFam4 `ncbiRefSeq` genePred, where `txStart` is 0-based by definition:
   the two existing BEDs in this repo agree with UCSC exactly (34,579/34,579 on
   both start and end), while `Unique.cf4genes2.txt` is +1 on start for all
   17,212 genes it shares with UCSC. Pass --no-coord-shift for a source that is
   already 0-based.

2. **Column order.** `get_gene_tuples_with_ensg()` reads positionally, and the
   input puts the gene name first.

3. **Chromosome names.** PheWeb codes chromosomes numerically with X last (39 for
   dog, 19 for cat) and *asserts* every chromosome is in the species chrom order,
   so an unmapped name aborts the load. The mapping inverts the species profile's
   `ucsc_chrom_map`, so it stays correct for cat too (cytogenetic names).

4. **Scaffolds.** Unplaced contigs (`Un_*`, `Y_unplaced_*`, anything not in the
   profile) are dropped -- PheWeb has no chromosome to place them on.

The `ensg` column is filled with the gene symbol. The format wants an Ensembl gene
id, but this is a RefSeq annotation with no ENSG; duplicating the symbol is what
the cat gene BED already does (see README) and it only feeds the gene-alias db.

--drop-loc
----------
`LOC<number>` is NCBI's placeholder symbol for a gene with no official symbol.
Dropping them gives cleaner labels but a *worse* nearest-gene call: on this
project's data it relabels ~47% of variants and moves the 90th-percentile distance
from ~67kb to ~409kb, because only the two flanking genes are ever considered.
It also discards 3,823 protein-coding genes that merely lack a symbol. Off by
default; `cf4_genes_noLOC.bed` was built with it on, deliberately.

USAGE
-----
    # what produced cf4_genes_noLOC.bed:
    ./prep_gene_bed.py --drop-loc Unique.cf4genes2.txt cf4_genes_noLOC.bed

    ./prep_gene_bed.py --species cat <input.txt> <output.bed>
    ./prep_gene_bed.py --compare-to <existing.bed> <input.txt> <output.bed>

Needs only the standard library.
'''

import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pheweb import species as species_module  # noqa: E402

EXPECTED_HEADER = ['Gene', 'Chrom', 'Start', 'End', 'Coding']


def build_chrom_lookup(profile):
    '''input chromosome name -> PheWeb numeric chromosome.

    `ucsc_chrom_map` maps pheweb-chrom -> assembly name (dog: {'39': 'X'}); invert
    it. Bare, `chr`-prefixed and case variants are all accepted, since sources
    differ on all three.'''
    lookup = {}
    for pheweb_chrom in profile['chrom_order_list']:
        assembly_name = profile['ucsc_chrom_map'].get(pheweb_chrom, pheweb_chrom)
        for spelling in (assembly_name, 'chr' + assembly_name):
            for variant in (spelling, spelling.upper(), spelling.lower()):
                lookup[variant] = pheweb_chrom
    return lookup


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='the Gene/Chrom/Start/End/Coding table')
    parser.add_argument('output', help='BED file to write')
    parser.add_argument('--species', default=species_module.DEFAULT_SPECIES,
                        choices=species_module.known_species())
    parser.add_argument('--drop-loc', action='store_true',
                        help='drop LOC* genes (see the note above before using this)')
    parser.add_argument('--no-coord-shift', action='store_true',
                        help='input is already 0-based; do not decrement start')
    parser.add_argument('--compare-to', default=None, metavar='BED',
                        help='an existing BED to diff the resulting gene set against')
    args = parser.parse_args()

    profile = species_module.get_profile(args.species)
    lookup = build_chrom_lookup(profile)
    chrom_rank = {c: i for i, c in enumerate(profile['chrom_order_list'])}

    rows = []
    dropped_chrom = collections.Counter()
    dropped_loc = 0
    biotypes = collections.Counter()
    loc_biotypes = collections.Counter()

    with open(args.input, newline='') as f:
        header = f.readline().rstrip('\r\n').split('\t')
        if header != EXPECTED_HEADER:
            raise SystemExit('unexpected header {!r}\nexpected {!r}'.format(header, EXPECTED_HEADER))
        for line_num, line in enumerate(f, start=2):
            line = line.rstrip('\r\n')
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) != 5:
                raise SystemExit('{}:{}: expected 5 columns, got {}'.format(
                    args.input, line_num, len(fields)))
            gene, chrom_raw, start_s, end_s, coding = fields

            if args.drop_loc and gene.startswith('LOC'):
                dropped_loc += 1
                loc_biotypes[coding] += 1
                continue

            pheweb_chrom = lookup.get(chrom_raw)
            if pheweb_chrom is None:
                dropped_chrom[chrom_raw] += 1
                continue
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                raise SystemExit('{}:{}: non-numeric coordinates {!r} {!r}'.format(
                    args.input, line_num, start_s, end_s))
            if not args.no_coord_shift:
                start -= 1
            if start < 0 or end < start:
                raise SystemExit('{}:{}: bad interval {}-{} for {}'.format(
                    args.input, line_num, start, end, gene))
            rows.append((pheweb_chrom, start, end, gene))
            biotypes[coding] += 1

    if not rows:
        raise SystemExit('no usable rows -- are the chromosome names what you expect?')

    rows.sort(key=lambda r: (chrom_rank[r[0]], r[1], r[2]))
    tmp_path = args.output + '.part'
    with open(tmp_path, 'w') as out:
        for chrom, start, end, gene in rows:
            # ensg column duplicates the symbol: RefSeq annotation, no ENSG available
            out.write('{}\t{}\t{}\t{}\t{}\n'.format(chrom, start, end, gene, gene))
    os.rename(tmp_path, args.output)

    print('wrote %s' % args.output)
    print('  genes kept   : %d' % len(rows))
    print('  coord shift  : %s' % ('none (--no-coord-shift)' if args.no_coord_shift
                                   else '1-based inclusive -> 0-based half-open (start-1)'))
    chroms_present = sorted({r[0] for r in rows}, key=lambda c: chrom_rank[c])
    print('  chromosomes  : %d of %d' % (len(chroms_present), len(profile['chrom_order_list'])))
    missing = [c for c in profile['chrom_order_list'] if c not in chroms_present]
    if missing:
        print('  NO GENES ON  : %s' % ','.join(missing))
        print('                 -> variants there get a blank nearest gene, silently')

    if args.drop_loc:
        print('  LOC* dropped : %d' % dropped_loc)
        pc = loc_biotypes.get('protein_coding', 0)
        if pc:
            print('                 of which %d are protein_coding -- real genes with no symbol' % pc)
    if dropped_chrom:
        print('  dropped      : %d rows on %d sequences with no PheWeb chromosome'
              % (sum(dropped_chrom.values()), len(dropped_chrom)))
        for name, cnt in dropped_chrom.most_common(4):
            print('      %-34s %d' % (name, cnt))
        if len(dropped_chrom) > 4:
            print('      ... and %d more' % (len(dropped_chrom) - 4))

    dupes = [g for g, cnt in collections.Counter(r[3] for r in rows).items() if cnt > 1]
    if dupes:
        print('  duplicate gene names: %d (e.g. %s)' % (len(dupes), ', '.join(sorted(dupes)[:3])))
        print('      fine for nearest-gene; /gene/<name> resolves each to one locus only')
    print('  biotypes     : %s' % ', '.join('%s=%d' % kv for kv in biotypes.most_common(5)))

    if args.compare_to:
        old = set()
        with open(args.compare_to) as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 4:
                    old.add(parts[3])
        new = {r[3] for r in rows}
        print('\ncompared to %s' % args.compare_to)
        print('  genes in old : %d' % len(old))
        print('  genes in new : %d' % len(new))
        print('  added        : %d' % len(new - old))
        print('  removed      : %d' % len(old - new))


if __name__ == '__main__':
    main()
