#!/usr/bin/env python3
'''
Convert a directory of per-chromosome recombination maps into the single
bgzipped + tabix-indexed file that PheWeb's region view serves.

Input: one file per chromosome, named `*_chr<N>_map.txt`, each with the header
`#POS  rate(cM/Mb)  Map(cM)` and three tab-separated columns sorted ascending by
position.  This is the layout of the canFam4 dog map from Kidd JM, Mamm Genome
2025;37(1):12, doi:10.1007/s00335-025-10178-0.

Output: `<chrom>  <pos>  <rate>  <pos_cm>`, bgzipped and tabixed, with
chromosomes emitted in the active species' `chrom_order_list` order.

Usage:
    ./prep_recomb_map.py <input_dir> <output_dir>
    ./prep_recomb_map.py --species cat <input_dir> <output_dir>

Needs only Python and pysam, so it runs directly in a cluster conda env; on a Mac
without pysam, `./run_local.sh recomb <data_dir> <input_dir>` runs it in the
container instead.

This is independent of `pheweb process` -- the map annotates the assembly, not the
association data, so no load step reads it and building it invalidates nothing.
Run it once per output location, before or after processing.

The output filename comes from the species profile's `recomb_map` key, so writing
into a data dir puts the file exactly where `conf.get_recomb_map_filepath()` looks
for it. Because the same assembly's map serves every dataset built on it, you can
instead write it to one shared directory and point each data dir's config.py at it
with an absolute path:

    recomb_map = '/work/pi_x/shared/dog_average_canFam4_recomb.tsv.gz'

Chromosome naming: input files use assembly names (`chr1`..`chr38`, `chrX`)
while PheWeb codes chromosomes numerically with X last (39 for dog, 19 for cat).
The mapping is derived by inverting the profile's `ucsc_chrom_map`, so this stays
correct for any species whose profile is filled in -- including cat, whose
chromosomes are cytogenetic names (A1, B2, ...) rather than numbers.
'''

import argparse
import os
import re
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pysam  # noqa: E402
from pheweb import species as species_module  # noqa: E402


def build_chrom_lookup(profile):
    '''file-chromosome-name -> PheWeb numeric chromosome, for this species.

    `ucsc_chrom_map` holds only the chromosomes that need remapping (e.g. dog's
    {'39': 'X'}); everything else is its own name. Both bare and `chr`-prefixed
    spellings are accepted since map files vary.'''
    lookup = {}
    for pheweb_chrom in profile['chrom_order_list']:
        assembly_name = profile['ucsc_chrom_map'].get(pheweb_chrom, pheweb_chrom)
        for spelling in (assembly_name, 'chr' + assembly_name):
            lookup[spelling] = pheweb_chrom
            lookup[spelling.upper()] = pheweb_chrom
            lookup[spelling.lower()] = pheweb_chrom
    return lookup


def find_input_files(input_dir, chrom_lookup):
    '''chrom -> filepath, for every `*_chr<N>_map.txt` we can place.'''
    found = {}
    unrecognized = []
    for filename in sorted(os.listdir(input_dir)):
        m = re.search(r'_(chr[^_]+)_map\.txt$', filename)
        if not m:
            continue
        raw_chrom = m.group(1)
        pheweb_chrom = chrom_lookup.get(raw_chrom)
        if pheweb_chrom is None:
            unrecognized.append(filename)
            continue
        if pheweb_chrom in found:
            raise SystemExit('two input files map to chromosome {}: {} and {}'.format(
                pheweb_chrom, os.path.basename(found[pheweb_chrom]), filename))
        found[pheweb_chrom] = os.path.join(input_dir, filename)
    if unrecognized:
        print('warning: ignoring {} file(s) with unrecognized chromosomes: {}'.format(
            len(unrecognized), ', '.join(unrecognized[:5])), file=sys.stderr)
    return found


def convert_one_chrom(pheweb_chrom, filepath, out):
    '''Stream one input file into `out`. Returns the number of rows written.

    Positions must come out ascending for tabix to index the file, so an
    out-of-order input is a hard error rather than something to silently sort --
    a map that isn't sorted is more likely malformed than merely unsorted.'''
    n_written = 0
    prev_pos = -1
    with open(filepath) as f:
        for line_num, line in enumerate(f, start=1):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 3:
                raise SystemExit('{}:{}: expected 3 tab-separated columns, got {}'.format(
                    filepath, line_num, len(fields)))
            pos_str, rate_str, cm_str = fields
            try:
                pos = int(pos_str)
                float(rate_str)  # validate only; write the original text
                float(cm_str)
            except ValueError:
                raise SystemExit('{}:{}: non-numeric value in {!r}'.format(
                    filepath, line_num, line.rstrip('\n')))
            if pos <= prev_pos:
                raise SystemExit('{}:{}: position {} is not greater than the previous {}; '
                                 'input must be sorted ascending'.format(
                                     filepath, line_num, pos, prev_pos))
            prev_pos = pos
            out.write('{}\t{}\t{}\t{}\n'.format(pheweb_chrom, pos, rate_str, cm_str))
            n_written += 1
    return n_written


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_dir', help='directory of per-chromosome *_chr<N>_map.txt files')
    parser.add_argument('output_dir', help='directory to write the map into (usually a PheWeb data dir)')
    parser.add_argument('--species', default=species_module.DEFAULT_SPECIES,
                        choices=species_module.known_species(),
                        help='species profile to use (default: %(default)s)')
    parser.add_argument('--output-name', default=None,
                        help="output filename; default comes from the species profile's recomb_map")
    args = parser.parse_args()

    profile = species_module.get_profile(args.species)
    output_name = args.output_name or profile.get('recomb_map')
    if not output_name:
        raise SystemExit('species {!r} has no recomb_map in its profile; pass --output-name'
                         .format(args.species))

    if not os.path.isdir(args.input_dir):
        raise SystemExit('no such directory: {}'.format(args.input_dir))
    if not os.path.isdir(args.output_dir):
        raise SystemExit('no such directory: {}'.format(args.output_dir))

    chrom_lookup = build_chrom_lookup(profile)
    inputs = find_input_files(args.input_dir, chrom_lookup)
    if not inputs:
        raise SystemExit('found no *_chr<N>_map.txt files in {}'.format(args.input_dir))

    expected = list(profile['chrom_order_list'])
    missing = [c for c in expected if c not in inputs]
    if missing:
        print('warning: no map for chromosome(s) {} -- the region view will show an '
              'empty recombination track there'.format(', '.join(missing)), file=sys.stderr)

    out_path = os.path.join(args.output_dir, output_name)
    # Write uncompressed to a temp file first: tabix_compress needs a complete,
    # position-sorted input, and this keeps a failed run from leaving a partial
    # map where get_recomb_map_filepath() would find it and trust it.
    tmp_fd, tmp_path = tempfile.mkstemp(suffix='.tsv', dir=args.output_dir)
    os.close(tmp_fd)
    total = 0
    try:
        with open(tmp_path, 'w') as out:
            out.write('#chrom\tpos\trate_cm_per_mb\tpos_cm\n')
            for pheweb_chrom in expected:
                filepath = inputs.get(pheweb_chrom)
                if filepath is None:
                    continue
                n = convert_one_chrom(pheweb_chrom, filepath, out)
                total += n
                print('chr{}: {:,} rows'.format(pheweb_chrom, n))

        print('compressing {:,} rows -> {}'.format(total, out_path))
        pysam.tabix_compress(tmp_path, out_path, force=True)
        pysam.tabix_index(
            filename=out_path, force=True,
            seq_col=0, start_col=1, end_col=1,  # pysam counts columns from 0
            line_skip=1,  # skip the header
        )
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    print('wrote {} and {}.tbi'.format(out_path, out_path))


if __name__ == '__main__':
    main()
