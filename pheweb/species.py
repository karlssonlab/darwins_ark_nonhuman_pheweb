
'''
Species profiles: the single source of truth for everything that differs
between the species this PheWeb fork can serve (dog, cat).

Upstream PheWeb is parametrized by human genome build; the DAP fork replaced
that with hardcoded dog values.  This module reintroduces an abstraction so one
codebase serves both dog and cat, selected by the `species` config key (see
`conf.get_species()` / `conf.get_species_profile()`).  Every former dog hardcode
becomes a lookup into the active profile.

Each profile provides:
  chrom_order_list  ordered list of accepted chromosomes as strings (X coded numerically)
  chrom_aliases     alias -> canonical chrom map (empty for both current species)
  genes_bed         default gene-annotation BED filename (5-col), resolved under the data dir
  ref_fasta_pattern reference-FASTA filename pattern, resolved under the data dir's resources/;
                    `{chrom}` is filled with the numeric chromosome
  recomb_map        bgzipped+tabixed recombination-map filename, resolved under the data
                    dir; None means no local map, and the region view falls back to
                    LocusZoom's remote human track (see prep_recomb_map.py)
  assembly_label    human-readable assembly name (e.g. 'canFam4')
  ucsc_db           UCSC database / GenArk hub id for the external browser link
  ucsc_chrom_map    numeric chrom -> UCSC chrom name; chroms absent from the map
                    fall through to `chr{chrom}`
  display_name      branding: project / dataset display name, as in "performed
                    against the <display_name> dataset"
  site_title        branding: browser-tab, navbar and header site title

`genes_bed`, `ref_fasta_pattern`, `recomb_map`, `site_title` and `display_name`
are defaults; a data dir's config.py may override any of them with a key of the
same name (applied in conf.get_species_profile()).

The two branding keys in particular are usually worth overriding: they describe
the *study*, not the species, and one codebase serves several dog datasets.
'''

DEFAULT_SPECIES = 'dog'

SPECIES = {
    'dog': {
        # 1-38 autosomes + 39 = X (unchanged from the DAP fork's dog behavior)
        'chrom_order_list': [str(c) for c in range(1, 39)] + ['39'],
        'chrom_aliases': {},
        'genes_bed': 'UU_CFAM_GSD_1.0_rosy.refseq.ensformat.bed',
        'ref_fasta_pattern': 'reference-canFam4-chrom-{chrom}.fa',
        # Fine-scale recombination map, natively on canFam4 coordinates -- no
        # liftover involved. Kidd JM, Mamm Genome 2025;37(1):12,
        # doi:10.1007/s00335-025-10178-0. See docs/recombination-track.md and
        # docs/recombination-track-implementation.md.
        'recomb_map': 'dog_average_canFam4_recomb.tsv.gz',
        'assembly_label': 'canFam4',
        'ucsc_db': 'canFam4',
        # autosomes map directly to chr1..chr38; only 39 -> X needs remapping
        'ucsc_chrom_map': {'39': 'X'},
        'display_name': 'Dog Aging Project',
        'site_title': 'DAP PheWeb',
    },
    'cat': {
        # 1-18 autosomes + 19 = X
        'chrom_order_list': [str(c) for c in range(1, 20)],
        'chrom_aliases': {},
        'genes_bed': 'cat_genes.bed',
        'ref_fasta_pattern': 'reference-F.catus_Fca126_mat1.0-chrom-{chrom}.fa',
        # no cat recombination map yet; region view falls back to the remote track
        'recomb_map': None,
        'assembly_label': 'F.catus_Fca126_mat1.0',
        # GenArk hub for GCF_018350175.1; UCSC displays cytogenetic chrom names
        'ucsc_db': 'hub_6476843_GCF_018350175.1',
        'ucsc_chrom_map': {
            '1': 'A1', '2': 'A2', '3': 'A3', '4': 'B1', '5': 'B2', '6': 'B3',
            '7': 'B4', '8': 'C1', '9': 'C2', '10': 'D1', '11': 'D2', '12': 'D3',
            '13': 'D4', '14': 'E1', '15': 'E2', '16': 'E3', '17': 'F1',
            '18': 'F2', '19': 'X',
        },
        'display_name': 'Cat GWAS',
        'site_title': 'Cat PheWeb',
    },
}


def known_species():
    return sorted(SPECIES)


def get_profile(species):
    '''Return a shallow copy of the profile for `species` (callers may mutate the copy).'''
    if species not in SPECIES:
        from .utils import PheWebError
        raise PheWebError(
            "unknown species {!r}; must be one of {}".format(species, known_species()))
    profile = dict(SPECIES[species])
    # copy nested dicts so overrides never leak back into the registry
    profile['chrom_aliases'] = dict(profile['chrom_aliases'])
    profile['ucsc_chrom_map'] = dict(profile['ucsc_chrom_map'])
    profile['chrom_order_list'] = list(profile['chrom_order_list'])
    return profile
