# Config for a dog GWAS dataset.
#
# `species = 'dog'` is the default, but it's stated explicitly here so this data
# directory is self-describing (the cat template does the same with 'cat'). It
# selects the dog profile in pheweb/species.py: chromosomes 1-38 plus 39 = X,
# the canFam4 assembly, the dog gene BED, and the canFam4 UCSC browser link.
species = 'dog'

# Optional per-dataset overrides (defaults come from the dog profile):
#   genes_bed = 'UU_CFAM_GSD_1.0_rosy.refseq.ensformat.bed'  # 5-col BED in this dir
#   ref_fasta_pattern = 'reference-canFam4-chrom-{chrom}.fa'  # in generated-by-pheweb/resources/

# hg_build_number only feeds the (non-functional-for-non-human) external
# LocusZoom links; the repo root uses 19, kept here for parity.
hg_build_number = 19
