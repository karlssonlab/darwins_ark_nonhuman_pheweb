# Config for a cat GWAS dataset.
#
# The presence of `species = 'cat'` is what makes this whole data directory a
# cat dataset: it selects the cat profile in pheweb/species.py (chromosome set
# 1-19 with 19=X, the F.catus_Fca126_mat1.0 assembly, the cat gene BED, and the
# GenArk UCSC hub). Everything else about running `pheweb process` / `pheweb
# serve` is identical to dog.
species = 'cat'

# Optional per-dataset overrides (defaults come from the cat profile):
#   genes_bed = 'cat_genes.bed'                                   # 5-col BED in this dir
#   ref_fasta_pattern = 'reference-F.catus_Fca126_mat1.0-chrom-{chrom}.fa'  # in generated-by-pheweb/resources/

# hg_build_number only feeds the (non-functional-for-non-human) external
# LocusZoom links; leave at the default.
