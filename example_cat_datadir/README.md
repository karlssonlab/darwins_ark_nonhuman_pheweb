# Example cat data directory

This folder is a template for a **cat** PheWeb dataset. The folder *is* the
dataset (see the deployment model in the main [README](../README.md)).

## What makes it "cat"

`config.py` declares `species = 'cat'`, which selects the cat profile in
`pheweb/species.py`:

- chromosomes `1-19` (numeric, `19` = X)
- assembly `F.catus_Fca126_mat1.0` (RefSeq GCF_018350175.1)
- gene BED `cat_genes.bed`
- UCSC GenArk hub `hub_6476843_GCF_018350175.1` (numeric chroms are mapped to
  cytogenetic names like `chrA2` server-side for the variant-page link)

## What you must add

```
example_cat_datadir/
├── config.py                 # provided (species = 'cat')
├── cat_genes.bed             # provided (5-col: chrom start end gene_name ensg)
├── pheno-list.json           # ADD: points at your cat summary-stat files
└── generated-by-pheweb/
    └── resources/
        ├── reference-F.catus_Fca126_mat1.0-chrom-1.fa   # ADD: per-chrom FASTAs
        ├── reference-F.catus_Fca126_mat1.0-chrom-2.fa   #      (headerless, newline-stripped,
        └── ...                                          #       chroms 1-19)
```

`cat_genes.bed` was built from the RefSeq annotation with:

```bash
awk -v OFS='\t' '{print $1,$2,$3,$4,$4}' \
    GCF_018350175.1_F.catus_Fca126_mat1.0_genes.bed > cat_genes.bed
```

(The 5th column duplicates the gene name into the `ensg` slot, matching what
`make_bed.sh` does for dog. No GTF conversion is needed.)

## Run

```bash
cd example_cat_datadir && pheweb process && pheweb serve --port 8001
# or, without cd:
PHEWEB_DATADIR=$PWD/example_cat_datadir pheweb process
PHEWEB_DATADIR=$PWD/example_cat_datadir pheweb serve --port 8001
```
