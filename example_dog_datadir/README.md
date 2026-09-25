# Example dog data directory

This folder is a template for a **dog** PheWeb dataset — the dog counterpart of
[`example_cat_datadir/`](../example_cat_datadir/). The folder *is* the dataset
(see the deployment model in the main [README](../README.md)), so you don't need
a second clone of the repo to run a second PheWeb: one installed `pheweb`, one
folder per dataset, one `serve` per folder on its own port.

## What makes it "dog"

`config.py` declares `species = 'dog'` (also the default), which selects the dog
profile in `pheweb/species.py`:

- chromosomes `1-38` plus `39` = X
- assembly `canFam4` (UU_Cfam_GSD_1.0)
- gene BED `UU_CFAM_GSD_1.0_rosy.refseq.ensformat.bed`
- reference FASTAs `reference-canFam4-chrom-{chrom}.fa`
- UCSC browser db `canFam4` (only `39` → `chrX` needs remapping)

## What you must add

```
example_dog_datadir/
├── config.py                                    # provided (species = 'dog')
├── README.md                                    # provided
├── UU_CFAM_GSD_1.0_rosy.refseq.ensformat.bed    # ADD: 5-col dog gene BED
├── pheno-list.json                              # ADD: points at your sumstats
├── data/                                        # ADD: your per-pheno csv files
│   └── <pheno>_N-####.loco.csv
└── generated-by-pheweb/
    └── resources/
        ├── reference-canFam4-chrom-1.fa         # ADD: per-chrom FASTAs
        ├── reference-canFam4-chrom-2.fa         #      (headerless, newline-stripped,
        └── ...                                  #       chroms 1-38 and 39 = X)
```

### 1. Gene BED

5 columns, tab-separated: `chrom start end gene_name ensg`, chroms numeric (no
`chr` prefix). [`make_bed.sh`](../make_bed.sh) builds it from the RefSeq GTF:

```bash
awk -F'\t' '$3 == "gene" { print $1,$4,$5,$9 }' UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.gtf |
    sed 's/gene_id "//' | sed 's/".*//' | sed 's/^chr//' |
    awk -v OFS='\t' '$1 ~ /^[0-9]+$/ { print $1,$2,$3,$4,$4 }' \
    > UU_CFAM_GSD_1.0_rosy.refseq.ensformat.bed
```

### 2. Association files

One csv per phenotype, as documented in the main README:

```
chrom,pos,ref,alt,pval
1,629,C,T,0.0856026
```

[`mlma_to_csv.py`](../mlma_to_csv.py) converts GCTA `.mlma` output to this.

### 3. `pheno-list.json`

Same format as the repo root's [`pheno-list.json`](../pheno-list.json) (the
313-phenotype DAP list — a good reference for `phenocode` / `category` /
`num_samples`), with paths relative to *this* folder. Generate one with
[`create_pheno_list.py`](../create_pheno_list.py) or by hand.

### 4. Reference FASTAs

One file per chromosome, **header stripped and newlines removed** so PheWeb can
`seek()` to a position:

```bash
for c in {1..38} 39; do
  gzip -cd chr${c}.fa.gz | tail -n +2 | tr -d '\n' \
    > generated-by-pheweb/resources/reference-canFam4-chrom-${c}.fa
done
```

(Chromosome `39` is X: use the `chrX.fa` sequence for it.)

## Run

Natively, with `pheweb` on your PATH:

```bash
cd example_dog_datadir && pheweb process && pheweb serve --port 8000
# or, without cd (e.g. sbatch):
PHEWEB_DATADIR=$PWD/example_dog_datadir pheweb process
PHEWEB_DATADIR=$PWD/example_dog_datadir pheweb serve --port 8000
```

Locally on a Mac (the pinned deps need Python 3.8 / linux wheels), use the
container runner from the repo root:

```bash
./run_local.sh build
./run_local.sh process example_dog_datadir
./run_local.sh serve   example_dog_datadir 8000
```

Cat runs identically — point either command at `example_cat_datadir` on a
different port to serve both at once.
