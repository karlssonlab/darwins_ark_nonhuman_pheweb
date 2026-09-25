
This is a fork of PheWeb for non-human GWAS. It serves both **dog** and **cat** data from one codebase, selected per dataset.

## Species support (dog & cat)

Upstream PheWeb is parametrized by human genome build. This fork replaces that
with a **species profile** (`pheweb/species.py`) selected by a `species` config
key. Each dataset directory has a `config.py` declaring `species = 'dog'` (the
default, for backward compatibility) or `species = 'cat'`, which selects:

| | dog | cat |
|---|---|---|
| chromosomes | `1-38` + `39`=X | `1-18` + `19`=X (numeric) |
| assembly | canFam4 | F.catus_Fca126_mat1.0 (GCF_018350175.1) |
| gene BED | `UU_CFAM_GSD_1.0_rosy.refseq.ensformat.bed` | `cat_genes.bed` |
| ref FASTAs | `reference-canFam4-chrom-{chrom}.fa` | `reference-F.catus_Fca126_mat1.0-chrom-{chrom}.fa` |
| UCSC browser | `canFam4` | GenArk hub `hub_6476843_GCF_018350175.1` (cytogenetic chrom names mapped server-side) |

**The folder is the dataset.** One installed `pheweb`; many data folders, each
declaring its own species:

```
/data/
├── dog_dap/       config.py: species='dog'   + pheno-list.json, dog genes.bed, resources/ (canFam4 FASTAs)
├── cat_study_A/   config.py: species='cat'   + pheno-list.json, cat_genes.bed, resources/ (F.catus FASTAs)
└── cat_study_B/   config.py: species='cat'   + pheno-list.json (different cat data); may share cat genes/FASTAs
```

Point `pheweb` at a folder via `cd` or the `PHEWEB_DATADIR` env var (a
`PHEWEB_SPECIES` env var can also override the config key). Nothing about the
commands changes between species:

```bash
cd /data/cat_study_A && pheweb process && pheweb serve --port 8001
# or, for sbatch/cluster use:
PHEWEB_DATADIR=/data/cat_study_A pheweb process
PHEWEB_DATADIR=/data/cat_study_A pheweb serve --port 8001
```

Run dog and cat simultaneously by starting one `serve` per folder on its own
port. Templates for both, each listing exactly what to drop in:

- [`example_dog_datadir/`](example_dog_datadir/) — dog (`config.py` + where to
  put the dog gene BED, sumstats, and canFam4 FASTAs)
- [`example_cat_datadir/`](example_cat_datadir/) — cat (`config.py` + prebuilt
  `cat_genes.bed` + where to drop FASTAs and `pheno-list.json`)

### Running locally (macOS / any machine without Python 3.8)

The pinned deps (Flask 1.1, numpy 1.19, pysam 0.16, gevent 21.1) have no wheels
for modern Python and none at all for macOS arm64, so local runs go through a
Python 3.8 / linux-amd64 container ([`docker/Dockerfile`](docker/Dockerfile),
pins in [`docker/requirements-py38.txt`](docker/requirements-py38.txt)). The
dataset folder is bind-mounted at `/data`, so the same image runs any dog or cat
dataset:

```bash
./run_local.sh build                              # once (and after editing pheweb/)
./run_local.sh process example_dog_datadir        # pheweb process
./run_local.sh serve   example_dog_datadir 8000   # dog -> http://localhost:8000
./run_local.sh serve   example_cat_datadir 8001   # cat -> http://localhost:8001
./run_local.sh test                               # pytest, no dataset needed
```

Why a container, what it pins and why, and how the dog dataset was built:
[`docs/local-dog-pheweb-implementation.md`](docs/local-dog-pheweb-implementation.md).

On the cluster, use the native install instead:
[`run_process_pheweb.sbatch`](run_process_pheweb.sbatch) and
[`run_serve_pheweb.sh`](run_serve_pheweb.sh), both of which take the data dir as
their first argument.

### Cat gene BED

Cat needs a 5-column BED (`chrom start end gene_name ensg`). The RefSeq
annotation is 4-column, so duplicate the gene name into the 5th column (no GTF
conversion; `make_bed.sh` stays dog-only):

```bash
awk -v OFS='\t' '{print $1,$2,$3,$4,$4}' \
    GCF_018350175.1_F.catus_Fca126_mat1.0_genes.bed > cat_genes.bed
```

### Genome-wide significance threshold

PheWeb's 5e-8 is a human GWAS convention — it encodes an assumption about the
number of independent tests in a *human* genome. Set your own per data directory
in `config.py`:

```python
significance_threshold = 4e-7    # default is 5e-8
```

One key moves four things:

| What | When it takes effect |
|---|---|
| Manhattan significance line + its tooltip | immediately on restart |
| Which peaks get gene labels (top 7 below the threshold) | immediately on restart |
| Region view's significance line | immediately on restart |
| `num_significant_in_peak` (Manhattan tooltips) | **only after reprocessing** |

The first three are read by the browser at page load, so a server restart is
enough. The fourth is computed during `pheweb process` and baked into the
Manhattan JSON, so it keeps the old value until you re-run the load steps. If you
only care about the line and the labels, no reprocessing is needed.

It must be stricter than `manhattan_peak_pval_threshold` (default 1e-6), since
peaks have to be found before their significant variants can be counted. Setting
it looser raises a clear error at startup rather than an unexplained assertion
mid-load. Set `manhattan_peak_variant_counting_pval_threshold` explicitly only if
you deliberately want the counting threshold to differ from the line.

### Out of scope

Local LD and GWAS-catalog tracks for non-human assemblies are not provided (the
region view's external human LocusZoom services return nothing useful for dog or
cat); and cat rsid/dbSNP annotation is stubbed out (as it already is for dog).

The recombination track *is* now local for dog — Jeffrey Kidd's fine-scale
canFam4 map, built by `prep_recomb_map.py`. Where no species-appropriate
map is configured (currently cat), the track is omitted rather than falling back
to the human one. See [docs/recombination-track.md](docs/recombination-track.md).


# How to Cite PheWeb
This is a slight modification of the original PheWeb which should be cited:
Gagliano Taliun, S.A., VandeHaar, P. et al. Exploring and visualizing large-scale genetic associations by using PheWeb. *Nat Genet* 52, 550–552 (2020).

The dog recombination track in the region view is not ours; cite it separately if
you use it:
Kidd JM. Fine-scale recombination rates inferred using the canFam4 assembly are
strongly correlated with previous maps of dog recombination. *Mamm Genome*. 2025
Dec 12;37(1):12. doi:10.1007/s00335-025-10178-0. PMID: 41387639. PMCID: PMC12701031.

# How to Build a PheWeb for DAP data

### 1. Install PheWeb

```bash
pip3 install pheweb
```
If on ASU SOL cluster then read the "sol_install_notes.txt"

### 2. Create a directory and `config.py` for your new dataset
Create a folder for the dataset and put a `config.py` in it declaring the
species, e.g. `species = 'dog'` (the default) or `species = 'cat'`. That key
selects the species profile (chromosomes, assembly, gene BED, reference FASTAs,
UCSC link, branding) — see the "Species support" section above. `pheweb` reads
this folder via `cd` into it or via `PHEWEB_DATADIR=/path/to/folder`.

### 3. Convert MLMA files into csv files
There are multiple options, but PheWeb expects files to look like:

    chrom,pos,ref,alt,pval
    1,629,C,T,0.0856026
    1,2076,G,T,0.835506
    1,2388,C,T,0.0574887

This conversion was done with a simple python script `mlma_to_csv.py`
and the files were stored in the `mlmas/` subdir.

### 4. Make a list of your phenotypes

Inside of your data directory, you need a file named `pheno-list.json` that looks like this:

```json
[
    {
        "assoc_files": [
            "mlmas/DogAgingProject_gp-0.70_biallelic-snps_N-6358_maf-0.01_geno-0.05_hwe-1.0E-20-midp-keep-fewhet_phe-dd_weight_lbs_N-6279_cov-dd_sex_N-6279_qcov-Estimated_Age_Years_at_HLES_N-6279_chr1.loco.csv"
        ],
        "phenocode": "Weight",
        "category": "Physical"
    },
    {
        "assoc_files": [
            "mlmas/DogAgingProject_gp-0.70_biallelic-snps_N-6358_maf-0.01_geno-0.05_hwe-1.0E-20-midp-keep-fewhet_phe-pa_activity_level_N-6279_cov-dd_sex_N-6279_qcov-Estimated_Age_Years_at_HLES-dd_weight_lbs_N-6279_chr1.loco.csv"
        ],
        "phenocode": "Activity level",
        "category": "Activity"
    },
    {
        "assoc_files": [
            "mlmas/DogAgingProject_gp-0.70_biallelic-snps_N-6358_maf-0.01_mp_dental_extraction_N-1524_cov-dd_sex_N-6279_Estimated_Age_Years_at_HLES-dd_weight_lbs_N-6279_chr1.loco.csv"
        ],
        "phenocode": "Dental extraction",
        "category": "Dental"
    }
]
```

Each phenotype needs `assoc_files` (a list of paths to association files) and `phenocode` (a string representing your phenotype that is used in filenames and URLs, comprised of `[A-Za-z0-9_~-]`).

If you want, you can also include:

- `phenostring` (string): a name for the phenotype. Shown in tables and tooltips and page headers.
- `category` (string): groups together phenotypes in the PheWAS plot. Shown in tables and tooltips.
- `num_cases`, `num_controls`, and/or `num_samples` (number): if your input data only has `AC` or `MAC`, this will be used to calculated `AF` or `MAF`.  Shown in tooltips.  If your input data has correctly-named columns for these, the command `pheweb phenolist read-info-from-association-files` will add them into your existing `pheno-list.json`.
- anything else you want, but you'll have to modify templates to use it.

### 5. Load your association files

Run `pheweb process`.

To distribute jobs across a cluster, follow [these instructions](etc/detailed-loading-instructions.md#distributing-jobs-across-a-cluster).

To include VEP annotations, follow [these instructions](etc/detailed-loading-instructions.md#annotating-with-vep).

If something breaks and you can't understand the error message or it's something that PheWeb should support by default, [open an issue on github](https://github.com/statgen/pheweb/issues/new) or email me.

### 6. Serve the website

Run `pheweb serve --open`.

That command should either open a browser to your new PheWeb, or it should give you a URL that you can open in your browser to access your new PheWeb.
If it doesn't, follow [the directions for hosting a PheWeb and accessing it from your browser](etc/detailed-webserver-instructions.md#hosting-a-pheweb-and-accessing-it-from-your-browser).

### More options:

To run pheweb through systemd, see sample file [here](etc/pheweb.service).
To use Apache2 or Nginx, see instructions [here](etc/detailed-webserver-instructions.md#using-apache2-or-nginx).
To require login via OAuth, see instructions [here](etc/detailed-webserver-instructions.md#using-oauth).
To track page views with Google Analytics, see instructions [here](etc/detailed-webserver-instructions.md#using-google-analytics).
To reduce storage use, see instructions [here](etc/detailed-webserver-instructions.md#reducing-storage-use).
To customize page contents, see instructions [here](etc/detailed-webserver-instructions.md#customizing-page-contents).

PheWeb can display phenotype correlations generated by [another tool](https://github.com/statgen/pheweb-rg-pipeline).
To use this feature, set `show_correlations = True`  in `config.py` and place the output of the rg pipeline as `pheno-correlations.txt` in the same folder as `pheno-list.json`.

To hide the button for downloading summary stats, add `download_pheno_sumstats = "secret"` and `SECRET_KEY = "your random string"` in `config.py`.  That will make a secret page (printed to the console when you start the server) to share summary stats.
To hide the button for downloading top hits and phenotypes, add `download_top_hits = "hide"` and `download_phenotypes = "hide"` respectively.

To allow dynamically filtering the manhattan plot, run `pheweb best-of-pheno` and set `show_manhattan_filter_button=True` in `config.py`.

# Modifying DAP PheWeb

Here are the steps I took to get a development installation of PheWeb running on my mac laptop:

1. Clone the DAP PheWeb repo

* `git clone https://github.com/karlssonlab/darwins_ark_nonhuman_pheweb.git`

2. cd into the repo and create a virtual environment. I used python 3.8 after I had some trouble with 3.11

* `cd DAP_pheweb`

* `python3.8 -m venv .venv`

* `source .venv/bin/activate`

3. Install wheel and then pheweb as editable

* `pip install wheel`

* `pip install -e .`

4. Check that pheweb is installed

* `which pheweb`

5. Run the included tests, hopefully they should all pass

* `pip install pytest`

* `python -m pytest`

6. Do a test run of the server. Run the following command, and while that is running, open a browser window and go to http://0.0.0.0:8000/ and there should be a small pheweb example application running.

* `./tests/run-all.sh`
