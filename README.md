# Darwin's Ark non-human PheWeb

A fork of [PheWeb](https://github.com/statgen/pheweb) for **non-human GWAS**. It
serves **dog** and **cat** data from one codebase, selected per dataset.

The live instance built from this repo is the **Darwin's Ark Dog Compulsive
Disorder PheWeb** (3 factor GWAS + 14 survey items, canFam4).

## Species support (dog & cat)

Upstream PheWeb is parametrized by human genome build. This fork replaces that
with a **species profile** ([`pheweb/species.py`](pheweb/species.py)) selected by
a `species` config key. Each dataset directory has a `config.py` declaring
`species = 'dog'` (the default, for backward compatibility) or `species = 'cat'`,
which selects:

| | dog | cat |
|---|---|---|
| chromosomes | `1-38` + `39`=X | `1-18` + `19`=X (numeric) |
| assembly | canFam4 | F.catus_Fca126_mat1.0 (GCF_018350175.1) |
| gene BED | `UU_CFAM_GSD_1.0_rosy.refseq.ensformat.bed` | `cat_genes.bed` |
| ref FASTAs | `reference-canFam4-chrom-{chrom}.fa` | `reference-F.catus_Fca126_mat1.0-chrom-{chrom}.fa` |
| recombination map | Kidd fine-scale canFam4 map | none (track omitted) |
| UCSC browser | `canFam4` | GenArk hub `hub_6476843_GCF_018350175.1` (cytogenetic chrom names mapped server-side) |

**The folder is the dataset.** One installed `pheweb`; many data folders, each
declaring its own species, branding, and thresholds:

```
/data/
├── dog_compulsive_disorder_pheweb/  config.py: species='dog' + pheno-list.json, gene BED, recomb map, custom_templates/
├── cat_study_A/                     config.py: species='cat' + pheno-list.json, cat_genes.bed
└── cat_study_B/                     config.py: species='cat' + pheno-list.json (may share cat genes/FASTAs)
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

The pheweb source is **baked into the image**, so after editing `pheweb/` you
must re-run `./run_local.sh build` before `serve` picks the change up (only the
last two layers rebuild — seconds).

On a cluster, use a native conda install instead (see "Building a PheWeb for a
new non-human dataset" below). One trap worth knowing there: `pip install -e .`
records the path of the directory it was run in *inside the environment*, so
unpacking the code somewhere new does not change which copy gets imported, and a
console script does not put the working directory on `sys.path`. Export
`PYTHONPATH=/path/to/repo` and check it:

```bash
python -c "import pheweb, os; print(os.path.dirname(pheweb.__file__))"
```

---

# Building a PheWeb for a new non-human dataset

### 1. Install

On a cluster, in a Python 3.8 environment:

```bash
pip install --no-deps -e .     # from a clone of this repo
```

On a Mac, skip this and use `./run_local.sh` (above).

### 2. Create the dataset directory and `config.py`

Copy [`example_dog_datadir/`](example_dog_datadir/) or
[`example_cat_datadir/`](example_cat_datadir/). The keys that matter:

```python
species = 'dog'                    # selects the species profile; 'cat' also supported
genes_bed = 'cf4_genes_noLOC.bed'  # overrides the species default (optional)
significance_threshold = 4e-7      # default is 5e-8 (see below)
site_title = "Darwin's Ark Dog Compulsive Disorder PheWeb"   # browser title / navbar
display_name = "Darwin's Ark"      # short label
num_procs = 8
```

Anything that differs **between species** belongs in `pheweb/species.py`;
anything that differs **between studies** (thresholds, branding, phenotype list)
belongs here.

### 3. Prepare the input files

Three standalone scripts do the preparation. Each needs only Python + pysam, so
they run fine in the cluster conda env.

**Summary statistics — [`normalize_sumstats.py`](normalize_sumstats.py).**
PheWeb expects `chrom,pos,ref,alt,pval,beta,sebeta,af`, sorted by chromosome (in
species order) then position, and it **never checks `ref` against the reference
genome** — upstream treats that as the user's obligation. GCTA MLMA output gives
`A1` (the *effect* allele) and `A2`, which is **not** the same as alt/ref:
assigning `ref=A2, alt=A1` mislabelled 16.2% of variants in an earlier build of
this site. This script anchors each variant against the reference FASTA and flips
the sign of `beta` whenever it has to swap the alleles:

```bash
python normalize_sumstats.py --species dog \
    UU_Cfam_GSD_1.0_ROSY.fa  CCDF1_sumstats.mlma  data/CCDF1_sumstats.norm.csv.gz
```

It auto-detects MLMA vs. the ITEM survey format, spools per chromosome so memory
stays bounded, and emits chromosomes in species order.

**Gene annotation — [`prep_gene_bed.py`](prep_gene_bed.py).** Gene downloads are
disabled in this fork (upstream fetches a human BED), so the gene BED ships in
the data dir. `--drop-loc` removes uncharacterized `LOC*` genes so the nearest-gene
labels name real genes:

```bash
python prep_gene_bed.py --species dog --drop-loc genes_table.txt cf4_genes_noLOC.bed
```

Two load steps ignore this file's timestamp. After changing it, delete
`generated-by-pheweb/best-phenos-by-gene.sqlite3` and
`generated-by-pheweb/resources/gene_aliases-v*.sqlite3` by hand.

**Recombination map — [`prep_recomb_map.py`](prep_recomb_map.py)** (optional,
dog only). Converts per-chromosome map files into one bgzipped, tabixed TSV in
the data dir:

```bash
python prep_recomb_map.py --species dog kidd_maps/ dog_compulsive_disorder_pheweb/
```

The region view shows the track only if both the `.gz` and its `.tbi` are
present; otherwise the right-hand axis is omitted entirely. The map file to use
for a species is named in `pheweb/species.py` (`recomb_map`).

### 4. Make `pheno-list.json`

In the data directory, one object per phenotype:

```json
[
    {
        "assoc_files": ["data/CCDF1_sumstats.norm.csv.gz"],
        "phenocode": "CCDF1",
        "phenostring": "Compulsive behaviour factor 1",
        "category": "Factors",
        "num_samples": 2414
    },
    {
        "assoc_files": ["data/ITEM145_sumstats.norm.csv.gz"],
        "phenocode": "item145",
        "phenostring": "item145: chases tail",
        "category": "Survey items",
        "num_samples": 2414
    }
]
```

`assoc_files` and `phenocode` (`[A-Za-z0-9_~-]`) are required. Optional:
`phenostring` (shown in tables, tooltips, page headers), `category` (groups
phenotypes in the PheWAS plot), and `num_cases` / `num_controls` /
`num_samples`. Beware: a GCTA-style `N` column is often an **allele** count (2N),
not a sample count — check it against the paper rather than importing it blindly.

Verify before spending a cluster job:

```bash
PHEWEB_DATADIR=/path/to/dataset pheweb phenolist verify
```

### 5. Process

```bash
cd /path/to/dataset && pheweb process
```

The log prints `TRYING` / `Downloading ... rsids-...`. **Nothing is fetched** —
this fork has the rsID download commented out and writes an empty file in its
place. (Do not set `disallow_downloads = True`; the permission check runs before
the no-op and raises.)

To distribute jobs across a cluster, see
[these instructions](etc/detailed-loading-instructions.md#distributing-jobs-across-a-cluster).
`run_process_pheweb_unity.sbatch` is a working Slurm example; set `num_procs` in
`config.py` to match the CPUs you request.

### 6. Serve

```bash
pheweb serve --open           # or ./run_local.sh serve <data_dir> 8000
```

### 7. Brand the site

Per-dataset, without touching the code:

- `site_title` / `display_name` in `config.py` — page titles and navbar.
- `custom_templates/about/content.html` in the data directory — the dataset's
  own About page (study description, citations, links). Each dataset gets its
  own; nothing study-specific lives in the repo templates.

`/` redirects to `/phenotypes`; there is no separate landing page.

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
Manhattan JSON, so it keeps the old value until you re-run the load steps.

It must be stricter than `manhattan_peak_pval_threshold` (default 1e-6), since
peaks have to be found before their significant variants can be counted. Setting
it looser raises a clear error at startup rather than an unexplained assertion
mid-load. Set `manhattan_peak_variant_counting_pval_threshold` explicitly only if
you deliberately want the counting threshold to differ from the line.

### Cat gene BED

Cat needs a 5-column BED (`chrom start end gene_name ensg`). The RefSeq
annotation is 4-column, so duplicate the gene name into the 5th column (no GTF
conversion; `make_bed.sh` stays dog-only):

```bash
awk -v OFS='\t' '{print $1,$2,$3,$4,$4}' \
    GCF_018350175.1_F.catus_Fca126_mat1.0_genes.bed > cat_genes.bed
```

### Out of scope

Local LD and GWAS-catalog tracks for non-human assemblies are not provided (the
region view's external human LocusZoom services return nothing useful for dog or
cat); and cat rsid/dbSNP annotation is stubbed out (as it already is for dog).

The recombination track *is* now local for dog — Jeffrey Kidd's fine-scale
canFam4 map, built by `prep_recomb_map.py`. Where no species-appropriate map is
configured (currently cat), the track is omitted rather than falling back to the
human one.

The governing rule: **never fall back to human data for a non-human species.** A
missing track is self-evidently missing; a human track silently mislabelled as
dog is not.

---

# Deploying a dataset as a container

[`docker/Dockerfile.serve`](docker/Dockerfile.serve) builds a self-contained,
serve-only image: the code plus one processed dataset, with no bind mounts. It
targets [SciLifeLab Serve](https://serve.scilifelab.se/), whose contract requires
the startup script to be named exactly `start-script.sh` and be the `ENTRYPOINT`.

```bash
docker build --platform linux/amd64 \
  -f docker/Dockerfile.serve \
  --build-arg DATASET=dog_compulsive_disorder_pheweb \
  -t <user>/dog-compulsive-disorder-pheweb:<tag> .

# check it locally before pushing (an arm64 Mac needs --platform)
docker run --rm --platform linux/amd64 -p 8000:8000 \
  <user>/dog-compulsive-disorder-pheweb:<tag>

docker push <user>/dog-compulsive-disorder-pheweb:<tag>
```

Notes worth knowing before you build:

- The image carries the dataset's `generated-by-pheweb/` output, so the host
  needs no data volume — but that makes it multi-GB. Keep
  [`docker/Dockerfile.serve.dockerignore`](docker/Dockerfile.serve.dockerignore)
  accurate; a missing `**/` prefix once turned a 250 MB image into 32 GB.
- `config.py`, `pheno-list.json`, and `custom_templates/` are copied in a **final
  layer**, after the big data layer, so editing the About page rebuilds in
  seconds rather than re-copying gigabytes.
- [`docker/start-script.sh`](docker/start-script.sh) validates `PORT` (3000–9999),
  requires `PHEWEB_DATADIR`, checks that `pheno-list.json` and
  `generated-by-pheweb/` are present, and then `exec`s `pheweb serve`.
- The image is amd64 by design. On an Apple-Silicon Mac, always pass
  `--platform linux/amd64` or the pull fails with "no matching manifest".

---

# Modifying this fork

```bash
git clone https://github.com/karlssonlab/darwins_ark_nonhuman_pheweb.git
cd darwins_ark_nonhuman_pheweb
python3.8 -m venv .venv && source .venv/bin/activate   # 3.11 does not work
pip install wheel && pip install -e .
pip install pytest && python -m pytest
```

On an Apple-Silicon Mac the native install will not work at all (no arm64 wheels
for the pinned deps) — use `./run_local.sh test` instead. Two test failures
(`test_all`, `test_detectref`) are **pre-existing** and download-related.

Read [`CLAUDE.md`](CLAUDE.md) first: it collects the gotchas that have each cost
real time (stale images, the `/app` mount that breaks the console script,
`PHEWEB_DATADIR` precedence, the gene-BED caches that never invalidate).

### Upstream options that still apply

To run pheweb through systemd, see the sample file [here](etc/pheweb.service).
To use Apache2 or Nginx, see [these instructions](etc/detailed-webserver-instructions.md#using-apache2-or-nginx).
To require login via OAuth, see [these instructions](etc/detailed-webserver-instructions.md#using-oauth).
To reduce storage use, see [these instructions](etc/detailed-webserver-instructions.md#reducing-storage-use).
To customize page contents, see [these instructions](etc/detailed-webserver-instructions.md#customizing-page-contents).

To hide the button for downloading summary stats, add `download_pheno_sumstats = "secret"` and `SECRET_KEY = "your random string"` in `config.py`. That makes a secret page (printed to the console when you start the server) to share summary stats.
To hide the buttons for downloading top hits and phenotypes, add `download_top_hits = "hide"` and `download_phenotypes = "hide"` respectively.

To allow dynamically filtering the manhattan plot, run `pheweb best-of-pheno` and set `show_manhattan_filter_button = True` in `config.py`.

---

# How to cite

This is a modification of the original PheWeb, which should be cited:

> Gagliano Taliun, S.A., VandeHaar, P. et al. Exploring and visualizing
> large-scale genetic associations by using PheWeb. *Nat Genet* 52, 550–552 (2020).

The dog recombination track in the region view is not ours; cite it separately if
you use it:

> Kidd JM. Fine-scale recombination rates inferred using the canFam4 assembly are
> strongly correlated with previous maps of dog recombination. *Mamm Genome*.
> 2025 Dec 12;37(1):12. doi:10.1007/s00335-025-10178-0. PMID: 41387639.
> PMCID: PMC12701031.

The 4e-7 threshold used by the dog instance comes from:

> Pettersson ME, Tengvall K, Meadows JRS. Towards a Standard Threshold for
> Genome Wide Significance in Dogs. *Animal Genetics*. 2026;57(5).
> doi:10.1002/age.70208.
