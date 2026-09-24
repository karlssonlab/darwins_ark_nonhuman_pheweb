# darwins_ark_nonhuman_pheweb

A fork of [statgen/pheweb](https://github.com/statgen/pheweb) (GWAS results
browser) serving **non-human** data — dog and cat. Upstream is parametrized by
human genome build; this fork replaced that with a species abstraction.

**The folder is the dataset.** One installed `pheweb`, many data directories, each
with its own `config.py` declaring its species and settings, each served on its own
port. Adding a dataset means adding a folder, not cloning the repo.

## Commands

Python 3.8 only (pysam pins). Nothing runs natively on an Apple-Silicon Mac —
everything goes through the container.

```bash
./run_local.sh build                      # build/rebuild the image (dap-pheweb:py38)
./run_local.sh process <data_dir>         # pheweb process
./run_local.sh serve   <data_dir> [port]  # pheweb serve (default 8000)
./run_local.sh run     <data_dir> <args>  # any pheweb subcommand
./run_local.sh shell   [data_dir]         # bash in the container
./run_local.sh test                       # pytest
./run_local.sh recomb  <data_dir> <map_dir>   # build the recombination map
```

On a cluster, use the conda env directly (see `docs/unity-*.md`) — the standalone
scripts (`normalize_sumstats.py`, `prep_recomb_map.py`, `prep_gene_bed.py`) need
only Python + pysam and are meant to run there.

## Gotchas that have each cost real time

- **`serve` does not pick up code edits.** The source is baked into the image, and
  `run_local.sh` only checks that *an* image exists. After editing `pheweb/`, run
  `./run_local.sh build` first.
- **Do not mount the repo at `/app`.** `-v "$PWD:/app"` hides the in-image
  egg-info; the `pheweb` console script then dies with
  `PackageNotFoundError: PheWeb`. It works for `python -m pytest` and direct
  `python` calls, but not the CLI.
- **`PHEWEB_DATADIR` is baked into the image** and beats the `data_dir` config key
  in `conf.get_data_dir()`. Tests pointing at a `tmp_path` must
  `monkeypatch.delenv('PHEWEB_DATADIR')` first.
- **`conf.overrides` is empty on a bare `import`.** `config.py` is loaded by the CLI
  entry point (`command_line.py:main`), not on import. To reproduce real behaviour
  in a script, call `conf.load_overrides_from_file(datadir + '/config.py')`.
- **Neither rsIDs nor genes are ever downloaded**, despite log lines that say
  otherwise. This fork commented out both `wget.download` calls (marked `#RB`):
  `download_rsids` writes a **0-byte file** in place of the dbSNP download, and
  `get_genes_for_build` is a no-op. So `generated-by-pheweb` is safe to delete —
  the empty `resources/rsids-v154-hg19.tsv.gz` is recreated automatically — but
  gene BEDs must be placed in the data dir by hand. Expect
  `TRYING` / `Downloading ... rsids-...` in the log; nothing is fetched.
  Setting `disallow_downloads = True` *would* break this: the permission check
  runs before the no-op and raises.
- **`generated-by-pheweb/sites/sites.tsv` is gzipped** despite the `.tsv` name.
- **Two load steps ignore the gene BED's timestamp** and will not rebuild when it
  changes. After changing `genes_bed`, delete both by hand:
  `generated-by-pheweb/best-phenos-by-gene.sqlite3` and
  `generated-by-pheweb/resources/gene_aliases-v*.sqlite3`.
- **PheWeb never validates `ref` against the reference genome.** Upstream states it
  as a user obligation. Assigning `ref=A2/alt=A1` from GCTA output is wrong for
  ~16% of variants — see `normalize_sumstats.py`.
- **Input must be sorted** by chromosome (in species order) then position. MLMA
  files are grouped by chromosome but scrambled (33, 36, 35, …).
- **`beta` means "effect of the alternate allele"** (upstream's definition). Any
  change to which allele is `alt` must flip its sign.
- **This machine:** disk sits near-full (~97%), Docker has ~7.7 GB, and Docker
  **cannot mount `~/Downloads`** (`operation not permitted`). Stage large inputs
  elsewhere and prefer gzipped intermediates — PheWeb reads gzip transparently.

## Where things live

- `pheweb/species.py` — the species registry: chromosomes, gene BED,
  reference-FASTA pattern, assembly label, UCSC link, branding, `recomb_map`.
  **Single source of truth for anything that differs between dog and cat.**
- `pheweb/conf.py` — config accessors. Per-data-dir overrides are applied in
  `get_species_profile()`.
- `pheweb/load/` — the `pheweb process` steps, in the order listed in
  `load/process_assoc_files.py`.
- `pheweb/serve/` — Flask app, templates, and static JS (LocusZoom 0.13.0, loaded
  from CDN, not vendored).
- `docs/` — dated design and implementation notes; read the relevant one before
  changing an area. `*-implementation.md` files record decisions and rationale,
  including rejected alternatives.

## Conventions

- **Species-specific values belong in `species.py`**; study-specific values belong
  in the data dir's `config.py`. A significance threshold is a property of the
  study, not the species; a chromosome list is the reverse.
- **Prefer a config key with a safe default** over editing code per dataset, and
  derive switches from artifacts where possible (e.g. the recombination track is
  driven by whether the map file exists, not by a boolean).
- **Never fall back to human data for a non-human species.** A missing track is
  self-evidently missing; a human track silently mislabelled as dog is not. Several
  LocusZoom sources in `serve/static/region.js` (`catalog`, `ld`, `gene`) still
  point at human services and are known-wrong for dog/cat — see
  `docs/recombination-track-implementation.md`.
- **Document decisions in `docs/`**, dated, with the alternatives that were
  rejected and why. Match the existing files' style.
- When adding a load-time behaviour, say explicitly whether it needs a reprocess.
  Serve-time changes need only a restart; load-time changes are baked into
  `generated-by-pheweb/`.

## Verification expectations

Claims about data should be measured, not assumed — this codebase has already
produced two silent, systematic errors (human recombination rates served as dog;
`ref`/`alt` swapped for 16% of variants) that looked fine on screen.

- `./run_local.sh test` — two failures (`test_all`, `test_detectref`) are
  **pre-existing** and network/reference-download related. Confirm against a
  stashed baseline before blaming a change.
- There is no browser automation here. Rendering claims cannot be verified from
  the CLI — say so rather than implying otherwise.
