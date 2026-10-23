---
name: m6a-metagene
description: Metagene (transcript-wide 5'UTR / CDS / 3'UTR) distribution plots of DART-seq / m6A sites with the Guitar R package on an SGE cluster, one overlay per site version (no-RAC and RAC filter) with all groups, built from the dart-diff-m6a group-level site tables. Use when the user asks for a Guitar plot, metagene / distribution of m6A sites along transcripts, stop-codon enrichment, after the dart-diff-m6a step.
---

# Guitar metagene plots (SGE)

Generalised from project scripts (GuitarPlot settings as in our earlier DART-seq distribution plots); group names now come from the site table. Re-running on the original inputs reproduced the original beds and plot data (repo `docs/validation.md`).

## What it does (`run_guitar.sh`)
1. **01_make_beds.py** (runs on the login node, seconds): per group, the sites of `group_sites.tsv` minus `snp_hits.tsv`, all sites (`<group>_norac.bed`) and RAC only (`<group>_rac.bed`), BED6 sorted. The set sizes are checked against `filter_counts.tsv` (`no_dbSNP`) of dart-diff-m6a, so the plot shows exactly the sites that were tested; a group outside every comparison is written but `not checked`.
2. **02_guitar.R** (job `guitar`, 4 slots, ~1 h): TxDb from the GTF (`makeTxDbFromGFF`), `GuitarPlot(headOrtail = TRUE, enableCI = FALSE, mapFilterTranscript = TRUE, pltTxType = "mrna")` with all groups, per version -> `results/guitar_<ver>_overlay.{png,pdf}`, Guitar's default pdf, `guitar_<ver>_plot_data.tsv` (densities; use it for subset plots), `guitar_<ver>_plot.rds`.

## Inputs
- `<diff_m6a>/sites/` (dart-diff-m6a steps 01-03: `group_sites.tsv`, `snp_hits.tsv`, `filter_counts.tsv`).
- **`GTF`** (env, required): the same annotation as the STAR index / refFlat (e.g. Ensembl 111). `ORGANISM` default "Homo sapiens".

**Software:** the job scripts activate the repo's conda environments when `ENV_SETUP` is set (see `docs/INSTALL.md`; pass it with the other env vars), otherwise they use the cluster's `module load` lines. Ask the user which applies before the first run.

## Rules for every run (from the user)
Copy `scripts/` into `<project>/guitar/scripts/` + `MANIFEST.md5`; `run_guitar_<YYYYMMDD>.sh` calls the COPY; keep logs (`error_file/`, `01_make_beds.log`); show the command and wait for approval; update STATUS.md.

## Typical command
```bash
GTF=<ref>/Homo_sapiens.GRCh38.111.gtf bash <project>/guitar/scripts/run_guitar.sh <project>/diff_m6a/sites <project>/guitar
```

## Checks
- `01_make_beds.log`: every checked set `OK`.
- DART-seq m6A sites should be enriched near the stop codon / start of the 3'UTR; a flat or 5'-biased profile suggests strand or coordinate problems (BED start must be pos-1).

## Gotchas
- `GuitarPlot` with `miscOutFilePrefix` only writes a pdf and returns NULL: not used.
- `txGuitarTxdb` expects a file, not the `makeGuitarTxdb` object: the TxDb is passed and Guitar builds its coordinates (slow but once per call).
- Guitar's ggplot layers refer to its own vectors: facets and `%+%` do not work; draw per-subset plots from `*_plot_data.tsv`.
