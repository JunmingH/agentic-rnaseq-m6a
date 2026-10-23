---
name: deseq2-expression
description: Differential expression (DESeq2) and standard RNA-seq PCA on an SGE cluster from an RSEM gene expected_count matrix, driven by a samplesheet (sample, group, subset) and a comparisons table (treated vs control, optionally within a subset such as a cell line); writes DE tables, volcano plots, an xlsx, PCA/scree/distance plots. Use when the user asks for differential expression, DESeq2, volcano plots, PCA or sample clustering of RNA-level counts, usually after the rsem-quantification step (e.g. the control samples of a DART-seq project).
---

# DESeq2 differential expression + PCA (SGE)

Generalised from project scripts that had been run and audited (size factors recomputed independently, fold-change directions checked). Groups and comparisons come from two tables; re-running on the original inputs reproduced the original results (repo `docs/validation.md`).

## What it does
1. **`diff_expr.R`** (job `diff_expr`, 4 slots): per comparison, `DESeqDataSetFromMatrix(design = ~ group)` on the samples of its subset, control group = reference level, `DESeq(betaPrior = FALSE)`, `results(name = "group_<treated>_vs_<control>", alpha = ALPHA)`, plus `lfcShrink(type = "apeglm")` as an extra column. One fit is shared by comparisons with the same subset and control (e.g. several treatments vs the same control). Gene filter per fit: >= `MIN_COUNT` counts in at least as many samples as the smallest group. Also a quick PCA of all samples (`qc/`).
   Output `results/tables/<name>.csv` (gene_id, symbol, baseMean, log2FoldChange, lfcSE, log2FC_shrunk_apeglm, pvalue, padj, mean normalised counts per group; sorted by padj), `summary_counts.csv`, `DE_all_comparisons.xlsx`, `volcano/volcano_<name>.{png,pdf}` (top 20 labelled), `sessionInfo.txt`.
2. **`pca_rnaseq.R`** (job `pca_rnaseq`): genes with >= `MIN_COUNT` in >= `MIN_SAMPLES` samples -> size factors -> `vst(blind = TRUE)` -> `prcomp` (centred, not scaled) on the `NTOP` most variable genes and on all genes; PC1-3 plots, scree, top loadings, sample-distance heatmap, `vst_matrix.csv`, `size_factors.csv`.

## Inputs
- **Matrix**: `matrix.genes.expected_count.tsv` from `rsem-quantification/make_matrix.py` (column `id` = `ENSG..._SYMBOL` with `--append-names`, one column per sample). Rounded to integers (RSEM expected counts are not integers; same as the original analysis; tximport is the alternative).
- **samplesheet.tsv** (tab): `sample  group  [subset]`. `sample` = matrix column; `subset` = e.g. cell line; group order = plot order.
- **comparisons.tsv** (tab): `name  treated  control  [subset  label]`. `subset` empty/`all` = all samples in one model. `name` is the output file name; dart-diff-m6a can join on it (`de_name`).
- `check_sheets.py <matrix> <samplesheet> <comparisons>` runs automatically before submitting and prints every comparison's samples: **show that list to the user** (who is control, who is treated, which samples).

## Decisions to ask the user (do not guess)
- **Which samples go in** (DART-seq: usually only the control samples have the transcriptome BAM and give the RNA level). Which group is the control of each comparison (the reference level: log2FC > 0 = higher in treated).
- **One model per subset or all samples together.** Per subset (e.g. cell line) = each subset has its own dispersion estimates; all together = more df but assumes equal dispersion; an interaction model (line x treatment) is not built in.
- Thresholds: `ALPHA` (0.05), `LFC` (1; only for the volcano colours and the up/down counts, not for padj), `MIN_COUNT` (10).

**Software:** the job scripts activate the repo's conda environments when `ENV_SETUP` is set (see `docs/INSTALL.md`; pass it with the other env vars), otherwise they use the cluster's `module load` lines. Ask the user which applies before the first run.

## Rules for every run (from the user)
1. **Save the code that actually runs in the project**: copy `scripts/` to `<project>/<de_dir>/scripts/`, `md5sum` into `MANIFEST.md5`, write `<de_dir>/run_deseq2_<YYYYMMDD>.sh` that calls the COPY with the sheets; keep the sheets next to it.
2. **Keep all logs** (`<out>/error_file/`); never delete logs of a successful step.
3. Failed job with a clear cause: delete only its outputs and logs and rerun; unclear: keep logs, report, wait.
4. **Show the sheets (check_sheets output), scripts and the submit command, and wait for approval.**
5. `command rm -f` (rm may be aliased to interactive). 6. Update the project STATUS.md.

## Typical commands
```bash
P=<project>/diff_expr; SK=$P/scripts
bash $SK/run_deseq2.sh <project>/rsem/matrix.genes.expected_count.tsv $P/samplesheet.tsv $P/comparisons.tsv $P
```
Env: `SGE_PROJECT` (qsub -P, optional), `ALPHA LFC MIN_COUNT MIN_SAMPLES NTOP`, `HOLD_JID`.

## Checks
- `error_file/diff_expr_output.txt`: genes kept per fit and the tested/significant counts per comparison; `sessionInfo.txt` (we used R 4.2.1, DESeq2 1.36.0).
- Direction: for a top gene, `mean_norm_<treated>` > `mean_norm_<control>` must go with log2FoldChange > 0.
- PCA: replicates of a group cluster; an outlier or a sample sitting with another group = possible label swap (check before trusting the DE).
- Very few genes tested or all padj NA = wrong matrix columns or a group with one sample.

## Gotchas
- Group names must be valid R names for the coefficient (`make.names`); avoid names starting with a digit or containing `-`.
- `lfcShrink(apeglm)` needs the `apeglm` package (Bioconductor).
- The PCA in `results/qc/` (>= 10 counts in >= 3 samples, all genes) and `pca/` (top 500) differ on purpose; report `pca/pca_PC1_PC2_top500` as the standard plot.
