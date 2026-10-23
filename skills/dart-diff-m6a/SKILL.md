---
name: dart-diff-m6a
description: Differential m6A from DART-seq on an SGE cluster - builds group-level site sets from the Bullseye per-sample beds (sites in >= 2 replicates, mean editing, optional RAC motif, minus dbSNP C>T/G>A SNVs), extracts edited/total read counts from the parseBAM matrices and tests every site and gene between two conditions with a quasi-binomial GLM (Wald) + IHW, Meyer-lab method. Use when the user asks for differential m6A / differential editing / which m6A sites change between conditions, group-level DART-seq site lists, or the dbSNP filter, after the bullseye-analysis step.
---

# Differential m6A (DART-seq) (SGE)

Generalised from project scripts that had been run and audited; comparisons, groups and replicate numbers now come from two tables, and IHW gets a fixed seed. Re-running on the original inputs reproduced the original results (repo `docs/validation.md`).
Method sources: site filtering Meyer 2019 / Tegowski 2024 (Nat Neurosci); differential test Flamand 2022 (GLM + Wald + IHW). Our own choices are marked below.

## Steps (`run_diff_m6a.sh` submits a hold_jid chain)
1. **01_candidates.py** (per group = condition): from the beds of the EDITED samples, keep sites called in >= `MIN_REPS` samples (2 of 3) with mean editing >= `MIN_MEAN` (5 %); flag RAC (R-A-C in transcript orientation, C = edited base, from the fasta). -> `sites/group_sites.tsv`, `per_sample_calls.tsv`, `snp_candidates.tsv`.
2. **02_dbsnp.job.sh**: a candidate is a SNP if the dbSNP bed has the same position with the matching change (`+` site <-> C>T, `-` site <-> G>A). -> `snp_hits.tsv`. Streams the whole bed (human: 8.7 GB, ~1 h).
3. **03_matrix.py** per comparison: per condition reps -> mean -> [RAC] -> no dbSNP, union of the two conditions, then `tabix -R` on each sample's parseBAM matrix: `mut` = T (+) / A (-), `cov` = total. A sample that did not call the site still gets its counts (0/0 if no row). -> `sites/<cmp>/sites_kept.tsv`, `counts.tsv`, `sites/filter_counts.tsv`.
4. **04_diff_m6a.R**, two versions `norac` / `rac`: per site, samples with cov >= `MIN_COV` (20), >= `MIN_PER_GROUP` (2) per group, `glm(mut/cov ~ condition, weights = cov, quasibinomial)`, control = reference, Wald P (t, as `summary.glm` reports; `wald_z_p` for comparison), **IHW with mean coverage as covariate (our choice)**, `set.seed(SEED)` before IHW. Per gene: sum mut/cov over the gene's kept sites, same GLM, editing score = sum(mut)/sum(cov) x n_sites, `GENE_FC` (1.25) + adj.P. Joined to the DESeq2 table (`expr_log2FC`, `expr_padj`) if given. -> `diff_m6a/<ver>/<cmp>/per_site_results.tsv`, `per_gene_results.tsv`, `volcano_{sites,genes}.{png,pdf}`, `diff_m6a/summary_counts.tsv`. Direction labels: `lost_in_treated` / `gained_in_treated` / `ns`.

## Inputs
- **Beds**: Bullseye `Find_edit_site` output per EDITED sample, `<bed_dir>/<sample>.bed` (bullseye-analysis, pairs edited vs control). Use the **stranded** beds + **stranded** matrices together (parseBAM `--stranded`); unstranded matrices work too (one row serves both strands).
- **Matrices**: `<matrix_dir>/<sample>.matrix.gz(.tbi)` of the same edited samples.
- **samplesheet.tsv**: `sample  group` (edited samples only; order inside a group = replicate order).
- **comparisons.tsv**: `name  treated  control  [de_name]` (`de_name` = file name of the deseq2-expression table if it differs from `name`).
- **`FASTA`** (with `.fai`) and **`DBSNP`** (env, required): the same genome/release as the BAMs. Make the dbSNP bed once per release: `OUT=<ref>/dbsnp_ensembl111 RELEASE=111 bash download_dbsnp_ensembl.sh` (Ensembl variation VCFs, C>T/G>A SNVs; needs internet; hours).

## Decisions to ask the user
- Which sample is edited vs control was settled in bullseye-analysis (DART-seq: active APOBEC1-YTH = edited; the YTH-mutant = control). Here only EDITED samples go in the sheet.
- Comparisons and which condition is the control (direction of `lost` / `gained`).
- Thresholds (`MIN_REPS`, `MIN_MEAN`, `MIN_COV`, `MIN_PER_GROUP`, `GENE_FC`), RAC version to report (both are always made), and the Bullseye thresholds used upstream.

## Known limitations (the Meyer-lab method is kept on purpose; report these with the results)
- **Wald fails when one group has 0 edited reads** (Hauck-Donner): for example, a hypothetical site edited ~50 % in one group and 0 % in the other gets Wald P ≈ 1, while an F / likelihood-ratio test would call it highly significant; such sites are essentially never called, even when the change is large. A likelihood-ratio / F test would catch them; not switched on purpose.
- Number of significant sites is sensitive to the multiple-testing method (IHW vs BH) and to the low quasi-binomial dispersion; no random effect for replicate.
- The input beds carry the Bullseye `Find_edit_site.pl` stranded-merge issue (a few percent of stranded sites lost, see repo `docs/bullseye_known_issues.md`).

**Software:** the job scripts activate the repo's conda environments when `ENV_SETUP` is set (see `docs/INSTALL.md`; pass it with the other env vars), otherwise they use the cluster's `module load` lines. Ask the user which applies before the first run.

## Rules for every run (from the user)
1. Copy `scripts/` into `<project>/diff_m6a/scripts/` + `MANIFEST.md5`; write `run_diff_m6a_<YYYYMMDD>.sh` calling the COPY with FASTA/DBSNP/thresholds; keep the sheets.
2. Keep all logs (`<out>/error_file/0{1..4}_*`). 3. Failed job: clear cause -> delete its outputs/logs and restart from that step with `START=<n>`; unclear -> report and wait.
4. Show the sheets (check_sheets output), scripts and the command, wait for approval. 5. `command rm -f`. 6. Update STATUS.md.

## Typical commands
```bash
P=<project>/diff_m6a; SK=$P/scripts; R=<ref>
FASTA=$R/grch38.fa DBSNP=$R/dbsnp_ensembl111/dbsnp_ensembl111.CT_GA.snv.bed \
  bash $SK/run_diff_m6a.sh <project>/bullseye/beds_stranded <project>/bullseye/matrix_stranded \
       $P/samplesheet.tsv $P/comparisons.tsv $P <project>/diff_expr/results/tables
START=4 FASTA=... DBSNP=... bash $SK/run_diff_m6a.sh ...    # rerun only the statistics
```

## Checks
- `error_file/01_output.txt`: per group union / >= reps / mean counts; `filter_counts.tsv` decreases monotonically along the steps; dbSNP hits are a few % of candidates.
- `sites_kept.tsv` rows = `union_of_conditions` (norac) of `filter_counts.tsv`; `counts.tsv` rows = sites x samples of the two groups.
- Spot-check one site: its `mut`/`cov` in `counts.tsv` equals the Bullseye bed (`mut=` in column 4, `dart_total` column 10) for a sample that called it .
- `summary_counts.tsv`: `adj_method` must be IHW (BH only if IHW failed); compare the numbers lost vs gained.

## Gotchas
- Compute-node Python may be old (3.6: no `capture_output`); the scripts avoid newer features; bcftools/tabix need `module load htslib/1.18` first.
- qsub runs a spooled copy: scripts are passed by absolute path, never found through `$0` inside a job.
- IHW needs >= ~1,000 tests per comparison to form covariate bins; with fewer it falls back to BH and the column says so.
- Give the DESeq2 tables dir even if it does not exist yet (chain waiting on `HOLD_JID`): `run_diff_m6a.sh` resolves it with `readlink -m` and prints `DE tables: <dir>`; an earlier version used `readlink -f`, which returned an empty path for a not-yet-existing dir and silently left the expression columns empty. Check that `expr_log2FC` is filled in `per_gene_results.tsv`.
- Output column names: `n_ctl` / `n_trt`, `f_ctl` / `f_trt`, `score_ctl` / `score_trt` in the gene table; per-site tables name them by group (`edit_<group>`, `n_samples_<group>`).
