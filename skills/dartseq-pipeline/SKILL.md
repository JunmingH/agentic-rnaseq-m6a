---
name: dartseq-pipeline
description: Run or resume the complete DART-seq m6A + RNA-seq analysis on an SGE cluster from raw paired-end fastq to differential expression, differential m6A sites/genes and metagene plots, by chaining the step skills (fastqc-qc, read-trimming, star-alignment, rsem-quantification, bullseye-analysis, deseq2-expression, dart-diff-m6a, m6a-metagene) with a pipeline-audit check after each step. Use when the user asks to run the whole DART-seq pipeline, start a new DART-seq project, asks "what is the next step", or wants to reproduce a previous DART-seq analysis.
---

# DART-seq pipeline (orchestration)

This skill only orchestrates. Each step is done by its own skill: load it, follow its rules, and use its scripts. The order and the hand-over files are fixed. Every decision is made by the user, and every step is checked before the next one starts.

```
fastq ─▶ 1 fastqc-qc ─▶ 2 read-trimming ─▶ 3 star-alignment ─┬─▶ 4 rsem-quantification ─▶ 6 deseq2-expression ──────────┐
                                                              │      (control/Mut samples: transcriptome BAM)          │ DE tables
                                                              └─▶ 5 bullseye-analysis ─▶ 7 dart-diff-m6a ◀─────────────┘
                                                                     (edited vs control)      └─▶ 8 m6a-metagene
              9 pipeline-audit: after every step (code, parameters, sample placement, outputs)
```

| step | skill | input | hand-over output (used by) |
|---|---|---|---|
| 1 | fastqc-qc | raw `*_1/_2.fq.gz` | MultiQC report (decide trimming) |
| 2 | read-trimming | raw fastq | `<s>_paired_R{1,2}.fastq.gz`, `trimm_summary.tsv` (3) |
| 3 | star-alignment | trimmed fastq | `<s>.markdup.bam(.bai)` (5), `<s>_Aligned.toTranscriptome.out.bam` (4), `align_summary.tsv` |
| 4 | rsem-quantification | transcriptome BAMs | `matrix.genes.expected_count.tsv` (6) |
| 5 | bullseye-analysis | genome BAMs | `matrix*/<s>.matrix.gz` + `beds*/<edited>.bed` (7) |
| 6 | deseq2-expression | count matrix + sheets | `results/tables/<cmp>.csv` (7), PCA |
| 7 | dart-diff-m6a | beds + matrices + sheets (+ DE tables) | `sites/` (8), per-site/per-gene results |
| 8 | m6a-metagene | `sites/` | Guitar plots |

## Before step 1: settle with the user (write the answers into STATUS.md)
1. **Sample meaning.** Decide which label is the active DART enzyme (edited) and which is the control (YTH-mutant or APOBEC1-only). Ask whether replicate numbers are matched pairs or only replicate numbers, and what each group/condition means. Labels such as WT/Mut are not self-explanatory ("WT" may be the active enzyme and "Mut" the control); confirm them with the bullseye-analysis motif check, not from the file names.
2. **Library strandedness.** Test it in step 4 on one sample (`none` / `reverse` / `forward`); do not take it from the user's wording. dUTP (ISR) = RSEM `reverse`, which matches `parseBAM.pl --stranded`.
3. **Reference.** Use one Ensembl release for everything: STAR index, RSEM reference, refFlat (patched, see `patches/`), dbSNP, and the GTF for Guitar.
4. **Comparisons.** Settle the treated vs control pairs, and whether the samples are modelled per subset (for example per cell line) or all together. Write `expr_samplesheet.tsv` / `expr_comparisons.tsv` (control samples, RNA level) and `m6a_samplesheet.tsv` / `m6a_comparisons.tsv` (edited samples); the formats are in `templates/` of the repo.

## Software
Either the cluster's modules (the `module load` fallback in the job scripts) or the three conda environments of the repo (`envs/`, installed once; `ENV_SETUP` = absolute path of a copy of `envs/env_setup.example.sh`, exported before every run script). Settle this, `SGE_PROJECT` and the reference paths before step 1; `docs/INSTALL.md` has the commands, including how to build every reference file.

## Project layout (one folder per step, each with its own `scripts/` copy, `MANIFEST.md5`, `run_<step>_<date>.sh`, `error_file/`)
```
<project>/STATUS.md  samples.txt  sheets/
  fastqc/  trimm_fastp/  bam/  rsem/  bullseye/  diff_expr/  diff_m6a/  guitar/
```

## For every step
1. Load the step skill and copy its `scripts/` into the step folder (do not run the live skill). Then write the dated run script, which sets every REQUIRED option explicitly.
2. Show the user the sample list or pairs (`check_sheets.py` / `make_pairs.py` output), the parameters and the command, then **wait for approval**.
3. Test ONE sample or pair first, check its output, runtime and memory (`qacct -j`), then submit the rest. Chain jobs with `HOLD_JID` only when the inputs will exist.
4. When a step finishes, run **pipeline-audit** on it (code + logs + outputs). Report problems with evidence and let the user decide. Do not silently change a method that the user has fixed.
5. Update STATUS.md: what ran, where, parameters, key numbers, open questions.

## Decision points that change results (always ask; recommend, do not decide)
- Trimming tool and settings: is there adapter or poly-G content in the MultiQC report?
- STAR `RM_DUP` (mark/remove) and `TRANSCRIPTOME`. DART-seq without UMIs: mark duplicates, and let Bullseye ignore them with `DEDUP=1`.
- RSEM `STRANDEDNESS`, after the one-sample test.
- Bullseye: pairing mode (`same_n` / `all`), `DEDUP`, stranded or unstranded matrices, thresholds (Meyer lab: 10–95 %, ≥ 1.5×, ≥ 2 C-to-U reads, coverage ≥ 10).
- Group-level sites: ≥ 2 of 3 replicates, mean editing ≥ 5 %, RAC version, dbSNP release.
- Statistics: DESeq2 per subset or all samples together; diff m6A cutoffs `MIN_COV` 20, `MIN_PER_GROUP` 2, `ALPHA` 0.05, `GENE_FC` 1.25.
- Known method limits, which have to be reported with the results: Bullseye stranded merge loses a few percent of sites; the Wald test is uninformative when one group has 0 edited reads. See `docs/bullseye_known_issues.md` and the dart-diff-m6a skill.

## Resuming ("continue" / "what next")
Read the project STATUS.md first, then check which hand-over files exist (table above). Next, run `qstat -u $USER` and look at the `error_file/` folders. Continue at the first step whose hand-over output is missing or not yet audited.

## Reproducing a previous analysis
Use the `scripts/` copy and `run_*.sh` saved in each step folder of that project; they record the exact code (md5 in `MANIFEST.md5`) and parameters. Do not use the current skill version. When the newer skills are run on the same inputs, compare the outputs with the old ones (see the repo `docs/validation.md`).
