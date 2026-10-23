# Validation of the downstream skills

## Scope and verification status (2026-10-23)

The project-recorded comparisons below cover three downstream skills and the specifically listed environment checks. They do not establish that all ten skills, the complete pipeline, or the agent's scientific judgments have been independently validated.

A review found that the original `docs/compare_tables.py` could report `IDENTICAL` when an output column was missing or a finite value changed to `NaN`. The comparator has been corrected, with regression tests on synthetic tables. **The original project output tables are not available in this repository, and comparisons that relied on the affected comparator have not yet been rerun with the correction.** The historical results below are retained as reported; they should not be described as newly reconfirmed by this fix. Independently obtained byte-level comparisons are not invalidated by a defect in this numeric table comparator.

To reconfirm an affected result, compare the saved original and skill-generated output tables using the corrected script, record its exit status and diagnostic output, and investigate every difference. An unchanged result alone does not establish scientific correctness; it establishes agreement with the selected reference output under the comparison criteria.

## Previously reported project comparisons

`deseq2-expression`, `dart-diff-m6a` and `m6a-metagene` were written by generalising project scripts that had already been run and audited on an unpublished DART-seq dataset. It has several conditions with edited and control libraries in replicate, and several treated-vs-control comparisons. In the original scripts, group names, comparisons and replicate numbers were hard-coded. In the skills they come from sample sheets.

To check that generalising changed nothing, the skills were re-run on the same inputs with the same parameters, and every output table was compared with the original one. Tables were matched row by row on their key (gene, site or sample). Numeric columns were compared as relative differences, and all other columns had to match exactly. No data from that project is included in this repository.

| step | outputs compared | result |
|---|---|---|
| `deseq2-expression` | all DE tables (all genes and columns), summary counts, size factors, PCA variance and scores | **identical** when the sample sheet lists samples in the matrix column order. With a different sample order, two-group fits were identical (≤ 1e-14) and a fit with more than two groups differed by ≤ 6e-5 relative in p-values (iterative fitting depends on column order); significance calls and top-gene order were unchanged. PCA: identical up to the arbitrary sign of some components. |
| `dart-diff-m6a` steps 01–03 | group site tables, dbSNP candidates and hits, kept sites and count matrices of all comparisons | **byte-identical**. The filter table is identical apart from one renamed step label. |
| `dart-diff-m6a` step 04 | per-site and per-gene tables (both RAC versions, all comparisons), summary | per-site tables and summary **identical**; per-gene statistics equal within 1e-14. The only differences are the joined expression columns, which inherit the DESeq2 sample-order effect above. |
| `m6a-metagene` | all per-group site BED files, Guitar plot data (both versions) | **byte-identical** |

Output naming changed on purpose: direction labels are `lost_in_treated` / `gained_in_treated`, and the internal column names are `*_ctl` / `*_trt`.

## Bug found by the validation

The first validation run of `dart-diff-m6a` left the expression columns empty. `run_diff_m6a.sh` resolved the DESeq2 tables directory with `readlink -f` at submit time. That directory did not exist yet, because the chain was waiting for the DESeq2 job, so the path came back empty and the step silently ran without it. The script now uses `readlink -m`, prints the directory it will use, and refuses a missing directory unless it is waiting on a job. The re-run gave the results in the table above.

## Conda environments (`envs/`)

The downstream steps and a slice of Bullseye were also run as real cluster jobs with the three conda environments, activated through `ENV_SETUP`. Their outputs were compared with the original module-based runs:

| step | software difference | result |
|---|---|---|
| Bullseye parseBAM `--stranded` + Find_edit_site (one pair, one chromosome) | Perl 5.32 / BioPerl 1.7.8 / samtools 1.19.2 vs Perl 5.26 / BioPerl 1.7.2 / samtools 1.12 | matrices and site bed **byte-identical** |
| `deseq2-expression` | DESeq2 1.38.0 / apeglm 1.20.0 vs 1.36.0 / 1.18.0 | identical gene sets and significance summary; p-values ≤ 6e-5 relative; apeglm-shrunk log2FC ≤ 0.04 absolute (r > 0.99999) |
| `dart-diff-m6a` step 04 | IHW 1.26.0 vs 1.24.0 (R 4.2.3 vs 4.2.1) | GLM coefficients and Wald p-values equal (≤ 1e-10); IHW-adjusted p-values differ slightly for part of the sites, but **no site or gene changes its significance call** and the summary is identical |
| `m6a-metagene` | Guitar 2.14.0 vs 2.12.0 | plot data **byte-identical** |
| `rsem-quantification` (one sample) | RSEM 1.3.3 from bioconda vs the cluster module, 26 vs 10 threads | gene and isoform results and the `.cnt` read counts **byte-identical** |

Every module-based tool also has a conda build of the same version (fastp, Trimmomatic, STAR, samtools, RSEM), except Picard (2.27.5 vs 2.25.2) and FastQC/MultiQC.

## How to repeat this on your own data

1. Run the original analysis and the skill on the same inputs, with every parameter set explicitly.
2. Compare each output table on its key with `python3 docs/compare_tables.py <original> <new> <key> [old:new column renames]`. Preserve the complete output and check the exit status; a reported schema, key, missing-value, or numeric mismatch must not be treated as a successful comparison.
3. For DESeq2, list the samples in the same order as before if exact identity is required.

The CLI retains its default relative tolerance of `1e-6` for finite numeric values. It is not a single acceptance rule for every historical experiment above: some reported environment or sample-order comparisons have larger numerical differences. Such differences require separate, explicitly recorded, step-specific assessment; do not relax a threshold solely to make a comparison pass. For a claim of byte-identical files, use a byte comparison such as `cmp` in addition to any table-level assessment.

Missing or extra columns/keys and unequal values return exit status 1. Duplicate headers or keys, rename collisions, empty or header-only tables, missing key fields, and malformed rows return exit status 2. Both-sided `NaN` and same-sign infinity are considered equal; one-sided `NaN` and other nonfinite mismatches fail. `IDENTICAL` means equality under the documented table-comparison rules, including the existing empty-string/`NA` convention, not byte-for-byte file identity.

## Comparator regression tests

From the repository root, run:

```bash
python3 -m unittest discover -s tests -p 'test_compare_tables.py' -v
```

These tests use synthetic CSV/TSV inputs to check the comparator. They do not rerun the genomics pipeline or verify the historical project results.
