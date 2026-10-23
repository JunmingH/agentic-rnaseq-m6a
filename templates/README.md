# Sample sheets (tab-separated, header required)

| file | used by | columns |
|---|---|---|
| `expr_samplesheet.tsv` | deseq2-expression | `sample` (= column of the RSEM matrix; control/Mut samples, which have the transcriptome BAM), `group`, optional `subset` (e.g. cell line) |
| `expr_comparisons.tsv` | deseq2-expression | `name` (output file), `treated`, `control` (group names), optional `subset` (model only these samples; empty/`all` = all), `label` (plot title) |
| `m6a_samplesheet.tsv` | dart-diff-m6a | `sample` (EDITED sample = Bullseye bed `<sample>.bed` and matrix `<sample>.matrix.gz`), `group`; order inside a group = replicate order |
| `m6a_comparisons.tsv` | dart-diff-m6a | `name`, `treated`, `control`, optional `de_name` (the deseq2-expression comparison of the same groups, joins expression to the m6A genes) |

Check them before running: `python3 skills/deseq2-expression/scripts/check_sheets.py <matrix.tsv | -> <samplesheet> <comparisons>`. It prints the samples of every comparison, control first.
