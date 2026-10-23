---
name: pipeline-audit
description: Audit (检查 / 复查 / 审查 / 核对 / code review) the code of a bioinformatics pipeline step on an SGE cluster and prove whether it is correct - read every script line by line for bugs, check that every parameter is the right option with the right value and really reaches the tool, that samples sit in the right place (R1/R2, sample lists, control vs treated/edited pairs, group parsing, matrix column order, reference level), that the right reference/annotation is used, then confirm with the tool's logs and outputs and write an audit report. Use this whenever the user asks to check, re-check, review or verify a step, a script or its results (FastQC, trimming, STAR, Picard, RSEM, DESeq2/PCA, Bullseye, differential m6A, Guitar, any qsub job script), asks whether a parameter is used or passed correctly, whether files or samples are in the right position, what a parameter means, or wants confidence before building on a result - even if they do not say "audit".
---

# Pipeline audit (code review + proof)

The main job is **reading the code for bugs**: wrong option or value, arguments in the wrong order, a sample in the wrong slot, a variable that never reaches the job, logic that only works for some records. Logs and outputs are then used to prove or clear every suspicion, and to show what actually ran. Lessons come from a full audit of a DART-seq project, where code reading found an unexported variable, a no-op option and an annotation converter that mislabelled many transcripts and coding regions (see `docs/bullseye_known_issues.md`).

## Ground rules
- **Read-only.** Do not edit, re-run, move or delete anything of the audited step. Report each problem with evidence, impact and a proposed fix, then ask; the user decides. Check code and its output go to `<step>/audit/` (or the scratchpad), never over results.
- **Data scope:** only the project folder the user named and the files its scripts point to; other folders only with explicit permission.
- Light checks on the login node; heavy ones (whole BAMs, GTF -> TxDb, many samples with samtools) via `qsub` with the job script saved in `audit/`.
- Report in the user's language, short and tabular.

## Workflow
1. **Map the code.** From `STATUS.md`/README find the run script of the step and follow every call to the line that runs the tool (`run_*.sh` -> `qsub ... *.job.sh` -> tool; R/Python/Perl scripts it calls). Note which copy ran: project copy (`MANIFEST.md5`, `md5sum -c`), live skill, or a version changed after the run (compare mtimes with the run time; if changed, the tool's own logs are the record of what ran).
2. **Read every script line by line** with `references/code_review.md` open. For each line ask: what does it assume, what if the input is missing/empty/duplicated/on the other strand, is the variable defined here. Write each suspicion down; do not decide yet.
3. **Parameters - right option, right value, really passed.** Table: parameter -> where set -> how passed -> received as -> value the tool got.
   - Caller arguments vs callee `$1..$n` / `sys.argv` / `commandArgs`: count, order, meaning (`scripts/check_qsub_vars.py <caller> <callee>` for every hand-over, including chain scripts that call run scripts).
   - Environment: `VAR=x bash`, `export`, `qsub -v`; variables only taking a `${X:-default}` because nobody set them.
   - Tool side: the option exists (`GetOptions`, `--help`, `formals()`), its meaning, units (percent vs fraction), default, and whether the code actually uses it (parsed-but-unused options exist).
4. **Samples - right place in the pipeline.** Trace samples from raw files to final columns (`references/code_review.md` section 1): R1/R2 at every step; the sample list is the intended set; pairing tables (which column is control, which is treated/edited, and how the worker uses them); group/condition parsed from names; matrix column order vs metadata; reference level and sign of fold changes; skip-if-exists logic mixing old and new outputs.
5. **References and annotation.** Same release and naming for FASTA, GTF, index, transcriptome, annotation, SNP set; `cmp` files that should be identical; derived annotation compared with its GTF **in every field** (names, strand, tx and CDS start/end, exons) plus a positive control (inject 2-3 errors, the check must catch them). Strandedness consistent across tools (RSEM `reverse` = dUTP = Bullseye `--stranded` labels = RNA strand), proven with data (`scripts/strand_vs_gene.py`; read 1 antisense on known genes).
6. **Prove or clear each suspicion with evidence, for all samples, with code** (`references/checklists.md` lists where each tool records what it did): aligner/tool logs, BAM `@PG`, fastp JSON, RSEM `.model`, Bullseye `Processing files`, R `sessionInfo`; outputs complete (n/n), thresholds satisfied (allow printed rounding), counts reconcile with the step before; `scripts/scan_logs.sh` to separate errors from progress lines. For a code bug: minimal reproduction on real input with the ORIGINAL code, a genome-/sample-wide count, and a test that the mechanism explains all affected records.
7. **Recompute one or two key numbers independently** (size factors, log2FC direction, coverage with/without duplicate-flagged reads, model forward probability, strand labels). State the difference.
8. **Method claims**: check the source of any "the paper does X" (docs, methods; Europe PMC / NCBI efetch if PMC is blocked); label the project's own choices.
9. **Audit the audit.** Before reporting a failure, re-check the check (`references/pitfalls.md`): wrong log line parsed, wrong argument order in your own function, `module load` in a pipe, a tool's default filters, rounding. Do not explain a failure away without evidence.
10. **Report**, save it, add one line to `STATUS.md`, and tell the user: verdict, bugs found (with proposed fixes, not applied), caveats/choices, what was not checked, your own corrected mistakes.

## Report
Save to `<step folder>/audit/audit_<step>_<YYYYMMDD>.md` with the check scripts and logs next to it, and add a one-line pointer to the project `STATUS.md`:

```markdown
# Audit: <step> (<project>), <date>
## Verdict
<one or two sentences: correct / correct with caveats / bug found>
## Code reviewed
<files, which version ran (md5 / mtime), call graph>
## Bugs and suspicions
| Where (file:line) | What | Evidence | Impact | Proposed fix (not applied) | Status: confirmed / cleared |
## Parameter map
| Parameter | Set in | Passed via | Received as | Value the tool got | OK? |
## Sample placement
| Hand-over | Rule | Checked on | Result |
## Other checks (reference, outputs, recomputation)
| Check | How | Result | Evidence |
## Caveats / choices (not errors)
## Not checked
## Corrections to my own checks
## Files
```

## Bundled scripts (`scripts/`, read-only; `python3 -I`)
- `check_qsub_vars.py <caller.sh> <callee.sh>` - positional-argument count per call vs `$N` used; env vars the callee reads that nobody passes (`-v`, including names inside a `-v "$VARS"` string, or `export`); silent `${X:-default}`s. Found the unexported `STRANDED` of the reference audit.
- `scan_logs.sh <dir> [glob]` - log summary with numbers masked, so repeated progress lines collapse.
- `check_annotation_vs_gtf.py <gtf> <refFlat>` - gene name per transcript vs GTF `gene_name`/`gene_id`; carry-over pattern.
- `compare_beds.py <old_dir> <new_dir> [names]` - Bullseye beds before/after a change: same sites? gene-name and region-label changes counted separately; anything else must be 0.
- `strand_vs_gene.py <refFlat> <gtf> <matrix_dir> <out.tsv> <tmp_dir> <samples...>` - stranded parseBAM labels vs known genes (needs `module load htslib`).

## References
- `references/code_review.md` - bug patterns by language (bash/SGE, Perl, Python, R) and the sample-placement checklist. Read before step 2.
- `references/checklists.md` - per step: what to check and where each tool records what it did.
- `references/pitfalls.md` - false alarms from the checking code itself, and known tool traps. Read before reporting a failure.
