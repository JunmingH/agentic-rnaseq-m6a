# Code review: bug patterns and sample placement

The core of an audit is reading the code line by line. Use this list while reading; every item has caused (or nearly caused) a wrong result in a real project. When a pattern matches, write it down as a suspicion and then prove or clear it (small reproduction, logs, outputs).

Contents: 1 Sample placement · 2 Argument and parameter hand-over · 3 bash / SGE · 4 Perl · 5 Python · 6 R / DESeq2 · 7 Annotation and coordinates · 8 How to prove a suspected bug

## 1. Sample placement (where each sample goes in the pipeline)
Trace every sample from its raw file to its final column/row. For each hand-over write down the rule and check it on ALL samples.
- **Paired files**: which file is R1 and which R2 at each step (raw `_1/_2` -> trimmed `_paired_R1/_R2` -> aligner `readFilesIn R1 R2`). Suffix variables (`R1S`, `R2S`) and their defaults; positional order in the tool call. Evidence: read headers `1:N`/`2:N`, aligner log, strand of read 1 on known genes.
- **Sample lists**: where the list comes from (`ls`, a file, a regex on names); is it the intended subset (e.g. control samples only for RSEM, or only one sub-study); duplicates or blanks; order. Compare against the files that exist.
- **Pairing / role tables** (control vs treated, control vs edited): column order in the table AND the order in which the worker uses them (`c e o` -> `--controlMatrix $c --EditedMatrix $e`). A swapped column order silently swaps biology. Check the same-number rule (`ctrl_1_Mut` with `ctrl_1_WT`), group consistency, output name = the right sample, no name collisions between runs that share an output dir.
- **Group / condition parsed from names**: the regex (`sub("_[0-9]+_Mut$", "")`, `rsplit("_", 2)`) on every sample name; names with extra underscores (`ctrl_typeA_2_WT`); prefixes that collide (`T1` vs `T10` with `startswith`).
- **Matrix columns vs metadata rows**: count matrix columns in the same order as `colData` rows (`stopifnot(identical(colnames(cts), rownames(coldata)))`); merged matrices keep the sample-list order; joins by name, not by position.
- **Reference level / direction**: `relevel(ref = "vehicle")`, `contrast = c(factor, numerator, denominator)`, sign of log2FC = treated / control; in editing tools "edited" vs "control" roles.
- **Replicates**: numbering means replicate (not pair) unless the user confirmed pairing; designs that should include a pair/batch term.
- **Skip-if-exists logic**: a re-run with new settings skips samples whose output exists from the old settings -> mixed results. Check mtimes and logs of every output against the intended run.

## 2. Argument and parameter hand-over
Build a table: parameter -> where set -> how passed -> name in the receiver -> value the tool got.
- Caller `cmd a b c` vs receiver `$1 $2 $3` (bash), `sys.argv[1:]`, `commandArgs(trailingOnly=TRUE)`: same count, same order, same meaning. Off-by-one after a `shift`.
- Options: does each option exist in the tool (`GetOptions`, `argparse`, `--help`, `formals()`)? Spelling/case (`--EditedMatrix` vs `--editedMatrix`); abbreviations Getopt accepts silently; options that are parsed but never used (Bullseye `Find_edit_site.pl --stranded`).
- Units and scales: percent vs fraction (`--minEdit 10` = 10 %, divided by 100 inside), 0-based vs 1-based, phred offsets, kb vs bp.
- Defaults: what happens when the parameter is NOT given (`${X:-default}`, `getopt` default, R default args). A default that hides a missing value is a bug source; prefer scripts that refuse to run without required values.
- Thresholds applied twice or in a different order than documented (filters before vs after merging replicates).

## 3. bash / SGE
- `VAR=x bash script` sets VAR for that command only; variables reach a qsub job only via `-v` (or `-V`). A chain script that sets a value but does not `export` it before calling the next script loses it (e.g. `STRANDED`). Use `scripts/check_qsub_vars.py`.
- `$0` / `dirname $0` inside a qsub job points to the spool copy.
- `cmd | tee log || exit 1` never exits on cmd failure without `set -o pipefail`; `module load X | head` loses the module (subshell).
- Missing `|| exit` after `cd`, `mkdir`, tool calls; `rm -rf "$DIR/"` with an empty variable.
- Unquoted variables with spaces/globs; `for f in $(ls ...)` ordering; `read -r` loops that drop the last line without newline.
- `-hold_jid` with an empty id list; job names truncated or shared between stages.
- Threads > slots (killed jobs), memory flags, `h_rt` too short (partial outputs that look complete).

## 4. Perl
- **Stale captures**: `$x =~ /re(...)/; my $v = $1;` - if the match fails `$1` keeps the previous value (gtf2genepred.pl gene names). Use `my ($v) = $x =~ /re(...)/;` and test `defined`.
- Variables named opposite to their content (`$first` holding the 0x80 bit in parseBAM `get_strand`): check bit masks against the SAM spec.
- Hash keys built from labels that are not unique (gene names) overwrite each other; keys from IDs are safe.
- Logic that only works when every record has an optional feature (CDS derived from start/stop codons only).
- `split` limits, `chomp` missing, numeric vs string comparison (`==` vs `eq`).

## 5. Python
- Function argument order vs call (`covering(c, s, p)` called as `covering(c, p, s)`): assert on one known case before running on all.
- Dict built from rows with duplicate keys keeps the last row silently; check `len(set(keys)) == len(rows)`.
- Coordinates: BED start is 0-based, end 1-based; matrix positions 1-based; mixing them shifts by one.
- Strand: keys must include strand where two strands can share a position.
- Python version on compute nodes (3.6 here: no `capture_output`, no walrus).

## 6. R / DESeq2
- `as.integer()` truncates, `round()` rounds; input must be raw counts (or tximport), not TPM / normalised counts.
- `colData` row order vs count columns; factor levels and reference; `design` variables of the right type (numeric vs factor).
- `results(name=)` / `contrast=` matches the intended comparison; `stopifnot(coef %in% resultsNames(dds))`.
- `prcomp(x)` treats ROWS as observations: a genes x samples matrix must be transposed; `scale.` choice; vst before PCA.
- Merges by symbol: duplicated symbols, missing symbols; `!duplicated()` keeps the first row only.
- IHW / random procedures without a seed give slightly different numbers between runs (only last digits here).

## 7. Annotation and coordinates
- Derived annotation (refFlat/genePred/BED from GTF) must agree with the GTF in EVERY field: names, chrom, strand, tx start/end, CDS start/end, exons. Compare field by field (`scripts/check_annotation_vs_gtf.py` for names; a structure check for the rest) and run a positive control (inject 2-3 errors, confirm they are caught).
- Same release and naming everywhere (FASTA, GTF, index, transcriptome, annotation, SNP database); `cmp` files that should be identical.
- Overlapping genes with identical transcripts: tools may pick either name per run.

## 8. How to prove a suspected bug
1. Minimal reproduction with real input: extract the few records that trigger it (e.g. two genes from the GTF) and run the ORIGINAL code on them.
2. Genome-/sample-wide count of affected records, and a test that the mechanism explains all of them (31,617 of 31,617).
3. Impact on the project's results (what changes, by how much, which downstream steps use it).
4. Proposed fix in a copy (never the shared original), its own validation (field-by-field vs source + positive control), and a one-sample test before re-running everything. Ask the user before applying.
