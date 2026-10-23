# Checklists per step

Contents: 0 SGE / qsub (every step) · 1 FastQC/MultiQC · 2 Trimming (fastp/Trimmomatic) · 3 STAR + Picard · 4 RSEM · 5 DESeq2 / PCA · 6 Bullseye (parseBAM, Find_edit_site, summarize) · 7 Differential m6A (site/gene GLM) · 8 Guitar metagene

Every check: do it for ALL samples with a script, print `n checked / n wrong`, keep the script in `audit/`.

## 0. SGE / qsub (applies to every step)
- Run script saved in the step folder? Calls a project copy (with `MANIFEST.md5`) or the live skill? `md5sum -c` from the directory the manifest paths are relative to; `cmp` copy vs skill.
- Variables: `VAR=x bash run.sh` sets VAR only for that command; jobs get only what `qsub -v` passes (or `-V`). A variable used in the worker but not in `-v` (or not `export`ed in a chain script) silently takes the worker's default. Grep the worker for every `$VAR` and check it is passed.
- `$0` inside a qsub job is the spooled copy (`/var/spool/sge/...`), so `dirname $0` fails; scripts must take their own directory as an argument.
- Threads vs slots: tools that use more cores than `-pe omp N` get killed (fastp `--thread N` uses ~N+1; STAR sorting threads); workers should use `NSLOTS-1/-2`. Check `qacct -j <id>` (`exit_status`, `maxvmem`, `ru_wallclock`; `timeout 60 qacct` if it hangs) and the e-mail for killed jobs (the error file may be empty).
- Job names are truncated (`${s:0:20}`) and different stages can share a name (`fe_ctrl_1_WT` stranded vs unstranded): identify jobs by `qstat -j <id>` -> `job_args`, `jid_predecessor_list`.
- `| tee log` hides the exit status of the command before it: chains that should stop need `set -o pipefail`.
- Re-runs that skip existing outputs: check that skipped outputs are from the right run (mtime, logs).

## 1. FastQC / MultiQC
- One `_fastqc.zip` per fastq (2 x samples); MultiQC job held on all FastQC jobs (`-hold_jid`); `multiqc_report.html` newer than the last zip.
- Total reads per file R1 = R2; compare with the trimming input counts later.

## 2. Trimming
- fastp: the JSON of every sample has the full `command`: parse `-i/-I` (R1 = `_1`/`R1`, R2 = `_2`), `-o/-O`, and every option; group samples by option set (only `--thread` may differ). Versions (`summary.fastp_version`). `read1_before == read2_before` and `read1_after == read2_after` per sample.
- Input symlinks point to the merged files of the SAME sample name (`readlink`).
- Adapter auto-detection may differ between samples (R2 often `unspecified`); for PE data fastp also trims by overlap: check `adapter_trimmed_reads` and the post-trim FastQC adapter module instead of the detected sequence.
- Error logs of fastp contain lines like "reads failed due to too short" - statistics, not errors.
- Trimmomatic: `ILLUMINACLIP` file and `keepBothReads`; summary file per sample.

## 3. STAR (+ index, Picard, samtools)
- Index: `genomeParameters.txt` -> `genomeFastaFiles`, `sjdbGTFfile`, `sjdbOverhang` (= max read length - 1), `versionGenome` (format version, e.g. 2.7.4a for STAR 2.7.11a is normal); `chrNameLength.txt` count = BAM `@SQ` count; naming `1` vs `chr1` same as the GTF/annotation used later.
- Per sample `<s>_Log.out` (`readFilesIn` line: fields 2 and 3 are R1, R2; the line ends with `~RE-DEFINED`), `genomeDir`, `quantMode`; BAM `@PG CL:` repeats the command. Check R1 first, R2 second, the sample's own files, for all samples.
- Read names in R1/R2: `1:N:0` vs `2:N:0`, same coordinates (pairing not shifted).
- R1/R2 from the data, independent of logs: share of read 1 (`-f 0x40 -F 0x904`) on the strand opposite to known genes (GAPDH, MALAT1 +; ACTB, EEF1A1 -). dUTP library: ~99 %; swapped R1/R2 would give ~1 %.
- Script versions: if `align.job.sh` (or the manifest) is newer than the BAMs, the record of what ran is `Log.out` / `@PG`, not the current script.
- Counts: `Number of input reads` in `Log.final.out` = trimmed pairs; `samtools flagstat` primary = 2 x pairs, "paired in sequencing" = primary.
- `TranscriptomeSAM` only where intended (e.g. control samples only), and transcriptome BAMs exist exactly for those.
- Picard MarkDuplicates: `@PG` shows `REMOVE_DUPLICATES`; flagstat `duplicates` > 0 for markdup BAMs.
- Default STAR filters (`--outFilterScoreMinOverLread 0.66` etc.) explain "too short" unmapped reads; a soft-masked (`dna_sm`) FASTA is fine.
- samtools 1.19.2 in some module systems needs `module load htslib/1.19.1` first.

## 4. RSEM
- Inputs: transcriptome BAMs of the intended samples; `--paired-end --bam`; `--strandedness` must match the library. The command line printed by RSEM does NOT show the strandedness: read line 3 of `<s>.stat/<s>.model` (forward probability: 1 = forward, 0 = reverse, 0.5 = none) for all samples.
- Strand test (one sample, none/forward/reverse): the right setting keeps 95-99 % of the `none` assigned reads; the wrong one a few %. MT-ND6 TPM far above other MT genes = wrong strand.
- Reference: built from the same GTF/FASTA as the STAR index; transcript IDs in the BAM `@SQ` == RSEM transcript IDs (`isoforms.results` first column, strip `_name`); compare sorted lists, expect 0 differences.
- `.cnt` line 1: N0 N1 N2 N_tot; N0 = STAR's unmapped reads (they are in the transcriptome BAM with `--outSAMunmapped Within`), so N0 should equal `Number of reads unmapped` summed from `Log.final.out`; sum of `expected_count` / N_tot (a large share with the right strand setting, a few % with the wrong one); merged matrix: row counts, column order = sample list, column sums = per-sample sums.
- Old/invalid run scripts (e.g. a forward-stranded run) still in the folder with misleading header comments; skip-if-exists in `run_rsem.sh` would keep old results after a settings change.

## 5. DESeq2 / PCA
- Input = un-normalised counts (RSEM `expected_count` rounded, or tximport); never TPM/normalised counts. Note `round()` vs `as.integer()` (truncation).
- Design and reference level (`relevel`); coefficient names asserted in code (`stopifnot(coef %in% resultsNames(dds))`); one model per subset or one global model - which, and why.
- Size factors: `DESeq()` calls `estimateSizeFactors` (median-of-ratios) unless size factors exist; recompute independently (geometric mean over genes with no zero, median ratio per sample) - should match to ~1e-12. Normalised counts are only for display; vst only for PCA/heatmaps.
- Direction: for baseMean >= 100, corr(log2FC, log2(mean_norm_num / mean_norm_den)) ~ 1; every significant gene's sign agrees.
- `lfcShrink` (apeglm) changes LFC only, not p/padj; which LFC the volcano/threshold uses.
- PCA: vst(blind) -> top 500 var genes -> `prcomp(t(X), center = TRUE, scale. = FALSE)` (DESeq2 convention). Old code with `prcomp(counts, scale.=TRUE)` on genes x samples does PCA of genes, not samples.
- Duplicate gene symbols (`id_symbol` split) - joins by symbol lose rows; prefer gene IDs.

## 6. Bullseye
- parseBAM: per sample the log's last line `done <bam>` = the intended BAM (markdup BAM if `--removeDuplicates`); options exist in `GetOptions`; matrix columns 8 (unstranded) or 9 (stranded). Dedup really applied: matrix coverage at a few positions ~ `samtools mpileup -x -B -Q 0 -q 0 --ff 0x404` (dups excluded, secondary kept) and well below `--ff 0x4`. Multi-mapped (secondary) reads are counted unless `--removeMultiMapped`.
- Stranded matrix: label = RNA strand? `scripts/strand_vs_gene.py` (>= 99 % of exonic coverage on the gene's own strand) and C->T high in `+` rows / G->A high in `-` rows. Must agree with RSEM's strandedness.
- Find_edit_site: per bed, the error log's `Processing files <refFlat>, <control matrix> and <edited matrix>` = pair table (control = first column); `Detecting bases with 0.1 to 0.95 edit` etc. per bed; `minEdit/maxEdit` are percent (`/= 100` in code). `--stranded` in Find_edit_site.pl is parsed but never used: strand comes from the matrix's 9th column. Bed rows satisfy thresholds (fold computed from 4-decimal ratios can show 1.4986 for an exact 1.5).
- Annotation: refFlat vs GTF in EVERY field - names (`scripts/check_annotation_vs_gtf.py`; genes without `gene_name` are mislabelled by `gtf2genepred.pl` v0.2), chrom/strand/tx start-end/exons, and cdsStart/cdsEnd vs the span of the GTF CDS + start_codon + stop_codon features (v0.2 gets 8,618 coding transcripts wrong). Rebuild the refFlat from the GTF in Python and compare field by field; inject 2-3 errors into a copy as a positive control. Gene names in beds are labels only (annotation is stored per transcript ID), so a fixed refFlat changes names, not sites - confirm with `scripts/compare_beds.py`. Overlapping/bicistronic genes (PIGY/PYURF, CHTF8/DERPC, FAM220A/SAGSIN1) get an arbitrary one of the names in each run.
- summarize: input = the 3 replicate beds of the group from the intended bed set; `--mut/--coverage` filter each per-sample record before counting replicates.

## 7. Differential m6A (site/gene GLM, Meyer-lab style)
- Site sets reproduce from upstream (reps >= 2 of 3, mean >= 5 %, RAC flag, dbSNP removal) - compare counts with the filter table.
- RAC check uses the same genome as the BAMs (`cmp` the FASTA), 1-based position = edited C, reverse complement on `-`.
- Counts from the matrix: `+` site mut = T, `-` site mut = A, cov = total; one strand row per site in stranded matrices.
- Test: quasibinomial GLM, Wald t vs z, IHW covariate (own choice?), coverage filter, min samples per group; gene level = sum over kept sites per gene (depends on correct gene names); expression join by symbol vs ID.
- Literature parameters (Meyer 2019, Tegowski 2024, Flamand 2022 = TRIBE-STAMP) - note deviations.

## 8. Guitar (metagene distribution)
- BED6 with strand; site sets reproduce the intended filtered sets (counts vs upstream table).
- TxDb from the same GTF; `pltTxType`, `headOrtail`, `enableCI`, `mapFilterTranscript` as intended; log reports how many mRNAs survive the length filter (Ensembl 111: 10,835 of 252,989 transcripts).
- `GuitarPlot` returns the ggplot only without `miscOutFilePrefix` (with it: writes a pdf, returns NULL); `txGuitarTxdb` must be a file read by `read.table`, not the `makeGuitarTxdb` object; the plot's layers use internal vectors, so `facet_wrap` and `%+%` fail - draw subsets from `p$data` (x, density, group).
