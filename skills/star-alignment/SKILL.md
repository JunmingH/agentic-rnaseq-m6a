---
name: star-alignment
description: Align paired-end (trimmed) fastq files to a genome with STAR on an SGE cluster and produce sorted, indexed BAMs, optionally duplicate-marked/removed and with a transcriptome BAM (one SGE job per sample); first checks the STAR index and builds it (separate step) if one with the right sjdbOverhang does not exist. Use when the user asks to align / map fastq to bam, build a STAR index, or prepare BAMs for Bullseye / DART-seq analysis, usually after the read-trimming step.
---

# STAR alignment (+ index build if needed, SGE)

Defaults suit a genome + Ensembl GTF and the Bullseye pipeline. Bullseye only needs coordinate-sorted, indexed BAMs; it reads the sequence + CIGAR, not the MD tag.

## What it does
1. **Index** (a separate step, done first). `READ_LEN` (default 150 = the maximum read length of the data) gives `--sjdbOverhang READ_LEN-1` and the index folder `<REF_DIR>/index_<READ_LEN>` (or set `INDEX_DIR`).
   - Check an existing index: `SA` and `genomeParameters.txt` must exist, with `sjdbOverhang` = READ_LEN-1.
   - If it is missing, build it with `scripts/build_star_index.sh <fasta> <gtf> <index_dir> <read_length>` (16 slots). It builds into `<index>.building`, renames at the end, and refuses to touch an existing dir.
   - `run_alignment.sh` never builds an index; it only waits for a running job named `star_index`.
2. **Per sample** (`align.job.sh`, 16 slots, `h_rt=24:00:00`):
   - `STAR --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within`, plus `--quantMode TranscriptomeSAM` only if the sample gets a transcriptome BAM;
   - `samtools index`;
   - Picard 2.25.2 `MarkDuplicates` + index, depending on `RM_DUP`: `mark` only sets the 0x400 flag (`*.markdup.bam`, all reads kept), `remove` deletes the duplicates (`*.rm_dup.bam`), `0` skips Picard.

   The plain sorted BAM is always kept.
3. **Summary** (`align_summary` job, after all samples): `align_summary.tsv` from STAR's `*_Log.final.out` (input reads, unique %, multi %, too-short %, mismatch rate, `status` = BAM index present).

Outputs in `<out_dir>`:
- `<s>_Aligned.sortedByCoord.out.bam(.bai)`;
- depending on `RM_DUP`, `.markdup.bam(.bai)` or `.rm_dup.bam(.bai)`, plus a `.csv` of duplicate metrics;
- `<s>_Log.final.out`, and logs in `error_file/`.

## STAR version
`STAR_BIN` = a STAR binary (we used a static 2.7.11a build, because the dynamic build needed a newer glibc than the nodes had). If it is empty or not executable, the scripts use `module load star/2.7.10b`. Index and alignment must use the same STAR; the scripts use one binary for both.

## Inputs to confirm with the user
- **in_dir** (trimmed fastq), **out_dir** (e.g. `<project>/bam`), **sample list** (`ls *_paired_R1.fastq.gz | sed 's/_paired_R1.fastq.gz//' > samples.txt`), suffixes (default `_paired_R1/R2.fastq.gz`).
- **Reference**: `REF_DIR` or `INDEX_DIR` (e.g. `/path/to/ref/GRCh38_ens111/index_150`), plus the genome fasta + GTF for an index build. Use the same Ensembl release for every later step (RSEM, refFlat, dbSNP, Guitar).
- **READ_LEN**: the maximum read length after trimming (e.g. 150), not the shortest. sjdbOverhang = max length - 1 is the recommended value; a value that is too small hurts reads that span junctions.
- **Two options the user must decide.** The script refuses to run unless `RM_DUP` (0|mark|remove) and `TRANSCRIPTOME` (0|1|file) are set. Always ask, explain in plain words, and give the recommendation:
  - **`RM_DUP`, duplicates (Picard MarkDuplicates): `0` / `mark` / `remove`?**
    - What it is: read pairs whose two mates start at exactly the same positions count as PCR copies.
    - The trade-off: without UMIs, independent reads on very highly expressed genes are also called duplicates, so deleting them lowers coverage there. Keeping them lets PCR copies inflate the C→T counts at low input.
    - `mark` only sets the 0x400 flag. Bullseye `parseBAM.pl --removeDuplicates` ignores flagged reads, so the same BAM can be analysed with and without dedup.
    - Recommend `mark`. Dedup can be added later: rerunning with the same out_dir skips STAR and only runs Picard.
  - **`TRANSCRIPTOME`, also write the transcript-coordinate BAM (`--quantMode TranscriptomeSAM`)?**
    - It is only an input for RSEM; Bullseye needs only the genome BAM.
    - It can be a file of sample names, so that only those samples get it (for example only the DART-seq control samples, which give the RNA level).
    - It cannot be added later without re-running STAR for those samples.

**Software:** the job scripts activate the repo's conda environments when `ENV_SETUP` is set (see `docs/INSTALL.md`; pass it with the other env vars), otherwise they use the cluster's `module load` lines. Ask the user which applies before the first run.

## Run rules
1. **Save the code that actually runs, in the right folder.**
   - Alignment: copy `scripts/` (without `build_star_index.sh`) into `<project>/bam/scripts/` and write the `md5sum` of the scripts + the STAR binary into `<project>/bam/MANIFEST.md5`. Then write `<project>/bam/run_align_<YYYYMMDD>.sh`, which calls the COPY, never the live skill.
   - Index: put `build_star_index.sh`, the exact qsub command (`build_index_<READ_LEN>_<YYYYMMDD>.sh`) and an md5 file next to the index.
2. **Keep all logs** (`<out_dir>/error_file/`).
3. A failed sample with a clear cause: delete only its outputs + logs and rerun (finished samples are skipped). Unclear or many failures: keep the logs, report and wait.
4. **Never submit before the user approves.** The index build and the alignment are approved separately. After the build, check it (sjdbOverhang, STAR version, files, logs) and report before submitting the alignment.
5. Run ONE sample first (the whole worker) and check output, time and memory. Then run the rest.
6. Use `command rm -f` in commands: `rm` may be aliased to interactive mode and hangs on the prompt.

## Steps
1. Confirm the inputs and the two options. Check the index and the disk space: a sorted BAM is roughly 0.5–1× the size of the trimmed fastq pair, and the marked and transcriptome BAMs add more.
2. If the index is missing: write the build script (rule 1), get approval and submit (about 1–2 h for human). Check the result and report.
3. Build the sample list and check that both files exist for every sample. Copy the code, write `run_align_<date>.sh`, get approval, test one sample, then submit the rest:
   ```bash
   REF_DIR=/path/to/ref/GRCh38_ens111 READ_LEN=150 RM_DUP=mark TRANSCRIPTOME=control_samples.txt SGE_PROJECT=<project> \
     bash <project>/bam/scripts/run_alignment.sh <project>/trimm_fastp <project>/bam <project>/samples.txt
   ```
4. Check `align_summary.tsv`:
   - every sample `ok`;
   - unique % (human RNA-seq is usually 80–95 %; much lower means investigate, e.g. contamination);
   - too-short %;
   - no sample far from the rest.

   Grep `error_file/*` for `terminated` / `exception`.
5. Update the project STATUS.md.

## Gotchas
- The cluster may kill jobs that use more CPU than their slots. STAR's sorting adds threads, so `--runThreadN` is `NSLOTS-2`.
- genomeGenerate for human needs about 30–40 GB RAM (`--limitGenomeGenerateRAM 60e9`). If the job dies out of memory, ask for more slots.
- qsub runs a spooled copy of job scripts: paths are passed as arguments, never found through `$0`.
- Index folder names can lie. Read `genomeParameters.txt`; an index named for 51 bp reads can turn out to have the default `sjdbOverhang 100`.
- Some `samtools` modules refuse to load unless a matching `htslib` module is loaded first. One full run lost its index/Picard step to this, although STAR had finished. Test every `module load` in `bash -lc` before submitting, and do a one-sample test whenever a script is new.
- The worker skips STAR when `<s>_Aligned.sortedByCoord.out.bam` and `<s>_Log.final.out` exist (STAR writes the log last), so a rerun only does the missing index/Picard steps.
