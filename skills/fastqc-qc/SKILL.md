---
name: fastqc-qc
description: Run FastQC (one SGE job per fastq file) plus a MultiQC summary on a directory of fastq/fq.gz files on an SGE cluster. Use when the user asks to do quality check / QC / fastqc on raw sequencing reads, e.g. after downloading or merging sequencing-provider data.
---

# FastQC quality check (SGE)

## Inputs to confirm with the user
- **fastq dir**: the directory holding the `*.fq.gz` files (symlinks are fine). Output goes to `<fastq dir>/fastqc`.
- **file list**: a text file with one fastq file name per line. By default, build it with `ls *.fq.gz` and show the user the count before submitting.
- **SGE project** (`SGE_PROJECT`, the `qsub -P` value), if the cluster requires one.

**Software:** the job scripts activate the repo's conda environments when `ENV_SETUP` is set (see `docs/INSTALL.md`; pass it with the other env vars), otherwise they use the cluster's `module load` lines. Ask the user which applies before the first run.

## Run rules
1. **Save the code first.** Write the exact commands into a script inside the target folder, e.g. `<fastq dir>/run_fastqc_<YYYYMMDD>.sh` (the `ls > list` step plus the `run_fastqc.sh` call). Show it to the user, then run that file, not ad-hoc commands.
2. **Keep all logs**, including qsub logs (`fastqc/error_file/`). Never delete the logs of a successful step.
3. **A failure with a clear cause** (killed job, node problem): delete only the logs and partial outputs of that file, then rerun it.
4. **A failure that needs the user's eyes** (unclear error, many files failing): leave the error logs untouched, report which files failed, and wait for the user.

## Steps
1. Build the list (file names relative to the fastq dir):
   ```bash
   cd <fastq dir> && ls *.fq.gz > list && wc -l list
   ```
2. Show the user the count and the scripts, and get approval before submitting (one job per file).
3. Submit, from any directory:
   ```bash
   SGE_PROJECT=<project> bash /path/to/skills/fastqc-qc/scripts/run_fastqc.sh <fastq dir> <fastq dir>/list
   ```
   - Each file gets one job (`-pe omp 4 -l h_rt=24:00:00`, `module load fastqc`) that runs `fastqc <file> -o <fastq dir>/fastqc -t $NSLOTS`.
   - One final job, `multiqc_all`, waits with `-hold_jid <all fastqc job ids>` and starts after every fastqc job has finished. It runs `multiqc . -f` inside `<fastq dir>/fastqc` and writes `multiqc_report.html` + `multiqc_data/`.
   - Logs: `<fastq dir>/fastqc/error_file/<file>_{output,error}_fastqc.txt` and `multiqc_{output,error}.txt`.
4. Watch `qstat -u $USER`. `multiqc_all` sits in state `hqw` while it waits; that is normal. When everything is done, check that every input has a `<name>_fastqc.zip`:
   ```bash
   cd <fastq dir>; for f in $(cat list); do b=${f%.fq.gz}; b=${b%.fastq.gz}; [ -s fastqc/${b}_fastqc.zip ] || echo MISSING $f; done
   ```
   - Look for real errors in `*_error_fastqc.txt`. FastQC prints progress to stdout, not stderr.
   - If a fastqc job failed, `multiqc_all` stays held forever. `qdel` it, rerun the missing files, then submit `scripts/multiqc.job.sh <fastq dir>` by hand.
5. Confirm that `<fastq dir>/fastqc/multiqc_report.html` exists and is newer than the fastqc outputs.
6. Report to the user:
   - any MISSING files;
   - the headline numbers from `multiqc_data/multiqc_fastqc.txt`: total sequences, %GC, and the per-base quality / adapter / overrepresented-sequence flags.

   Do not call a sample "pass" without looking at those flags.

## Gotchas
- `fastqc.job.sh` and `multiqc.job.sh` are the qsub workers. Keep them in the same `scripts/` dir as `run_fastqc.sh`; the submit script uses absolute paths, so it can be run from any directory.
- FastQC reports each file separately, so R1 and R2 get separate reports. Judge a sample on both.
- Overrepresented sequences that are not adapters, or a low mapping rate later, can mean contamination (for example a host bacterium in a library). BLAST a few of them before trimming hides the evidence.
- Update the project's STATUS.md (if it has one) with the job IDs, date and result.
