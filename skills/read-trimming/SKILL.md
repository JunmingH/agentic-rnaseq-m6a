---
name: read-trimming
description: Adapter trimming and quality filtering of paired-end fastq files on an SGE cluster with a choice of Trimmomatic or fastp (one SGE job per sample), plus a trimming summary table. Use when the user asks to trim, remove adapters, remove poly-G, or quality-filter raw reads, usually after the fastqc-qc step showed adapter content.
---

# Read trimming + quality filter (Trimmomatic or fastp, SGE)

Choose the tool with the env var `TOOL` (default `trimmomatic`).

## What each option does (per sample)
| | Trimmomatic 0.39 (`TOOL=trimmomatic`) | fastp 0.23.4 (`TOOL=fastp`) |
|---|---|---|
| adapters | `ILLUMINACLIP:<ADAPTER_DIR>/TruSeq3-PE-2.fa:2:30:10:8:true` (keeps both reads) | auto-detected (`--detect_adapter_for_pe`), or `ADAPTER_R1`/`ADAPTER_R2` |
| read ends | `LEADING:3 TRAILING:3` | `--cut_front/--cut_tail`, window 1, mean Q3 |
| sliding window | `SLIDINGWINDOW:4:20` | `--cut_right`, window 4, mean Q20 |
| min length | `MINLEN:50` | `--length_required 50` |
| poly-G | **not removed** | `--trim_poly_g` |
| extra | none | fastp defaults stay on (drop reads with > 40 % of bases < Q15 or > 5 N); HTML + JSON report per sample |

Outputs per sample in `<out_dir>` (same names for both tools):
- `<s>_paired_R1.fastq.gz` and `<s>_paired_R2.fastq.gz`: use these downstream;
- `<s>_unpaired_R1/R2.fastq.gz`;
- `<s>_trimm_summary.txt` (Trimmomatic) or `<s>_fastp.json/.html` (fastp).

After all jobs, `trimm_summary.tsv` lists per sample: tool, input pairs, pairs where both reads survived, dropped, %, and Q30 before/after for fastp.

**How to choose.** In a side-by-side test on a subsample of a two-colour-chemistry library, both tools removed adapters well and kept about the same number of pairs. Only fastp cleared the poly-G tails that two-colour chemistry leaves. Choose Trimmomatic to stay consistent with older data processed the same way; choose fastp when poly-G matters or for speed. Ask the user if they have not said.

## Inputs to confirm with the user
- **tool**: trimmomatic or fastp.
- **in_dir**: the raw fastq dir. **out_dir**: a new dir (e.g. `<project>/trimm_fastp`); it is never deleted or emptied.
- **sample list**: sample names (not file names), one per line: `ls *_1.fq.gz | sed 's/_1.fq.gz//' > samples.txt`.
- **file suffixes**: default `_1.fq.gz` / `_2.fq.gz`; pass others as arguments 4 and 5 (e.g. `_R1.fastq.gz` / `_R2.fastq.gz`).
- **MINLEN** (default 50). For Trimmomatic also the adapter file (`ADAPTER`) and `ADAPTER_DIR` (the folder of Trimmomatic's adapter fasta files). Ask before changing these.

**Software:** the job scripts activate the repo's conda environments when `ENV_SETUP` is set (see `docs/INSTALL.md`; pass it with the other env vars), otherwise they use the cluster's `module load` lines. Ask the user which applies before the first run.

## Run rules
1. **Save the code first.** Write the exact commands (sample-list step + `run_trimm.sh` call + follow-up QC) into a script inside the target folder, e.g. `<project>/trimm/run_trimm_<YYYYMMDD>.sh`. Show it to the user, then run that file.
2. **Keep all logs** (`<out_dir>/error_file/`). Never delete the logs of a successful step.
3. **A sample fails with a clear cause:** delete only that sample's outputs and logs, then rerun. Rerunning `run_trimm.sh` skips finished samples.
4. **Unclear or many failures:** keep the logs, report, and wait.
5. Never submit before the user approves the script.

## Steps
1. Build the sample list, check that every sample has both files, and show the count.
2. Write the run script (rule 1) and get approval.
3. Submit:
   ```bash
   TOOL=fastp SGE_PROJECT=<project> bash <project>/trimm_fastp/scripts/run_trimm.sh <in_dir> <out_dir> <sample_list> [r1_suffix r2_suffix]
   ```
   Each sample gets one job (`-pe omp 4 -l h_rt=10:00:00`). A final `trimm_summary` job waits until all samples finish. Samples whose outputs exist are skipped.
4. When the jobs are done, read `<out_dir>/trimm_summary.tsv`. Every sample should be `ok`; compare `both_pct` across samples and flag clearly lower ones.
5. Check the `*_paired_R1/R2.fastq.gz` files with the `fastqc-qc` skill: adapter (and, for fastp, poly-G) content should be gone.
6. Report the surviving percentages (min / median / max), odd samples and any residual adapter or poly-G. Update the project's STATUS.md.

## Gotchas
- **The cluster may kill jobs that use more CPU than their slots.** Example message: "using 5.1 processors but was allocated only 4". fastp `--thread 4` really used about 5.1 cores, so short jobs finished but long ones were killed with exit 137. The workers therefore use `NSLOTS-2` threads for fastp and `NSLOTS-1` for Trimmomatic. After such a kill, read the error file or the e-mail, not only `qacct`.
- Trimmomatic's `TruSeq3-PE.fa` holds only the PrefixPE sequences (palindrome mode) and left clearly more adapter reads than `TruSeq3-PE-2.fa`. That is why `-PE-2` is the default.
- qsub runs a spooled copy of each job script: never find sibling files through `$0`; pass paths as arguments.
- Always use the `_paired_` files downstream; `_unpaired_` files hold reads whose mate was dropped. fastp's summary counts whole pairs, so `forward_only` / `reverse_only` are `-` for fastp.
- Trimmed reads no longer have a fixed length, so downstream tools must not assume one. For STAR, use the maximum length.
