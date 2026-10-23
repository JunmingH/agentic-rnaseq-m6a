---
name: bullseye-analysis
description: Run the Bullseye DART-seq C-to-U (C>T) editing-site analysis on an SGE cluster from STAR BAMs - parseBAM (per-sample base-count matrix), Find_edit_site (edited vs control, one job per pair) and summarize_sites (merge replicates), with strand, replicate and motif checks. Use when the user asks to run Bullseye / find DART-seq or m6A editing sites / compare edited vs control C>T, usually after the star-alignment step.
---

# Bullseye (DART-seq editing sites, SGE)

Bullseye is a set of Perl scripts by M. Flamand (<https://github.com/mflamand/Bullseye>, MIT, `scripts/LICENSE.Bullseye`). The unmodified copies we ran are in `scripts/`: parseBAM 1.4.4, Find_edit_site 1.5.2, summarize_sites 1.21. The environment is a conda env `Bullseye` (`module load miniconda; conda activate Bullseye`) with perl and its modules, bgzip, tabix and samtools.

## Steps (3 stages, one SGE job per sample / per pair)
1. **parseBAM** (`run_parsebam.sh`, `parsebam.job.sh`, 16 slots, `h_rt=24:00:00`). Per sample, base counts at every covered position → `<s>.matrix.gz` + `.tbi`.
   - `--minCoverage 10`.
   - `--removeDuplicates` only if `DEDUP=1`: it ignores reads flagged 0x400 by Picard, so it needs the `markdup`/`rm_dup` BAM.
   - Input BAMs must be coordinate sorted, with `.bai`, and the same chromosome names as the refFlat (Ensembl: `1`, `2`, …).
2. **Find_edit_site** (`make_pairs.py`, then `run_find_edit.sh` / `find_edit.job.sh`, 8 slots). Per pair, the edited matrix is compared with the control matrix → `<name>.bed` (C>T sites; each row has gene|type|…|n_mut|fold, plus control/edited ratio and totals).
   - Annotation (`REFFLAT`, REQUIRED): a refFlat made by Bullseye's `gtf2genepred.pl` **with `patches/gtf2genepred_refflat_fix.patch` applied**, from the same Ensembl GTF as the STAR index.
   - The unpatched converter gives wrong gene names to unnamed genes and wrong CDS bounds to transcripts without start/stop codon features.
3. **summarize** (`summarize.sh`, light): merges beds and keeps sites found in ≥ `MinRep` samples.

## Inputs to confirm with the user (the scripts refuse to run without the REQUIRED ones)
- **Which samples are the control and which are the edited (active DART enzyme) samples.**
  - Labels differ between labs: "WT" can mean the active APOBEC1-YTH, and "Mut" the YTH-mutant control. Do not guess: ask.
  - Confirm with the checks below (the motif check is the decisive one).
  - Pass the tags to `make_pairs.py <samples> <control_tag> <edit_tag> <mode>`.
- **Pairing mode**: `same_n` (`ctrl_1_Mut` vs `ctrl_1_WT`; only if the number means a matched pair, not just a replicate) or `all` (every edited sample vs every control of the group). Show the pair list and the count.
- **`DEDUP` 0/1** (REQUIRED, parseBAM). Without UMIs, real reads on very highly expressed genes are also called duplicates. If unsure, run both, in separate out dirs.
- **Find_edit_site thresholds** (REQUIRED env: `MIN_EDIT`, `MAX_EDIT` in percent, `FOLD`, `MIN_EDIT_SITES`). Optional `MIN_COV` (default 10, control and edited). Meyer-lab bulk values: 10 / 95 / 1.5 / 2 with coverage 10. Thresholds change the site list a lot, so always record which were used.
- **MinRep** for summarize (in how many replicates a site must be found).
- Low-depth samples (check `align_summary.tsv`) give fewer sites; mention it when a sample is far below the rest.

## Strandedness (`--stranded`) and direction/quality checks
- **`parseBAM.pl --stranded` is made for dUTP (ISR) libraries**, i.e. RSEM `reverse` (read 1 antisense to the RNA).
  - Its `get_strand` has swapped flag variable names (`$first` really holds 0x80 = read 2), but the labels come out right for dUTP. Do NOT swap `+`/`-`.
  - `Find_edit_site.pl` uses the strand column (column 9) of the matrix automatically; `--stranded` there is a no-op.
  - For a library that is NOT dUTP the labels are wrong: first determine the library type (rsem-quantification strand test).
- **Test the labels on one sample before the batch** (`scripts/strand_check.sh matrix.gz 1`). In `+` rows C→T must be clearly above background and G→A at background; in `-` rows the reverse. Unstranded rows show both at an intermediate level.
- **Single-pair Find_edit_site beds are mostly noise.** A single pair gives very many sites, most supported by only a few mutated reads and often with no edits in the control; only a minority reproduce across replicates. Never treat them as m6A sites. Filter by replicate agreement and motif (or use the dart-diff-m6a skill).
- **Direction check** (`scripts/validate_rep3.py`):
  1. Choose sites from replicates 1+2: present in both beds, mut ≥ 5, coverage ≥ 20.
  2. Read the editing ratio at the same positions in the untouched replicate 3: T/(C+T) on `+` genes, A/(A+G) on `-` genes; with stranded matrices use the row of the site's strand.
  3. Do it BOTH ways, with the roles swapped in a second run.

  Both directions can reproduce, because the control construct (active APOBEC1 + an m6A-blind YTH) edits non-specifically. Counting sites after swapping the roles does NOT discriminate, and neither does the global C→T rate.
- **Motif check that discriminates** (`scripts/motif_context.py <bed_dir> <validate_dir> <fasta> <out>`; gene-oriented 5-mer centred on the C).
  - The edited C follows the m6A, so true sites are R-A-C-H.
  - A at −1 is high for both roles, because APOBEC1 prefers AC. The R-A-C-H share is what separates edited-specific sites (candidate m6A) from the control background and from single-replicate noise.
  - The centre must be C in ~100 % of sites; otherwise there is a coordinate or strand error.

## Known Bullseye issues (details: repo `docs/bullseye_known_issues.md`)
- **Stranded merge bug in `Find_edit_site.pl`** (not fixed). With stranded matrices, positions where a matrix has both `+` and `-` rows lose sites (a few percent of all sites; none are added). Unstranded runs are not affected. Say so whenever stranded site counts are reported or compared with unstranded ones.
- `Find_edit_site.pl --stranded` is a no-op. `parseBAM.pl get_strand` has swapped flag names but correct dUTP labels.
- Upstream has moved on (`Find_edit_sites.pl` 1.6.1, 2026-05) and keeps the same merge logic.

**Software:** the job scripts activate the repo's conda environments when `ENV_SETUP` is set (see `docs/INSTALL.md`; pass it with the other env vars), otherwise they use the cluster's `module load` lines. Ask the user which applies before the first run.

## Run rules
1. **Save the code that actually runs in the project.**
   - Copy `scripts/` into `<project>/bullseye/scripts/` and `md5sum` them into `MANIFEST.md5`.
   - Write `run_<step>_<YYYYMMDD>.sh`, which sets the env vars (DEDUP, thresholds, REFFLAT…) and calls the COPY.
   - Keep `pairs.tsv` and the sample list.
2. **Keep all logs** (`<out_dir>/error_file/`).
3. A failed job with a clear cause: delete only its outputs and logs and rerun (finished ones are skipped). Unclear or many failures: keep the logs, report, wait.
4. **Show scripts, pairs and the submit command, and wait for approval.** Run ONE sample/pair first and check its output, runtime and memory.
5. `command rm -f` (rm may be aliased to interactive). 6. Update the project STATUS.md.

## Typical commands
```bash
P=<project>/bullseye; SK=$P/scripts          # copy of the skill scripts
DEDUP=1 SGE_PROJECT=<project> bash $SK/run_parsebam.sh <project>/bam $P/matrix $P/samples.txt
python3 -I $SK/make_pairs.py $P/samples.txt Mut WT same_n > $P/pairs.tsv     # control tag, edited tag
REFFLAT=/path/to/ref/GRCh38_ens111/Homo_sapiens.GRCh38.111.fixed.refFlat MIN_EDIT=10 MAX_EDIT=95 FOLD=1.5 MIN_EDIT_SITES=2 \
  bash $SK/run_find_edit.sh $P/matrix $P/beds $P/pairs.tsv
bash $SK/summarize.sh 2 $P/beds/groupA_sum.bed $P/beds/groupA_*_WT.bed
```

## Checks
- Site-level QC, in this order:
  1. strand labels (`strand_check.sh`);
  2. replicate agreement (`summarize.sh`);
  3. independent-replicate validation in both directions;
  4. motif context (centre = C ~100 %, R-A-C-H enrichment of edited-specific sites).

  Only then call sites candidates.
- parseBAM: every sample has `.matrix.gz.tbi`; the error files have no `die` / `failed`; `zcat | head` of a matrix shows per-position A/C/G/T counts.
- Find_edit: every pair has a `.bed`. An empty or header-only bed means no site passed; check thresholds and depth before concluding anything.
- Chromosome names in the matrix must match the refFlat.

## Gotchas
- **qsub runs a spooled copy of the job script.** Never find helper scripts through `$0` / `dirname $0` inside a worker; pass the scripts dir as an argument (`validate_generic.job.sh … $S`, `motif_generic.job.sh … $S <fasta>`).
- `run_find_edit.sh` skips pairs whose matrix is missing when it submits. Submit Find_edit_site after parseBAM has finished, or use `HOLD_JID` only when the matrices already exist.
- The cluster may kill jobs that use more CPU than their slots: `--cpu` = `NSLOTS-2`. Check the first job's memory (`qacct` maxvmem) and raise the slots if it is killed.
- `parseBAM.pl` bgzips and tabix-indexes its output itself: pass `--output <name>.matrix`, not `.gz`.
- `Find_edit_site.pl` needs both matrices indexed. The batch script skips an existing `.bed`, so delete a bad one before rerunning.
- Pairs come from an explicit `pairs.tsv`, never from parsing names inside the batch loop.
- Thresholds tuned on one read length may not transfer to data with a different read length.
