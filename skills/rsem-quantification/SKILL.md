---
name: rsem-quantification
description: Quantify gene and transcript expression (linear RNA level) with RSEM on an SGE cluster from the STAR transcriptome BAMs (Aligned.toTranscriptome.out.bam), one SGE job per sample, then merge into gene/isoform count and TPM matrices. Use when the user asks to run RSEM / quantify expression / get counts or TPM from the BAMs (e.g. the control group of a DART-seq project), after the star-alignment step with TRANSCRIPTOME enabled.
---

# RSEM quantification (SGE)

## What it does
1. **Per sample** (`run_rsem.sh` → `rsem.job.sh`, 12 slots, `h_rt=24:00:00`). Runs `rsem-calculate-expression --paired-end --bam --append-names --strandedness <S> (--output-genome-bam | --no-bam-output) <s>_Aligned.toTranscriptome.out.bam <ref> <out>/<s>`, which writes `<s>.genes.results`, `<s>.isoforms.results` and `<s>.stat/`.
   - Modules: `R/4.2.1 rsem/1.3.3 bowtie2/2.3.4.1 python2/2.7.16`, the same ones that built the reference.
   - `-p` = `NSLOTS-2`, because the cluster may kill jobs that use more CPU than their slots.
2. **Matrices** (`make_matrix.py <out_dir> <samples> <prefix>`): gene and isoform × `expected_count` / `TPM`, columns in sample-list order. `expected_count` is not an integer, so DESeq2 needs `round()` (as deseq2-expression does) or tximport.

## Inputs
- **Reference** (`REF`, required): an RSEM reference prefix (`rsem-prepare-reference` output, e.g. `/path/to/ref/GRCh38_ens111/rsem_ref/human_ensembl`, with `.ti .grp .transcripts.fa`). It must be the SAME annotation as the STAR index that made the BAMs.
  - The BAM `@SQ` names are Ensembl transcript ids without version, exactly as in the reference.
  - Check: `samtools view -H <bam> | grep -c '^@SQ'` = the number of transcripts in the reference.
- **BAMs**: only samples aligned with `--quantMode TranscriptomeSAM` have `<s>_Aligned.toTranscriptome.out.bam`. A missing BAM is reported as MISSING and not submitted.
- **Sample list**: one name per line.

## Options the user must decide (the script refuses to run unless set; ask, explain, recommend)
- **`STRANDEDNESS` none|forward|reverse.** It tells RSEM which strand read 1 comes from:
  - `forward` = read 1 is on the RNA's sense strand;
  - `reverse` = read 1 is antisense, as in the usual dUTP / Illumina stranded kits.

  Do NOT take the user's word "forward stranded" as the RSEM setting, because tools and people use the words differently. **Always test it on ONE sample before the batch** (see "Strandedness test"). If untested, use `none`.
- **`GENOME_BAM` 0|1.** `1` also writes a genome-coordinate BAM (large), which is only needed for visualising the RSEM result. Recommend `0`.

## Strandedness test (before every new library type or project)
Run the SAME transcriptome BAM of one sample with `none`, `reverse` and `forward` (separate out dirs, about 20 min each at 12 slots), then compare:
- **Sum of `expected_count` over genes.** The right setting keeps about 95–99 % of the `none` sum. The wrong one keeps only a few percent, because it discards reads of the wrong orientation. `none` keeps 100 % but also counts antisense reads.
- **MT-ND6** (the only mitochondrial protein-coding gene on the L strand, antisense-overlapping MT-ND5). With the wrong setting its TPM jumps far above all other MT genes. With the right one, the MT-ND6/MT-ND5 TPM ratio is far below 1.

Lesson: a library described as "forward stranded" turned out to be `reverse` in RSEM terms. `forward` kept only a few percent of the reads and made MT-ND6 the top gene. Choose the setting from the test, then run the batch.

**Software:** the job scripts activate the repo's conda environments when `ENV_SETUP` is set (see `docs/INSTALL.md`; pass it with the other env vars), otherwise they use the cluster's `module load` lines. Ask the user which applies before the first run.

## Run rules
1. **Save the code that actually runs in the project.** Copy `scripts/` to `<project>/rsem/scripts/`, write `md5sum` into `MANIFEST.md5`, and write `<project>/rsem/run_rsem_<YYYYMMDD>.sh`, which sets REF, STRANDEDNESS and GENOME_BAM and calls the COPY. Keep the sample list.
2. **Keep all logs** (`<out_dir>/error_file/`).
3. A failed job with a clear cause: delete only its outputs and logs, then rerun (finished samples are skipped). Unclear or many failures: keep the logs, report, wait.
4. **Show the scripts and the submit command, and wait for approval.** Run ONE sample first and check its results, time and memory (`qacct -j`).
5. `command rm -f` (rm may be aliased to interactive). 6. Update the project STATUS.md.

## Typical commands
```bash
P=<project>/rsem; SK=$P/scripts
REF=/path/to/ref/GRCh38_ens111/rsem_ref/human_ensembl STRANDEDNESS=reverse GENOME_BAM=0 SGE_PROJECT=<project> \
  bash $SK/run_rsem.sh <project>/bam $P/results $P/samples_control.txt
python3 -I $SK/make_matrix.py $P/results $P/samples_control.txt $P/matrix
```

## Checks
- **Strand sanity.** Each sample's summed expected_count should be a large share of `.cnt` N_tot; a few percent means the wrong `STRANDEDNESS`. The MT-ND6/MT-ND5 TPM ratio should be far below 1.
- Every sample has `.genes.results` + `.isoforms.results`, and the error files contain no `Error` / `failed`.
- `<s>.stat/<s>.cnt` looks sane. Its first line is N0 (unalignable), N1 (unique), N2 (multi), N_tot. With STAR `--outSAMunmapped Within`, N0 = the reads STAR did not align, which are also written into the transcriptome BAM.
- If the library contains non-host reads (contamination), RSEM only sees host-aligned reads and TPM is normalised to them. Keep that in mind when sequencing depth varies a lot.
- Known marker genes should be among the top expressed. A sample with a very different top list is suspect (possible label swap).

## Gotchas
- RSEM `--bam` with paired-end data needs the two mates of a pair adjacent in the file. STAR's transcriptome BAM is written that way; do not sort it.
- Do not mix transcriptome BAMs from a different STAR index or annotation with this reference; mismatching `@SQ` names fail immediately.
- A killed job leaves `<s>.temp/` and a partial `.stat`; the worker removes them before rerunning that sample.
- `--append-names` gives ids like `ENSG..._TSPAN6`; strip the suffix before joining with other tables.
