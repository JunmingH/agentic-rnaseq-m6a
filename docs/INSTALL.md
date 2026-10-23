# Installing the software and references

The job scripts work in two ways:

- **With `ENV_SETUP` (recommended outside the original cluster).** Every job runs `. "$ENV_SETUP" <tools|r|bullseye>`, which activates one of the three conda environments in `envs/`.
- **Without `ENV_SETUP`.** The scripts fall back to the `module load ...` lines of the cluster they were developed on. Those module names will not exist elsewhere.

The scripts submit jobs with SGE `qsub` (`-pe omp <slots>`, `-hold_jid`, `-v`). On SLURM or another scheduler, the `qsub` lines in `scripts/run_*.sh` must be adapted; the job scripts themselves do not depend on the scheduler, apart from `$NSLOTS`.

## 1. Software (conda)

```bash
# Miniforge (conda + mamba): https://github.com/conda-forge/miniforge
mamba env create -f envs/dartseq-tools.yml      # FastQC, MultiQC, fastp, Trimmomatic, STAR, samtools, bcftools, Picard, RSEM, python
mamba env create -f envs/dartseq-r.yml          # R 4.2 + DESeq2, apeglm, IHW, Guitar, GenomicFeatures, ggplot2, ggrepel, openxlsx, data.table
mamba env create -f envs/dartseq-bullseye.yml   # Perl + BioPerl, MCE, Array::IntSpan, Math::CDF, samtools/tabix for Bullseye

cp envs/env_setup.example.sh ~/dartseq_env.sh   # set CONDA_BASE inside if conda is not in ~/miniforge3
export ENV_SETUP=$HOME/dartseq_env.sh           # absolute path; it is passed to every job (qsub -v)
```

The environments and the `ENV_SETUP` file must be on storage that the compute nodes can see (home or project space, not a node-local `/scratch` or `/tmp`). Otherwise the jobs stop at once with `EnvironmentLocationNotFound`. Installing into a shared prefix works too: `mamba env create -p /path/to/shared/envs/dartseq-tools -f envs/dartseq-tools.yml`, then activate by path in your `ENV_SETUP` copy.

`nodefaults` in the channel list matters. If Anaconda's `defaults` channel is mixed in, the solver can pick an `r-stringi` built against another ICU, and `library(GenomicFeatures)` then fails with `libicui18n.so.58: cannot open shared object file`.

Check the setup:
```bash
. $ENV_SETUP tools    && STAR --version && samtools --version | head -1 && picard MarkDuplicates --version 2>&1 | tail -1
. $ENV_SETUP r        && Rscript -e 'library(DESeq2); library(IHW); library(Guitar); cat("R ok\n")'
. $ENV_SETUP bullseye && perl -MMCE::Loop -MArray::IntSpan -MBio::DB::Fasta -MMath::CDF -e 'print "perl ok\n"'
bash docs/repro_find_edit_stranded_merge/run_repro.sh skills/bullseye-analysis/scripts/Find_edit_site.pl   # Bullseye runs (expected: A 1 of 4, B 4 of 4)
```

**Versions.** `envs/*.yml` pin the versions the pipeline was run with, with these exceptions:

| package | used for validation | conda env | why |
|---|---|---|---|
| R / Bioconductor | 4.2.1 / 3.15 (DESeq2 1.36.0, IHW 1.24.0, Guitar 2.12.0, apeglm 1.18.0) | 4.2 / 3.16 (DESeq2 1.38.0, IHW 1.26.0, Guitar 2.14.0, apeglm 1.20.0) | bioconda has no Bioconductor 3.15 builds of these packages |
| Perl / BioPerl | 5.26.2 / 1.7.2 | 5.32 / 1.7.8 | bioconda's Array::IntSpan needs Perl ≥ 5.32 |
| Picard | 2.25.2 | 2.27.5 | nearest bioconda build; MarkDuplicates options unchanged |
| FastQC / MultiQC | cluster default | 0.12.1 / 1.21 | QC only |

RSEM 1.3.3 prints `Current version: RSEM v1.3.1` (upstream version string, the same with the cluster module).

[`validation.md`](validation.md) lists what was re-checked with the conda environments.

## 2. References (one Ensembl release for everything)

Example: human GRCh38, Ensembl 111. Use the same release for every file.

```bash
R=/path/to/ref/GRCh38_ens111; mkdir -p $R && cd $R
E=https://ftp.ensembl.org/pub/release-111
curl -O $E/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
curl -O $E/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
gunzip *.gz; mv Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa grch38.fa
. $ENV_SETUP tools; samtools faidx grch38.fa
```

| file | how | used by |
|---|---|---|
| STAR index | `bash skills/star-alignment/scripts/build_star_index.sh grch38.fa Homo_sapiens.GRCh38.111.gtf $R/index_150 150` (as a qsub job, 16 slots, ~40 GB RAM) | star-alignment (`INDEX_DIR`) |
| RSEM reference | `rsem-prepare-reference --gtf Homo_sapiens.GRCh38.111.gtf grch38.fa $R/rsem_ref/human_ensembl` | rsem-quantification (`REF`) |
| refFlat (patched) | `git clone https://github.com/mflamand/Bullseye`; `cp Bullseye/Code/gtf2genepred.pl gtf2genepred_fixed.pl`; `patch gtf2genepred_fixed.pl < patches/gtf2genepred_refflat_fix.patch`; `. $ENV_SETUP bullseye; perl gtf2genepred_fixed.pl --gtf Homo_sapiens.GRCh38.111.gtf --out Homo_sapiens.GRCh38.111.fixed.refFlat` | bullseye-analysis (`REFFLAT`) |
| dbSNP C>T/G>A bed | `OUT=$R/dbsnp_ensembl111 RELEASE=111 bash skills/dart-diff-m6a/scripts/download_dbsnp_ensembl.sh` (needs internet, several hours) | dart-diff-m6a (`DBSNP`) |
| genome fasta + `.fai` | above | dart-diff-m6a (`FASTA`, RAC motif), bullseye motif check |
| GTF | above | m6a-metagene (`GTF`) |
| Trimmomatic adapters | inside the env: `$CONDA_PREFIX/share/trimmomatic/adapters` (the job script uses it when `ADAPTER_DIR` is unset; `run_trimm.sh` asks for `ADAPTER_DIR` explicitly) | read-trimming |

Check the refFlat against the GTF: `python3 skills/pipeline-audit/scripts/check_annotation_vs_gtf.py Homo_sapiens.GRCh38.111.gtf Homo_sapiens.GRCh38.111.fixed.refFlat`.

## 3. Environment variables per step

All of them are explicit. A dated run script in the project folder should set them (see each `SKILL.md`).

```bash
export ENV_SETUP=$HOME/dartseq_env.sh SGE_PROJECT=<your qsub project, if the cluster needs one>
R=/path/to/ref/GRCh38_ens111
# star-alignment:      INDEX_DIR=$R/index_150 READ_LEN=150 RM_DUP=mark TRANSCRIPTOME=<0|1|sample file>
# rsem-quantification: REF=$R/rsem_ref/human_ensembl STRANDEDNESS=<after the one-sample test> GENOME_BAM=0
# bullseye-analysis:   REFFLAT=$R/Homo_sapiens.GRCh38.111.fixed.refFlat DEDUP=1 MIN_EDIT=10 MAX_EDIT=95 FOLD=1.5 MIN_EDIT_SITES=2
# dart-diff-m6a:       FASTA=$R/grch38.fa DBSNP=$R/dbsnp_ensembl111/dbsnp_ensembl111.CT_GA.snv.bed
# m6a-metagene:        GTF=$R/Homo_sapiens.GRCh38.111.gtf
# read-trimming:       TOOL=fastp   (or TOOL=trimmomatic ADAPTER_DIR=<env>/share/trimmomatic/adapters)
```
