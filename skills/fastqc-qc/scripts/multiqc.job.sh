#!/bin/bash -l
# qsub worker. $1 = fastq dir; summarises everything in $1/fastqc into $1/fastqc/multiqc_report.html
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" tools || exit 1; else module load multiqc; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
cd "$1/fastqc" || exit 1
multiqc . -f
