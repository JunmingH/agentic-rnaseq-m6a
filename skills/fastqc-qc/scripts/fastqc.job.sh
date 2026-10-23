#!/bin/bash -l
# qsub worker. $1 = fastq file (full path), $2 = fastq dir (output -> $2/fastqc)
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" tools || exit 1; else module load fastqc; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
fastqc "$1" -o "$2/fastqc" -t "${NSLOTS:-1}"
