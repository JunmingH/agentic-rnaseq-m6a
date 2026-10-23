#!/bin/bash -l
# qsub worker: $1 02_guitar.R (absolute path), $2 gtf, $3 bed dir, $4 out dir, $5 organism. R 4.2.1 + Guitar 2.12.0 / GenomicFeatures.
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" r || exit 1; else module load R/4.2.1; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
Rscript "$1" "$2" "$3" "$4" "$5" || { echo "guitar failed"; exit 1; }
