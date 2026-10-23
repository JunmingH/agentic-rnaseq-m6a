#!/bin/bash -l
# qsub worker: $1 R script (absolute path), then its arguments. R 4.2.1 + DESeq2 1.36.0 / apeglm / ggrepel / openxlsx (cluster module).
# env ALPHA LFC MIN_COUNT MIN_SAMPLES NTOP are passed through by run_deseq2.sh (-v)
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" r || exit 1; else module load R/4.2.1; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
R=$1; shift
Rscript "$R" "$@" || { echo "$(basename "$R") failed"; exit 1; }
