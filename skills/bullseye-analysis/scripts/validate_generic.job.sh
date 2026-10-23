#!/bin/bash -l
# qsub worker: $1 group, $2 edited tag (WT|Mut), $3 bed dir, $4 matrix dir, $5 out tsv, $6 scripts dir
# (qsub runs a spooled copy of this script: the scripts dir is passed as an argument, never derived from $0)
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" bullseye || exit 1; else module load miniconda; conda activate Bullseye || exit 1; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
python3 -I $6/validate_rep3.py $1 $3 $4 $5 5 20 $2
