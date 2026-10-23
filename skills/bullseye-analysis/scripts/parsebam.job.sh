#!/bin/bash -l
# qsub worker: $1 bam, $2 output prefix (-> $2.matrix.gz/.tbi), $3 dir with parseBAM.pl
# env: DEDUP (1 = --removeDuplicates), MIN_COV, STRANDED (1 = --stranded)
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" bullseye || exit 1; else module load miniconda; conda activate Bullseye || exit 1; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
# parseBAM runs bgzip with -@ cpu and sort threads; The cluster kills jobs above their slots, so keep 2 free.
T=$(( ${NSLOTS:-16} - 2 )); [ "$T" -lt 1 ] && T=1
OPT=""; [ "$DEDUP" = 1 ] && OPT="$OPT --removeDuplicates"; [ "$STRANDED" = 1 ] && OPT="$OPT --stranded"
rm -f "$2.matrix" "$2.matrix.gz" "$2.matrix.gz.tbi"     # leftovers of a killed run of this same sample
perl "$3/parseBAM.pl" --input "$1" --output "$2.matrix" --cpu "$T" --minCoverage "$MIN_COV" $OPT --verbose \
    || { echo "parseBAM failed for $1"; exit 1; }
[ -s "$2.matrix.gz.tbi" ] || { echo "no $2.matrix.gz.tbi"; exit 1; }
echo "done $1"
