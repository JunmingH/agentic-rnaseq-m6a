#!/bin/bash -l
# qsub worker: $1 control matrix.gz, $2 edited matrix.gz, $3 output bed, $4 refFlat, $5 dir with Find_edit_site.pl
# env: MIN_EDIT MAX_EDIT FOLD MIN_EDIT_SITES MIN_COV STRANDED
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" bullseye || exit 1; else module load miniconda; conda activate Bullseye || exit 1; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
T=$(( ${NSLOTS:-8} - 2 )); [ "$T" -lt 1 ] && T=1
OPT=""; [ "$STRANDED" = 1 ] && OPT="--stranded"
rm -f "$3"
perl "$5/Find_edit_site.pl" --annotationFile "$4" --controlMatrix "$1" --EditedMatrix "$2" --outfile "$3" \
    --minEdit "$MIN_EDIT" --maxEdit "$MAX_EDIT" --editFoldThreshold "$FOLD" --MinEditSites "$MIN_EDIT_SITES" \
    --ControlMinCoverage "$MIN_COV" --EditedMinCoverage "$MIN_COV" --cpu "$T" $OPT --verbose \
    || { echo "Find_edit_site failed for $3"; rm -f "$3"; exit 1; }
echo "done $3"
