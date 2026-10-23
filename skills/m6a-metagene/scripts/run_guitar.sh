#!/bin/bash -l
# Usage: bash run_guitar.sh <diff_m6a_sites_dir> <out_dir>
#   1. per-group BED files (no-RAC / RAC) from the dart-diff-m6a site tables, checked against its filter_counts.tsv (runs here, light) -> <out_dir>/beds
#   2. GuitarPlot job -> <out_dir>/results/guitar_{norac,rac}_{overlay.png/pdf,plot_data.tsv,plot.rds}; logs <out_dir>/error_file
# env REQUIRED: GTF (the same annotation as the STAR index / refFlat)   optional: ORGANISM ("Homo sapiens"), SGE_PROJECT (qsub -P)
SD=$(readlink -f "$1"); OUT=$(readlink -f -m "$2")
[ -s "$SD/group_sites.tsv" ] && [ -n "$2" ] || { echo "Usage: $0 <diff_m6a_sites_dir> <out_dir>   (sites dir with group_sites.tsv)"; exit 1; }
[ -s "$GTF" ] || { echo "Set GTF= the annotation gtf used for STAR"; exit 1; }
HERE=$(dirname "$(readlink -f "$0")")
mkdir -p "$OUT/error_file" "$OUT/results"
set -o pipefail   # so that a failure of step 1 is not hidden by tee
python3 -I "$HERE/01_make_beds.py" "$SD" "$OUT/beds" | tee "$OUT/error_file/01_make_beds.log" || exit 1
qsub -terse -N guitar ${HOLD_JID:+-hold_jid $HOLD_JID} ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} -pe omp 4 -l h_rt=04:00:00 \
    -e "$OUT/error_file/02_error.txt" -o "$OUT/error_file/02_output.txt" \
    "$HERE/guitar.job.sh" "$HERE/02_guitar.R" "$GTF" "$OUT/beds" "$OUT/results" "${ORGANISM:-Homo sapiens}"
