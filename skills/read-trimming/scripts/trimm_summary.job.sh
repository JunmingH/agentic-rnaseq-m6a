#!/bin/bash -l
# qsub worker: $1 out_dir, $2 sample list, $3 scripts dir (qsub runs a spooled copy of this file, so $0 cannot be used to find it). Collects the per-sample Trimmomatic stats into trimm_summary.tsv
OUT=$1; LIST=$2; HERE=$3
python3 -I "$HERE/trimm_summary.py" "$OUT" "$LIST" > "$OUT/trimm_summary.tsv"
