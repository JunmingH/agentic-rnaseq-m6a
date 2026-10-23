#!/bin/bash -l
# qsub worker: $1 out_dir, $2 sample_list, $3 scripts dir
python3 -I "$3/align_summary.py" "$1" "$2"
