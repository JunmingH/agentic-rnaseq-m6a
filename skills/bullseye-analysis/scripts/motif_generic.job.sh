#!/bin/bash -l
# qsub worker: $1 bed dir (main direction beds), $2 validate dir (<g>.tsv, <g>_rev.tsv), $3 out tsv, $4 scripts dir (passed as argument, not derived from $0), $5 genome fasta (+ .fai)
python3 -I $4/motif_context.py $1 $2 $5 $3
