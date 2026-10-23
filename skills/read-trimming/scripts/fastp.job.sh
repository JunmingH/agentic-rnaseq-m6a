#!/bin/bash -l
# qsub worker (fastp): $1 sample, $2 in_dir, $3 out_dir, $4 r1 suffix, $5 r2 suffix
# env: MINLEN, ADAPTER_R1 / ADAPTER_R2 (optional; default = auto-detect)
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" tools || exit 1; else module load fastp/0.23.4; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
s=$1; IN=$2; OUT=$3; R1S=$4; R2S=$5
# The cluster kills a job that uses more CPU than its slots. fastp --thread N is only the worker count; reader/writer threads come on top
# (--thread 4 used 5.1 cores on a 4-slot job and was killed after several minutes), so keep 2 slots free.
T=$(( ${NSLOTS:-4} - 2 )); [ "$T" -lt 1 ] && T=1; [ "$T" -gt 14 ] && T=14
if [ -n "$ADAPTER_R1" ]; then AD="--adapter_sequence $ADAPTER_R1 --adapter_sequence_r2 ${ADAPTER_R2:-$ADAPTER_R1}"; else AD="--detect_adapter_for_pe"; fi
# Settings mirror the Trimmomatic run (LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN):
#   --cut_front  window 1 mean Q3   (LEADING:3)      --cut_tail  window 1 mean Q3  (TRAILING:3)
#   --cut_right  window 4 mean Q20  (SLIDINGWINDOW)  --length_required MINLEN
# plus poly-G tail trimming (what Trimmomatic cannot do). fastp's own defaults stay on: reads with >40% bases below Q15 and
# reads with >5 N are dropped.
fastp -i "$IN/$s$R1S" -I "$IN/$s$R2S" \
    -o "$OUT/${s}_paired_R1.fastq.gz" -O "$OUT/${s}_paired_R2.fastq.gz" \
    --unpaired1 "$OUT/${s}_unpaired_R1.fastq.gz" --unpaired2 "$OUT/${s}_unpaired_R2.fastq.gz" \
    $AD --trim_poly_g \
    --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 \
    --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 \
    --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
    --length_required "${MINLEN:-50}" --thread "$T" \
    -j "$OUT/${s}_fastp.json" -h "$OUT/${s}_fastp.html" -R "$s"
