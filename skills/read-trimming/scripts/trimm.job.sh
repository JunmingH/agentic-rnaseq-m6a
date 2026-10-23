#!/bin/bash -l
# qsub worker: $1 sample, $2 in_dir, $3 out_dir, $4 r1 suffix, $5 r2 suffix
# env: ADAPTER (adapter file name), ADAPTER_DIR (folder with the Trimmomatic adapter fasta files, e.g. <trimmomatic install>/adapters), MINLEN, ILLUMINACLIP (optional full string)
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" tools || exit 1; else module load trimmomatic/0.39; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
[ -z "$ADAPTER_DIR" ] && [ -n "$CONDA_PREFIX" ] && ADAPTER_DIR=$CONDA_PREFIX/share/trimmomatic/adapters   # bioconda location
s=$1; IN=$2; OUT=$3; R1S=$4; R2S=$5
CLIP=${ILLUMINACLIP:-ILLUMINACLIP:$ADAPTER_DIR/${ADAPTER:-TruSeq3-PE-2.fa}:2:30:10:8:true}
# ILLUMINACLIP seedMismatches:palindromeClipThreshold:simpleClipThreshold:minAdapterLength:keepBothReads
# keep 1 slot free: reader/writer/GC threads come on top of -threads (the cluster kills jobs that use more CPU than their slots)
T=$(( ${NSLOTS:-4} - 1 )); [ "$T" -lt 1 ] && T=1
trimmomatic PE -threads "$T" -summary "$OUT/${s}_trimm_summary.txt" \
    "$IN/$s$R1S" "$IN/$s$R2S" \
    "$OUT/${s}_paired_R1.fastq.gz" "$OUT/${s}_unpaired_R1.fastq.gz" \
    "$OUT/${s}_paired_R2.fastq.gz" "$OUT/${s}_unpaired_R2.fastq.gz" \
    $CLIP LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:${MINLEN:-50}
