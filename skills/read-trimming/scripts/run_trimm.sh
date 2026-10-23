#!/bin/bash -l
# Usage: bash run_trimm.sh <in_dir> <out_dir> <sample_list> [r1_suffix] [r2_suffix]
#   <in_dir>       directory with the raw fastq files (symlinks are fine)
#   <out_dir>      where trimmed files and logs go (created if missing, never deleted)
#   <sample_list>  one sample name per line, e.g. ctrl_1_Mut   (file = <sample><r1_suffix> / <sample><r2_suffix>)
#   r1/r2 suffix   default _1.fq.gz / _2.fq.gz   (the old pipeline used _R1.fastq.gz / _R2.fastq.gz)
# Optional env vars: TOOL        trimmomatic (default) | fastp
#                    MINLEN      minimum read length after trimming (default 50)
#                    ADAPTER     [trimmomatic] adapter file in the module's adapters dir (default TruSeq3-PE-2.fa)
#                    ILLUMINACLIP [trimmomatic] full clip string; overrides ADAPTER
#                    ADAPTER_R1 / ADAPTER_R2  [fastp] adapter sequences; default = fastp auto-detection
# Submits one qsub job per sample (job name trm_*), then one job `trimm_summary` held until all of them finish.
# Samples whose paired outputs already exist are skipped (nothing is ever rm -rf'd).

IN=$(readlink -f "$1"); OUT=$(readlink -f -m "$2"); LIST=$(readlink -f "$3")
R1S=${4:-_1.fq.gz}; R2S=${5:-_2.fq.gz}
[ -d "$IN" ] && [ -f "$LIST" ] || { echo "Usage: $0 <in_dir> <out_dir> <sample_list> [r1_suffix] [r2_suffix]"; exit 1; }
HERE=$(dirname "$(readlink -f "$0")")
mkdir -p "$OUT/error_file"

TOOL=${TOOL:-trimmomatic}
case "$TOOL" in trimmomatic) WORKER=trimm.job.sh;; fastp) WORKER=fastp.job.sh;; *) echo "TOOL must be trimmomatic or fastp"; exit 1;; esac
VARS="MINLEN=${MINLEN:-50}"
if [ "$TOOL" = trimmomatic ] && [ -z "$ILLUMINACLIP" ] && [ ! -s "$ADAPTER_DIR/${ADAPTER:-TruSeq3-PE-2.fa}" ]; then echo "Set ADAPTER_DIR= the folder with the Trimmomatic adapter fasta (${ADAPTER:-TruSeq3-PE-2.fa})"; exit 1; fi
if [ "$TOOL" = trimmomatic ]; then VARS="$VARS,ADAPTER=${ADAPTER:-TruSeq3-PE-2.fa}"; [ -n "$ADAPTER_DIR" ] && VARS="$VARS,ADAPTER_DIR=$ADAPTER_DIR"; [ -n "$ILLUMINACLIP" ] && VARS="$VARS,ILLUMINACLIP=$ILLUMINACLIP"
else [ -n "$ADAPTER_R1" ] && VARS="$VARS,ADAPTER_R1=$ADAPTER_R1"; [ -n "$ADAPTER_R2" ] && VARS="$VARS,ADAPTER_R2=$ADAPTER_R2"; fi
echo "tool: $TOOL"
ids=""; n=0; skipped=0; missing=0
while read -r s; do
    [ -z "$s" ] && continue
    if [ ! -e "$IN/$s$R1S" ] || [ ! -e "$IN/$s$R2S" ]; then echo "MISSING input for $s"; missing=$((missing+1)); continue; fi
    if [ -s "$OUT/${s}_paired_R1.fastq.gz" ] && [ -s "$OUT/${s}_paired_R2.fastq.gz" ]; then echo "SKIP $s (outputs exist)"; skipped=$((skipped+1)); continue; fi
    id=$(qsub -terse -N "trm_${s:0:20}" -l h_rt=10:00:00 ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} -pe omp 4 \
        -e "$OUT/error_file/${s}_error_trimm.txt" -o "$OUT/error_file/${s}_output_trimm.txt" \
        -v "$VARS" \
        "$HERE/$WORKER" "$s" "$IN" "$OUT" "$R1S" "$R2S")
    echo "$s -> job $id"; ids="${ids:+$ids,}$id"; n=$((n+1))
done < "$LIST"
echo "submitted $n, skipped $skipped, missing $missing"
[ -z "$ids" ] && exit 0
mid=$(qsub -terse -N trimm_summary -hold_jid "$ids" -l h_rt=00:30:00 ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} \
    -e "$OUT/error_file/trimm_summary_error.txt" -o "$OUT/error_file/trimm_summary_output.txt" \
    "$HERE/trimm_summary.job.sh" "$OUT" "$LIST" "$HERE")
echo "summary job $mid -> $OUT/trimm_summary.tsv (starts after all trimming jobs finish)"
