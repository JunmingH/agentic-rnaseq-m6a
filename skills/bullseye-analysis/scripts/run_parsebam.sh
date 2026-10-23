#!/bin/bash -l
# Usage: bash run_parsebam.sh <bam_dir> <out_dir> <sample_list> [bam_suffix]
#   bam_dir      dir with <sample><bam_suffix> (+ .bai); bam_suffix default _Aligned.sortedByCoord.out.markdup.bam
#   out_dir      matrix output dir (created if missing, never deleted); <sample>.matrix.gz + .tbi
#   sample_list  one sample name per line
# env (REQUIRED, the user decides): DEDUP=0|1   1 = parseBAM --removeDuplicates (ignore reads flagged 0x400; needs a markdup/rm_dup BAM)
# env (optional): MIN_COV (default 10), STRANDED=1 (adds --stranded), SLOTS (default 16)
# Samples whose <sample>.matrix.gz.tbi exists are skipped.
BAMD=$(readlink -f "$1"); OUT=$(readlink -f -m "$2"); LIST=$(readlink -f "$3")
SUF=${4:-_Aligned.sortedByCoord.out.markdup.bam}
[ -d "$BAMD" ] && [ -f "$LIST" ] || { echo "Usage: $0 <bam_dir> <out_dir> <sample_list> [bam_suffix]"; exit 1; }
case "$DEDUP" in 0|1) ;; *) echo "Set DEDUP=0 or 1 explicitly; ask the user"; exit 1;; esac
HERE=$(dirname "$(readlink -f "$0")")
mkdir -p "$OUT/error_file"
n=0; skipped=0; missing=0
while read -r s; do
    [ -z "$s" ] && continue
    if [ ! -s "$BAMD/$s$SUF" ] || [ ! -s "$BAMD/$s$SUF.bai" ]; then echo "MISSING bam/bai for $s"; missing=$((missing+1)); continue; fi
    if [ -s "$OUT/$s.matrix.gz.tbi" ]; then echo "SKIP $s (matrix exists)"; skipped=$((skipped+1)); continue; fi
    id=$(qsub -terse -N "pb_${s:0:20}" -l h_rt=24:00:00 ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} -pe omp ${SLOTS:-16} \
        -e "$OUT/error_file/${s}_error_parsebam.txt" -o "$OUT/error_file/${s}_output_parsebam.txt" \
        -v "DEDUP=$DEDUP,MIN_COV=${MIN_COV:-10},STRANDED=${STRANDED:-0}" \
        "$HERE/parsebam.job.sh" "$BAMD/$s$SUF" "$OUT/$s" "$HERE")
    echo "$s -> job $id"; n=$((n+1))
done < "$LIST"
echo "submitted $n, skipped $skipped, missing $missing"
