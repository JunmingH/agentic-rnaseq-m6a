#!/bin/bash -l
# Usage: bash run_alignment.sh <in_dir> <out_dir> <sample_list> [r1_suffix] [r2_suffix]
#   in_dir       trimmed fastq dir      out_dir  BAM output dir (created if missing, never deleted)
#   sample_list  one sample name per line;  files = <sample><r1_suffix> / <sample><r2_suffix>
#   suffixes     default _paired_R1.fastq.gz / _paired_R2.fastq.gz  (fastp / Trimmomatic output names)
# env: RM_DUP (0|mark|remove), TRANSCRIPTOME (0|1|<file with sample names, only those get the transcriptome BAM>)  REQUIRED, set explicitly - the user decides, no silent default
#      READ_LEN   max read length of the data (default 150) -> index $INDEX_DIR or $REF_DIR/index_<READ_LEN>, sjdbOverhang READ_LEN-1
#      INDEX_DIR  use this index instead
#      STAR_BIN
# The index must exist, or a build job named star_index must be queued/running (then every alignment job is held until it finishes).
# This script never builds an index itself.
# Samples whose final BAM (+.bai) already exists are skipped.
IN=$(readlink -f "$1"); OUT=$(readlink -f -m "$2"); LIST=$(readlink -f "$3")
R1S=${4:-_paired_R1.fastq.gz}; R2S=${5:-_paired_R2.fastq.gz}
[ -d "$IN" ] && [ -f "$LIST" ] || { echo "Usage: $0 <in_dir> <out_dir> <sample_list> [r1_suffix] [r2_suffix]"; exit 1; }
HERE=$(dirname "$(readlink -f "$0")")
[ -n "$REF_DIR" ] || [ -n "$INDEX_DIR" ] || { echo "Set INDEX_DIR= the STAR index (or REF_DIR= the folder holding index_<READ_LEN>)"; exit 1; }
REF=$REF_DIR
READ_LEN=${READ_LEN:-150}
IDX=${INDEX_DIR:-$REF/index_$READ_LEN}
case "$RM_DUP" in 0|mark|remove) ;; *) echo "Set RM_DUP=0, mark or remove explicitly; ask the user"; exit 1;; esac
case "$TRANSCRIPTOME" in 0|1) ;; *) [ -f "$TRANSCRIPTOME" ] || { echo "Set TRANSCRIPTOME=1, 0, or the path of a file listing the samples that need Aligned.toTranscriptome.out.bam; ask the user"; exit 1; };; esac
mkdir -p "$OUT/error_file"

hold=""
if [ -s "$IDX/SA" ] && [ -s "$IDX/genomeParameters.txt" ]; then
    echo "index: $IDX (sjdbOverhang $(grep -P '^sjdbOverhang' "$IDX/genomeParameters.txt" | cut -f2)), version $(grep -P '^versionGenome' "$IDX/genomeParameters.txt" | cut -f2)"
else
    [ -e "$IDX" ] && { echo "$IDX exists but is incomplete (no SA); remove or fix it first"; exit 1; }
    # index build is a separate step (build_star_index.sh, saved next to the index). A build job named star_index that is still queued/running: wait for it.
    bid=$(qstat -u "$USER" 2>/dev/null | awk '$3=="star_index" {print $1; exit}')
    if [ -n "$bid" ]; then echo "index build job $bid still queued/running -> alignment jobs wait for it"; hold="-hold_jid $bid"
    else echo "index $IDX does not exist and no star_index job is running: build it first (see the star-alignment skill), then rerun"; exit 1; fi
fi
VARS="RM_DUP=$RM_DUP"; [ -n "$STAR_BIN" ] && VARS="$VARS,STAR_BIN=$STAR_BIN"
case "$RM_DUP" in mark) FINAL=_Aligned.sortedByCoord.out.markdup.bam;; remove) FINAL=_Aligned.sortedByCoord.out.rm_dup.bam;; *) FINAL=_Aligned.sortedByCoord.out.bam;; esac
ids=""; n=0; skipped=0; missing=0
while read -r s; do
    [ -z "$s" ] && continue
    if [ ! -e "$IN/$s$R1S" ] || [ ! -e "$IN/$s$R2S" ]; then echo "MISSING input for $s"; missing=$((missing+1)); continue; fi
    if [ -s "$OUT/$s$FINAL.bai" ]; then echo "SKIP $s (BAM exists)"; skipped=$((skipped+1)); continue; fi
    case "$TRANSCRIPTOME" in 0|1) TQ=$TRANSCRIPTOME;; *) grep -qx "$s" "$TRANSCRIPTOME" && TQ=1 || TQ=0;; esac
    id=$(qsub -terse -N "aln_${s:0:20}" $hold -l h_rt=24:00:00 ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} -pe omp 16 \
        -e "$OUT/error_file/${s}_error_align.txt" -o "$OUT/error_file/${s}_output_align.txt" \
        -v "$VARS,TRANSCRIPTOME=$TQ" "$HERE/align.job.sh" "$s" "$IN" "$OUT" "$R1S" "$R2S" "$IDX")
    echo "$s -> job $id"; ids="${ids:+$ids,}$id"; n=$((n+1))
done < "$LIST"
echo "submitted $n, skipped $skipped, missing $missing"
[ -z "$ids" ] && exit 0
mid=$(qsub -terse -N align_summary -hold_jid "$ids" -l h_rt=00:30:00 ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} \
    -e "$OUT/error_file/align_summary_error.txt" -o "$OUT/error_file/align_summary_output.txt" \
    "$HERE/align_summary.job.sh" "$OUT" "$LIST" "$HERE")
echo "summary job $mid -> $OUT/align_summary.tsv (after all alignment jobs)"
