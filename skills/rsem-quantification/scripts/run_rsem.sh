#!/bin/bash -l
# Usage: bash run_rsem.sh <bam_dir> <out_dir> <sample_list> [bam_suffix]
#   bam_dir      dir with <sample><bam_suffix>; suffix default _Aligned.toTranscriptome.out.bam (STAR --quantMode TranscriptomeSAM)
#   out_dir      results dir (created if missing, never deleted): <s>.genes.results, <s>.isoforms.results, <s>.stat/
#   sample_list  one sample name per line
# env (REQUIRED, the user decides): STRANDEDNESS=none|forward|reverse   GENOME_BAM=0|1 (1 = also write the genome-coordinate BAM, large)
# env (REQUIRED): REF = RSEM reference prefix (rsem-prepare-reference output, e.g. /path/to/ref/rsem_ref/human_ensembl)
# env (optional): SLOTS (default 12), SGE_PROJECT (qsub -P)
# Samples whose <s>.genes.results and <s>.isoforms.results exist are skipped.
BAMD=$(readlink -f "$1"); OUT=$(readlink -f -m "$2"); LIST=$(readlink -f "$3")
SUF=${4:-_Aligned.toTranscriptome.out.bam}
[ -d "$BAMD" ] && [ -f "$LIST" ] || { echo "Usage: $0 <bam_dir> <out_dir> <sample_list> [bam_suffix]"; exit 1; }
case "$STRANDEDNESS" in none|forward|reverse) ;; *) echo "Set STRANDEDNESS=none, forward or reverse explicitly; ask the user"; exit 1;; esac
case "$GENOME_BAM" in 0|1) ;; *) echo "Set GENOME_BAM=0 or 1 explicitly; ask the user"; exit 1;; esac
HERE=$(dirname "$(readlink -f "$0")")
[ -n "$REF" ] || { echo "Set REF= the RSEM reference prefix"; exit 1; }
[ -s "$REF.ti" ] && [ -s "$REF.grp" ] && [ -s "$REF.transcripts.fa" ] || { echo "RSEM reference $REF incomplete"; exit 1; }
mkdir -p "$OUT/error_file"
n=0; skipped=0; missing=0
while read -r s; do
    [ -z "$s" ] && continue
    if [ ! -s "$BAMD/$s$SUF" ]; then echo "MISSING $BAMD/$s$SUF"; missing=$((missing+1)); continue; fi
    if [ -s "$OUT/$s.genes.results" ] && [ -s "$OUT/$s.isoforms.results" ]; then echo "SKIP $s (results exist)"; skipped=$((skipped+1)); continue; fi
    id=$(qsub -terse -N "rsem_${s:0:20}" -l h_rt=24:00:00 ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} -pe omp ${SLOTS:-12} \
        -e "$OUT/error_file/${s}_error_rsem.txt" -o "$OUT/error_file/${s}_output_rsem.txt" \
        -v "STRANDEDNESS=$STRANDEDNESS,GENOME_BAM=$GENOME_BAM" \
        "$HERE/rsem.job.sh" "$BAMD/$s$SUF" "$OUT/$s" "$REF")
    echo "$s -> job $id"; n=$((n+1))
done < "$LIST"
echo "submitted $n, skipped $skipped, missing $missing"
