#!/bin/bash -l
# Usage: bash run_deseq2.sh <gene_expected_count.tsv> <samplesheet.tsv> <comparisons.tsv> <out_dir>
#   submits two SGE jobs: DESeq2 differential expression (-> <out_dir>/results) and the standard PCA (-> <out_dir>/pca); logs in <out_dir>/error_file
# env (optional): SGE_PROJECT (qsub -P, optional), ALPHA (0.05), LFC (1), MIN_COUNT (10), MIN_SAMPLES (3, PCA), NTOP (500, PCA), HOLD_JID
MAT=$(readlink -f "$1"); SS=$(readlink -f "$2"); CMP=$(readlink -f "$3"); OUT=$(readlink -f -m "$4")
[ -s "$MAT" ] && [ -s "$SS" ] && [ -s "$CMP" ] && [ -n "$4" ] || { echo "Usage: $0 <matrix.tsv> <samplesheet.tsv> <comparisons.tsv> <out_dir>"; exit 1; }
HERE=$(dirname "$(readlink -f "$0")")
# check the sheets before submitting: every samplesheet sample is a matrix column, every comparison group exists
python3 -I "$HERE/check_sheets.py" "$MAT" "$SS" "$CMP" || exit 1
mkdir -p "$OUT/error_file"
V="ALPHA=${ALPHA:-0.05},LFC=${LFC:-1},MIN_COUNT=${MIN_COUNT:-10},MIN_SAMPLES=${MIN_SAMPLES:-3},NTOP=${NTOP:-500}"
H=${HOLD_JID:+-hold_jid $HOLD_JID}
j1=$(qsub -terse -N diff_expr $H ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} -pe omp 4 -l h_rt=04:00:00 -v "$V" \
    -e "$OUT/error_file/diff_expr_error.txt" -o "$OUT/error_file/diff_expr_output.txt" \
    "$HERE/deseq2.job.sh" "$HERE/diff_expr.R" "$MAT" "$SS" "$CMP" "$OUT/results") || exit 1
j2=$(qsub -terse -N pca_rnaseq $H ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} -pe omp 4 -l h_rt=02:00:00 -v "$V" \
    -e "$OUT/error_file/pca_error.txt" -o "$OUT/error_file/pca_output.txt" \
    "$HERE/deseq2.job.sh" "$HERE/pca_rnaseq.R" "$MAT" "$SS" "$OUT/pca") || exit 1
echo "diff_expr job $j1; pca job $j2"
