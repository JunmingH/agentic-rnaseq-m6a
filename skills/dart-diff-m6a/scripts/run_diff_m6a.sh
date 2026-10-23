#!/bin/bash -l
# Usage: bash run_diff_m6a.sh <bed_dir> <matrix_dir> <samplesheet.tsv> <comparisons.tsv> <out_dir> [<diff_expr_tables_dir>]
#   bed_dir      Bullseye Find_edit_site beds of the EDITED samples, <sample>.bed (edited vs its control)
#   matrix_dir   parseBAM matrices <sample>.matrix.gz(.tbi) of the same edited samples (use the --stranded ones if the beds are stranded)
#   samplesheet  columns sample, group (edited samples only)        comparisons  columns name, treated, control [, de_name]
#   out_dir      -> sites/ (01-03), diff_m6a/ (04), error_file/ (logs)
# Submits a job chain 01 candidates -> 02 dbSNP -> 03 matrix -> 04 GLM + IHW (hold_jid).
# env REQUIRED: FASTA (genome fasta + .fai, same chromosome names as the beds), DBSNP (C>T/G>A SNV bed from download_dbsnp_ensembl.sh)
# env optional: MIN_REPS (2) MIN_MEAN (0.05) MIN_COV (20) MIN_PER_GROUP (2) ALPHA (0.05) GENE_FC (1.25) EPS (0.01) SEED (1)
#               SGE_PROJECT (qsub -P, optional) HOLD_JID (job ids to wait for, e.g. the Find_edit_site jobs) START (1..4, default 1: restart from a later step)
BED=$(readlink -f "$1"); MAT=$(readlink -f "$2"); SS=$(readlink -f "$3"); CMP=$(readlink -f "$4"); OUT=$(readlink -f -m "$5")
DE=${6:+$(readlink -m "$6")}; DE=${DE:--}      # -m: the DE tables may not exist yet when the chain waits for the DESeq2 job (HOLD_JID)
[ "$DE" = - ] || [ -d "$DE" ] || [ -n "$HOLD_JID" ] || { echo "DE tables dir $DE does not exist (and no HOLD_JID to wait for it)"; exit 1; }
[ "$DE" = - ] && echo "no diff_expr tables given: expression columns will be empty" || echo "DE tables: $DE"
[ -d "$BED" ] && [ -d "$MAT" ] && [ -s "$SS" ] && [ -s "$CMP" ] && [ -n "$5" ] || { echo "Usage: $0 <bed_dir> <matrix_dir> <samplesheet.tsv> <comparisons.tsv> <out_dir> [<diff_expr_tables_dir>]"; exit 1; }
[ -s "$FASTA" ] && [ -s "$FASTA.fai" ] || { echo "Set FASTA= genome fasta with .fai (same chromosome names as the beds)"; exit 1; }
[ -s "$DBSNP" ] || { echo "Set DBSNP= C>T/G>A SNV bed (download_dbsnp_ensembl.sh)"; exit 1; }
HERE=$(dirname "$(readlink -f "$0")")
python3 -I "$HERE/check_sheets.py" - "$SS" "$CMP" || exit 1
START=${START:-1}
if [ "$START" -le 1 ]; then         # inputs of step 01 must exist now unless we wait for the jobs that make them
  miss=0
  while IFS=$'\t' read -r s g _; do
    [ "$s" = sample ] && continue
    [ -n "$HOLD_JID" ] || [ -f "$BED/$s.bed" ] || { echo "MISSING $BED/$s.bed"; miss=1; }
    [ -n "$HOLD_JID" ] || [ -s "$MAT/$s.matrix.gz.tbi" ] || { echo "MISSING $MAT/$s.matrix.gz.tbi"; miss=1; }
  done < "$SS"
  [ $miss = 0 ] || exit 1
fi
mkdir -p "$OUT/sites" "$OUT/error_file"; O=$OUT/sites
V="MIN_REPS=${MIN_REPS:-2},MIN_MEAN=${MIN_MEAN:-0.05},MIN_COV=${MIN_COV:-20},MIN_PER_GROUP=${MIN_PER_GROUP:-2},ALPHA=${ALPHA:-0.05},GENE_FC=${GENE_FC:-1.25},EPS=${EPS:-0.01},SEED=${SEED:-1}"
J=$HOLD_JID
sub() {   # name step h_rt slots args...   (prints the job id; holds on $J if set)
  local n=$1 st=$2 rt=$3 sl=$4; shift 4
  qsub -terse -N "$n" ${J:+-hold_jid $J} ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} -pe omp "$sl" -l h_rt=$rt -v "$V" \
     -e "$OUT/error_file/${st}_error.txt" -o "$OUT/error_file/${st}_output.txt" "$@"
}
[ "$START" -le 1 ] && { J=$(sub m6a_cand 01 04:00:00 1 "$HERE/m6a_step.job.sh" "$HERE/01_candidates.py" "$BED" "$FASTA" "$SS" "$O") || exit 1; echo "01 candidates: $J"; }
[ "$START" -le 2 ] && { J=$(sub m6a_snp 02 06:00:00 1 "$HERE/02_dbsnp.job.sh" "$O/snp_candidates.tsv" "$DBSNP" "$O/snp_hits.tsv") || exit 1; echo "02 dbSNP: $J"; }
[ "$START" -le 3 ] && { J=$(sub m6a_mat 03 04:00:00 1 "$HERE/m6a_step.job.sh" "$HERE/03_matrix.py" "$O/group_sites.tsv" "$O/per_sample_calls.tsv" "$O/snp_hits.tsv" "$MAT" "$CMP" "$O") || exit 1; echo "03 matrix: $J"; }
J=$(sub diff_m6a 04 06:00:00 4 "$HERE/m6a_step.job.sh" "$HERE/04_diff_m6a.R" "$O" "$CMP" "$OUT/diff_m6a" "$DE") || exit 1; echo "04 diff_m6a: $J"
