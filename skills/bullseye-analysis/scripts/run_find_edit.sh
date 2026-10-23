#!/bin/bash -l
# Usage: bash run_find_edit.sh <matrix_dir> <out_dir> <pairs.tsv>
#   pairs.tsv  columns control_sample, edit_sample, output_name (make_pairs.py); matrices = <matrix_dir>/<sample>.matrix.gz
# env (all REQUIRED, the user decides): MIN_EDIT, MAX_EDIT (percent), FOLD, MIN_EDIT_SITES, REFFLAT
#   REFFLAT = refFlat made by gtf2genepred.pl WITH patches/gtf2genepred_refflat_fix.patch (the unpatched one has wrong gene names/CDS),
#   same Ensembl release as the STAR index
# env (optional): HOLD_JID (comma separated job ids the jobs wait for), MIN_COV (default 10, both control and edited), STRANDED=1, SLOTS (default 8)
# Pairs whose <out_dir>/<output_name>.bed exists are skipped.
MD=$(readlink -f "$1"); OUT=$(readlink -f -m "$2"); PAIRS=$(readlink -f "$3")
[ -d "$MD" ] && [ -f "$PAIRS" ] || { echo "Usage: $0 <matrix_dir> <out_dir> <pairs.tsv>"; exit 1; }
for v in MIN_EDIT MAX_EDIT FOLD MIN_EDIT_SITES; do [ -n "${!v}" ] || { echo "Set $v explicitly; ask the user"; exit 1; }; done
HERE=$(dirname "$(readlink -f "$0")")
[ -s "$REFFLAT" ] || { echo "Set REFFLAT= the PATCHED refFlat (see docs/bullseye_known_issues.md); got: '$REFFLAT'"; exit 1; }
mkdir -p "$OUT/error_file"
n=0; skipped=0; missing=0
while IFS=$'\t' read -r c e o; do
    [ -z "$c" ] && continue
    if [ ! -s "$MD/$c.matrix.gz.tbi" ] || [ ! -s "$MD/$e.matrix.gz.tbi" ]; then echo "MISSING matrix for $c or $e"; missing=$((missing+1)); continue; fi
    if [ -s "$OUT/$o.bed" ]; then echo "SKIP $o (bed exists)"; skipped=$((skipped+1)); continue; fi
    id=$(qsub -terse -N "fe_${o:0:20}" ${HOLD_JID:+-hold_jid $HOLD_JID} -l h_rt=24:00:00 ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} -pe omp ${SLOTS:-8} \
        -e "$OUT/error_file/${o}_error_findedits.txt" -o "$OUT/error_file/${o}_output_findedits.txt" \
        -v "MIN_EDIT=$MIN_EDIT,MAX_EDIT=$MAX_EDIT,FOLD=$FOLD,MIN_EDIT_SITES=$MIN_EDIT_SITES,MIN_COV=${MIN_COV:-10},STRANDED=${STRANDED:-0}" \
        "$HERE/find_edit.job.sh" "$MD/$c.matrix.gz" "$MD/$e.matrix.gz" "$OUT/$o.bed" "$REFFLAT" "$HERE")
    echo "$o (control $c) -> job $id"; n=$((n+1))
done < "$PAIRS"
echo "submitted $n, skipped $skipped, missing $missing"
