#!/bin/bash -l
# Usage: bash run_fastqc.sh <fastq_dir> <list_file>
#   <fastq_dir>  directory with the fastq files (output -> <fastq_dir>/fastqc)
#   <list_file>  one fastq file name per line, relative to <fastq_dir>
# Submits one qsub job per fastq file, then one MultiQC job that is held until
# ALL fastqc jobs have finished (-hold_jid), so the summary is made automatically.

FQDIR=$(readlink -f "$1")
LIST=$(readlink -f "$2")
[ -d "$FQDIR" ] && [ -f "$LIST" ] || { echo "Usage: $0 <fastq_dir> <list_file>"; exit 1; }
HERE=$(dirname "$(readlink -f "$0")")

mkdir -p "$FQDIR/fastqc/error_file"
echo "Start qc: $(wc -l < "$LIST") files in $FQDIR"
ids=""
while read -r f; do
    [ -z "$f" ] && continue
    id=$(qsub -terse -N "fqc_${f:0:20}" -l h_rt=24:00:00 ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} -pe omp 4 \
        -e "$FQDIR/fastqc/error_file/${f}_error_fastqc.txt" \
        -o "$FQDIR/fastqc/error_file/${f}_output_fastqc.txt" \
        "$HERE/fastqc.job.sh" "$FQDIR/$f" "$FQDIR")
    echo "$f -> job $id"
    ids="${ids:+$ids,}$id"
done < "$LIST"

mid=$(qsub -terse -N multiqc_all -hold_jid "$ids" -l h_rt=02:00:00 ${SGE_PROJECT:+-P $SGE_PROJECT} ${ENV_SETUP:+-v ENV_SETUP=$ENV_SETUP} \
    -e "$FQDIR/fastqc/error_file/multiqc_error.txt" \
    -o "$FQDIR/fastqc/error_file/multiqc_output.txt" \
    "$HERE/multiqc.job.sh" "$FQDIR")
echo "MultiQC job $mid (starts after all fastqc jobs finish) -> $FQDIR/fastqc/multiqc_report.html"
