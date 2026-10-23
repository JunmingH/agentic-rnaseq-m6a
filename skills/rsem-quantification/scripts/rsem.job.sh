#!/bin/bash -l
# qsub worker: $1 transcriptome bam, $2 output prefix, $3 RSEM reference prefix
# env: STRANDEDNESS (none|forward|reverse), GENOME_BAM (1 = --output-genome-bam, else --no-bam-output)
# Use the same RSEM version as the one that built the reference (we used rsem 1.3.3).
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" tools || exit 1; else module load R/4.2.1 rsem/1.3.3 bowtie2/2.3.4.1 python2/2.7.16; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
# The cluster kills jobs above their slots: keep 2 free.
T=$(( ${NSLOTS:-12} - 2 )); [ "$T" -lt 1 ] && T=1
if [ "$GENOME_BAM" = 1 ]; then BOPT="--output-genome-bam"; else BOPT="--no-bam-output"; fi
rm -rf "$2.genes.results" "$2.isoforms.results" "$2.stat" "$2.temp"     # leftovers of a killed run of this same sample
rsem-calculate-expression -p "$T" --paired-end --bam --append-names --strandedness "$STRANDEDNESS" $BOPT \
    "$1" "$3" "$2" || { echo "RSEM failed for $1"; exit 1; }
[ -s "$2.genes.results" ] && [ -s "$2.isoforms.results" ] || { echo "no results for $1"; exit 1; }
echo "done $1"
