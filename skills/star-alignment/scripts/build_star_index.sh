#!/bin/bash -l
# Build a STAR genome index (qsub worker or run directly).
# Usage: build_star_index.sh <fasta> <gtf> <index_dir> <read_length>
#   sjdbOverhang = read_length - 1 (use the MAXIMUM read length; trimmed reads of 50-150 bp -> read_length 150).
# Skips if <index_dir> already holds a finished index with the same sjdbOverhang. Builds into <index_dir>.building and renames at
# the end, so a killed job never leaves a half-built <index_dir>.
# env: STAR_BIN (optional STAR binary; default module star/2.7.10b)
FA=$1; GTF=$2; IDX=$3; RL=$4
[ -n "$RL" ] || { echo "Usage: $0 <fasta> <gtf> <index_dir> <read_length>"; exit 1; }
OH=$((RL - 1))
if [ -s "$IDX/SA" ] && [ -s "$IDX/SAindex" ] && [ -s "$IDX/genomeParameters.txt" ] && grep -q -P "^sjdbOverhang\t$OH\$" "$IDX/genomeParameters.txt"; then
    echo "index $IDX already exists (sjdbOverhang $OH), skipping"; exit 0
fi
[ -e "$IDX" ] && { echo "$IDX exists but is not a finished index with sjdbOverhang $OH; not touching it"; exit 1; }
STAR=${STAR_BIN:-}             # a STAR binary (e.g. 2.7.11a); otherwise STAR from ENV_SETUP (conda) or module star/2.7.10b
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" tools || exit 1; [ -x "$STAR" ] || STAR=STAR
elif [ ! -x "$STAR" ]; then module load star/2.7.10b; STAR=STAR; fi
T=$(( ${NSLOTS:-16} - 2 )); [ "$T" -lt 1 ] && T=1
rm -rf "$IDX.building"; mkdir -p "$IDX.building"
echo "STAR $($STAR --version), sjdbOverhang $OH, threads $T"
$STAR --runMode genomeGenerate --runThreadN "$T" --genomeDir "$IDX.building" \
    --genomeFastaFiles "$FA" --sjdbGTFfile "$GTF" --sjdbOverhang "$OH" \
    --limitGenomeGenerateRAM 60000000000 --outFileNamePrefix "$IDX.building/" || { echo "STAR genomeGenerate failed"; exit 1; }
mv "$IDX.building" "$IDX" && echo "index ready: $IDX"
