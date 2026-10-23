#!/bin/bash -l
# qsub worker: $1 sample, $2 in_dir, $3 out_dir, $4 r1 suffix, $5 r2 suffix, $6 index_dir
# env: RM_DUP (0 = no Picard step | mark = Picard MarkDuplicates, flag only | remove = Picard, duplicates deleted), STAR_BIN, TRANSCRIPTOME (1 = add --quantMode TranscriptomeSAM)
s=$1; IN=$2; OUT=$3; R1S=$4; R2S=$5; IDX=$6
STAR=${STAR_BIN:-}             # a STAR binary (e.g. 2.7.11a); otherwise STAR from ENV_SETUP (conda) or module star/2.7.10b
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" tools || exit 1; [ -x "$STAR" ] || STAR=STAR
elif [ ! -x "$STAR" ]; then module load star/2.7.10b; STAR=STAR; fi
[ -n "$ENV_SETUP" ] || { module load htslib/1.19.1; module load samtools/1.19.2; }   # samtools/1.19.2 refuses to load without htslib
# The cluster kills jobs that use more CPU than their slots; STAR spawns extra sorting threads, so keep 2 slots free.
T=$(( ${NSLOTS:-16} - 2 )); [ "$T" -lt 1 ] && T=1
BAM=$OUT/${s}_Aligned.sortedByCoord.out.bam
EXTRA=""; [ "$TRANSCRIPTOME" = 1 ] && EXTRA="--quantMode TranscriptomeSAM"
if [ -s "$BAM" ] && [ -s "$OUT/${s}_Log.final.out" ]; then
    echo "STAR for $s already finished (Log.final.out is written last), skipping alignment; only the missing steps run"
    samtools quickcheck -v "$BAM" || { echo "$BAM is truncated"; exit 1; }
    [ -s "$BAM.bai" ] || samtools index -@ "$T" "$BAM" || exit 1
else
    rm -rf "$OUT/${s}__STARtmp"      # leftover of a killed run of this same sample
    $STAR --genomeDir "$IDX" --readFilesCommand zcat \
        --readFilesIn "$IN/$s$R1S" "$IN/$s$R2S" \
        --outFileNamePrefix "$OUT/${s}_" \
        --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
        --runThreadN "$T" $EXTRA || { echo "STAR failed for $s"; exit 1; }
    samtools index -@ "$T" "$BAM" || exit 1
fi
case "$RM_DUP" in
  mark|remove)
    if [ -n "$ENV_SETUP" ]; then PIC=(picard -Xmx16g -Djava.io.tmpdir="${TMPDIR:-/tmp}")          # bioconda wrapper
    else module load picard/2.25.2; PIC=(java -Xmx16g -Djava.io.tmpdir="${TMPDIR:-/tmp}" -jar "$PICARD_JAR"); fi
    if [ "$RM_DUP" = mark ]; then SUF=markdup; RMD=false; else SUF=rm_dup; RMD=true; fi
    # mark: duplicates only get flag 0x400 (every read stays; Bullseye parseBAM --removeDuplicates then ignores them, or not);
    # remove: duplicates are deleted from the BAM.
    "${PIC[@]}" MarkDuplicates \
        I="$BAM" O="$OUT/${s}_Aligned.sortedByCoord.out.$SUF.bam" M="$OUT/${s}_Aligned.sortedByCoord.out.$SUF.csv" \
        REMOVE_DUPLICATES=$RMD REMOVE_SEQUENCING_DUPLICATES=$RMD || { echo "MarkDuplicates failed for $s"; exit 1; }
    samtools index -@ "$T" "$OUT/${s}_Aligned.sortedByCoord.out.$SUF.bam" || exit 1;;
esac
rm -rf "$OUT/${s}__STARtmp"      # STAR leaves empty BAMsort dirs behind
echo "done $s"
