#!/bin/bash -l
# Download the Ensembl variation VCFs of one release (per chromosome, includes dbSNP; chromosome names 1..22,X,Y,MT = Ensembl, the same as an
# Ensembl STAR / refFlat reference) and convert them to a bed of C>T / G>A SNVs for removing SNP positions from the DART-seq candidates.
# Run on a node with internet (a login node or a qsub job, if compute nodes have internet), e.g.
#   OUT=<ref>/dbsnp_ensembl111 qsub -N dbsnp_dl -P <project> -l h_rt=24:00:00 -V -e <log> -o <log> download_dbsnp_ensembl.sh
# env: OUT (REQUIRED, output dir), RELEASE (111; must match the genome/annotation release), SPECIES (homo_sapiens), CHROMS (default "1..22 X Y MT")
# Re-running is safe: finished downloads (verified by bcftools) and finished beds are skipped. Human release 111: ~8.7 GB bed, takes hours.
# Output:
#   vcf/<species>-chr*.vcf.gz(.csi)        raw downloads, kept
#   bed/chr*.snv.bed                       chrom, start0, end, ref, alt   (SNVs only: single-base REF, single-base ALT, multi-allelic split)
#   dbsnp_ensembl<REL>.CT_GA.snv.bed       only REF>ALT = C>T or G>A (the DART-seq C-to-U and its minus-strand mirror); the DBSNP file of run_diff_m6a.sh
#   n_snv_per_chr.tsv                      counts for the report
OUT=${OUT:?set OUT= output dir}; REL=${RELEASE:-111}; SP=${SPECIES:-homo_sapiens}
BASE=https://ftp.ensembl.org/pub/release-$REL/variation/vcf/$SP
CHRS=${CHROMS:-"$(seq 1 22) X Y MT"}
mkdir -p "$OUT/vcf" "$OUT/bed" "$OUT/error_file"
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" tools || exit 1; else module load htslib/1.18 bcftools/1.18; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
# (htslib before bcftools; before set -u, because module scripts use unset variables)
set -uo pipefail
command -v bcftools >/dev/null || { echo "bcftools not found"; exit 1; }
command -v bgzip >/dev/null || { echo "bgzip not found"; exit 1; }
echo "bcftools: $(bcftools --version | head -1); start $(date)"

fail=0
for c in $CHRS; do
  f=$SP-chr$c.vcf.gz
  # 1) download (resume with -C -), then verify; a bad file is deleted and retried up to 3 times
  ok=0
  for try in 1 2 3; do
    if [ -s "$OUT/vcf/$f" ] && [ -s "$OUT/vcf/$f.csi" ] && bcftools index -n "$OUT/vcf/$f" >/dev/null 2>&1 \
       && gzip -t "$OUT/vcf/$f" 2>/dev/null; then ok=1; break; fi
    echo "[chr$c] download try $try"
    curl -sSL --retry 3 -C - -o "$OUT/vcf/$f" "$BASE/$f" && curl -sSL --retry 3 -o "$OUT/vcf/$f.csi" "$BASE/$f.csi" \
      || { rm -f "$OUT/vcf/$f" "$OUT/vcf/$f.csi"; }
  done
  [ $ok = 1 ] || { echo "[chr$c] FAILED to download/verify"; fail=1; continue; }

  # 2) SNVs -> bed (-m- splits multi-allelic sites; -v snps keeps SNVs; the awk guard keeps single-base REF/ALT only)
  b=$OUT/bed/chr$c.snv.bed
  if [ -s "$b" ] && [ "$b" -nt "$OUT/vcf/$f" ]; then echo "[chr$c] bed exists, skip"; continue; fi
  bcftools view -v snps -G "$OUT/vcf/$f" 2>/dev/null | bcftools norm -m- -Ov 2>/dev/null \
    | awk 'BEGIN{OFS="\t"} !/^#/ && length($4)==1 && length($5)==1 && $4!="." && $5!="." {print $1,$2-1,$2,$4,$5}' > "$b.tmp" \
    && mv "$b.tmp" "$b" || { echo "[chr$c] conversion FAILED"; rm -f "$b.tmp"; fail=1; continue; }
  echo "[chr$c] $(wc -l < "$b") SNVs"
done
[ $fail = 0 ] || { echo "some chromosomes failed; merged files NOT written. Re-run to retry."; exit 1; }

# 3) merge the C>T / G>A SNVs only (the full SNV set is ~1e9 lines; per-chromosome beds are kept). Each chr bed is already
#    position-sorted (VCF order), so concatenating in chromosome order needs no global sort; uniq drops repeated records.
cd "$OUT"
for c in $CHRS; do awk '($4=="C"&&$5=="T")||($4=="G"&&$5=="A")' bed/chr$c.snv.bed | uniq; done > dbsnp_ensembl$REL.CT_GA.snv.bed
( echo -e "chrom\tn_snv"; for c in $CHRS; do echo -e "$c\t$(wc -l < bed/chr$c.snv.bed)"; done ) > n_snv_per_chr.tsv
md5sum dbsnp_ensembl$REL.CT_GA.snv.bed > MANIFEST.md5
echo "C>T/G>A SNVs: $(wc -l < dbsnp_ensembl$REL.CT_GA.snv.bed); all SNVs: $(awk 'NR>1{s+=$2}END{print s}' n_snv_per_chr.tsv); done $(date)"
