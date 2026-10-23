#!/bin/bash -l
# Strand-label check of a STRANDED parseBAM matrix (columns: chr pos A T C G N total strand). $1 = matrix.gz, $2 = chromosome (default 1).
# Per strand label (+ / -): ref-majority mismatch rates (C>T, G>A, ... share of reads at positions whose majority base is the reference). If the label = RNA strand (dUTP, as Bullseye assumes):
#   '+' rows: C->T high, G->A at background;  '-' rows: G->A high, C->T at background.
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" bullseye || exit 1; else module load miniconda; conda activate Bullseye || exit 1; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
tabix $1 ${2:-1} | awk 'BEGIN{FS="\t"} {t=$8; s=$9; if(t<20||(s!="+"&&s!="-")) next; A=$3;T=$4;C=$5;G=$6;
  if(C>=0.8*t){rc[s]+=t; ct[s]+=T; ca[s]+=A}
  else if(G>=0.8*t){rg[s]+=t; ga[s]+=A; gt[s]+=T}}
  END{printf "strand\tC_tot\tC2T%%\tC2A%%\tG_tot\tG2A%%\tG2T%%\n"; for(k in rc){printf "%s\t%d\t%.3f\t%.3f\t%d\t%.3f\t%.3f\n",k,rc[k],ct[k]/rc[k]*100,ca[k]/rc[k]*100,rg[k],ga[k]/rg[k]*100,gt[k]/rg[k]*100}}'
