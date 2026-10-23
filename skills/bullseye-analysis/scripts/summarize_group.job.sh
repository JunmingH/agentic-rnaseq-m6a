#!/bin/bash -l
# qsub worker: $1 group, $2 bed dir, $3 out dir, $4 scripts dir
g=$1; B=$2; O=$3; S=$4
for mr in 2 3; do
  bash $S/summarize.sh $mr $O/${g}_minrep${mr}.bed $B/${g}_1_WT.bed $B/${g}_2_WT.bed $B/${g}_3_WT.bed || exit 1
  EXTRA="--mut 5 --coverage 20" bash $S/summarize.sh $mr $O/${g}_minrep${mr}_mut5_cov20.bed $B/${g}_1_WT.bed $B/${g}_2_WT.bed $B/${g}_3_WT.bed || exit 1
done
