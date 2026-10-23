#!/bin/bash -l
# qsub worker for steps 01, 03 (python3, needs tabix: htslib/1.18) and 04 (R 4.2.1 + IHW, data.table, ggplot2, ggrepel).
# $1 step script (absolute path), then its arguments. Env thresholds are passed by run_diff_m6a.sh (-v).
S=$1; shift
case "$S" in *.R) E=r;; *) E=tools;; esac      # python steps need tabix, step 04 needs R + IHW
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" $E || exit 1; else module load htslib/1.18 R/4.2.1; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
case "$S" in
  *.py) python3 -I "$S" "$@" ;;
  *.R)  Rscript "$S" "$@" ;;
esac || { echo "$(basename "$S") failed"; exit 1; }
echo "done $(basename "$S")"
