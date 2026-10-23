#!/bin/bash -l
# Usage: bash summarize.sh <min_rep> <out_bed> <bed1> <bed2> ...     (run on a compute node / qrsh, it is light)
# Merge the per-sample .bed files: keep sites found in >= min_rep samples (summarize_sites.pl --MinRep). Extra summarize_sites.pl options via env EXTRA (e.g. "--edit 0.05 --coverage 20").
HERE=$(dirname "$(readlink -f "$0")")
MINREP=$1; OUTB=$2; shift 2
[ $# -ge 2 ] || { echo "need at least 2 bed files"; exit 1; }
if [ -n "$ENV_SETUP" ]; then . "$ENV_SETUP" bullseye || exit 1; else module load miniconda; conda activate Bullseye || exit 1; fi   # software: ENV_SETUP (conda envs, see envs/) or the cluster modules
perl "$HERE/summarize_sites.pl" --MinRep "$MINREP" $EXTRA "$@" > "$OUTB" || { echo "summarize failed"; exit 1; }
echo "wrote $OUTB: $(grep -vc '^#' "$OUTB") sites"
