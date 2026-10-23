#!/bin/bash
# Usage: bash scan_logs.sh <log_dir> [file_glob]      (default glob: *error* and *output* files; read-only)
# Summarises SGE logs of a step: how many logs, how many empty, and every distinct line that mentions error/fail/kill/warn/halt/exception,
# with numbers masked (N) and counted, so the same progress or statistics line repeated in many logs appears once with its count.
# A line being listed is not yet an error: read it (fastp prints "reads failed due to too short", R prints package warnings).
D=${1:?log dir}; G=${2:-}
if [ -n "$G" ]; then files=$(ls $D/$G 2>/dev/null); else files=$(ls $D/*error* $D/*output* 2>/dev/null); fi
[ -z "$files" ] && { echo "no log files in $D"; exit 1; }
n=$(echo "$files" | wc -l); e=0; for f in $files; do [ -s "$f" ] || e=$((e+1)); done
echo "log files: $n   empty: $e   non-empty: $((n-e))"
echo "== distinct lines with error|fail|kill|warn|halt|exception|abort|cannot|denied|no such (numbers masked):"
cat $files | grep -iE "error|fail|kill|warn|halt|exception|abort|cannot|denied|no such" | sed -E 's/[0-9]+(\.[0-9]+)?/N/g' | sort | uniq -c | sort -rn | head -40
echo "== last line of each non-empty *output* log (distinct, numbers masked):"
for f in $files; do case $f in *output*) [ -s "$f" ] && tail -1 "$f";; esac; done | sed -E 's/[0-9]+(\.[0-9]+)?/N/g; s#/[^ ]*/#<path>/#g' | sort | uniq -c | sort -rn | head -10
