#!/bin/bash
# Minimal reproduction of the Find_edit_site.pl stranded-merge bug (docs/bullseye_known_issues.md, issue 2).
# Usage: bash run_repro.sh /path/to/Find_edit_site.pl     (needs perl with the Bullseye modules, bgzip, tabix; ~10 s)
# Toy data: one chromosome "1", an antisense gene pair GPLUS (+) / GMINUS (-) over 1000-2000, four true edited sites:
#   1100 '-' (edited has only a '-' row; control has '+' then '-')   1200 '+'   1300 '+' (edited has '+' then '-'; control '+')   1400 '-'
# Run A uses the matrices as parseBAM --stranded writes them; run B removes only the two extra opposite-strand rows, which carry no edits.
# Expected with a correct merge: 4 sites in both runs.
set -e
FE=${1:?usage: bash run_repro.sh /path/to/Find_edit_site.pl}
W=$(mktemp -d); cp "$(dirname "$0")"/{toy.refFlat,edited.matrix,control.matrix} "$W/"; cd "$W"
awk '!($2==1100 && $9=="+")' control.matrix > control_B.matrix; awk '!($2==1300 && $9=="-")' edited.matrix > edited_B.matrix
for m in edited control edited_B control_B; do bgzip -f $m.matrix; tabix -f -s1 -b2 -e2 $m.matrix.gz; done
run() { perl "$FE" --annotationFile toy.refFlat --controlMatrix $1 --EditedMatrix $2 --outfile $3 --minEdit 10 --maxEdit 95 \
          --editFoldThreshold 1.5 --MinEditSites 2 --ControlMinCoverage 10 --EditedMinCoverage 10 --cpu 1 >/dev/null 2>&1; }
run control.matrix.gz edited.matrix.gz A.bed; run control_B.matrix.gz edited_B.matrix.gz B.bed
echo "A (as written by parseBAM --stranded): $(grep -vc '^#' A.bed) of 4 sites: $(grep -v '^#' A.bed | awk '{printf "%s%s ", $3, $6}')"
echo "B (opposite-strand rows removed):      $(grep -vc '^#' B.bed) of 4 sites: $(grep -v '^#' B.bed | awk '{printf "%s%s ", $3, $6}')"
rm -rf "$W"
