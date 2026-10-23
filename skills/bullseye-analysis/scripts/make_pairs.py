#!/usr/bin/env python3
"""Make the control/edit pair list for Find_edit_site.
Usage: make_pairs.py <samples.txt> <control_tag> <edit_tag> <mode> > pairs.tsv
  sample names are <group>_<n>_<tag>  (e.g. ctrl_1_WT, ctrl_1_Mut); tags = the last field.
  mode  same_n : edit <group>_<n> is paired with the control <group>_<n>
        all    : each edit sample is paired with every control of the same group
Output columns: control_sample, edit_sample, output_name (= edit sample).  The user must say which tag is the
control (background) and which is the edited (DART enzyme) sample.
"""
import sys, collections
lst, ctag, etag, mode = sys.argv[1:5]
if mode not in ("same_n", "all"): sys.exit("mode must be same_n or all")
S = [l.strip() for l in open(lst) if l.strip()]
info = {}
for s in S:
    g, n, t = s.rsplit("_", 2)
    info[s] = (g, n, t)
ctrl = collections.defaultdict(list)
for s, (g, n, t) in info.items():
    if t == ctag: ctrl[g].append((n, s))
npairs = 0
for s, (g, n, t) in info.items():
    if t != etag: continue
    for cn, cs in ctrl.get(g, []):
        if mode == "same_n" and cn != n: continue
        print(f"{cs}\t{s}\t{s}"); npairs += 1
    if not any(cn == n or mode == "all" for cn, _ in ctrl.get(g, [])):
        print(f"WARNING: no control for {s}", file=sys.stderr)
print(f"{npairs} pairs", file=sys.stderr)
