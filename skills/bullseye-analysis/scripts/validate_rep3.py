#!/usr/bin/env python3
"""Cross-validation of the direction (WT = edited, Mut = control) with an independent replicate.
Usage: validate_rep3.py <group> <bed_dir> <matrix_dir> <out_tsv> [min_mut=5] [min_cov=20] [edited_tag=WT]
edited_tag = WT (main direction: beds <g>_n_WT.bed of the edited-vs-control run) or Mut (reverse direction, negative control: beds <g>_n_Mut.bed of the run with the roles swapped).
Sites = C->T sites present in the Find_edit_site bed of pair 1 AND pair 2 (WT_1, WT_2 as edited, Mut_1, Mut_2 as control) with
mut reads >= min_mut and WT coverage >= min_cov (selected without replicate 3). In replicate 3 the editing ratio is read from the
matrices of WT_3 and Mut_3 at the same positions: T/(C+T) on '+' genes, A/(A+G) on '-' genes (unstranded matrix; columns chr pos A T C G N total).
If the sites are WT-specific edits, WT_3 > Mut_3 at these sites. Needs tabix on PATH."""
import sys, subprocess, tempfile, os, statistics as st
g, bd, md, out = sys.argv[1:5]
mn = int(sys.argv[5]) if len(sys.argv) > 5 else 5; mc = int(sys.argv[6]) if len(sys.argv) > 6 else 20; et = sys.argv[7] if len(sys.argv) > 7 else 'WT'
ct = 'Mut' if et == 'WT' else 'WT'
def load(f):
    d = {}
    for l in open(f):
        if l[0] == '#': continue
        p = l.rstrip('\n').split('\t'); n = int(p[3].split('|')[3].split('=')[1])
        if n >= mn and int(p[9]) >= mc: d[(p[0], int(p[2]), p[5])] = float(p[8])
    return d
a = load(f'{bd}/{g}_1_{et}.bed'); b = load(f'{bd}/{g}_2_{et}.bed')
sites = sorted(set(a) & set(b))
tmp = tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False)
for c, e, s in sites: tmp.write(f'{c}\t{e-1}\t{e}\n')   # tabix -R wants BED (0-based start); the matrix pos is 1-based = bed 'end'
tmp.close()
def fetch(sample):
    r = subprocess.run(['tabix', '-R', tmp.name, f'{md}/{sample}.matrix.gz'], capture_output=True, text=True, check=True).stdout
    d = {}
    for l in r.split('\n'):
        if not l: continue
        p = l.split('\t')
        key = (p[0], int(p[1]), p[8] if len(p) > 8 else None)      # stranded matrix: one row per strand label (+ / - / NS), use the row of the site's strand
        d[key] = (int(p[2]), int(p[3]), int(p[4]), int(p[5]), int(p[7]))   # A T C G tot
    return d
W3 = fetch(f'{g}_3_{et}'); M3 = fetch(f'{g}_3_{ct}'); os.unlink(tmp.name)   # W3 = edited sample rep 3, M3 = control sample rep 3
rows = []
def ratio(x, strand):
    A, T, C, G, t = x
    den = (C + T) if strand == '+' else (A + G)
    return (T if strand == '+' else A) / den if den else None
for c, e, s in sites:
    w = W3.get((c, e, s)) or W3.get((c, e, None)); m = M3.get((c, e, s)) or M3.get((c, e, None))
    if not w or not m: continue
    rw, rm = ratio(w, s), ratio(m, s)
    if rw is None or rm is None or w[4] < 10 or m[4] < 10: continue
    rows.append((c, e, s, rw, rm, w[4], m[4]))
with open(out, 'w') as o:
    o.write('chr\tpos\tstrand\tedited3_ratio\tcontrol3_ratio\tedited3_cov\tcontrol3_cov\n')
    for r in rows: o.write('\t'.join(str(x) for x in r) + '\n')
dw = [r[3] - r[4] for r in rows]
print(f'{g} [edited={et}]: sites in {et}_1 and {et}_2 beds {len(sites)}; usable in rep 3 {len(rows)}; median {et}3 {st.median(r[3] for r in rows):.3f} {ct}3 {st.median(r[4] for r in rows):.3f}; {et}3>{ct}3 in {sum(1 for x in dw if x>0)/len(dw)*100:.0f}%')
