#!/usr/bin/env python3
"""Compare Find_edit_site beds made with two annotations (e.g. original vs fixed refFlat).
Usage: compare_beds.py <old_bed_dir> <new_bed_dir> [name ...]      (names = bed base names; default: every *.bed in new_bed_dir)
Column 4 of a Bullseye bed is gene|region|C2T|mut=n|score. Per bed:
  same set of sites (chr, end, strand)? for shared sites: changes of the gene name, of the region label (5UTR/CDS/3UTR/ncRNA/...),
  and of anything else (positions, ratios, coverages, mut, score) - the latter must be 0.
Gene-name changes to a non-ENSG name are counted separately (overlapping genes with identical transcripts can swap names between runs).
At the end: the region transitions (old -> new) summed over all beds. Exit 1 if sites or 'other' columns differ."""
import sys, os, glob, collections
old, new = sys.argv[1:3]
names = sys.argv[3:] or sorted(os.path.basename(f)[:-4] for f in glob.glob(os.path.join(new, '*.bed')))
def load(p):
    d = {}
    for l in open(p):
        if l[0] == '#': continue
        f = l.rstrip('\n').split('\t'); d[(f[0], f[2], f[5])] = f
    return d
bad = 0; trans = collections.Counter()
print('bed\tsites_old\tsites_new\tonly_old\tonly_new\tother_columns_differ\tgene_name_changed\tchanged_to_non_ENSG\tregion_changed')
for n in names:
    o, w = load(os.path.join(old, n + '.bed')), load(os.path.join(new, n + '.bed'))
    so, sn = set(o), set(w); oth = chg = nonid = reg = 0
    for k in so & sn:
        a, b = o[k], w[k]
        ga, ra, rest_a = (a[3].split('|', 2) + ['', ''])[:3]; gb, rb, rest_b = (b[3].split('|', 2) + ['', ''])[:3]
        if a[:3] + a[4:] != b[:3] + b[4:] or rest_a != rest_b: oth += 1
        if ga != gb:
            chg += 1
            if not gb.startswith('ENSG'): nonid += 1
        if ra != rb: reg += 1; trans[(ra, rb)] += 1
    if len(so - sn) or len(sn - so) or oth: bad += 1
    print(f'{n}\t{len(o)}\t{len(w)}\t{len(so - sn)}\t{len(sn - so)}\t{oth}\t{chg}\t{nonid}\t{reg}')
if trans: print('region transitions (old -> new), all beds:', ', '.join(f'{a}->{b} {c}' for (a, b), c in trans.most_common()))
print('PASS: same sites, only gene names / region labels differ' if not bad else f'FAIL: {bad} bed(s) differ beyond gene names and region labels'); sys.exit(1 if bad else 0)
