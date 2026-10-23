#!/usr/bin/env python3
"""Sequence context of the C->T sites (DART edits the C right after an m6A: the m6A is the A in the DRACH motif, so the edited C sits in an 'A[C]H' context: base -1 = A).
Usage: motif_context.py <bed_dir_main> <validate_dir> <fasta> <out_tsv>
Sets per group: WT-chosen (validate_rep3/<g>.tsv, sites from WT_1+WT_2), Mut-chosen (<g>_rev.tsv), and 'single' = sites only in the WT_1 bed but in neither WT_2 nor WT_3 beds (mostly noise, baseline).
Context in gene orientation (reverse complement for '-' genes). Reports the share of A at -1, of the DRACH-like R,A,C,H pattern, and the sanity check that the centre is C."""
import sys, os, glob, collections
bd, vd, fa, out = sys.argv[1:5]
G = sorted(os.path.basename(f)[:-4] for f in glob.glob(os.path.join(vd, '*.tsv')) if not f.endswith('_rev.tsv'))   # groups = <g>.tsv of validate_rep3.py
fai = {l.split('\t')[0]: tuple(int(x) for x in l.split('\t')[1:5]) for l in open(fa + '.fai')}
fh = open(fa, 'rb')
def seq(chrom, start, end):          # 1-based inclusive
    L, off, lb, lw = fai[chrom]; start = max(start, 1)
    s0 = off + (start - 1) // lb * lw + (start - 1) % lb; n = end - start + 1
    fh.seek(s0); raw = fh.read(n + n // lb + 2).replace(b'\n', b'')[:n]
    return raw.decode().upper()
comp = str.maketrans('ACGTN', 'TGCAN')
def ctx(chrom, pos, strand):
    if chrom not in fai: return None
    s = seq(chrom, pos - 2, pos + 2)
    if len(s) != 5: return None
    return s if strand == '+' else s.translate(comp)[::-1]
def keys(f):
    d = set()
    for l in open(f):
        if l[0] == '#': continue
        p = l.split('\t', 6); d.add((p[0], int(p[2]), p[5]))
    return d
def stats(sites):
    n = c = a1 = dr = 0
    comp5 = collections.Counter()
    for chrom, pos, strand in sites:
        x = ctx(chrom, pos, strand)
        if not x: continue
        n += 1; c += x[2] == 'C'
        if x[2] != 'C': continue
        a1 += x[1] == 'A'; dr += (x[0] in 'AG' and x[1] == 'A' and x[3] in 'ACT')
        comp5[x[1]] += 1
    return n, c / max(n, 1) * 100, a1 / max(c, 1) * 100 if n else 0, dr / max(c, 1) * 100 if n else 0
rows = []
for g in G:
    sets = {}
    for tag, f in (('WT-chosen', f'{vd}/{g}.tsv'), ('Mut-chosen', f'{vd}/{g}_rev.tsv')):
        r = [l.rstrip('\n').split('\t') for l in open(f)][1:]
        sets[tag] = [(x[0], int(x[1]), x[2]) for x in r]
        e = [(x[0], int(x[1]), x[2]) for x in r if float(x[3]) - float(x[4]) >= 0.1]
        sets[tag + ' (strong: ratio diff >= 0.1)'] = e
    w1 = keys(f'{bd}/{g}_1_WT.bed'); w2 = keys(f'{bd}/{g}_2_WT.bed'); w3 = keys(f'{bd}/{g}_3_WT.bed')
    sets['single (baseline, noise)'] = sorted(w1 - w2 - w3)[:30000]
    for tag, s in sets.items():
        n, cc, a1, dr = stats(s); rows.append((g, tag, n, cc, a1, dr)); print(g, tag, n, f'centre=C {cc:.0f}%  A at -1: {a1:.1f}%  R-A-C-H: {dr:.1f}%', flush=True)
with open(out, 'w') as o:
    o.write('group\tset\tn\tcentre_is_C_pct\tA_at_minus1_pct\tRACH_pct\n')
    for r in rows: o.write('\t'.join(str(x) if not isinstance(x, float) else f'{x:.2f}' for x in r) + '\n')
