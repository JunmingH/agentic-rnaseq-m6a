#!/usr/bin/env python3
"""Check a samplesheet / comparisons pair before submitting (deseq2-expression, dart-diff-m6a).
Usage: check_sheets.py <matrix.tsv | -> <samplesheet.tsv> <comparisons.tsv>
  matrix.tsv : header = id + sample columns (rsem make_matrix.py); '-' skips the matrix check
Prints the groups, the samples of every comparison (control first) and exits 1 on any problem."""
import sys, csv, collections
mat, ss_f, cmp_f = sys.argv[1:4]
def rows(f):
    with open(f) as fh: return list(csv.DictReader(fh, delimiter='\t'))
ss, cmps = rows(ss_f), rows(cmp_f)
err = []
for col in ('sample', 'group'):
    if not ss or col not in ss[0]: err.append(f'samplesheet has no column {col}')
for col in ('name', 'treated', 'control'):
    if not cmps or col not in cmps[0]: err.append(f'comparisons has no column {col}')
if err: sys.exit('\n'.join(err))
dup = [s for s, n in collections.Counter(r['sample'] for r in ss).items() if n > 1]
if dup: err.append('duplicated samples: ' + ', '.join(dup))
if mat != '-':
    cols = open(mat).readline().rstrip('\n').split('\t')[1:]
    miss = [r['sample'] for r in ss if r['sample'] not in cols]
    if miss: err.append('samples not in the matrix: ' + ', '.join(miss))
    extra = [c for c in cols if c not in {r['sample'] for r in ss}]
    if extra: print('note: matrix columns not in the samplesheet (ignored):', ', '.join(extra))
grp = collections.defaultdict(list)
for r in ss: grp[r['group']].append(r)
print('groups:', ', '.join(f'{g} ({len(v)})' for g, v in grp.items()))
for c in cmps:
    sub = (c.get('subset') or '').strip() or 'all'
    for role in ('control', 'treated'):
        g = c[role]
        if g not in grp: err.append(f"{c['name']}: {role} group {g} not in the samplesheet"); continue
        smp = [r['sample'] for r in grp[g] if sub == 'all' or r.get('subset', 'all') == sub]
        if not smp: err.append(f"{c['name']}: no {role} samples of group {g} in subset {sub}")
        print(f"{c['name']:24s} {role:8s} {g:10s} subset={sub:8s} {', '.join(smp)}")
    if c['treated'] == c['control']: err.append(f"{c['name']}: treated == control")
if err: sys.exit('ERROR\n' + '\n'.join(err))
print('sheets OK')
