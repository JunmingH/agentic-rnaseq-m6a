#!/usr/bin/env python3
"""Compare two csv/tsv tables (used for docs/validation.md) row by row on a key: same key set, max relative difference per numeric column, identical text columns.
Usage: compare_tables.py <original> <new> <key_col[,key_col2]> [column renames old:new,...]
Exit 0 and prints 'IDENTICAL'/'EQUAL within 1e-6' or the differences."""
import sys, csv, math
a_f, b_f, key = sys.argv[1:4]; ren = dict(x.split(':') for x in sys.argv[4].split(',')) if len(sys.argv) > 4 and sys.argv[4] else {}
def load(f, rename):
    d = '\t' if f.endswith(('.tsv', '.txt')) else ','
    with open(f) as fh:
        r = csv.DictReader(fh, delimiter=d); rows = {}
        for x in r:
            x = {rename.get(k, k): v for k, v in x.items()}
            rows[tuple(x[k] for k in key.split(','))] = x
        return rows, [rename.get(c, c) for c in r.fieldnames]
A, ca = load(a_f, ren); B, cb = load(b_f, {})
print(f'{a_f}\n  vs {b_f}\n  rows {len(A)} / {len(B)}; keys only in original {len(set(A)-set(B))}, only in new {len(set(B)-set(A))}')
common = [c for c in ca if c in cb]
miss = [c for c in ca if c not in cb]
if miss: print('  columns of the original missing in new:', miss)
worst = 0.0; txt = 0
for c in common:
    md = 0.0; td = 0
    for k in set(A) & set(B):
        x, y = A[k][c], B[k][c]
        if x == y: continue
        try:
            fx, fy = float(x), float(y)
            if math.isnan(fx) and math.isnan(fy): continue
            md = max(md, abs(fx - fy) / max(abs(fx), abs(fy), 1e-300))
        except ValueError:
            if x in ('', 'NA') and y in ('', 'NA'): continue
            td += 1
    if md > 0 or td: print(f'  {c}: max rel diff {md:.3g}; text differences {td}')
    worst = max(worst, md); txt += td
ok = len(set(A) ^ set(B)) == 0 and txt == 0 and worst < 1e-6
print('  RESULT:', ('IDENTICAL' if worst == 0 else f'EQUAL within 1e-6 (max rel diff {worst:.2g})') if ok else 'DIFFERENT')
sys.exit(0 if ok else 1)
