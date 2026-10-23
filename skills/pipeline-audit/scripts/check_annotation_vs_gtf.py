#!/usr/bin/env python3
"""Check a refFlat/genePred-style annotation against the GTF it was made from (read-only).
Usage: python3 -I check_annotation_vs_gtf.py <gtf> <refFlat> [--name-col 0] [--tx-col 1]
For every annotation row (transcript ID in --tx-col, version suffix ignored) compare the gene name in --name-col with the GTF:
  - GTF has gene_name  -> annotation name should equal it
  - GTF has no gene_name -> annotation name should be the gene_id (or anything that is not another gene's name)
Also tests the carry-over pattern of Bullseye's gtf2genepred.pl v0.2: is a wrong name equal to the most recent gene_name seen before the
transcript in the GTF? Prints counts, biotypes of affected transcripts and examples. Exit 1 if mislabelled rows are found."""
import sys, re, collections, argparse
ap = argparse.ArgumentParser(); ap.add_argument('gtf'); ap.add_argument('annot')
ap.add_argument('--name-col', type=int, default=0); ap.add_argument('--tx-col', type=int, default=1)
a = ap.parse_args()
rx_n, rx_t, rx_g, rx_b = (re.compile(r'%s "([^"]+)"' % k) for k in ('gene_name', 'transcript_id', 'gene_id', 'gene_biotype'))
tr = {}; last = None
for l in open(a.gtf):
    if l[0] == '#': continue
    f = l.split('\t'); m = rx_n.search(f[8])
    if m: last = m.group(1)
    if f[2] == 'transcript':
        b = rx_b.search(f[8])
        tr[rx_t.search(f[8]).group(1)] = (m.group(1) if m else None, rx_g.search(f[8]).group(1), last, b.group(1) if b else '?')
c = collections.Counter(); bt = collections.Counter(); ex = []; genes = set()
for l in open(a.annot):
    if l[0] == '#': continue
    p = l.rstrip('\n').split('\t'); t = p[a.tx_col].split('.')[0]; nm = p[a.name_col]
    if t not in tr: c['transcript not in GTF'] += 1; continue
    gn, gid, carried, biot = tr[t]; c['rows'] += 1
    if gn is not None:
        c['named gene: name == gene_name' if nm == gn else 'named gene: DIFFERENT name'] += 1
        if nm != gn and len(ex) < 5: ex.append((p[a.tx_col], nm, 'GTF gene_name', gn))
    elif nm == gid: c['unnamed gene: name == gene_id'] += 1
    else:
        c['unnamed gene: name is something else (mislabelled)'] += 1; bt[biot] += 1; genes.add(gid)
        if nm == carried: c['  ... equal to the previous gene_name in the GTF (carry-over)'] += 1
        if len(ex) < 5: ex.append((p[a.tx_col], nm, 'GTF gene_id', gid))
for k in ['rows', 'transcript not in GTF', 'named gene: name == gene_name', 'named gene: DIFFERENT name', 'unnamed gene: name == gene_id',
          'unnamed gene: name is something else (mislabelled)', '  ... equal to the previous gene_name in the GTF (carry-over)']:
    print(f'{c[k]:>9d}  {k}')
bad = c['named gene: DIFFERENT name'] + c['unnamed gene: name is something else (mislabelled)']
if bad:
    print(f'mislabelled rows: {bad} ({100 * bad / max(c["rows"], 1):.1f}%), unnamed genes affected: {len(genes)}; biotypes: {bt.most_common(6)}')
    for e in ex: print('  e.g.', e)
print('PASS' if not bad else 'FAIL'); sys.exit(1 if bad else 0)
