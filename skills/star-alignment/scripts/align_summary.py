#!/usr/bin/env python3
"""Usage: align_summary.py <out_dir> <sample_list>  -> <out_dir>/align_summary.tsv (from STAR *_Log.final.out)"""
import sys, os
out, lst = sys.argv[1:3]
keys = [('Number of input reads', 'input_reads'), ('Uniquely mapped reads number', 'unique'), ('Uniquely mapped reads %', 'unique_pct'),
        ('% of reads mapped to multiple loci', 'multi_pct'), ('% of reads unmapped: too short', 'too_short_pct'),
        ('Mismatch rate per base, %', 'mismatch_pct'), ('Average input read length', 'avg_len')]
rows = ['\t'.join(['sample'] + [k[1] for k in keys] + ['status'])]
for s in [l.strip() for l in open(lst) if l.strip()]:
    f = f'{out}/{s}_Log.final.out'; bam = f'{out}/{s}_Aligned.sortedByCoord.out.bam'
    if not os.path.exists(f): rows.append(f'{s}\t' + '\t'.join('-' for _ in keys) + '\tNO_LOG'); continue
    d = {l.split('|')[0].strip(): l.split('|')[1].strip() for l in open(f) if '|' in l}
    ok = 'ok' if os.path.exists(bam + '.bai') else 'NO_INDEX'
    rows.append(f'{s}\t' + '\t'.join(d.get(k[0], '-') for k in keys) + f'\t{ok}')
open(f'{out}/align_summary.tsv', 'w').write('\n'.join(rows) + '\n'); print('\n'.join(rows[:3]), '...', len(rows) - 1, 'samples')
