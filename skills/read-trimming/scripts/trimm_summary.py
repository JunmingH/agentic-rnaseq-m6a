#!/usr/bin/env python3
"""Usage: trimm_summary.py <out_dir> <sample_list>  -> TSV on stdout
Works for both tools: <sample>_trimm_summary.txt (Trimmomatic) or <sample>_fastp.json (fastp). Missing -> MISSING."""
import sys, os, json
out, lst = sys.argv[1:3]
print('sample\ttool\tinput_pairs\tboth_surviving\tboth_pct\tforward_only\treverse_only\tdropped\tdropped_pct\tq30_before\tq30_after\tstatus')
for s in [l.strip() for l in open(lst) if l.strip()]:
    tp, fp = f'{out}/{s}_trimm_summary.txt', f'{out}/{s}_fastp.json'
    if os.path.exists(tp):
        d = dict(l.rstrip('\n').split(': ', 1) for l in open(tp) if ': ' in l)
        g = lambda k: d.get(k, '-')
        print('\t'.join([s, 'trimmomatic', g('Input Read Pairs'), g('Both Surviving Reads'), g('Both Surviving Read Percent'),
                         g('Forward Only Surviving Reads'), g('Reverse Only Surviving Reads'), g('Dropped Reads'), g('Dropped Read Percent'), '-', '-', 'ok']))
    elif os.path.exists(fp):
        j = json.load(open(fp)); sm = j['summary']
        inp = sm['before_filtering']['total_reads'] // 2; ok = j['filtering_result']['passed_filter_reads'] // 2
        # fastp reports whole pairs; reads whose mate failed are written to the unpaired files and are not counted separately
        print('\t'.join([s, 'fastp', str(inp), str(ok), f'{ok/inp*100:.2f}', '-', '-', str(inp - ok), f'{(inp-ok)/inp*100:.2f}',
                         f"{sm['before_filtering']['q30_rate']:.4f}", f"{sm['after_filtering']['q30_rate']:.4f}", 'ok']))
    else:
        print(f'{s}\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tMISSING')
