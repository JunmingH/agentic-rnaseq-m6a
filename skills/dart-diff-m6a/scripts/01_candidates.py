#!/usr/bin/env python3
"""Step 1a: per-group site tables from the Bullseye per-sample beds of the EDITED (active DART) samples.
Usage: 01_candidates.py <bed_dir> <fasta> <samplesheet.tsv> <out_dir>
  samplesheet.tsv : columns sample, group  (edited samples only; bed = <bed_dir>/<sample>.bed, Find_edit_site output, edited vs its control)
  fasta           : genome fasta with .fai, same chromosome names as the beds (used for the RAC flag)
env (optional): MIN_REPS (2) = a site must be called in >= MIN_REPS samples of the group; MIN_MEAN (0.05) = mean editing over those samples
Filters already applied INSIDE Bullseye (not repeated here): Find_edit_site thresholds (e.g. editing 10-95 %, >= 1.5x over control,
>= 2 C-to-U reads, coverage >= 10). Applied here, per group (= condition), in this order:
  reps  : site called in >= MIN_REPS samples of the group
  mean  : mean editing (dart_ratio, column 9) over the samples where it was called >= MIN_MEAN (Bullseye's 10 % minimum makes 5 % always true)
  rac   : flag only (not a filter here): R-A-C in transcript orientation, C = edited base. Both versions (with/without) are made downstream.
Output:
  group_sites.tsv      group, chr, pos(1-based), strand, gene, n_reps, mean_edit, rac   (gene = from the first sample, in samplesheet order, that called it)
  per_sample_calls.tsv sample, n_sites (Bullseye calls per sample, for the filter table)
  snp_candidates.tsv   unique chr, pos, strand (input for the dbSNP overlap step)
"""
import sys, os, csv, collections
bd, fa, ss_f, out = sys.argv[1:5]
MIN_REPS = int(os.environ.get('MIN_REPS', '2')); MIN_MEAN = float(os.environ.get('MIN_MEAN', '0.05'))
print(f'MIN_REPS {MIN_REPS}  MIN_MEAN {MIN_MEAN}', flush=True)
groups = collections.OrderedDict()
with open(ss_f) as fh:
    for r in csv.DictReader(fh, delimiter='\t'): groups.setdefault(r['group'], []).append(r['sample'])
miss = [s for v in groups.values() for s in v if not os.path.isfile(f'{bd}/{s}.bed')]
if miss: sys.exit('missing beds: ' + ', '.join(miss))
os.makedirs(out, exist_ok=True)
fai = {l.split('\t')[0]: tuple(int(x) for x in l.split('\t')[1:5]) for l in open(fa + '.fai')}
fh = open(fa, 'rb')
def seq(chrom, start, end):          # 1-based inclusive, upper case
    L, off, lb, lw = fai[chrom]; start = max(start, 1)
    s0 = off + (start - 1) // lb * lw + (start - 1) % lb; n = end - start + 1
    fh.seek(s0); raw = fh.read(n + n // lb + 2).replace(b'\n', b'')[:n]
    return raw.decode().upper()
comp = str.maketrans('ACGTN', 'TGCAN')
def is_rac(chrom, pos, strand):
    if chrom not in fai: return False
    if strand == '+': s = seq(chrom, pos - 2, pos)
    else: s = seq(chrom, pos, pos + 2).translate(comp)[::-1]
    return len(s) == 3 and s[0] in 'AG' and s[1] == 'A' and s[2] == 'C'
rac_cache = {}
calls = []; rows = []; cand = set()
for g, samples in groups.items():
    reps = []
    for smp in samples:
        d = {}
        n = 0
        for l in open(f'{bd}/{smp}.bed'):
            if l[0] == '#': continue
            p = l.rstrip('\n').split('\t')
            d[(p[0], int(p[2]), p[5])] = (p[3].split('|')[0], float(p[8])); n += 1
        calls.append((smp, n)); reps.append(d)
    allk = set().union(*[set(d) for d in reps])
    n_all = len(allk); n2 = n3 = 0
    for k in sorted(allk):
        v = [d[k] for d in reps if k in d]
        if len(v) < MIN_REPS: continue
        n2 += 1
        m = sum(x[1] for x in v) / len(v)
        if m < MIN_MEAN: continue
        n3 += 1
        if k not in rac_cache: rac_cache[k] = is_rac(*k)
        rows.append((g, k[0], k[1], k[2], v[0][0], len(v), round(m, 4), int(rac_cache[k])))
        cand.add(k)
    print(f'{g}: {len(samples)} samples; union of calls {n_all}; >={MIN_REPS} reps {n2}; mean>={MIN_MEAN} {n3}', flush=True)
with open(f'{out}/group_sites.tsv', 'w') as o:
    o.write('group\tchr\tpos\tstrand\tgene\tn_reps\tmean_edit\trac\n')
    for r in rows: o.write('\t'.join(map(str, r)) + '\n')
with open(f'{out}/per_sample_calls.tsv', 'w') as o:
    o.write('sample\tgroup\tn_sites\n')
    for (s, n), g in zip(calls, [g for g, v in groups.items() for _ in v]): o.write(f'{s}\t{g}\t{n}\n')
with open(f'{out}/snp_candidates.tsv', 'w') as o:
    for c, p, s in sorted(cand): o.write(f'{c}\t{p}\t{s}\n')
print('unique candidate sites for dbSNP check:', len(cand))
