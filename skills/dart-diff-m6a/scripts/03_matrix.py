#!/usr/bin/env python3
"""Step 1b + 2: apply the dbSNP filter, report site counts after every filter, build the union of kept sites per comparison and extract
mut/cov for ALL edited samples of the two groups (also where the site was not called).
Usage: 03_matrix.py <group_sites.tsv> <per_sample_calls.tsv> <snp_hits.tsv> <matrix_dir> <comparisons.tsv> <out_dir>
  comparisons.tsv : columns name, treated, control (group names of the samplesheet used in step 01)
  matrix_dir      : parseBAM matrices <sample>.matrix.gz(.tbi) of the edited samples; STRANDED (--stranded, 9 columns) or unstranded (8)
Filter order per condition: reps -> mean (step 01) -> [RAC, only in the RAC version] -> no dbSNP ; then union of the two conditions.
Matrix columns chr pos A T C G N total [strand]: cov = total (col 8, same definition as Bullseye; includes N, so at a few positions it can exceed A+C+G+T by 1-2 reads);
  mut = T count at '+' strand sites, A count at '-' strand sites. Unstranded matrix: the single row of the position is used for both strands.
Output per comparison: <out>/<cmp>/sites_kept.tsv (union, column in_rac = passes RAC), counts.tsv (mut/cov long format); one shared filter_counts.tsv.
"""
import sys, os, csv, subprocess, collections
gs, pc, sh, md, cmp_f, out = sys.argv[1:7]
with open(cmp_f) as fh: CMPS = [(r['name'], r['treated'], r['control']) for r in csv.DictReader(fh, delimiter='\t')]
snp = {tuple(l.split('\t')[:3]) for l in open(sh)}                      # (chr, pos, strand)
samples = collections.defaultdict(list); calls = {}
for l in list(open(pc))[1:]:
    s, g, n = l.rstrip('\n').split('\t'); samples[g].append(s); calls[s] = int(n)
for name, trt, ctl in CMPS:
    for g in (trt, ctl):
        if g not in samples: sys.exit(f'{name}: group {g} not in {pc}')
sites = collections.defaultdict(dict)                                     # group -> key -> (gene, nrep, mean, rac)
for l in list(open(gs))[1:]:
    g, c, p, s, gene, nr, m, rac = l.rstrip('\n').split('\t')
    sites[g][(c, p, s)] = (gene, int(nr), float(m), int(rac))
fc = open(f'{out}/filter_counts.tsv', 'w'); fc.write('comparison\tversion\tcondition\tstep\tn_sites\n')
def log(cmp_, ver, cond, step, n): fc.write(f'{cmp_}\t{ver}\t{cond}\t{step}\t{n}\n')
for name, trt, ctl in CMPS:
    d = f'{out}/{name}'; os.makedirs(d, exist_ok=True)
    keep = {}                                                             # key -> [gene, in_rac]
    for ver in ('norac', 'rac'):
        un = set()
        for cond in (trt, ctl):
            if ver == 'norac':
                for i, smp in enumerate(samples[cond], 1): log(name, 'both', cond, f'bullseye_calls_rep{i}', calls[smp])
            s0 = sites[cond]                                              # already reps and mean filtered (01_candidates.py)
            log(name, ver, cond, 'reps_and_mean', len(s0))
            s1 = {k: v for k, v in s0.items() if ver == 'norac' or v[3] == 1}
            if ver == 'rac': log(name, ver, cond, 'RAC', len(s1))
            s2 = {k: v for k, v in s1.items() if k not in snp}
            log(name, ver, cond, 'no_dbSNP', len(s2))
            un |= set(s2)
            for k, v in s2.items(): keep.setdefault(k, [v[0], v[3]])
        log(name, ver, 'union', 'union_of_conditions', len(un))
    with open(f'{d}/sites_kept.tsv', 'w') as o:
        o.write('chr\tpos\tstrand\tgene\tin_rac\n')
        for k in sorted(keep, key=lambda k: (k[0], int(k[1]))): o.write('\t'.join(k) + f'\t{keep[k][0]}\t{keep[k][1]}\n')
    # matrix extraction: tabix -R per sample (regions = kept positions); keep the row with the right strand; missing row -> mut=0, cov=0
    reg = f'{d}/regions.tmp'
    with open(reg, 'w') as o:
        for c, p in sorted({(k[0], int(k[1])) for k in keep}): o.write(f'{c}\t{p}\t{p}\n')
    with open(f'{d}/counts.tsv', 'w') as o:
        o.write('chr\tpos\tstrand\tsample\tcondition\tmut\tcov\n')
        for cond in (ctl, trt):
            for smp in samples[cond]:
                got = {}
                pr = subprocess.run(['tabix', '-R', reg, f'{md}/{smp}.matrix.gz'], stdout=subprocess.PIPE, universal_newlines=True, check=True)
                for l in pr.stdout.splitlines():
                    f = l.split('\t')
                    if len(f) < 8: continue
                    for st in ([f[8]] if len(f) >= 9 else ['+', '-']):     # unstranded: the same row serves both strands
                        k = (f[0], f[1], st)
                        if k in keep: got[k] = (int(f[3]) if st == '+' else int(f[2]), int(f[7]))   # + : T (col 4) ; - : A (col 3)
                for k in sorted(keep, key=lambda k: (k[0], int(k[1]))):
                    m, cv = got.get(k, (0, 0)); o.write(f'{k[0]}\t{k[1]}\t{k[2]}\t{smp}\t{cond}\t{m}\t{cv}\n')
    os.remove(reg)
    print(name, 'kept (union, no RAC filter):', len(keep), ' with RAC:', sum(v[1] for v in keep.values()), flush=True)
fc.close()
