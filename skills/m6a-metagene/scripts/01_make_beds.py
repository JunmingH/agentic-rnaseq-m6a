#!/usr/bin/env python3
"""Per-group site BED files (no-RAC and RAC versions) for the Guitar distribution plots, from the dart-diff-m6a site tables.
Usage: 01_make_beds.py <diff_m6a_sites_dir> <out_bed_dir>
Input  : <sites_dir>/group_sites.tsv (group chr pos strand gene n_reps mean_edit rac), snp_hits.tsv (chr pos strand change),
         filter_counts.tsv (used only to check that the rebuilt site sets have the same sizes as the ones used by dart-diff-m6a)
Filters, same order as dart-diff-m6a (03_matrix.py), per group (= condition):
  Bullseye -> called in >= MIN_REPS samples and mean editing >= MIN_MEAN (already applied inside group_sites.tsv)
  -> [RAC version only] R-A-C in transcript orientation (rac == 1) -> not in dbSNP (snp_hits.tsv)
Output : <out>/<group>_norac.bed, <out>/<group>_rac.bed   BED6 (chr, pos-1, pos, gene, mean_edit, strand), sorted; groups.txt (group order)
Groups = order in group_sites.tsv;
groups that are in no comparison have no filter_counts entry and are reported as 'not checked'."""
import sys, os, collections
sd, out = sys.argv[1:3]
os.makedirs(out, exist_ok=True)
snp = {tuple(l.rstrip('\n').split('\t')[:3]) for l in open(os.path.join(sd, 'snp_hits.tsv'))}
sites = collections.OrderedDict()
for l in list(open(os.path.join(sd, 'group_sites.tsv')))[1:]:
    g, c, p, s, gene, nr, m, rac = l.rstrip('\n').split('\t')
    sites.setdefault(g, [])
    if (c, p, s) in snp: continue
    sites[g].append((c, int(p), s, gene, m, int(rac)))
exp = {}                                     # expected sizes from dart-diff-m6a: step 'no_dbSNP' per condition and version
for l in list(open(os.path.join(sd, 'filter_counts.tsv')))[1:]:
    cmp_, ver, cond, step, n = l.rstrip('\n').split('\t')
    if step == 'no_dbSNP': exp[(cond, ver)] = int(n)
bad = 0
for g in sites:
    for ver in ('norac', 'rac'):
        rows = sorted((r for r in sites[g] if ver == 'norac' or r[5] == 1), key=lambda r: (r[0], r[1]))
        with open(os.path.join(out, f'{g}_{ver}.bed'), 'w') as fh:
            for c, p, s, gene, m, rac in rows: fh.write(f'{c}\t{p - 1}\t{p}\t{gene}\t{m}\t{s}\n')
        e = exp.get((g, ver))
        st = 'not checked' if e is None else ('OK' if e == len(rows) else 'MISMATCH')
        if st == 'MISMATCH': bad += 1
        print(f'{g:8s} {ver:6s} {len(rows):7d} sites   diff_m6a: {e}   {st}')
with open(os.path.join(out, 'groups.txt'), 'w') as fh: fh.write('\n'.join(sites) + '\n')
if bad: sys.exit(f'{bad} site sets differ from the sets used by dart-diff-m6a')
print(f'{len(sites)} groups written; checked sets match filter_counts.tsv')
