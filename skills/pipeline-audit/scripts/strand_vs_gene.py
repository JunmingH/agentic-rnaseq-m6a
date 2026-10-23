#!/usr/bin/env python3
"""Check that the strand label of parseBAM --stranded (column 9 of matrix_stranded/<sample>.matrix.gz) is the RNA (transcript) strand.
For well-known, highly expressed genes the coverage (matrix column 8, total) on exons is summed per strand label; if the label is the
transcript strand, ~100 % of an exon's coverage carries the gene's own strand. (RSEM --strandedness reverse = read 1 antisense = dUTP,
which is what parseBAM --stranded assumes; the two agree if this share is ~100 %.)
Usage: strand_vs_gene.py <refFlat> <gtf> <matrix_dir> <out.tsv> <tmp_dir> <sample> [<sample> ...]
Run with `python3 -I`; needs tabix on PATH (module load htslib).
The gene's chromosome, strand and span are taken from the Ensembl GTF 'gene' record that carries that gene_name (one record per gene here); of the refFlat
transcripts with that name only those inside that span on that strand are used, and per gene the one with the most exonic bases (refFlat columns: gene, tx,
chrom, strand, txStart, txEnd, cdsStart, cdsEnd, nExons, exonStarts, exonEnds). This matters: in this refFlat (made by gtf2genepred.pl) the transcripts of genes
that have no gene_name in the GTF carry the name of another gene (e.g. 'GAPDH' transcripts on chr12 '-', 'FTL' on chr19 '-'), so the name alone is not enough.
Output columns: sample, gene, chrom, gene_strand, cov_plus, cov_minus, pct_gene_strand (empty if no covered position, i.e. all positions < minCoverage 10)."""
import sys, os, re, subprocess
ref, gtf, md, out, tmp = sys.argv[1:6]; samples = sys.argv[6:]
GENES = ["GAPDH", "ACTB", "B2M", "EEF1A1", "TPT1", "RPL13A", "FTL", "VIM", "ACTG1", "HSP90AB1", "MALAT1", "NEAT1"]
span = {}
for l in open(gtf):
    if "\tgene\t" not in l: continue
    f = l.split("\t"); m = re.search(r'gene_name "([^"]+)"', f[8])
    if f[2] == "gene" and m and m.group(1) in GENES and m.group(1) not in span: span[m.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
best = {}
for l in open(ref):
    p = l.rstrip("\n").split("\t")
    g, chrom, strand = p[0], p[2], p[3]
    if g not in span: continue
    c0, s0, e0, st0 = span[g]
    if chrom != c0 or strand != st0 or int(p[4]) + 1 < s0 or int(p[5]) > e0: continue
    S = [int(x) for x in p[9].strip(",").split(",")]; E = [int(x) for x in p[10].strip(",").split(",")]
    size = sum(e - s for s, e in zip(S, E))
    if g not in best or size > best[g][0]: best[g] = (size, chrom.replace("chr", ""), strand, S, E)
os.makedirs(tmp, exist_ok=True)
rows = []
for smp in samples:
    m = f"{md}/{smp}.matrix.gz"
    for g in GENES:
        if g not in best: rows.append((smp, g, "", "", "", "", "")); continue
        _, chrom, strand, S, E = best[g]
        reg = f"{tmp}/{smp}_{g}.bed"
        open(reg, "w").write("".join(f"{chrom}\t{s + 1}\t{e}\n" for s, e in zip(S, E)))
        r = subprocess.run(["tabix", "-R", reg, m], stdout=subprocess.PIPE, universal_newlines=True, check=True)
        cov = {"+": 0, "-": 0}
        for l in r.stdout.splitlines():
            f = l.split("\t")
            if len(f) >= 9 and f[8] in cov: cov[f[8]] += int(f[7])
        a = cov[strand]; o = cov["-" if strand == "+" else "+"]
        rows.append((smp, g, chrom, strand, cov["+"], cov["-"], f"{100 * a / (a + o):.1f}" if a + o else ""))
with open(out, "w") as fh:
    fh.write("sample\tgene\tchrom\tgene_strand\tcov_plus\tcov_minus\tpct_gene_strand\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
tot = {}
for r in rows:
    if r[6] == "": continue
    a = r[4] if r[3] == "+" else r[5]; o = r[5] if r[3] == "+" else r[4]
    t = tot.setdefault(r[0], [0, 0]); t[0] += a; t[1] += o
for s, (a, o) in tot.items(): print(f"{s}: {100 * a / (a + o):.2f}% of exonic coverage carries the gene's own strand label ({a + o} total)")
print("wrote", out)
