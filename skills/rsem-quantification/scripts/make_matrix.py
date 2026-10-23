#!/usr/bin/env python3
"""Merge RSEM per-sample results into gene/isoform matrices.
Usage: make_matrix.py <rsem_out_dir> <sample_list> <out_prefix>
Writes <out_prefix>.{genes,isoforms}.{expected_count,TPM}.tsv (rows = gene_id / transcript_id with --append-names suffix,
columns = samples in list order). expected_count is not an integer; DESeq2 needs round() or tximport.
Fails if a sample is missing."""
import sys, os
d, lst, pre = sys.argv[1:4]
S = [l.strip() for l in open(lst) if l.strip()]
for kind in ("genes", "isoforms"):
    for col in ("expected_count", "TPM"):
        rows, order = {}, []
        for s in S:
            f = os.path.join(d, f"{s}.{kind}.results")
            if not os.path.exists(f): sys.exit(f"missing {f}")
            with open(f) as fh:
                h = fh.readline().rstrip("\n").split("\t"); ci = h.index(col)
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if p[0] not in rows: rows[p[0]] = {}; order.append(p[0])
                    rows[p[0]][s] = p[ci]
        with open(f"{pre}.{kind}.{col}.tsv", "w") as o:
            o.write("id\t" + "\t".join(S) + "\n")
            for k in order: o.write(k + "\t" + "\t".join(rows[k][s] for s in S) + "\n")
        print(f"{pre}.{kind}.{col}.tsv: {len(order)} rows x {len(S)} samples")
