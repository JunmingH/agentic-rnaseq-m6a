#!/bin/bash -l
# qsub worker: $1 snp_candidates.tsv (chr pos strand), $2 dbsnp C>T/G>A bed (chr start0 end ref alt), $3 out snp_hits.tsv
# A candidate is a dbSNP hit if a SNV sits on the same position with the matching change on the + strand of the genome:
#   '+' strand site (C>T)  <-> SNV C>T ;  '-' strand site (C>T on the minus strand = G>A on the genome) <-> SNV G>A
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR { want[$1":"$2":"($3=="+"?"C>T":"G>A")]=$3; next }
     { k=$1":"$3":"$4">"$5; if (k in want) print $1,$3,want[k],$4">"$5 }' "$1" "$2" | sort -u > "$3".tmp \
  && mv "$3".tmp "$3" && echo "dbSNP hits: $(wc -l < "$3")"
