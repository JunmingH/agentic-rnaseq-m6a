# Pitfalls

Read before reporting any failure. Part A: bugs in the checking code that produced false alarms in the reference audit. Part B: real traps in tools and data.

## A. The audit itself (false alarms)
| Symptom | Cause | How to avoid |
|---|---|---|
| "N/N parseBAM inputs wrong" (every sample) | `grep -m1 '^done '` matched `done with chromosome: 18` instead of the final `done /path.bam` | match the exact final marker (`grep '^done /' | tail -1`); print what was parsed for 1-2 samples before trusting a count |
| "N PATH MISMATCH" (every sample) in STAR logs | `readFilesIn` line ends with `~RE-DEFINED`; wrong field taken | print the raw line once; use `awk '{print $2, $3}'` on the first `^readFilesIn` line |
| samtools "command not found", every ID "different" | `module load ... | head` runs in a subshell, the environment is lost; an empty list was compared | never pipe `module load`; check that inputs of a comparison are non-empty (`wc -l`) |
| matrix coverage "between" dedup and non-dedup depth | `samtools mpileup` excludes UNMAP, SECONDARY, QCFAIL, DUP by default | isolate one flag: `--ff 0x4` vs `--ff 0x404` |
| "0 transcripts cover the site" | function called with (chr, pos, strand) but defined as (chr, strand, pos) | assert on a known positive case before running on all |
| `.ti` "1,002,754 differences" | RSEM `.ti` is a multi-line record format, not one ID per line | compare IDs from `isoforms.results` |
| relative paths fail | two parallel shell calls; one `cd` changed the cwd of the other | absolute paths in parallel calls |
| a job was submitted although step 1 failed | `python ... | tee log || exit 1` - tee's status wins | `set -o pipefail` |
| gene picked by name lands on another locus | refFlat names were wrong (see B); `GAPDH` also named an unnamed gene on chr12 '-' | select genes by GTF gene record (chrom, strand, span), not by name alone |
| hundreds of rows "fold < 1.5" | bed ratios printed with 4 decimals (2/29 vs 3/29 shows 1.4986) | test whether the threshold lies in the rounding interval |
| "FAIL: name changed to another gene" | overlapping identical transcripts; Bullseye picks either name | check both names cover the site; then accept as tie |

## B. Tools and data
- **`gtf2genepred.pl` v0.2 (Bullseye repo)**: `$group =~ /gene_name "(.+?)"/; my $name = $1;` - a failed match keeps `$1`, so every transcript of a gene without `gene_name` (Ensembl: 20,601 genes / 31,617 transcripts in release 111, mostly lncRNA) gets the previous named gene's name. Sites unchanged, gene names and any gene-level summary wrong. Fix: `my ($name) = /gene_name "(.+?)"/; ($name) = /gene_id "(.+?)"/ unless defined $name;`
- **Find_edit_site.pl `--stranded`**: parsed, never used. Strand comes from the matrix.
- **parseBAM `--stranded`**: designed for dUTP (= RSEM reverse); label = RNA strand. For other library types the labels are wrong - test, don't trust the word "forward".
- **RSEM**: printed commands omit strandedness; `.model` line 3 has it. Users' "forward stranded" may mean the opposite convention.
- **DESeq2**: test uses raw counts + size factors; `counts(dds, normalized=TRUE)` only for display; padj unaffected by `lfcShrink`.
- **STAR**: `versionGenome 2.7.4a` with STAR 2.7.11a is normal. With `--outSAMunmapped Within` STAR 2.7.11a also writes the unmapped reads into `Aligned.toTranscriptome.out.bam`, so RSEM's N0 (`.cnt` line 1) equals STAR's unmapped count (checked on every sample of a project, difference 0) - it is not "mapped to the genome but not to a transcript". Exclude flag 0x4 (`-F 0x4`) when computing strand shares from that BAM.
- **Scripts changed after the run**: e.g. `align.job.sh` edited (htslib fix) after STAR had run and `MANIFEST.md5` rewritten for the new version; the manifest then describes the current code, not what ran. Compare script mtimes with output/log times and use the tool's recorded command (Log.out, `@PG`) as the record.
- **gtf2genepred.pl v0.2, coding region**: cdsStart/cdsEnd are derived from `start_codon`/`stop_codon` only (default cdsStart = cdsEnd = txEnd); the `CDS` features are never read. Ensembl 111: 8,618 of 111,163 coding transcripts wrong (no start codon on '+' -> cdsStart > cdsEnd, 5,562 invalid rows; neither codon -> looks non-coding; split codons -> 1-2 nt off). Bullseye uses these columns to label sites 5UTR/CDS/3UTR/ncRNA; sites themselves unchanged. Fix: span of CDS + start_codon + stop_codon features.
- **fastp**: `--thread N` uses > N cores (killed on N slots); adapter "unspecified" for R2 is common.
- **samtools 1.19.2** (some module systems) needs `module load htslib/1.19.1`; **bcftools 1.18** needs `htslib/1.18`; compute-node `python3` is 3.6 (no `capture_output`).
- **Guitar 2.12**: see checklists section 8.
- **Literature**: Flamand 2022 (PMC9732904) is the TRIBE-STAMP paper, not DART-seq m6A; DART-seq site rules: Meyer 2019 (PMC6884681), Tegowski 2024 (PMC11614689). PMC pages may be blocked - use `eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=<id>&rettype=xml`.
