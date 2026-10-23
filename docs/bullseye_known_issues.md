# Bullseye: issues found while running it on DART-seq data

Bullseye (Flamand & Meyer; <https://github.com/mflamand/Bullseye>, MIT license) is the C-to-U calling tool used by `skills/bullseye-analysis`.
This repository bundles three unmodified scripts at the versions we ran:

| script | bundled version | upstream `main` (checked 2026-10-23) |
|---|---|---|
| `parseBAM.pl` | 1.4.4 (2022-12-06) | one line differs (exit code of the usage message) |
| `Find_edit_site.pl` | 1.5.2 (2022-12-12) | renamed `Find_edit_sites.pl`, v1.6.1 (2026-05-23), uses `Bio::DB::HTS::Faidx` |
| `summarize_sites.pl` | 1.21 | identical |

The bundled copies are kept as they are, so results stay reproducible. The table below lists what we found and what was done.

| # | Issue | File | Effect | Status |
|---|---|---|---|---|
| 1 | refFlat annotation built wrongly: genes without `gene_name` take the previous gene's name; CDS taken from start/stop codons only | `gtf2genepred.pl` (makes the refFlat; not bundled) | wrong gene names on 31,617 transcripts and wrong CDS boundaries on 8,618 coding transcripts (Ensembl 111) | **Fixed**: `patches/gtf2genepred_refflat_fix.patch` |
| 2 | Merging stranded matrices skips sites | `Find_edit_site.pl` | a few percent of stranded sites lost (unstranded runs not affected); reproduced on toy data | **Documented, not fixed** |
| 3 | `Find_edit_site.pl --stranded` does nothing | `Find_edit_site.pl` | none; strand comes from column 9 of the matrix | Documented |
| 4 | Flag variable names in `get_strand` are swapped | `parseBAM.pl` | none; the labels are correct for dUTP (ISR) libraries | Documented |

---

## 1. refFlat from `gtf2genepred.pl` (fixed)

**Problem A, gene names.** `$1` from the previous regex match is not cleared. A gene without `gene_name` therefore gets the name of the gene parsed before it. In Ensembl 111 that affects 20,601 genes (31,617 transcripts). For example, `ENSG00000285238` came out as `GAPDH`.

**Problem B, coding region.** `cdsStart`/`cdsEnd` are derived from `start_codon`/`stop_codon` only, which goes wrong for 8,618 of 111,163 coding transcripts:
- If the start codon is missing on a `+` strand transcript, `cdsStart > cdsEnd` (5,562 invalid rows).
- If both codons are missing, the whole transcript looks non-coding.
- If a codon is split across exons, the boundary is off by 1–2 nt.

Bullseye uses these boundaries to assign 5'UTR/CDS/3'UTR to each site.

**Fix.** `patches/gtf2genepred_refflat_fix.patch` does two things:
- Genes without a name are named by their `gene_id`.
- `cdsStart`/`cdsEnd` become the span of the GTF `CDS` + `start_codon` + `stop_codon` features.

How the fix was checked:
- Every field of the output was compared against the GTF: 0 mismatches.
- A positive control with 3 injected errors was detected.
- The patch applies cleanly to the v0.2 script we used, and the result equals our fixed script (apart from comment lines).
- It also applies to upstream `main` (2026-10-23). On a small test GTF that patched upstream script gives the same refFlat as ours, while unpatched upstream still shows Problem A.

Effect on results: the refFlat does not change which positions are called, only the gene names and region labels attached to them (and, through the names, gene-level summaries).

```bash
cp /path/to/Bullseye/Code/gtf2genepred.pl gtf2genepred_fixed.pl
patch gtf2genepred_fixed.pl < patches/gtf2genepred_refflat_fix.patch
perl gtf2genepred_fixed.pl --gtf Homo_sapiens.GRCh38.111.gtf --out Homo_sapiens.GRCh38.111.fixed.refFlat
python3 skills/pipeline-audit/scripts/check_annotation_vs_gtf.py Homo_sapiens.GRCh38.111.gtf Homo_sapiens.GRCh38.111.fixed.refFlat
```

## 2. Stranded matrix merge in `Find_edit_site.pl` (documented, not fixed)

`parseBAM.pl --stranded` writes up to two rows per position, `+` then `-`. `Find_edit_site.pl` walks the edited and the control matrix together. When the strands of the current rows differ, the code assumes `-` sorts before `+`; the comment reads `# case when dart is ahead of control # - before + #`. The matrix order is the opposite. In addition, after a strand mismatch the edited line is read once more. So whenever one of the two matrices has rows on both strands at a position, the site there can be skipped. The walk can then also fall out of step and skip the following sites.

**Minimal reproduction** ([`repro_find_edit_stranded_merge/`](repro_find_edit_stranded_merge/), toy data, about 10 s):

```bash
bash docs/repro_find_edit_stranded_merge/run_repro.sh skills/bullseye-analysis/scripts/Find_edit_site.pl
```

The toy data is one antisense gene pair and four true edited sites (`1100 -`, `1200 +`, `1300 +`, `1400 -`). Each site is called when it is tested alone.

| run | matrices | sites called |
|---|---|---|
| A | as `parseBAM --stranded` writes them (control has `+` and `-` rows at 1100; edited has `+` and `-` at 1300) | **1 of 4** (`1300 +`) |
| B | same, with only those two extra opposite-strand rows removed (they carry no edits) | 4 of 4 |

The bundled 1.5.2 and upstream `Find_edit_sites.pl` 1.6.1 (2026-05) give the same result. On real stranded data the loss is a few percent of all sites, concentrated where genes overlap on opposite strands; no sites are added. Unstranded matrices (8 columns) go through another code path and are not affected.

**Suggested fix (not applied).** Advance the edited and control rows together by (position, strand order `+` < `-`), move only the side that is behind, and drop the extra read after a strand mismatch. Re-sorting the matrices to put `-` first is not enough: it only moves the loss to the `+` strand.

**Consequence.** Stranded runs lose sites, so differences between stranded and unstranded site counts are mostly caused by this bug, not by strand information. The bundled script is kept unchanged so that results stay reproducible, and every stranded result made with this repository (including `dart-diff-m6a`) carries the loss. It has not been reported upstream yet.

## 3. `Find_edit_site.pl --stranded` is a no-op

The option is parsed (`my $stranded`, `GetOptions`) and never used again. Strand handling switches on automatically when the matrix has a 9th (strand) column. Passing or omitting `--stranded` makes no difference; use stranded matrices from `parseBAM.pl --stranded` instead.

## 4. `parseBAM.pl get_strand`: swapped variable names

In `get_strand($first, $second, $reverse, $m_reverse)`, `$first` actually holds flag 0x80 (read 2) and `$second` holds 0x40 (read 1). With that mapping the result is right for dUTP / Illumina stranded (ISR) libraries, where read 1 is antisense; this is the same as RSEM `--strandedness reverse`. Do not "fix" the names by swapping `+`/`-`. For other library types the labels are wrong, so test the labels on one sample first (`skills/bullseye-analysis/scripts/strand_check.sh`): C→T should be high on `+` rows, and G→A high on `-` rows.
