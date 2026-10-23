# Patches to third-party code

`gtf2genepred_refflat_fix.patch` fixes Bullseye's `gtf2genepred.pl`, the script that makes the refFlat annotation for `Find_edit_site.pl`. It fixes two problems:
- gene names carried over from the previous gene;
- CDS boundaries taken from the start/stop codons only.

Details and checks are in `docs/bullseye_known_issues.md` (issue 1). Only the diff is distributed here. The original script is MIT-licensed, (c) 2021 mflamand, <https://github.com/mflamand/Bullseye>.

Base versions it was tested on:
- v0.2, as used in our analyses (md5 `173a2c9622a397c645b311913a06895c`): the patch applies exactly; the result is identical to our fixed script except for comment lines.
- upstream `main` of 2026-10-23 (md5 `a74a0ccfe7b53d31e5a482b5e78c7de2`): the patch applies with line offsets, and the output on a test GTF is identical.

```bash
cp Bullseye/Code/gtf2genepred.pl gtf2genepred_fixed.pl
patch gtf2genepred_fixed.pl < gtf2genepred_refflat_fix.patch
```
