import sys
out = sys.argv[1]
W, H = 1280, 930
S = []
def esc(s): return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
def rect(x, y, w, h, fill, stroke, rx=10, dash=None, sw=1.6):
    d = f' stroke-dasharray="{dash}"' if dash else ''
    S.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')
def text(x, y, s, size=13, weight='normal', fill='#1f2933', anchor='middle', style='normal'):
    S.append(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}" font-style="{style}">{esc(s)}</text>')
def arrow(pts, color='#52606d', dash=None, w=1.6):
    d = f' stroke-dasharray="{dash}"' if dash else ''
    p = ' '.join(f'{x},{y}' for x, y in pts)
    S.append(f'<polyline points="{p}" fill="none" stroke="{color}" stroke-width="{w}"{d} marker-end="url(#ah)"/>')
BLUE = ('#eef4fb', '#4a78b0'); GREEN = ('#edf7ee', '#4f9a5a'); PURP = ('#f4effa', '#7b5ea7'); ORANGE = '#b35c00'; GREY = '#616e7c'
def step(x, y, num, name, tool, decide, col, w=170, h=92):
    rect(x, y, w, h, col[0], col[1])
    text(x + w/2, y + 22, f'{num}  {name}', 13.5, 'bold', col[1])
    for i, t in enumerate(tool): text(x + w/2, y + 41 + 15*i, t, 12, fill='#1f2933')
    if decide: text(x + w/2, y + h - 10, 'decide: ' + decide, 11, fill=ORANGE, style='italic')
    return (x, y, w, h)
S.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif">')
S.append('<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#52606d"/></marker></defs>')
rect(0, 0, W, H, '#ffffff', '#ffffff', rx=0)
text(W/2, 34, 'DART-seq m6A + RNA-seq pipeline: AI agent skills on an SGE cluster', 20, 'bold')
# frame: orchestrator
rect(8, 50, W-16, 495, 'none', '#9aa5b1', rx=14, dash='6 5', sw=1.2)
text(30, 70, 'dartseq-pipeline (orchestrator skill): order, decisions, hand-over files, resume / reproduce', 12.5, 'bold', GREY, 'start')
# inputs
rect(40, 86, 400, 44, '#f5f7fa', '#9aa5b1', rx=8); text(240, 113, 'Raw paired-end fastq (_1 / _2.fq.gz)', 13, 'bold')
rect(470, 86, 440, 44, '#f5f7fa', '#9aa5b1', rx=8); text(690, 104, 'Reference: ONE Ensembl release', 13, 'bold'); text(690, 121, 'fasta + GTF, STAR index, RSEM ref, patched refFlat, dbSNP bed', 11, fill=GREY)
rect(940, 86, 300, 44, '#f5f7fa', '#9aa5b1', rx=8); text(1090, 104, 'Sample sheets (templates/)', 13, 'bold'); text(1090, 121, 'groups, control vs treated, comparisons', 11, fill=GREY)
# steps
step(40, 250, '1', 'fastqc-qc', ['FastQC per file', '+ MultiQC'], '', BLUE, w=180)
step(260, 250, '2', 'read-trimming', ['fastp or Trimmomatic', 'adapters, poly-G, quality'], 'tool, MINLEN', BLUE, w=180)
step(480, 250, '3', 'star-alignment', ['STAR index + align', 'sort, Picard duplicates'], 'RM_DUP, TRANSCRIPTOME', BLUE, w=205)
step(750, 165, '4', 'rsem-quantification', ['RSEM, control samples', 'strand test on 1 sample'], 'STRANDEDNESS', GREEN, w=190)
step(750, 345, '5', 'bullseye-analysis', ['parseBAM, Find_edit_site', 'edited vs control pairs'], 'roles, DEDUP, thresholds', PURP, w=190)
step(985, 165, '6', 'deseq2-expression', ['DESeq2 per comparison', '+ PCA'], 'comparisons, subsets', GREEN, w=180)
step(985, 345, '7', 'dart-diff-m6a', ['sites in >=2 reps, -dbSNP', 'GLM + Wald + IHW'], 'cutoffs, RAC', PURP, w=180)
step(1085, 462, '8', 'm6a-metagene', ["Guitar 5'UTR / CDS / 3'UTR"], '', PURP, w=170, h=62)
# arrows from inputs
arrow([(130, 130), (130, 250)]); arrow([(350, 130), (350, 250)])
text(138, 200, 'raw reads', 11, fill=GREY, anchor='start')
arrow([(600, 130), (600, 250)], dash='5 4'); arrow([(845, 130), (845, 165)], dash='5 4')
text(608, 158, 'reference (also Bullseye refFlat,', 11, fill=GREY, anchor='start'); text(608, 172, 'dbSNP filter, Guitar GTF)', 11, fill=GREY, anchor='start')
arrow([(1075, 130), (1075, 165)], dash='5 4')
# main flow
arrow([(220, 296), (260, 296)]); text(240, 288, 'QC', 10.5, fill=GREY)
arrow([(440, 296), (480, 296)]); text(460, 288, 'fastq', 10.5, fill=GREY)
arrow([(685, 280), (715, 280), (715, 211), (750, 211)]); text(710, 230, 'transcriptome', 11, fill=GREEN[1], anchor='end'); text(710, 244, 'BAM', 11, fill=GREEN[1], anchor='end')
arrow([(685, 312), (715, 312), (715, 391), (750, 391)]); text(710, 372, 'genome BAM', 11, fill=PURP[1], anchor='end'); text(710, 386, '(dups marked)', 11, fill=PURP[1], anchor='end')
arrow([(940, 211), (985, 211)]); text(962, 203, 'counts', 10.5, fill=GREY)
arrow([(940, 391), (985, 391)]); text(962, 383, 'beds +', 10.5, fill=GREY); text(962, 410, 'matrices', 10.5, fill=GREY)
arrow([(1075, 257), (1075, 345)]); text(1083, 305, 'DE tables', 11, fill=GREEN[1], anchor='start')
arrow([(1165, 391), (1185, 391), (1185, 462)]); text(1192, 430, 'sites', 10.5, fill=GREY, anchor='start')
# results
rect(1185, 165, 72, 150, '#fbfbfb', '#9aa5b1', rx=8, dash='4 3'); text(1221, 185, 'Results', 12.5, 'bold')
for i, s in enumerate(['QC report', 'DE tables,', 'volcano,', 'PCA', 'diff m6A', 'sites, genes', 'metagene']): text(1221, 205 + 15*i, s, 10.5, fill=GREY)
arrow([(1165, 211), (1185, 211)])
# audit bar
rect(40, 490, 1030, 46, '#fff6e6', '#d08a2e', rx=10)
text(555, 509, '9  pipeline-audit: after every step', 13.5, 'bold', '#8a4b00')
text(555, 527, 'read every script, trace each parameter to the tool, check sample placement (R1/R2, control vs edited, reference level), prove with logs + outputs', 11.5, fill='#5c3d00')
for x in (130, 350, 582): S.append(f'<line x1="{x}" y1="342" x2="{x}" y2="490" stroke="#d08a2e" stroke-width="1" stroke-dasharray="2 4"/>')
for x in (845, 1030): S.append(f'<line x1="{x}" y1="437" x2="{x}" y2="490" stroke="#d08a2e" stroke-width="1" stroke-dasharray="2 4"/>')
# legend (below the frame)
yl = 560
for i, (col, lab) in enumerate([(BLUE, 'shared steps'), (GREEN, 'RNA level (control samples)'), (PURP, 'm6A level (edited vs control)')]):
    rect(40 + i*200, yl, 14, 14, col[0], col[1], rx=3); text(60 + i*200, yl + 12, lab, 11.5, fill=GREY, anchor='start')
text(650, yl + 12, 'decide: …  = a choice the user makes (the scripts refuse to run without it)', 11.5, fill=ORANGE, anchor='start', style='italic')
# panel B: how every step runs
text(30, 620, 'How every step runs (the same loop for steps 1–8)', 16, 'bold', anchor='start')
boxes = [('1  Load the skill', ['copy scripts/ into the', 'step folder + MANIFEST.md5']),
         ('2  Write run script', ['run_<step>_<date>.sh sets', 'every required option']),
         ('3  User approves', ['sample sheet / pairs,', 'parameters, command']),
         ('4  Test 1 sample', ['check output, runtime,', 'memory (qacct)']),
         ('5  Submit the batch', ['qsub: one job per sample,', 'logs kept in error_file/']),
         ('6  Audit the step', ['pipeline-audit: code,', 'parameters, outputs']),
         ('7  Record', ['update STATUS.md,', 'then the next step'])]
bw, gap, x0, y0 = 162, 18, 22, 640
for i, (h, ls) in enumerate(boxes):
    x = x0 + i*(bw + gap)
    col = ('#fff6e6', '#d08a2e') if i in (2, 5) else ('#f5f7fa', '#7b8794')
    rect(x, y0, bw, 78, col[0], col[1], rx=9)
    text(x + bw/2, y0 + 22, h, 12.5, 'bold')
    for j, l in enumerate(ls): text(x + bw/2, y0 + 44 + 15*j, l, 11, fill=GREY)
    if i < len(boxes) - 1: arrow([(x + bw, y0 + 39), (x + bw + gap, y0 + 39)])
xe = x0 + 6*(bw + gap) + bw/2
arrow([(xe, y0 + 78), (xe, y0 + 108), (x0 + bw/2, y0 + 108), (x0 + bw/2, y0 + 80)], dash='5 4')
text((xe + x0 + bw/2)/2, y0 + 124, 'next step   ·   a failed job with a clear cause: delete only its outputs + logs and rerun; unclear: keep logs, report, wait', 11.5, fill=GREY)
# panel C: reproducibility
text(30, 818, 'What makes a run reproducible', 16, 'bold', anchor='start')
items = ['code frozen per run (script copy + md5)', 'no silent defaults for result-changing options', 'explicit sample sheets, printed before submit', 'all logs kept; seeds fixed (IHW)', 'patched Bullseye refFlat; known issues documented']
for i, s in enumerate(items):
    x = 30 + (i % 3) * 420; y = 845 + (i // 3) * 26
    S.append(f'<circle cx="{x + 6}" cy="{y - 4}" r="4" fill="#4a78b0"/>'); text(x + 18, y, s, 12.5, anchor='start')
S.append('</svg>')
open(out, 'w').write('\n'.join(S) + '\n')
