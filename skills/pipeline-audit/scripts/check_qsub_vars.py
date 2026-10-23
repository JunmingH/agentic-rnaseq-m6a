#!/usr/bin/env python3
"""Hand-over check between a caller script and a script it runs (qsub job or `bash script`) - read-only, heuristic.
Usage: python3 -I check_qsub_vars.py <caller.sh> <callee.sh>
Reports:
  - positional arguments the callee uses ($1..$9) vs the number of arguments in each call of it in the caller
  - environment variables the callee reads but neither sets itself nor receives (qsub -v, including names inside a -v "$VARS"
    string built earlier in the caller, or `export` in the caller) -> they are empty or take their default in the callee
  - variables read with a silent default ${X:-...} in callee and caller -> check that they are really set when they should be
SGE/shell built-ins (NSLOTS, TMPDIR, JOB_ID, HOME, USER, PATH, ...) are ignored. A pointer for reading the code, not a proof.
Note: for `bash callee.sh` calls the callee also inherits variables of the caller's own environment (e.g. `VAR=x bash caller.sh`):
the report lists them as missing; check whether the person running the caller sets them."""
import re, sys, os
caller, callee = sys.argv[1:3]
BUILTIN = {'NSLOTS', 'TMPDIR', 'JOB_ID', 'JOB_NAME', 'HOME', 'USER', 'PATH', 'PWD', 'SGE_TASK_ID', 'HOSTNAME', 'LANG', 'SHELL', 'OLDPWD',
           'PERL5LIB', 'LD_LIBRARY_PATH', 'CONDA_PREFIX', 'IFS', 'REPLY', 'RANDOM', 'LINENO', 'BASH_SOURCE', 'PICARD_JAR'}
NAME = r'[A-Za-z_][A-Za-z0-9_]*'
def code(path):  # drop comments (not inside quotes; good enough for these scripts) and the shebang
    out = []
    for l in open(path).read().splitlines():
        if l.lstrip().startswith('#'): continue
        out.append(re.sub(r'\s#.*$', '', l))
    return '\n'.join(out)
W, C = code(callee), code(caller)
used = set(re.findall(r'\$\{?(' + NAME + ')', W)) - BUILTIN
assigned = set(re.findall(r'(?:^|[\s;(&|])(?:local\s+|export\s+|readonly\s+)?(' + NAME + r')=', W, re.M))
assigned |= set(re.findall(r'\bfor\s+(' + NAME + r')\s+in', W))
for m in re.findall(r'\bread\s+(?:-r\s+)?((?:' + NAME + r'\s*)+)', W): assigned |= set(m.split())
pos = sorted({int(x) for x in re.findall(r'\$\{?([1-9])', W)}); shifts = len(re.findall(r'\bshift\b', W))
# caller: -v names (expand -v "$VARS,..." through VARS="..." assignments), exports
assign_c = {}
for k, v in re.findall(r'(?:^|[\s;])(' + NAME + r')=("[^"]*"|\S+)', C, re.M): assign_c.setdefault(k, []).append(v.strip('"'))
def v_names(s, depth=0):
    names = set()
    for part in s.split(','):
        part = part.strip()
        m = re.fullmatch(r'\$\{?(' + NAME + r')\}?', part)
        if m and depth < 3:
            for val in assign_c.get(m.group(1), []): names |= v_names(val, depth + 1)
        elif '=' in part: names.add(part.split('=')[0].lstrip('$'))
    return names
passed_v = set()
for m in re.finditer(r'\s-v\s+("[^"]*"|\S+)', C): passed_v |= v_names(m.group(1).strip('"'))
exported = set()
for m in re.findall(r'\bexport\s+([^;\n]+)', C):
    exported |= {t.split('=')[0] for t in m.split() if re.fullmatch(NAME + r'(=.*)?', t)}
name = os.path.basename(callee)
print(f'callee {name}: positional args used: {pos or "none"}' + (f' (with {shifts} shift)' if shifts else ''))
for l in open(caller).read().splitlines():
    if name not in l or l.lstrip().startswith('#'): continue
    after = l.split(name, 1)[1].lstrip('"\'')
    after = re.split(r'\s*(?:\||\)\s*$|;)', after)[0]
    n = len(re.findall(r'"[^"]*"|\'[^\']*\'|\S+', after))
    flag = '' if not pos or n == max(pos) else '   <-- differs from the highest $N used'
    print(f'  call passes {n} args{flag}: {l.strip()[:170]}')
missing = sorted(v for v in used - assigned - passed_v - exported if not v.isdigit())
print('env vars read by the callee, not set there, not passed (-v) and not exported by the caller:', missing or 'none')
print('callee silent defaults ${X:-...}:', sorted(set(re.findall(r'\$\{(' + NAME + r'):?-', W))) or 'none')
print('caller passes with -v:', sorted(passed_v) or 'none', '| caller exports:', sorted(exported) or 'none')
print('caller silent defaults ${X:-...}:', sorted(set(re.findall(r'\$\{(' + NAME + r'):?-', C))) or 'none')
