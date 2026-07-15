#!/usr/bin/env python3
"""
relabel.py - correct a mischaracterisation.

M-D was labelled "Plan4RES deterministic". That is wrong and a reviewer would
catch it: Plan4RES's SSV *is* SDDP (the binary is literally `sddp_solver`, and
their bellmanvalues.csv contains 55 cuts over 19 timesteps). Plan4RES is a
STOCHASTIC method. Our M-D is our own SDDP fitted on a single mean-inflow
scenario - a deterministic-equivalent benchmark, not their method.

Until their solver runs (blocked by the HiGHS bug, see
plan4res_fixes/PLAN4RES_SDDP_BUG.md) we have NO Plan4RES output to compare
against, and we must not imply otherwise.
"""
import pathlib, re

H = pathlib.Path.home() / 'mosshoos-plan4res'
edits = 0

def sub(path, pairs):
    global edits
    p = H / path
    if not p.exists():
        print(f"  skip (missing): {path}"); return
    t = p.read_text(); n0 = t
    for a, b in pairs:
        t = t.replace(a, b)
    if t != n0:
        p.write_text(t); edits += 1
        print(f"  updated: {path}")
    else:
        print(f"  no change: {path}")

PAIRS = [
    ('M-D Plan4RES', 'M-D det.-equiv.'),
    ('M-D  Plan4RES', 'M-D  det.-equiv.'),
    ('M-D</span> Plan4RES', 'M-D</span> deterministic-equivalent'),
    ('vs deterministic (M-D)', 'vs deterministic-equivalent (M-D)'),
    ('M-D**  | Plan4RES | deterministic mean-inflow Benders cuts',
     '**M-D** | benchmark | our SDDP on a single mean-inflow scenario (deterministic-equivalent)'),
    ('| **M-D**  | Plan4RES | deterministic mean-inflow Benders cuts |',
     '| **M-D**  | benchmark | our SDDP on one mean-inflow scenario - a deterministic-equivalent, NOT Plan4RES output |'),
    ('deterministic mean-inflow Benders cuts',
     'our SDDP on a single mean-inflow scenario'),
    ('M-D  Plan4RES deterministic', 'M-D  deterministic-equivalent benchmark'),
    ("MD:'M-D Plan4RES'", "MD:'M-D det.-equiv.'"),
    ('VSS (M-D vs M-S)', 'M-S vs det.-equiv. benchmark'),
    ('M-S vs deterministic (M-D)** - the additional value of modelling uncertainty',
     'M-S vs the deterministic-equivalent benchmark (M-D)** - the value of representing uncertainty at all'),
]

print("Relabelling M-D across the repository:")
for f in ['README.md', 'RESULTS.md', 'dashboard_app.html', 'build_dashboard.py',
          'write_results.py', 'src/sddp_clean.py', 'compare_plan4res.py']:
    sub(f, PAIRS)

# add an explicit caveat to the README's comparison table
p = H / 'README.md'
if p.exists():
    t = p.read_text()
    marker = 'VSS = C(M-D) - C(M-S)'
    if marker not in t:
        marker = 'VSS = C(M-D) − C(M-S)'
    caveat = '''
> **Important.** M-D is **not** Plan4RES output. Plan4RES's seasonal storage module (SSV)
> is itself SDDP - the binary is `sddp_solver` and its `bellmanvalues.csv` holds Benders
> cuts over scenarios. It is a *stochastic* method, not a deterministic one. M-D here is
> **our own SDDP fitted on a single mean-inflow scenario**: a deterministic-equivalent
> benchmark that isolates the value of representing uncertainty at all.
>
> A true MOSSHOOS-vs-Plan4RES comparison requires running their solver, which is currently
> blocked by a bug in their SDDP on open-source-solver builds
> (see [`plan4res_fixes/PLAN4RES_SDDP_BUG.md`](plan4res_fixes/PLAN4RES_SDDP_BUG.md)).
> Until that is resolved, **no claim of superiority over Plan4RES is made here.**
>
> The comparison that does test the MOSSHOOS contribution is **M-S vs M-S0** - our
> wavelet + block-bootstrap scenarios against plain scenarios, same solver. On the current
> runs that difference is ~0%, which is an open question we are investigating, not a result
> we are claiming.
'''
    if marker in t and 'not Plan4RES output' not in t:
        t = t.replace(marker, marker + '\n' + caveat, 1)
        p.write_text(t); edits += 1
        print("  README: added the M-D caveat")

print(f"\n{edits} files changed.")
print("\nWhat this fixes: we were labelling our own deterministic-equivalent")
print("benchmark as 'the Plan4RES method'. Plan4RES is stochastic (SDDP).")
print("No Plan4RES output has been compared against yet - their solver does not run.")
