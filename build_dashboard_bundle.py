#!/usr/bin/env python3
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, 'src')
R = 'results'

def csv_text(p):
    return open(p).read() if os.path.exists(p) else ''

b = {}
# water-value comparison (Plan4RES vs MOSSHOOS)
b['wv2'] = json.load(open(os.path.join(R, 'water_value_comparison_v2.json')))
# storage sweep
sw = json.load(open('storage_sweep_results.json'))
b['storage_sweep'] = sw
# policy costs from the MC runs
pol = {}
for m in ['MS', 'MS0', 'MD', 'MD0']:
    f = os.path.join(R, f'WP4_MC_{m}.csv')
    if os.path.exists(f):
        d = pd.read_csv(f)
        pol[m] = dict(mean=round(float(d['annual_cost'].mean()), 1),
                      p90=round(float(np.percentile(d['annual_cost'], 90)), 1),
                      deficit=round(float(d['deficit_gwh'].mean()), 2))
b['policies'] = pol
# summary (wavelet/scenarios/convergence)
summ = {}
sf = os.path.join(R, 'MOSSHOOS_summary.csv')
if os.path.exists(sf):
    for line in open(sf):
        if ',' in line:
            k, v = line.strip().split(',', 1); summ[k] = v
b['summary'] = summ
# scenarios
sc = []
scf = os.path.join(R, 'WP3_scenarios.csv')
if os.path.exists(scf):
    d = pd.read_csv(scf, index_col=0)
    stage_cols = [c for c in d.columns if c not in ('weight', 'type')]
    for i, (_, row) in enumerate(d.iterrows()):
        sc.append(dict(id=i+1, type=str(row.get('type', 'median')),
                       weight=round(float(row.get('weight', 0)), 4),
                       path=[round(float(row[c]), 2) for c in stage_cols]))
b['scenarios'] = sc
# validation numbers (from summary)
b['validation'] = [
    ['SDDP', 'M-S vs M-D, converged', f"{summ.get('cost_saving_pct','?')}% , conv {summ.get('sddp_converged','?')}"],
    ['SDE convergence', 'EM 0.5 / Milstein 1.0', 'path-coupled (see run log)'],
    ['Wavelet-Galerkin', 'O(M^-1) for H1', 'see WP4_WG_convergence.csv'],
]
# embedded files for download
b['files'] = {}
for name, path in [('OMVS_inflows.csv', 'data/OMVS_inflows.csv'),
                   ('WP3_scenarios.csv', os.path.join(R, 'WP3_scenarios.csv')),
                   ('MOSSHOOS_summary.csv', os.path.join(R, 'MOSSHOOS_summary.csv')),
                   ('water_value_comparison.json', os.path.join(R, 'water_value_comparison_v2.json'))]:
    t = csv_text(path)
    if t: b['files'][name] = t
# storage sweep as CSV
rows = ['frac,S_max_GWh,vss_vs_greedy_pct,vss_vs_determ_pct']
for r in sw['sweep']:
    rows.append(f"{r['frac']},{round(r['S_max'],1)},{round(r['vss_vs_greedy'],2)},{round(r['vss'],2)}")
b['files']['storage_sensitivity.csv'] = '\n'.join(rows)

json.dump(b, open(os.path.join(R, 'dashboard_bundle.json'), 'w'))
print('bundle written:', os.path.join(R, 'dashboard_bundle.json'))
print('  scenarios:', len(sc), '| policies:', list(pol.keys()),
      '| wv agreement:', round(b['wv2'].get('agreement_pct', 0), 1), '%')
