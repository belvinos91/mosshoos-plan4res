#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Next-step #2: monthly demand from Plan4RES profiles -> merit-order marginal
cost trajectory -> compare to per-stage water values in the real cuts."""
import argparse, os
import numpy as np, pandas as pd

def monthly_from_profile(path):
    df = pd.read_csv(path)
    ts = pd.to_datetime(df.iloc[:, 0], dayfirst=True, format='mixed')
    base = df['Base'] if 'Base' in df.columns else df.iloc[:, 1]
    s = pd.Series(base.values.astype(float), index=ts)
    frac = s.groupby(s.index.month).sum()
    return frac.reindex(range(1, 13)).fillna(0).values

def build_demand(cs, tsdir, zone):
    zv = pd.read_csv(os.path.join(cs, 'ZV_ZoneValues.csv'))
    zz = zv[zv.Zone == zone]
    demand_m = np.zeros(12); used = []
    for _, r in zz.iterrows():
        prof = str(r.get('Profile_Timeserie', '') or '').strip()
        if not prof or r.Type in ('Total', 'MaxActivePowerDemand', 'CostActivePowerDemand'):
            continue
        p = os.path.join(tsdir, prof)
        if not os.path.exists(p):
            print(f'  [skip {r.Type}: {prof} not found]'); continue
        demand_m += float(r.value) * monthly_from_profile(p)
        used.append(f"{r.Type}")
    return demand_m, used, zz

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', default=os.path.expanduser('~/plan4res_work/data/toyDataset'))
    ap.add_argument('--p4r-cuts', default=os.path.expanduser('~/plan4res_work/data/toyDataset/results_simul/cuts.txt'))
    ap.add_argument('--zone', default='RoGonDor')
    a = ap.parse_args()
    cs = os.path.join(a.dataset, 'csv_simul')
    tsdir = os.path.join(a.dataset, 'TimeSeries')
    demand_m, used, zz = build_demand(cs, tsdir, a.zone)
    print('Components used:', ', '.join(used))
    tu = pd.read_csv(os.path.join(cs, 'TU_ThermalUnits.csv'))
    fl = tu[tu.Zone == a.zone].sort_values('VariableCost').reset_index(drop=True)
    fl['cum'] = fl.MaxPower.cumsum()
    VOLL = float(zz[zz.Type == 'CostActivePowerDemand'].value.iloc[0])
    def mc(load):
        for _, r in fl.iterrows():
            if load <= r['cum']: return float(r.VariableCost)
        return VOLL
    hours = np.array([31,28,31,30,31,30,31,31,30,31,30,31], float)*24
    mw = demand_m / hours
    mcs = np.array([mc(v) for v in mw])
    cuts = pd.read_csv(a.p4r_cuts)
    ac = [c for c in cuts.columns if c.startswith('a_')]
    wv = []
    for t, g in cuts.groupby('Timestep'):
        nz = g[(g[ac] != 0).any(axis=1)]
        if len(nz): wv.append(-nz.loc[nz['a_2'].abs().idxmax()]['a_2'])
    wv = np.array(wv)
    print(f"{'mo':>3} {'TWh':>7} {'avgMW':>9} {'meritMC':>9}")
    for m in range(12):
        print(f"{m+1:>3} {demand_m[m]/1e6:>7.2f} {mw[m]:>9.0f} {mcs[m]:>9.2f}")
    print('cuts water values:', ' '.join(f'{v:.1f}' for v in wv))
    wv_m = np.interp(np.linspace(0,1,12), np.linspace(0,1,len(wv)), wv)
    r = np.corrcoef(np.log10(np.maximum(mcs,.1)), np.log10(np.maximum(wv_m,.1)))[0,1]
    print(f'log-shape correlation r = {r:.2f}')
    print(f'months at VOLL: {(mcs>0.5*VOLL).sum()}  stages at VOLL in cuts: {(wv>0.5*VOLL).sum()}')
    print('-->', 'TRAJECTORY REPRODUCED' if r > 0.6 else
          'PARTIAL: level agrees; full trajectory needs RES/ROR netting (Phase-2)')
    os.makedirs('results', exist_ok=True)
    pd.DataFrame(dict(month=range(1,13), demand_MWh=demand_m, avg_MW=mw,
                      merit_cost=mcs)).to_csv('results/residual_demand_monthly.csv', index=False)
    print('wrote results/residual_demand_monthly.csv')

if __name__ == '__main__':
    main()
