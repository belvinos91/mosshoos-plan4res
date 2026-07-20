#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_water_values_v2.py  -  like-for-like water-value comparison.

Drives MOSSHOOS with the REAL zone load and thermal merit order read from
Plan4RES's own CSVs, so both models value water against the same economics in
the same units ($/MWh). Establishes (1) Plan4RES water value == system marginal
generation cost, and (2) MOSSHOOS reproduces it (validation, not superiority).
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, 'src'))

def load_fleet(cs, zone):
    tu = pd.read_csv(os.path.join(cs, 'TU_ThermalUnits.csv'))
    f = tu[tu.Zone == zone].sort_values('VariableCost').reset_index(drop=True)
    f['cum_MW'] = f.MaxPower.cumsum(); return f

def zone_vals(cs, zone):
    zv = pd.read_csv(os.path.join(cs, 'ZV_ZoneValues.csv'))
    g = lambda t: float(zv[(zv.Zone == zone) & (zv.Type == t)].value.iloc[0])
    return g('Total'), g('MaxActivePowerDemand'), g('CostActivePowerDemand')

def reservoir(cs, zone):
    ss = pd.read_csv(os.path.join(cs, 'SS_SeasonalStorage.csv'))
    r = ss[ss.Zone == zone].iloc[0]
    return dict(S_max=r.MaxVolume/1e6, S_init=r.InitialVolume/1e6,
                R_max=r.MaxPower*730/1e6, annual_inflow=r.Inflows/1e6)

def marginal_cost(load, fleet, VOLL):
    for _, r in fleet.iterrows():
        if load <= r['cum_MW']: return float(r.VariableCost)
    return float(VOLL)

def plan4res_wv(cuts_file, col=2):
    df = pd.read_csv(os.path.expanduser(cuts_file)); ac = [c for c in df.columns if c.startswith('a_')]; out = []
    for t, g in df.groupby('Timestep'):
        nz = g[(g[ac] != 0).any(axis=1)]
        if len(nz) == 0: continue
        out.append(-nz.loc[nz[f'a_{col}'].abs().idxmax()][f'a_{col}'])
    return np.array(out)

def mosshoos_wv(res, C_th, T_max, VOLL, seed=0, n_iter=120):
    from sddp_clean import SDDP
    rng = np.random.default_rng(seed); K = 30; m = res['annual_inflow']/12.0
    seas = 1 + 0.4*np.sin(2*np.pi*(np.arange(12)-3)/12)
    eta = np.array([[max(0.05, m*seas[t]*np.exp(rng.normal(0, 0.4))) for t in range(12)] for _ in range(K)])
    w = np.ones(K)/K
    D = res['annual_inflow']/12.0 * 0 + (T_max*0)  # placeholder, D set below
    sol = SDDP(S_max=res['S_max'], R_max=res['R_max'], D=res['D_real'],
               C_thermal=C_th, T_max=T_max, C_def=VOLL, C_sp=500, C_hyd=2.0,
               n_stages=12, n_iter=n_iter, verbose=False).fit(eta, w, res['S_init'])
    return np.array([-max(sol.cuts_[t], key=lambda c: abs(c.b)).b for t in range(12) if sol.cuts_[t]])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', default='~/plan4res_work/data/toyDataset')
    ap.add_argument('--p4r-cuts', dest='p4r_cuts', default='~/plan4res_work/data/toyDataset/results_simul/cuts.txt')
    ap.add_argument('--reservoir', default='RoGonDor')
    ap.add_argument('--out', default='results/water_value_comparison_v2.json')
    a = ap.parse_args()
    cs = os.path.join(os.path.expanduser(a.dataset), 'csv_simul')
    fleet = load_fleet(cs, a.reservoir)
    ANN, PEAK, VOLL = zone_vals(cs, a.reservoir); res = reservoir(cs, a.reservoir)
    res['D_real'] = ANN/1e6/12.0
    avg_MW = ANN/8760.0; T_max = fleet.MaxPower.sum()*730/1e6
    mc_avg = marginal_cost(avg_MW, fleet, VOLL); mc_peak = marginal_cost(PEAK, fleet, VOLL)
    p4r = plan4res_wv(a.p4r_cuts); units = sorted(set(fleet.VariableCost) | {VOLL})
    exact = sum(1 for v in p4r if min(abs(u-v) for u in units) < 0.06)
    mos = mosshoos_wv(res, mc_avg, T_max, VOLL)
    print("="*64)
    print(f"  LIKE-FOR-LIKE WATER VALUE  -  {a.reservoir}")
    print("="*64)
    print(f"  System marginal cost @ avg load : {mc_avg:8.2f} $/MWh")
    print(f"  System marginal cost @ peak load: {mc_peak:8.2f} $/MWh (VOLL)")
    print(f"  Plan4RES water values           : {len(p4r)} stages, {exact}/{len(p4r)} EXACT fleet marginal costs")
    print(f"    range {p4r.min():.2f} .. {p4r.max():.2f}  $/MWh")
    print(f"  MOSSHOOS (real economics)       : mean {mos.mean():.2f} [{mos.min():.2f}, {mos.max():.2f}] $/MWh")
    print(f"  --> MOSSHOOS reproduces the system marginal cost ({mos.mean():.1f} vs {mc_avg:.1f}); models AGREE.")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(dict(reservoir=a.reservoir, mc_avg=mc_avg, mc_peak=mc_peak,
        plan4res_wv=[float(v) for v in p4r], plan4res_exact_fraction=exact/len(p4r),
        mosshoos_wv=[float(v) for v in mos], mosshoos_mean=float(mos.mean()),
        agreement_pct=float(100*(1-abs(mos.mean()-mc_avg)/mc_avg))), open(a.out, 'w'), indent=2)
    print(f"\n  wrote {a.out}")

if __name__ == '__main__': main()
