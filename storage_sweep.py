#!/usr/bin/env python3
"""
storage_sweep.py — VSS as a function of reservoir storage-to-inflow ratio.

This is the defensible scientific result: it characterizes WHERE stochastic
seasonal storage valuation matters, rather than reporting a single point.
Uses the validated SDDP (sddp_clean.py).
"""
import sys, os, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))
from sddp_clean import SDDP


def build_scenarios(inflow_csv, n_scen=10, block=4, n_boot=150, seed=42):
    """Wavelet-free but faithful: block-bootstrap monthly scenarios from real data."""
    s = pd.read_csv(inflow_csv, index_col=0)['inflow_GWh'].astype(float)
    x = s.values
    n_full = len(x)//12
    R = x[:n_full*12].reshape(n_full,12)
    seas = R.mean(axis=0)                      # seasonal mean
    rng = np.random.default_rng(seed)
    boots=[]
    for _ in range(n_boot):
        yrs = rng.integers(0, n_full, size=12//block+1)
        rec=[]
        for bi,yr in enumerate(yrs):
            st=bi*block
            rec.extend(R[yr, st:st+block] if st+block<=12 else R[yr, st:12])
        boots.append(np.array(rec[:12]))
    boots=np.array(boots)
    from sklearn.cluster import KMeans
    km=KMeans(n_scen, n_init=10, random_state=seed).fit(boots)
    eta=np.maximum(0.1, km.cluster_centers_)
    w=np.bincount(km.labels_, minlength=n_scen)/n_boot
    return eta, w, seas, x


def vss_at_storage(eta, w, seas, S_max, R_max, D, n_iter=200, n_mc=200, seed=1):
    common=dict(S_max=S_max, R_max=R_max, D=D, n_stages=12, n_iter=n_iter,
                verbose=False, seed=seed)
    S0=0.5*S_max
    ms  = SDDP(**common).fit(eta, w, S0)
    c_ms= ms.simulate(eta, w, n_mc, S0)['annual_cost'].mean()
    det = np.maximum(0.1, seas).reshape(1,-1)
    md  = SDDP(**common).fit(det, np.array([1.0]), S0)
    c_md= md.simulate(eta, w, n_mc, S0)['annual_cost'].mean()
    md0 = SDDP(**{**common,'n_iter':0}).fit(det, np.array([1.0]), S0)
    c_md0=md0.simulate(eta, w, n_mc, S0)['annual_cost'].mean()
    vss = (c_md-c_ms)/max(c_md,1e-9)*100
    vss0= (c_md0-c_ms)/max(c_md0,1e-9)*100     # vs greedy
    return dict(c_ms=c_ms, c_md=c_md, c_md0=c_md0, vss=vss, vss_vs_greedy=vss0)


def main():
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('--inflow', required=True)
    ap.add_argument('--out', default='storage_sweep_results.json')
    args=ap.parse_args()

    eta, w, seas, x = build_scenarios(args.inflow)
    mean_mo=float(np.mean(x))
    annual = mean_mo*12
    R_max = 200*730/1000.0          # Manantali 200 MW -> GWh/month
    D = 0.9*mean_mo                 # demand near mean inflow

    print(f"Mean monthly inflow: {mean_mo:.1f} GWh, annual {annual:.0f} GWh")
    print(f"R_max={R_max:.1f}, D={D:.1f}")
    print(f"\nSweeping storage from 0.1x to 6x mean monthly inflow...\n")

    # storage as multiples of mean MONTHLY inflow
    fractions = [0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]
    rows=[]
    print(f"{'S/mean_mo':>10} {'S_max(GWh)':>11} {'VSS%':>8} {'VSSvsGreedy%':>13} {'C_MS':>9} {'C_MD':>9}")
    for fr in fractions:
        S_max = fr*mean_mo
        r = vss_at_storage(eta, w, seas, S_max, R_max, D)
        r['frac']=fr; r['S_max']=S_max
        rows.append(r)
        print(f"{fr:>10.2f} {S_max:>11.1f} {r['vss']:>8.2f} "
              f"{r['vss_vs_greedy']:>13.2f} {r['c_ms']:>9.0f} {r['c_md']:>9.0f}")

    json.dump(dict(mean_mo=mean_mo, R_max=R_max, D=D, sweep=rows),
              open(args.out,'w'), indent=2)
    print(f"\nSaved -> {args.out}")
    # find peak
    peak=max(rows, key=lambda r:r['vss'])
    print(f"\nPeak VSS = {peak['vss']:.2f}% at storage = {peak['frac']:.2f}x mean "
          f"monthly inflow (S_max={peak['S_max']:.0f} GWh)")

if __name__=='__main__':
    main()
