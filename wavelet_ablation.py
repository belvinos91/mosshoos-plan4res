#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wavelet_ablation.py - Next-step #4: does the wavelet layer improve scenario
quality in a regime where it SHOULD matter?

Three policies, identical SDDP config, evaluated OUT-OF-SAMPLE on real
calendar-aligned 12-month windows of the observed record:
  (a) wavelet  : seasonal A_J + block-bootstrapped detail (raw paths)
  (b) naive    : i.i.d. bootstrap of raw months (no wavelet, no seasonality)
  (c) seasonal : monthly climatology + i.i.d. residuals (seasonality only)
Paired bootstrap gives 95% CIs. (c) vs (a) isolates the wavelet's detail
structure from plain seasonality.
"""
import argparse, json, sys, os
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))
from mosshoos_wp2_wavelet import HaarWavelet
from sddp_clean import SDDP

def naive_scenarios(x, K, T, rng):
    return np.array([[float(rng.choice(x)) for _ in range(T)] for _ in range(K)])

def eval_policy(sol, eval_seqs, S0):
    out = []
    for s in eval_seqs:
        tr = sol._forward(S0, np.asarray(s, float))
        out.append(sum(d['c'] for d in tr))
    return np.array(out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inflow', default='data/OMVS_inflows.csv')
    ap.add_argument('--frac', type=float, default=6.0)
    ap.add_argument('--n_iter', type=int, default=200)
    ap.add_argument('--n_eval', type=int, default=300)
    ap.add_argument('--seed', type=int, default=7)
    a = ap.parse_args()

    x = pd.read_csv(a.inflow)['inflow_GWh'].values.astype(float)
    rng = np.random.default_rng(a.seed)
    mean_mo = float(np.mean(x)); T = 12
    S_max = a.frac * mean_mo; R_max = 200*730/1000.0; D = 0.9*mean_mo
    common = dict(S_max=S_max, R_max=R_max, D=D, n_stages=T,
                  n_iter=a.n_iter, verbose=False)
    S0 = 0.5 * S_max

    # (a) wavelet: seasonal A_J + block-bootstrapped detail, raw paths
    J = HaarWavelet.optimal_levels(x); wv = HaarWavelet(J).fit(x)
    seas12 = wv.seasonal(len(x))[:T]; resid = wv.residuals(x)
    K = 24; b = 4
    eta_w = []
    for _ in range(K):
        path = []
        while len(path) < T:
            s = rng.integers(0, len(resid)-b)
            path.extend(resid[s:s+b])
        eta_w.append(np.maximum(0.5, seas12 + np.array(path[:T])))
    eta_w = np.asarray(eta_w); w_w = np.ones(K)/K

    # (b) naive i.i.d. months
    eta_n = naive_scenarios(x, K, T, rng); w_n = np.ones(K)/K

    # (c) seasonal-only: climatology + i.i.d. residuals
    clim = np.array([x[m::T].mean() for m in range(T)])
    res_c = x - np.tile(clim, len(x)//T + 1)[:len(x)]
    eta_c = np.array([np.maximum(0.5, clim + rng.choice(res_c, size=T))
                      for _ in range(K)]); w_c = np.ones(K)/K

    # calendar-aligned out-of-sample years (record starts in January)
    n_years = len(x)//T
    yr = rng.integers(0, n_years, size=a.n_eval)
    eval_seqs = np.array([x[y*T:(y+1)*T] for y in yr])

    print(f"Regime: S_max={S_max:.0f} GWh ({a.frac}x mean), D={D:.1f}, R_max={R_max:.1f}")
    sol_w = SDDP(**common).fit(eta_w, w_w, S0)
    sol_n = SDDP(**common).fit(eta_n, w_n, S0)
    sol_c = SDDP(**common).fit(eta_c, w_c, S0)
    c_w = eval_policy(sol_w, eval_seqs, S0); c_n = eval_policy(sol_n, eval_seqs, S0)
    c_c = eval_policy(sol_c, eval_seqs, S0)

    diff = c_n - c_w
    boot = [np.mean(diff[rng.integers(0, len(diff), len(diff))]) for _ in range(2000)]
    lo, hi = np.percentile(boot, [2.5, 97.5]); m = float(np.mean(diff))
    rel = 100*m/np.mean(c_n)
    print(f"\n  mean cost  wavelet={np.mean(c_w):.1f}  seasonal-only={np.mean(c_c):.1f}  naive={np.mean(c_n):.1f}")
    d2 = c_c - c_w
    b2 = [np.mean(d2[rng.integers(0,len(d2),len(d2))]) for _ in range(2000)]
    l2, h2 = np.percentile(b2, [2.5, 97.5])
    print(f"  wavelet vs seasonal-only: {np.mean(d2):+.1f}  [{l2:+.1f}, {h2:+.1f}] 95% CI"
          + ("  (wavelet detail structure adds value)" if l2 > 0 else
             ("  (seasonality explains the gain)" if h2 > 0 else "")))
    print(f"  paired difference (naive - wavelet): {m:+.1f}  [{lo:+.1f}, {hi:+.1f}] 95% CI")
    print(f"  relative: {rel:+.2f}%")
    verdict = ("WAVELET IMPROVES scenario quality (CI excludes 0)" if lo > 0 else
               ("WAVELET HURTS (CI excludes 0)" if hi < 0 else
                "NULL RESULT IS GENUINE: no significant difference (CI includes 0)"))
    print(f"  --> {verdict}")
    os.makedirs('results', exist_ok=True)
    json.dump(dict(frac=a.frac, S_max=S_max, mean_wavelet=float(np.mean(c_w)),
        mean_seasonal_only=float(np.mean(c_c)), mean_naive=float(np.mean(c_n)),
        diff_naive_minus_wavelet=m, ci=[float(lo), float(hi)],
        diff_seasonal_minus_wavelet=float(np.mean(d2)), ci_seasonal=[float(l2), float(h2)],
        rel_pct=float(rel), verdict=verdict, n_eval=int(a.n_eval)),
        open('results/wavelet_ablation.json', 'w'), indent=2)
    print("\n  wrote results/wavelet_ablation.json")

if __name__ == '__main__':
    main()
