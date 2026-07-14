#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
compare_plan4res.py
===================
Runs the full MOSSHOOS pipeline (WP1-WP4) on Plan4RES toyDataset inputs,
compares M-S (MOSSHOOS) vs M-D (Plan4RES deterministic method), and writes
results/comparison_results.json for the dashboard.

This proves the donor/supervisor requirement: SAME Plan4RES inputs ->
comparable MOSSHOOS vs Plan4RES outputs.

Usage
-----
    python3 compare_plan4res.py --dataset ~/plan4res_work/data/toyDataset --reservoir RoGonDor
    python3 compare_plan4res.py --inflow OMVS_inflows.csv     # your real OMVS data
"""
import argparse, json, os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))

from mosshoos_wp2_wavelet import HaarWavelet, fourier_fit_seasonal, acf
from mosshoos_wp3_scenarios import ScenarioGenerator, validate_scenarios
from mosshoos_wp4_sddp import compare_MS_MD


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', default='', help='Plan4RES toyDataset path')
    ap.add_argument('--reservoir', default='RoGonDor')
    ap.add_argument('--inflow', default='', help='Real OMVS CSV (overrides --dataset)')
    ap.add_argument('--output', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results'))
    ap.add_argument('--n_iter', type=int, default=100)
    ap.add_argument('--n_mc', type=int, default=300)
    ap.add_argument('--n_bootstrap', type=int, default=200)
    ap.add_argument('--n_clusters', type=int, default=10)
    args = ap.parse_args()
    os.makedirs(args.output, exist_ok=True)
    t0 = time.time()

    print("="*64)
    print("  MOSSHOOS (M-S) vs Plan4RES (M-D) — same inputs, comparable outputs")
    print("="*64)

    # ── 1. Load data — either Plan4RES inputs or real OMVS ────────────────────
    source = ""
    if args.inflow:
        from mosshoos_wp1_data import load_inflows, validate_inflows
        inflows = load_inflows(args.inflow)
        stats = validate_inflows(inflows)
        params = dict(reservoir='OMVS', S_max=1200.0, S_init=600.0, R_max=200.0,
                      C_def=300.0, annual_inflow=float(inflows.resample('YE').sum().mean()))
        source = f"Real OMVS data: {args.inflow}"
        bellman = None
    else:
        from plan4res_io import (load_plan4res_inflows, load_plan4res_params,
                                 load_plan4res_bellman)
        ds = args.dataset or os.path.expanduser('~/plan4res_work/data/toyDataset')
        inflows = load_plan4res_inflows(ds, args.reservoir)
        params  = load_plan4res_params(ds, args.reservoir)
        try:
            bellman = load_plan4res_bellman(ds)
        except Exception:
            bellman = None
        source = f"Plan4RES toyDataset ({args.reservoir})"

    x = inflows.values.astype(float)

    # ── 2. WP2 wavelet ────────────────────────────────────────────────────────
    J = HaarWavelet.optimal_levels(x)
    wv = HaarWavelet(J).fit(x)
    seasonal = wv.seasonal(len(x))
    residuals = wv.residuals(x)
    fourier = fourier_fit_seasonal(seasonal, dt=1/12., K=5)

    # ── 3. WP3 scenarios ──────────────────────────────────────────────────────
    gen = ScenarioGenerator(n_bootstrap=args.n_bootstrap, n_clusters=args.n_clusters,
                            block_size=4)
    eta, weights, types = gen.fit(seasonal, residuals, n_stages=12)
    val = validate_scenarios(gen, x, n_stages=12)

    # ── 4. WP4 SDDP: M-S vs M-D on identical inputs ───────────────────────────
    results = compare_MS_MD(eta, weights, seasonal, n_stages=12,
                            n_iter=args.n_iter, n_mc=args.n_mc,
                            S0=params['S_init'] if 'S_init' in params else 600.0,
                            C_def=params.get('C_def', 300.0))

    def cost(m): return float(results[m][1]['annual_cost'].mean())
    def defc(m): return float(results[m][1]['deficit_gwh'].mean())
    def p90(m):  return float(np.percentile(results[m][1]['annual_cost'], 90))

    c_ms, c_md, c_md0 = cost('MS'), cost('MD'), cost('MD0')
    VSS  = (c_md - c_ms)/max(c_md,1e-9)*100
    P90sav = (p90('MD')-p90('MS'))/max(p90('MD'),1e-9)*100
    defsav = (defc('MD')-defc('MS'))/max(defc('MD'),1e-9)*100

    # ── 5. Assemble results JSON (dashboard reads THIS) ───────────────────────
    out = dict(
        source=source, reservoir=params.get('reservoir','?'),
        runtime_s=round(time.time()-t0,1),
        params={k:(round(v,3) if isinstance(v,float) else v) for k,v in params.items()},
        wp2=dict(J_optimal=int(J), AJ_energy_pct=round(float(wv.energy_[0]),2),
                 fourier_R2=round(float(fourier['R2']),4),
                 residual_std=round(float(residuals.std()),2),
                 energy_spectrum=[round(float(e),2) for e in wv.energy_]),
        wp3=dict(n_scenarios=int(len(eta)),
                 scenarios=eta.round(3).tolist(),
                 weights=[round(float(w),4) for w in weights],
                 types=list(types),
                 annual_volumes=[round(float(v),1) for v in eta.sum(1)],
                 ks_stat=round(float(val['ks_stat']),4),
                 es_improvement=round(float(val['es_improvement']),2),
                 acf_max_err=round(float(val['acf_max_err']),4)),
        wp4=dict(
            models={m:dict(mean_cost=round(cost(m),1), p90_cost=round(p90(m),1),
                           deficit=round(defc(m),2)) for m in ['MS','MS0','MD','MD0']},
            VSS_pct=round(VSS,2), P90_saving_pct=round(P90sav,2),
            deficit_saving_pct=round(defsav,2)),
        monthly_climatology=[round(float(v),2) for v in
                             pd.Series(x, index=inflows.index).groupby(inflows.index.month).mean().values],
        seasonal_profile=[round(float(v),3) for v in seasonal[:12]],
    )
    if bellman is not None:
        out['plan4res_reference'] = dict(
            n_timesteps=int(bellman.Timestep.nunique()),
            n_cuts=int(len(bellman)),
            sample=bellman.head(6).round(3).to_dict('records'))

    json.dump(out, open(os.path.join(args.output,'comparison_results.json'),'w'), indent=2)
    # flat CSV too
    pd.DataFrame(out['wp4']['models']).T.to_csv(
        os.path.join(args.output,'comparison_results.csv'))

    print("\n"+"="*64)
    print(f"  SOURCE: {source}")
    print(f"  VSS (M-D vs M-S) = {VSS:.1f}%   P90 saving = {P90sav:.1f}%   "
          f"deficit saving = {defsav:.1f}%")
    print(f"  Results -> {args.output}/comparison_results.json")
    print("="*64)

if __name__ == '__main__':
    main()
