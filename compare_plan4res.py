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
from sddp_clean import compare_MS_MD   # corrected, validated SDDP


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
        # Real Manantali-scale parameters (from hydro_parameters_dam.csv).
        # Storage/turbine sized in GWh/month energy terms to match the inflow.
        mean_mo = float(inflows.mean())               # ~102 GWh/month
        annual  = float(inflows.resample('Y').sum().mean())
        # Manantali: p_max=200MW -> ~144 GWh/month turbine cap; active storage
        # ~7.9e9 m3 * 13.77 MWh/day/m3s equiv -> large seasonal buffer.
        params = dict(reservoir='OMVS-Manantali',
                      S_max=3.0*mean_mo,              # ~3 months of mean inflow as usable storage
                      S_init=1.5*mean_mo,
                      R_max=200.0*730/1000.0,         # 200 MW -> ~146 GWh/month
                      C_def=300.0, annual_inflow=annual)
        # Demand set so dry months (inflow near 0) create real deficits that
        # storage management must cover -> this is where stochastic value appears.
        params['demand_override'] = 0.9*mean_mo       # ~92 GWh/month
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

    # Size the reservoir demand to what it can actually deliver.
    # The zone total in ZV_ZoneValues.csv is met by the whole fleet (solar, gas,
    # ...), not by this reservoir alone, so it must not be used as D here.
    # Without this the Plan4RES branch fell back to the solver default D=95
    # while RoGonDor receives only 6.55 GWh/month -> demand 14x inflow -> every
    # policy forced to release everything -> all models identical, VSS void.
    # Same convention on both branches so the comparison is like-for-like.
    if 'demand_override' not in params:
        params['demand_override'] = 0.9 * float(inflows.mean())
    print(f"[SETUP] mean inflow={float(inflows.mean()):.2f} GWh/mo  ->  "
          f"D={params['demand_override']:.2f} GWh/mo  "
          f"(S_max={params.get('S_max', float('nan')):.1f}, "
          f"R_max={params.get('R_max', float('nan')):.1f})")

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
    # Pass the reservoir parameters through to the solver. Without this the
    # SDDP silently runs on its defaults (S_max=1500, R_max=180, D=95), which
    # have nothing to do with the reservoir being studied.
    _sddp_kw = dict(n_iter=args.n_iter, n_mc=args.n_mc,
                    S0=params.get('S_init', 600.0),
                    C_def=params.get('C_def', 300.0))
    if 'S_max' in params:           _sddp_kw['S_max'] = params['S_max']
    if 'R_max' in params:           _sddp_kw['R_max'] = params['R_max']
    if 'demand_override' in params: _sddp_kw['D']     = params['demand_override']
    print('[SETUP] solver kwargs: ' + ', '.join(f'{k}={v}' for k, v in _sddp_kw.items()))
    results = compare_MS_MD(eta, weights, seasonal, n_stages=12, **_sddp_kw)

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
    # WP4-B: validated continuous-time strand (Ito SDE calibration + EM convergence)
    try:
        from hjb_integrate import run_hjb
        hjb_res = run_hjb(seasonal, residuals, params, fast=True, verbose=True)
        out['wp4b_continuous'] = hjb_res
    except Exception as e:
        print(f"[WP4-B] skipped: {e}")

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
