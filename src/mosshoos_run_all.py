#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mosshoos_run_all.py
====================
Master runner — executes WP1 → WP2 → WP3 → WP4-A → WP4-B/C in sequence.
Produces all CSVs and figures. Prints a final alignment summary.

Usage
-----
# With synthetic data (default, no files needed):
    python mosshoos_run_all.py

# With real OMVS data:
    python mosshoos_run_all.py --inflow OMVS_inflows.csv

# Fast test (reduced resolution):
    python mosshoos_run_all.py --fast

# Full production run:
    python mosshoos_run_all.py --inflow OMVS_inflows.csv --n_iter 200 --n_mc 500
"""

import argparse, os, time
import numpy as np
import pandas as pd

OUT = 'results'  


def print_header(title):
    w = 62
    print('\n' + '='*w)
    print(f'  {title}')
    print('='*w)


def print_section(wp, title):
    print(f'\n{"─"*62}')
    print(f'  {wp} — {title}')
    print(f'{"─"*62}')


if __name__ == '__main__':
    t_start = time.time()

    parser = argparse.ArgumentParser(
        description='MOSSHOOS master runner — all work packages')
    parser.add_argument('--inflow',     default='',
                        help='Path to OMVS monthly inflow CSV')
    parser.add_argument('--output',     default=OUT)
    parser.add_argument('--n_stages',   type=int, default=12,
                        help='12=monthly  52=weekly')
    parser.add_argument('--n_clusters', type=int, default=10)
    parser.add_argument('--n_bootstrap',type=int, default=200)
    parser.add_argument('--n_iter',     type=int, default=100,
                        help='SDDP iterations (200 for production)')
    parser.add_argument('--n_mc',       type=int, default=200,
                        help='Monte Carlo runs (500 for production)')
    parser.add_argument('--n_paths',    type=int, default=300,
                        help='Itô SDE paths (1000 for production)')
    parser.add_argument('--N_S',        type=int, default=40,
                        help='HJB storage grid (80 for production)')
    parser.add_argument('--N_Q',        type=int, default=15,
                        help='HJB inflow grid (30 for production)')
    parser.add_argument('--N_t',        type=int, default=24,
                        help='HJB time steps (52 for production)')
    parser.add_argument('--J_wg',       type=int, default=3,
                        help='Wavelet-Galerkin levels (4 for production)')
    parser.add_argument('--fast',       action='store_true',
                        help='Fast test: reduced resolution')
    parser.add_argument('--skip_hjb',   action='store_true',
                        help='Skip HJB solver (saves time)')
    args = parser.parse_args()

    # fast mode overrides
    if args.fast:
        args.n_iter  = 30
        args.n_mc    = 50
        args.n_paths = 100
        args.N_S     = 25
        args.N_Q     = 10
        args.N_t     = 12
        args.J_wg    = 2
        args.n_bootstrap = 50
        print('[INFO] Fast mode: reduced resolution for testing')

    os.makedirs(args.output, exist_ok=True)

    print_header('MOSSHOOS — Full Pipeline')
    print(f'  Inflow CSV  : {args.inflow or "synthetic (West African)"}')
    print(f'  Stages      : {args.n_stages}')
    print(f'  SDDP iter   : {args.n_iter}')
    print(f'  MC runs     : {args.n_mc}')
    print(f'  Itô paths   : {args.n_paths}')
    print(f'  Output dir  : {args.output}')

    # ══════════════════════════════════════════════════════
    # WP1 — DATA
    # ══════════════════════════════════════════════════════
    print_section('WP1', 'Data preparation')
    from mosshoos_wp1_data import load_inflows, validate_inflows, plot_inflows
    inflows = load_inflows(args.inflow)
    stats   = validate_inflows(inflows)
    plot_inflows(inflows, stats,
                 save_path=os.path.join(args.output, 'WP1_inflows.png'))
    x = inflows.values.astype(float)

    # ══════════════════════════════════════════════════════
    # WP2 — WAVELET
    # ══════════════════════════════════════════════════════
    print_section('WP2', 'Haar DWT decomposition')
    from mosshoos_wp2_wavelet import (HaarWavelet, fourier_fit_seasonal,
                                       plot_decomposition, acf)
    J         = HaarWavelet.optimal_levels(x)
    wavelet   = HaarWavelet(J).fit(x)
    wavelet.print_summary()
    seasonal  = wavelet.seasonal(len(x))
    residuals = wavelet.residuals(x)
    fourier   = fourier_fit_seasonal(seasonal, dt=1.0/12.0, K=5)

    # save decomposition
    dec_df = pd.DataFrame({
        'inflow': x, 'seasonal': seasonal, 'residual': residuals,
        'fourier': fourier['theta'](np.arange(len(x))/12.0),
    }, index=inflows.index)
    dec_df.to_csv(os.path.join(args.output, 'WP2_decomposition.csv'),
                  float_format='%.4f')

    plot_decomposition(inflows, wavelet, seasonal, residuals, fourier,
                       save_path=os.path.join(args.output, 'WP2_wavelet.png'))

    # ══════════════════════════════════════════════════════
    # WP3 — SCENARIOS
    # ══════════════════════════════════════════════════════
    print_section('WP3', 'Block-bootstrap scenario generation')
    from mosshoos_wp3_scenarios import (ScenarioGenerator,
                                          validate_scenarios, plot_scenarios)
    gen = ScenarioGenerator(
        n_bootstrap  = args.n_bootstrap,
        n_clusters   = args.n_clusters,
        block_size   = 4,
    )
    eta, weights, types = gen.fit(seasonal, residuals,
                                   n_stages=args.n_stages)
    val = validate_scenarios(gen, x, n_stages=args.n_stages)

    # save
    T = args.n_stages
    eta_df = pd.DataFrame(
        eta,
        index=[f'scen_{k+1}' for k in range(len(eta))],
        columns=[f's{t+1}' for t in range(T)])
    eta_df['weight'] = weights
    eta_df['type']   = types
    eta_df.to_csv(os.path.join(args.output, 'WP3_scenarios.csv'),
                  float_format='%.4f')

    val_s = {k: v for k, v in val.items() if not isinstance(v, np.ndarray)}
    pd.Series(val_s).to_csv(os.path.join(args.output, 'WP3_validation.csv'))

    plot_scenarios(gen, x, val, n_stages=args.n_stages,
                   save_path=os.path.join(args.output, 'WP3_scenarios.png'))

    # ══════════════════════════════════════════════════════
    # WP4-A — SDDP
    # ══════════════════════════════════════════════════════
    print_section('WP4-A', 'SSV-SDDP (M-S vs M-D comparison)')
    from sddp_clean import compare_MS_MD          # validated, corrected SDDP core

    # Real Manantali regime (matches storage_sweep.py / RESULTS.md) so the
    # comparison is on the reservoir actually studied, not solver defaults.
    mean_mo = float(np.mean(x))
    S_max   = 3.0 * mean_mo                        # ~3 months usable storage
    R_max   = 200.0 * 730 / 1000.0                 # 200 MW  -> ~146 GWh/month
    D       = 0.9 * mean_mo                         # demand near mean inflow

    results = compare_MS_MD(
        eta, weights, seasonal,
        n_stages = args.n_stages,
        n_iter   = args.n_iter,
        n_mc     = args.n_mc,
        S0       = 0.5 * S_max,
        S_max    = S_max, R_max = R_max, D = D,
    )

    for name, (solver, df) in results.items():
        df.to_csv(os.path.join(args.output, f'WP4_MC_{name}.csv'),
                  index=False, float_format='%.4f')

    try:
        from mosshoos_wp4_sddp import plot_results
        plot_results(results, eta, weights, n_stages=args.n_stages,
                     save_path=os.path.join(args.output, 'WP4_SDDP_results.png'))
    except Exception as _e:
        print(f'  [WP4-A plot skipped: {_e}]')

    # ══════════════════════════════════════════════════════
    # WP4-B/C — ITÔ SDE + HJB + WAVELET–GALERKIN
    # ══════════════════════════════════════════════════════
    if not args.skip_hjb:
        print_section('WP4-B/C', 'Itô SDE + HJB + Wavelet–Galerkin')
        from mosshoos_wp4_ito_hjb import (ItoSDEModel, HJBWaveletSolver,
                                           plot_ito_hjb)

        ito = ItoSDEModel(sigma_type='sqrt', dt=1.0/52.0)
        ito.calibrate(seasonal, residuals, dt_obs=1.0/12.0)

        Q0    = float(ito.theta(np.array([0.0]))[0])
        Q_em  = ito.simulate_em(Q0, args.N_t, args.n_paths)
        Q_mil = ito.simulate_milstein(Q0, args.N_t, args.n_paths)
        # Path-coupled strong-convergence test (validated). The old inline test
        # re-seeded per discretisation, so both slopes capped near 0.3 and hid
        # Milstein's first order; this drives every level with one Brownian path.
        from sde_convergence import strong_convergence
        conv  = strong_convergence(Q0, ito.kappa_, ito.sigma0_, ito.theta,
                                   kind='sqrt', n_paths=min(args.n_paths, 500),
                                   verbose=False)

        Q_max = float(np.percentile(x, 97))
        # Validated wavelet-Galerkin convergence. The legacy in-class WG solver
        # returned a constant L2 error (spurious O(M^0.00)); this uses the
        # corrected FD reference and recovers the H^1-optimal O(M^-1) rate.
        from wavelet_galerkin import wg_convergence
        # hjb_clean expects theta_fn(t)->scalar; ito.theta is vectorised.
        theta_scalar = lambda t: float(np.ravel(ito.theta(np.atleast_1d(float(t))))[0])
        wg_params = dict(S_max=300.0, Q_max=Q_max,
                         R_max=float(np.percentile(x, 95)),
                         D=0.9 * float(np.mean(x)))
        # N_S>=64 is required for the O(M^-1) regime to be visible; a coarser
        # grid leaves the WG error dominated by FD-reference error and the fitted
        # slope is meaningless. This diagnostic is cheap, so fix the resolution.
        wg_rows, wg_slope = wg_convergence(
            theta_scalar, ito.kappa_, ito.sigma0_, wg_params,
            N_S=64, N_Q=12, N_t=8, verbose=False)
        pd.DataFrame(wg_rows).to_csv(
            os.path.join(args.output, 'WP4_WG_convergence.csv'),
            index=False, float_format='%.6e')

        # HJB surface figure (guarded: a plotting mismatch never aborts the run).
        try:
            conv_arr = {k: (np.asarray(v) if isinstance(v, list) else v)
                        for k, v in conv.items()}
            hjb = HJBWaveletSolver(
                ito=ito, S_max=300.0, Q_max=Q_max,
                R_max=float(np.percentile(x, 95)),
                N_S=args.N_S, N_Q=args.N_Q, N_t=args.N_t, J_wg=args.J_wg)
            hjb.solve_fd()
            hjb.solve_wavelet_galerkin()
            wg_table, _ = hjb.wg_convergence_table()
            plot_ito_hjb(ito, hjb, Q_em, Q_mil, conv_arr, wg_table, wg_slope,
                         save_path=os.path.join(args.output, 'WP4_ito_hjb.png'))
        except Exception as _e:
            print(f'  [WP4-C plot skipped: {_e}]')
    else:
        print('[INFO] HJB solver skipped (--skip_hjb)')

    # ══════════════════════════════════════════════════════
    # FINAL ALIGNMENT SUMMARY
    # ══════════════════════════════════════════════════════
    elapsed = time.time() - t_start
    ms_solver, mc_ms = results['MS']
    _,         mc_md = results['MD']
    c_ms = mc_ms['annual_cost'].mean()
    c_md = mc_md['annual_cost'].mean()
    d_ms = mc_ms['deficit_gwh'].mean()
    d_md = mc_md['deficit_gwh'].mean()

    print_header('FINAL RESULTS SUMMARY')
    print(f'\n  Data source         : {"REAL CSV: "+args.inflow if args.inflow else "SYNTHETIC (replace with OMVS_inflows.csv)"}')
    print(f'\n  WP2 — Wavelet:')
    print(f'    J* (optimal)      : {J}  levels')
    print(f'    A_J energy        : {wavelet.energy_[0]:.1f}%  '
          f'(target ≥75%)')
    print(f'    Residual std      : {residuals.std():.2f} GWh/stage')
    print(f'    Fourier R²        : {fourier["R2"]:.4f}')
    print(f'\n  WP3 — Scenarios:')
    print(f'    |η|               : {len(eta)}  scenarios')
    print(f'    KS stat / p       : {val["ks_stat"]:.4f} / {val["ks_p"]:.4f}')
    print(f'    ACF max error     : {val["acf_max_err"]:.4f}')
    print(f'    Energy score impr.: {val["es_improvement"]:.1f}%  vs i.i.d.')
    print(f'\n  WP4-A — SDDP:')
    print(f'    M-S  C^opex       : {c_ms:.1f}  $/cycle')
    print(f'    M-D  C^opex       : {c_md:.1f}  $/cycle')
    print(f'    Saving M-S vs M-D : {(c_md-c_ms)/c_md*100:.1f}%')
    print(f'    M-S  E[deficit]   : {d_ms:.1f} GWh/yr')
    print(f'    M-D  E[deficit]   : {d_md:.1f} GWh/yr')
    print(f'    Deficit reduction : {(d_md-d_ms)/max(d_md,1)*100:.1f}%')
    print(f'    SDDP converged    : {ms_solver.converged_}')
    if not args.skip_hjb:
        print(f'\n  WP4-B — Itô SDE:')
        print(f'    κ               : {ito.kappa_:.4f} /yr')
        print(f'    σ₀              : {ito.sigma0_:.4f}')
        print(f'    KS p-value      : {ito.cal_stats_["ks_p"]:.4f}')
        print(f'    EM slope        : {conv["slope_em"]:.2f}  (theory 0.50)')
        print(f'    Milstein slope  : {conv["slope_mil"]:.2f}  (theory 1.00)')
        print(f'\n  WP4-C — HJB Wavelet–Galerkin:')
        print(f'    W-G slope       : {wg_slope:.2f}  (theory -1.00)')

    print(f'\n  Output files in     : {args.output}')
    print(f'  Total runtime       : {elapsed:.1f}s')
    print('\n' + '='*62)

    # save summary CSV
    summary = {
        'data_source'       : args.inflow or 'synthetic',
        'wavelet_J'         : J,
        'wavelet_AJ_pct'    : wavelet.energy_[0],
        'fourier_R2'        : fourier['R2'],
        'n_scenarios'       : len(eta),
        'ks_stat'           : val['ks_stat'],
        'ks_p'              : val['ks_p'],
        'acf_max_err'       : val['acf_max_err'],
        'es_improvement_pct': val['es_improvement'],
        'C_opex_MS'         : c_ms,
        'C_opex_MD'         : c_md,
        'cost_saving_pct'   : (c_md - c_ms) / c_md * 100,
        'deficit_MS_gwh'    : d_ms,
        'deficit_MD_gwh'    : d_md,
        'deficit_reduction_pct': (d_md - d_ms) / max(d_md, 1) * 100,
        'sddp_converged'    : ms_solver.converged_,
        'runtime_s'         : elapsed,
    }
    pd.Series(summary).to_csv(
        os.path.join(args.output, 'MOSSHOOS_summary.csv'), header=False)
    print(f'  Summary saved → MOSSHOOS_summary.csv')
