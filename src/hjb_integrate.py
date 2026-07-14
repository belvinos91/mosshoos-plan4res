#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hjb_integrate.py
================
Runs the WP4-B/C continuous-time strand (Ito SDE calibration + HJB PDE solver
+ wavelet-Galerkin convergence) from the SAME inputs the SDDP comparison uses,
and returns a compact results dict for the dashboard JSON.

This wires mosshoos_wp4_ito_hjb.py into the integrated pipeline so the
continuous-time MOSSHOOS strand is reproducible from Plan4RES or OMVS inputs,
alongside the discrete-time SDDP comparison.

Interface
---------
run_hjb(seasonal, residuals, params, fast=True) -> dict
    seasonal, residuals : from WP2 Haar wavelet (mosshoos_wp2_wavelet)
    params              : dict with S_max, R_max, D, C_def, nu (reservoir setup)
    returns             : dict of calibration + convergence results (JSON-safe)
"""
import numpy as np


def run_hjb(seasonal, residuals, params, fast=True, verbose=True):
    """
    Calibrate the Ito SDE to the wavelet output and solve the HJB PDE.

    Returns a JSON-serializable dict capturing:
      - SDE calibration (kappa, sigma0, Feller, KS p, Fourier R2)
      - EM / Milstein strong-convergence slopes (theory 0.5 / 1.0)
      - wavelet-Galerkin L2 convergence rate (theory O(M^-1))
    """
    from mosshoos_wp4_ito_hjb import ItoSDEModel, HJBWaveletSolver

    seasonal = np.asarray(seasonal, float)
    residuals = np.asarray(residuals, float)

    # grid resolution: modest in fast mode so it runs quickly in the pipeline
    if fast:
        N_S, N_Q, N_t, J_wg, n_paths = 25, 10, 12, 3, 120
    else:
        N_S, N_Q, N_t, J_wg, n_paths = 60, 20, 52, 4, 400

    # ---- 1. Calibrate the Ito SDE from the wavelet seasonal + residuals ----
    ito = ItoSDEModel(sigma_type='sqrt', dt=1.0/52.0, K_fourier=5)
    ito.calibrate(seasonal, residuals, dt_obs=1.0/12.0)

    # ---- 2. Strong convergence test (EM vs Milstein) ----
    Q0 = float(ito.theta(np.array([0.0]))[0])
    conv = ito.convergence_test(Q0, n_paths=n_paths)

    # ---- 3. HJB PDE solver: FLAGGED as under validation ----
    # The explicit finite-difference HJB solver currently exhibits a numerical
    # instability (value function diverges) and the wavelet-Galerkin convergence
    # depends on it, so we DO NOT ship those numbers. Only the validated
    # continuous-time results (SDE calibration + Euler-Maruyama convergence)
    # are reported here. See RESULTS.md open items.
    out = dict(
        sde=dict(
            sigma_type=ito.sigma_type,
            kappa=round(float(ito.kappa_), 4),
            sigma0=round(float(ito.sigma0_), 4),
            fourier_R2=round(float(ito.cal_stats_['R2']), 4),
            ks_p=round(float(ito.cal_stats_['ks_p']), 4),
            feller=round(float(ito.cal_stats_['feller']), 3),
            feller_ok=bool(ito.cal_stats_['feller'] >= 1.0),
        ),
        convergence=dict(
            em_slope=round(float(conv['slope_em']), 3),
            em_theory=0.5,
            milstein_slope=round(float(conv['slope_mil']), 3),
            milstein_theory=1.0,
            milstein_status="under_validation",
        ),
        hjb_pde=dict(status="under_validation",
                     note="explicit FD solver numerically unstable; "
                          "wavelet-Galerkin convergence pending stable solver"),
        grid=dict(N_S=N_S, N_Q=N_Q, N_t=N_t, n_paths=n_paths, fast=fast),
    )
    if verbose:
        print("\n[HJB-INTEGRATE] summary:")
        print(f"  kappa={out['sde']['kappa']}, sigma0={out['sde']['sigma0']}, "
              f"Feller={out['sde']['feller']} ({'ok' if out['sde']['feller_ok'] else 'may hit 0'})")
        print(f"  EM slope={out['convergence']['em_slope']} (theory 0.5), "
              f"Milstein slope={out['convergence']['milstein_slope']} (theory 1.0)")
        print("  HJB PDE solver: under validation (see RESULTS.md)")
    return out


if __name__ == '__main__':
    # smoke test with synthetic wavelet-like inputs
    import argparse, sys, os, json
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    rng = np.random.default_rng(0)
    N = 360
    t = np.arange(N) / 12.0
    seasonal = 100 + 80*np.sin(2*np.pi*t - 1.3)      # seasonal-ish
    seasonal = np.maximum(seasonal, 1.0)
    residuals = rng.normal(0, 40, N)
    params = dict(S_max=300.0, R_max=146.0, D=92.0, C_def=300.0, nu=0.15)
    res = run_hjb(seasonal, residuals, params, fast=True)
    print("\nJSON:")
    print(json.dumps(res, indent=2))
