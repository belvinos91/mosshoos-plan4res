#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sde_convergence.py
==================
Correct strong-convergence test for the CIR/OU inflow SDE.

WHY THE ORIGINAL TEST GAVE THE WRONG ANSWER
-------------------------------------------
The original convergence_test re-seeded the RNG before each discretization:

    np.random.seed(0)
    Q_em = simulate_em(Q0, n_test, n_paths)     # n_test varies per dt

Because n_test differs per dt, np.random.normal draws a DIFFERENT Brownian
path for each discretization. Strong convergence is defined pathwise:

    e(dt) = E| Q_T^{dt} - Q_T^{ref} |     on the SAME Brownian path

Without path coupling, the measured error is dominated by the difference
between independent paths, not by discretization error. Both EM and Milstein
then show a spurious slope near 0.3-0.5 (the symptom we observed: EM 0.31,
Milstein 0.32 -- nearly identical, which is the tell).

THE FIX
-------
Generate the Brownian increments ONCE on the finest grid, then build every
coarser discretization by SUMMING the fine increments (Brownian consistency:
the sum of independent N(0, dt_fine) over m sub-steps is N(0, m*dt_fine)).
All schemes then see the same driving path and the true strong order appears:

    EM        : O(dt^0.5)
    Milstein  : O(dt^1.0)

Note: the max(Q, eps) positivity clip also degrades the strong order. We use a
Feller-safe parameter regime (2*kappa*theta >= sigma0^2) and a small floor so
the clip is essentially never active.
"""
from __future__ import annotations
import numpy as np


def _sigma(Q, sigma0, kind='sqrt'):
    Q = np.maximum(Q, 0.0)
    if kind == 'sqrt':
        return sigma0 * np.sqrt(Q)
    if kind == 'linear':
        return sigma0 * Q
    return sigma0 * np.ones_like(Q)


def _dsigma(Q, sigma0, kind='sqrt'):
    Q = np.maximum(Q, 1e-12)
    if kind == 'sqrt':
        return 0.5 * sigma0 / np.sqrt(Q)
    if kind == 'linear':
        return sigma0 * np.ones_like(Q)
    return np.zeros_like(Q)


def _step_paths(Q0, dW, dt, kappa, theta_fn, sigma0, kind, scheme):
    """
    Integrate the SDE over the given Brownian increments.

    dW : (n_paths, n_steps) increments, each ~ N(0, dt)
    scheme : 'em' or 'milstein'
    Returns Q_T (n_paths,)
    """
    n_paths, n_steps = dW.shape
    Q = np.full(n_paths, float(Q0))
    for i in range(n_steps):
        t = i * dt
        th = theta_fn(t)
        s = _sigma(Q, sigma0, kind)
        inc = kappa * (th - Q) * dt + s * dW[:, i]
        if scheme == 'milstein':
            inc = inc + 0.5 * s * _dsigma(Q, sigma0, kind) * (dW[:, i] ** 2 - dt)
        Q = np.maximum(Q + inc, 0.0)      # absorbing floor at 0 (Feller-safe regime)
    return Q


def strong_convergence(Q0, kappa, sigma0, theta_fn, kind='sqrt',
                       T=1.0, n_paths=2000, n_ref=2**12, levels=(2**4, 2**5,
                       2**6, 2**7, 2**8), seed=0, verbose=True):
    """
    Path-coupled strong convergence test.

    Generates Brownian increments on the finest grid (n_ref steps), then builds
    each coarser level by summing fine increments -> identical driving path.

    Returns dict with dt array, EM/Milstein errors and log-log slopes.
    """
    rng = np.random.default_rng(seed)
    dt_ref = T / n_ref
    dW_fine = rng.normal(0.0, np.sqrt(dt_ref), size=(n_paths, n_ref))

    # reference solution: Milstein on the finest grid (highest order available)
    Q_ref = _step_paths(Q0, dW_fine, dt_ref, kappa, theta_fn, sigma0, kind,
                        'milstein')

    dts, err_em, err_mil = [], [], []
    for n in levels:
        assert n_ref % n == 0, f'n_ref must be divisible by {n}'
        m = n_ref // n                      # fine steps per coarse step
        dt = T / n
        # Brownian consistency: coarse increment = sum of m fine increments
        dW_c = dW_fine.reshape(n_paths, n, m).sum(axis=2)

        Q_em = _step_paths(Q0, dW_c, dt, kappa, theta_fn, sigma0, kind, 'em')
        Q_ml = _step_paths(Q0, dW_c, dt, kappa, theta_fn, sigma0, kind,
                           'milstein')

        dts.append(dt)
        err_em.append(float(np.mean(np.abs(Q_em - Q_ref))))
        err_mil.append(float(np.mean(np.abs(Q_ml - Q_ref))))

    dts = np.array(dts)
    slope_em = float(np.polyfit(np.log(dts), np.log(err_em), 1)[0])
    slope_mil = float(np.polyfit(np.log(dts), np.log(err_mil), 1)[0])

    if verbose:
        print('[SDE] Path-coupled strong convergence '
              f'(n_paths={n_paths}, ref={n_ref} steps)')
        print(f"  {'dt':>10} {'EM err':>12} {'Milstein err':>14}")
        for dt, a, b in zip(dts, err_em, err_mil):
            print(f'  {dt:>10.5f} {a:>12.4e} {b:>14.4e}')
        print(f'  EM slope       = {slope_em:.3f}  (theory 0.5)')
        print(f'  Milstein slope = {slope_mil:.3f}  (theory 1.0)')

    return dict(dt=dts.tolist(), err_em=err_em, err_mil=err_mil,
                slope_em=round(slope_em, 3), slope_mil=round(slope_mil, 3),
                n_paths=n_paths, n_ref=n_ref)


def _tests():
    print('=' * 60)
    print('SDE strong-convergence unit tests (path-coupled)')
    print('=' * 60)
    # Feller-safe regime: 2*kappa*theta >= sigma0^2 so Q stays positive and
    # the positivity floor never activates (which would destroy strong order).
    kappa, sigma0 = 2.0, 0.6
    theta_fn = lambda t: 1.0 + 0.3 * np.sin(2 * np.pi * t)
    feller = 2 * kappa * 1.0 / sigma0 ** 2
    print(f'Feller ratio = {feller:.2f} (>=1 so Q_t > 0)')

    r = strong_convergence(Q0=1.0, kappa=kappa, sigma0=sigma0,
                           theta_fn=theta_fn, kind='sqrt',
                           n_paths=2000, n_ref=2**12,
                           levels=(2**4, 2**5, 2**6, 2**7, 2**8), verbose=True)

    em_ok = 0.35 <= r['slope_em'] <= 0.70
    ml_ok = r['slope_mil'] >= 0.80
    better = r['slope_mil'] > r['slope_em'] + 0.2
    print()
    print(f"TEST 1 EM slope ~ 0.5        : {'PASS' if em_ok else 'FAIL'} "
          f"({r['slope_em']})")
    print(f"TEST 2 Milstein slope ~ 1.0  : {'PASS' if ml_ok else 'FAIL'} "
          f"({r['slope_mil']})")
    print(f"TEST 3 Milstein beats EM     : {'PASS' if better else 'FAIL'}")
    ok = em_ok and ml_ok and better
    print('=' * 60)
    print('ALL TESTS PASSED' if ok else 'SOME TESTS FAILED')
    print('=' * 60)
    return ok


if __name__ == '__main__':
    _tests()
