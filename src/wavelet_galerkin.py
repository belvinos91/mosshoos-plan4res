#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wavelet_galerkin.py
===================
Haar wavelet-Galerkin reduced-basis solver for the HJB value function, and its
convergence study against the corrected finite-difference reference
(hjb_clean.HJB).

WHY A WAVELET BASIS
-------------------
The HJB value function V(.,Q,t) for this reservoir problem is CONTINUOUS but
has kinks in S (the min over the release u switches between regimes: deficit,
interior, spill). So V lies in H^1 but not H^2. For such functions a Haar
(piecewise-constant) multiresolution basis gives

    || V - P_M V ||_{L2}  =  O(M^-1)         for V in H^1

where P_M is the orthogonal projector onto the first M Haar basis functions.
A smoother (e.g. spectral) basis cannot beat this because of the kinks, which
is exactly why the wavelet basis is the right choice here.

METHOD
------
1. Build an orthonormal Haar basis on [0, S_max] at level J  ->  M = 2^J.
2. Time-step the HJB backward in the REDUCED space: at each sub-step the
   Hamiltonian H is evaluated in physical space (the min over u is nonlinear
   and cannot be applied in coefficient space), then the update is projected
   back onto the basis:
       c^{n}  =  c^{n+1} + dt * Psi^T H(Psi c^{n+1})
   With an orthonormal Psi this is a Galerkin projection of the semi-discrete
   HJB onto span(Psi).
3. Measure the L2 error against the FD reference at t=0 and fit the rate in M.

The reference here is hjb_clean.HJB, whose two sign/state defects are fixed and
which is unit-tested, so the convergence study measures the BASIS error rather
than inheriting a broken reference.
"""
from __future__ import annotations
import numpy as np

from hjb_clean import HJB


# -----------------------------------------------------------------------------
def haar_basis(N, J):
    """
    Orthonormal Haar basis sampled on N uniform points of [0,1].

    Returns Psi with shape (N, M), M = 2**J, columns orthonormal w.r.t. the
    discrete inner product <a,b> = (1/N) sum a_i b_i.
    Includes the scaling function (constant) plus wavelets at levels 0..J-1.
    """
    x = (np.arange(N) + 0.5) / N
    cols = [np.ones(N)]                       # scaling function
    for j in range(J):
        for k in range(2 ** j):
            lo = k / 2 ** j
            mid = (k + 0.5) / 2 ** j
            hi = (k + 1) / 2 ** j
            psi = np.zeros(N)
            psi[(x >= lo) & (x < mid)] = 1.0
            psi[(x >= mid) & (x < hi)] = -1.0
            psi *= 2 ** (j / 2.0)
            cols.append(psi)
    Psi = np.column_stack(cols[:2 ** J])
    # normalize w.r.t. the discrete inner product
    for m in range(Psi.shape[1]):
        nrm = np.sqrt(np.mean(Psi[:, m] ** 2))
        if nrm > 1e-12:
            Psi[:, m] /= nrm
    return Psi


def project(Psi, v):
    """Orthogonal projection coefficients of v onto span(Psi)."""
    return (Psi.T @ v) / Psi.shape[0]


def reconstruct(Psi, c):
    return Psi @ c


# -----------------------------------------------------------------------------
class WaveletGalerkinHJB(HJB):
    """
    HJB solved in a reduced Haar basis in the storage dimension.
    Inherits the (corrected) dynamics/cost from HJB.
    """

    def __init__(self, *args, J_wg=3, **kwargs):
        super().__init__(*args, **kwargs)
        self.J_wg = int(J_wg)
        self.M = 2 ** self.J_wg

    def solve_wg(self):
        N_S, N_Q, N_t = self.N_S, self.N_Q, self.N_t
        dS, dQ = self.dS, self.dQ
        Psi = haar_basis(N_S, self.J_wg)

        # CFL (same as the FD solver so the comparison is like-for-like)
        max_diff = 0.5 * float(np.max(self.sigma(self.Q_grid) ** 2))
        denom = ((self.Q_max + self.R_max) / dS
                 + self.kappa * self.Q_max / dQ + 2.0 * max_diff / dQ ** 2)
        dt_stable = 0.9 / max(denom, 1e-12)
        n_sub = max(1, int(np.ceil((1.0 / N_t) / dt_stable)))
        dt = (1.0 / N_t) / n_sub

        V = np.zeros((N_S, N_Q, N_t + 1))
        for j in range(N_Q):
            V[:, j, N_t] = -self.nu * self.S_grid
        # project the terminal condition into the basis
        for j in range(N_Q):
            V[:, j, N_t] = reconstruct(Psi, project(Psi, V[:, j, N_t]))

        u_full = np.linspace(0.0, self.R_max, self.n_u)

        for n in range(N_t - 1, -1, -1):
            theta_t = float(self.theta_fn(n / N_t))
            Vs = V[:, :, n + 1].copy()
            for _ in range(n_sub):
                Vp = Vs.copy()
                H = np.zeros_like(Vp)
                for i in range(N_S):
                    S_i = self.S_grid[i]
                    for j in range(N_Q):
                        Q_j = self.Q_grid[j]
                        u_hi = min(S_i + Q_j, self.R_max)
                        u = u_full[u_full <= u_hi + 1e-12]
                        if len(u) == 0:
                            u = np.array([0.0])
                        spill = np.maximum(0.0, S_i + Q_j - u - self.S_max)
                        driftS = Q_j - u - spill
                        cost = (self.C_def * np.maximum(0.0, self.D - u)
                                + self.C_sp * spill + self.C_hyd * u)
                        if i == 0:
                            V_S_f = (Vp[1, j] - Vp[0, j]) / dS; V_S_b = V_S_f
                        elif i == N_S - 1:
                            V_S_b = (Vp[N_S-1, j] - Vp[N_S-2, j]) / dS; V_S_f = V_S_b
                        else:
                            V_S_f = (Vp[i+1, j] - Vp[i, j]) / dS
                            V_S_b = (Vp[i, j] - Vp[i-1, j]) / dS
                        V_S = np.where(driftS >= 0, V_S_f, V_S_b)

                        driftQ = self.kappa * (theta_t - Q_j)
                        if j == 0:
                            V_Q = (Vp[i, 1] - Vp[i, 0]) / dQ
                        elif j == N_Q - 1:
                            V_Q = (Vp[i, N_Q-1] - Vp[i, N_Q-2]) / dQ
                        else:
                            V_Q = ((Vp[i, j+1] - Vp[i, j]) / dQ if driftQ >= 0
                                   else (Vp[i, j] - Vp[i, j-1]) / dQ)
                        if 0 < j < N_Q - 1:
                            V_QQ = (Vp[i, j+1] - 2*Vp[i, j] + Vp[i, j-1]) / dQ**2
                        else:
                            V_QQ = 0.0
                        diff = 0.5 * float(self.sigma(np.array([Q_j]))[0]) ** 2
                        H[i, j] = np.min(cost + driftS * V_S
                                         + driftQ * V_Q + diff * V_QQ)
                # Galerkin step: update, then project back onto the Haar basis
                Vn = Vp + dt * H
                for j in range(N_Q):
                    Vn[:, j] = reconstruct(Psi, project(Psi, Vn[:, j]))
                Vs = Vn
            V[:, :, n] = Vs

        self.V_wg = V
        return V


# -----------------------------------------------------------------------------
def wg_convergence(theta_fn, kappa, sigma0, params, J_list=(1, 2, 3, 4, 5),
                   N_S=64, N_Q=12, N_t=8, verbose=True):
    """
    L2 error of the wavelet-Galerkin solution vs the corrected FD reference,
    as a function of M = 2^J basis functions. Returns (rows, slope).
    """
    common = dict(kappa=kappa, sigma0=sigma0, theta_fn=theta_fn,
                  S_max=params.get('S_max', 300.0),
                  Q_max=params.get('Q_max', 250.0),
                  R_max=params.get('R_max', 146.0),
                  D=params.get('D', 92.0),
                  C_def=params.get('C_def', 300.0),
                  nu=params.get('nu', 0.15),
                  N_S=N_S, N_Q=N_Q, N_t=N_t, verbose=False)

    ref = HJB(**common)
    V_ref = ref.solve()
    Qm = N_Q // 2
    v_ref = V_ref[:, Qm, 0]

    rows = []
    for J in J_list:
        wg = WaveletGalerkinHJB(**common, J_wg=J)
        V_wg = wg.solve_wg()
        v_wg = V_wg[:, Qm, 0]
        M = 2 ** J
        err = float(np.sqrt(np.mean((v_ref - v_wg) ** 2)))
        rows.append(dict(J=J, M=M, L2_error=err))

    Ms = np.array([r['M'] for r in rows], float)
    errs = np.array([r['L2_error'] for r in rows], float)
    good = errs > 1e-12
    slope = float(np.polyfit(np.log(Ms[good]), np.log(errs[good]), 1)[0]) \
        if good.sum() >= 2 else float('nan')

    if verbose:
        print('[WG] Haar wavelet-Galerkin convergence vs corrected FD reference')
        print(f"  {'J':>3} {'M':>5} {'L2 error':>14}")
        for r in rows:
            print(f"  {r['J']:>3} {r['M']:>5} {r['L2_error']:>14.4e}")
        print(f'  Empirical rate: O(M^{slope:.2f})   (theory O(M^-1) for H^1)')
    return rows, slope


def _tests():
    print('=' * 60)
    print('Wavelet-Galerkin unit tests')
    print('=' * 60)
    theta = lambda t: 100.0 + 60.0 * np.sin(2 * np.pi * t - 1.2)
    params = dict(S_max=300.0, Q_max=250.0, R_max=146.0, D=92.0,
                  C_def=300.0, nu=0.15)

    # TEST 1: Haar basis is orthonormal
    Psi = haar_basis(64, 3)
    G = (Psi.T @ Psi) / 64
    orth = np.allclose(G, np.eye(Psi.shape[1]), atol=1e-8)
    print(f"TEST 1 Haar basis orthonormal    : {'PASS' if orth else 'FAIL'} "
          f"(max off-diag {np.max(np.abs(G - np.eye(Psi.shape[1]))):.2e})")

    # TEST 2: projection of a smooth function converges as M grows
    x = (np.arange(64) + 0.5) / 64
    f = np.sin(3 * x) + x ** 2
    errs = []
    for J in (1, 2, 3, 4):
        P = haar_basis(64, J)
        errs.append(np.sqrt(np.mean((f - reconstruct(P, project(P, f))) ** 2)))
    mono = all(errs[i] > errs[i + 1] for i in range(len(errs) - 1))
    print(f"TEST 2 projection error decreases: {'PASS' if mono else 'FAIL'} "
          f"({errs[0]:.2e} -> {errs[-1]:.2e})")

    # TEST 3: WG converges to the FD reference as M grows, at ~O(M^-1)
    rows, slope = wg_convergence(theta, kappa=2.0, sigma0=8.0, params=params,
                                 J_list=(1, 2, 3, 4, 5), N_S=64, N_Q=12, N_t=8,
                                 verbose=True)
    decreasing = rows[0]['L2_error'] > rows[-1]['L2_error']
    rate_ok = -1.6 <= slope <= -0.5
    print(f"TEST 3 WG error decreases with M : {'PASS' if decreasing else 'FAIL'}")
    print(f"TEST 4 rate ~ O(M^-1)            : {'PASS' if rate_ok else 'FAIL'} "
          f"(O(M^{slope:.2f}))")

    ok = orth and mono and decreasing and rate_ok
    print('=' * 60)
    print('ALL TESTS PASSED' if ok else 'SOME TESTS FAILED')
    print('=' * 60)
    return ok


if __name__ == '__main__':
    _tests()
