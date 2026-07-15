#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hjb_clean.py
============
Correct, unit-tested Hamilton-Jacobi-Bellman solver for the seasonal reservoir
control problem under a mean-reverting stochastic inflow (CIR/OU).

Fixes two real defects found in the original mosshoos_wp4_ito_hjb.py:

  BUG 1 (sign).  The backward HJB is  -dV/dt = H,  so stepping backward from
  the terminal time gives  V(t) = V(t+dt) + dt*H.  The original used
  `V - dt*H`, which flips the diffusion sign and turns the equation into a
  BACKWARD heat equation (anti-diffusion) -> exponential blow-up. Observed.

  BUG 2 (wrong state for diffusion).  The randomness lives in the INFLOW Q,
  not in storage S: storage evolves deterministically once (u, Q) are known.
  The original applied (1/2)sigma^2 * V_SS. The correct generator is
      kappa*(theta(t) - Q) * V_Q  +  (1/2)*sigma^2(Q) * V_QQ
  i.e. drift and diffusion act on Q, while S carries only the drift
      (Q - u - spill) * V_S.

FORMULATION (cost minimization, consistent with sddp_clean.py)
--------------------------------------------------------------
State   S in [0, S_max]   storage        [GWh]
        Q in [0, Q_max]   inflow rate    [GWh/period]
Control u in [0, min(R_max, S+Q)]  release

    -dV/dt = min_u { C(S,Q,u)
                     + (Q - u - spill) * dV/dS
                     + kappa*(theta(t) - Q) * dV/dQ
                     + 0.5*sigma(Q)^2 * d2V/dQ2 }
    V(S,Q,T) = -nu * S            (terminal water value: stored water is worth
                                   money, so it REDUCES cost -> negative)

Discretization: explicit upwind in S and Q, central second difference in Q,
with automatic CFL sub-stepping for stability.
"""
from __future__ import annotations
import numpy as np


class HJB:
    """Correct HJB solver. Cost-minimization convention."""

    def __init__(self, kappa, sigma0, theta_fn, sigma_type='sqrt',
                 S_max=300.0, Q_max=None, R_max=146.0, D=92.0,
                 C_def=300.0, C_sp=500.0, C_hyd=2.0, nu=0.15,
                 N_S=40, N_Q=20, N_t=12, n_u=40, verbose=True):
        self.kappa = float(kappa)
        self.sigma0 = float(sigma0)
        self.theta_fn = theta_fn          # callable: t (in [0,1]) -> theta
        self.sigma_type = sigma_type
        self.S_max = float(S_max)
        self.Q_max = float(Q_max) if Q_max else 3.0 * float(theta_fn(0.5))
        self.R_max = float(R_max)
        self.D = float(D)
        self.C_def = float(C_def)
        self.C_sp = float(C_sp)
        self.C_hyd = float(C_hyd)
        self.nu = float(nu)
        self.N_S, self.N_Q, self.N_t, self.n_u = N_S, N_Q, N_t, n_u
        self.verbose = verbose

        self.S_grid = np.linspace(0.0, self.S_max, N_S)
        self.Q_grid = np.linspace(0.0, self.Q_max, N_Q)
        self.dS = self.S_grid[1] - self.S_grid[0]
        self.dQ = self.Q_grid[1] - self.Q_grid[0]
        self.V = None
        self.policy = None

    def sigma(self, Q):
        Q = np.maximum(Q, 1e-6)
        if self.sigma_type == 'sqrt':
            return self.sigma0 * np.sqrt(Q)
        if self.sigma_type == 'linear':
            return self.sigma0 * Q
        return self.sigma0 * np.ones_like(np.asarray(Q, float))

    def solve(self):
        N_S, N_Q, N_t = self.N_S, self.N_Q, self.N_t
        dS, dQ = self.dS, self.dQ

        # ---- CFL: explicit scheme needs dt <= 1 / (drift/dS + drift/dQ + 2*diff/dQ^2)
        max_diff = 0.5 * float(np.max(self.sigma(self.Q_grid) ** 2))
        max_driftS = self.Q_max + self.R_max
        max_driftQ = self.kappa * self.Q_max
        denom = max_driftS / dS + max_driftQ / dQ + 2.0 * max_diff / dQ ** 2
        dt_stable = 0.9 / max(denom, 1e-12)
        dt_stage = 1.0 / N_t
        n_sub = max(1, int(np.ceil(dt_stage / dt_stable)))
        dt = dt_stage / n_sub
        if self.verbose:
            print(f'[HJB] grid {N_S}x{N_Q}, {N_t} stages, CFL sub-steps={n_sub}, '
                  f'dt={dt:.2e}')

        V = np.zeros((N_S, N_Q, N_t + 1))
        pol = np.zeros((N_S, N_Q, N_t))

        # terminal condition: stored water reduces future cost -> V = -nu*S
        for j in range(N_Q):
            V[:, j, N_t] = -self.nu * self.S_grid

        u_grid_full = np.linspace(0.0, self.R_max, self.n_u)

        for n in range(N_t - 1, -1, -1):
            t_now = n / N_t
            theta_t = float(self.theta_fn(t_now))
            Vs = V[:, :, n + 1].copy()

            for _ in range(n_sub):
                Vp = Vs.copy()
                Vn = Vp.copy()
                for i in range(N_S):
                    S_i = self.S_grid[i]
                    for j in range(N_Q):
                        Q_j = self.Q_grid[j]
                        # feasible releases
                        u_hi = min(S_i + Q_j, self.R_max)
                        u = u_grid_full[u_grid_full <= u_hi + 1e-12]
                        if len(u) == 0:
                            u = np.array([0.0])

                        spill = np.maximum(0.0, S_i + Q_j - u - self.S_max)
                        driftS = Q_j - u - spill
                        cost = (self.C_def * np.maximum(0.0, self.D - u)
                                + self.C_sp * spill + self.C_hyd * u)

                        # upwind dV/dS
                        if i == 0:
                            V_S_f = (Vp[1, j] - Vp[0, j]) / dS
                            V_S_b = V_S_f
                        elif i == N_S - 1:
                            V_S_b = (Vp[N_S-1, j] - Vp[N_S-2, j]) / dS
                            V_S_f = V_S_b
                        else:
                            V_S_f = (Vp[i+1, j] - Vp[i, j]) / dS
                            V_S_b = (Vp[i, j] - Vp[i-1, j]) / dS
                        V_S = np.where(driftS >= 0, V_S_f, V_S_b)

                        # inflow drift (mean reversion) — upwind in Q
                        driftQ = self.kappa * (theta_t - Q_j)
                        if j == 0:
                            V_Q = (Vp[i, 1] - Vp[i, 0]) / dQ
                        elif j == N_Q - 1:
                            V_Q = (Vp[i, N_Q-1] - Vp[i, N_Q-2]) / dQ
                        else:
                            V_Q = ((Vp[i, j+1] - Vp[i, j]) / dQ if driftQ >= 0
                                   else (Vp[i, j] - Vp[i, j-1]) / dQ)

                        # diffusion in Q (central second difference)
                        if 0 < j < N_Q - 1:
                            V_QQ = (Vp[i, j+1] - 2*Vp[i, j] + Vp[i, j-1]) / dQ**2
                        else:
                            V_QQ = 0.0
                        diff = 0.5 * float(self.sigma(np.array([Q_j]))[0]) ** 2

                        H_vec = (cost + driftS * V_S
                                 + driftQ * V_Q + diff * V_QQ)
                        k = int(np.argmin(H_vec))
                        # BACKWARD step: V(t) = V(t+dt) + dt*H   (correct sign)
                        Vn[i, j] = Vp[i, j] + dt * H_vec[k]
                        pol[i, j, n] = u[k]
                Vs = Vn
            V[:, :, n] = Vs

            if self.verbose and (n % max(1, N_t // 4) == 0):
                print(f'  t={1-n/N_t:.2f}  V[mid,mid]={V[N_S//2, N_Q//2, n]:12.2f}')

        self.V = V
        self.policy = pol
        if self.verbose:
            print('[HJB] solve complete.')
        return V


# ─────────────────────────────────────────────────────────────────────────────
def _tests():
    print('=' * 60)
    print('HJB unit tests')
    print('=' * 60)
    theta = lambda t: 100.0 + 60.0 * np.sin(2 * np.pi * t - 1.2)
    h = HJB(kappa=2.0, sigma0=8.0, theta_fn=theta, S_max=300.0, Q_max=250.0,
            R_max=146.0, D=92.0, N_S=30, N_Q=15, N_t=12, verbose=False)
    V = h.solve()

    # TEST 1: finite (no blow-up)
    finite = np.all(np.isfinite(V)) and np.max(np.abs(V)) < 1e8
    print(f'TEST 1 solution finite (no blow-up): {"PASS" if finite else "FAIL"} '
          f'(max|V|={np.max(np.abs(V)):.3g})')

    # TEST 2: V is U-shaped in storage (decreases from empty as water becomes
    # available, then rises near full where spill is penalized). The physically
    # correct property is that V(S=0) is high (deficit risk) and there is an
    # interior minimum -- NOT a monotone decrease.
    V0 = V[:, h.N_Q // 2, 0]
    argmin_S = int(np.argmin(V0))
    interior_min = 0 < argmin_S < h.N_S - 1
    drop_from_empty = V0[0] > V0[argmin_S]
    print(f'TEST 2 V U-shaped in storage       : '
          f'{"PASS" if (interior_min and drop_from_empty) else "FAIL"} '
          f'(V(0)={V0[0]:.0f}, min={V0[argmin_S]:.0f} at S={h.S_grid[argmin_S]:.0f}, '
          f'V(max)={V0[-1]:.0f})')
    test2 = interior_min and drop_from_empty

    # TEST 3: at LOW storage, more inflow reduces cost (water is scarce there,
    # so extra inflow strictly helps -- no spill risk at S=0).
    Vq_low = V[0, :, 0]
    dec_q = Vq_low[-1] <= Vq_low[0] + 1e-6 * (abs(Vq_low[0]) + 1)
    print(f'TEST 3 more inflow helps when empty: {"PASS" if dec_q else "FAIL"} '
          f'(V(Q=0)={Vq_low[0]:.0f} -> V(Q=max)={Vq_low[-1]:.0f})')

    # TEST 4: policy releases more when storage is high
    pol_low = h.policy[2, h.N_Q // 2, 0]
    pol_high = h.policy[-3, h.N_Q // 2, 0]
    mono_u = pol_high >= pol_low - 1e-9
    print(f'TEST 4 release increases with storage: {"PASS" if mono_u else "FAIL"} '
          f'(u(low S)={pol_low:.1f}, u(high S)={pol_high:.1f})')

    ok = finite and test2 and dec_q and mono_u
    print('=' * 60)
    print('ALL TESTS PASSED' if ok else 'SOME TESTS FAILED')
    print('=' * 60)
    return ok


if __name__ == '__main__':
    _tests()
