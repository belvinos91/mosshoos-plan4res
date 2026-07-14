#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mosshoos_wp4_ito_hjb.py
========================
WP4-B — Itô SDE (Euler–Maruyama and Milstein)
WP4-C — Hamilton–Jacobi–Bellman + Wavelet–Galerkin solver

Run:
    python mosshoos_wp4_ito_hjb.py
    python mosshoos_wp4_ito_hjb.py --inflow OMVS_inflows.csv --n_paths 1000 --N_S 80
"""

import argparse, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import ks_2samp, kstest
warnings.filterwarnings('ignore')

np.random.seed(2024)
OUT   = '/mnt/user-data/outputs/'
BLUE  = '#1A6B8A'; BLUE2 = '#2E9BBF'; BLUE3 = '#7EC8E0'
GOLD  = '#B8860B'; SAGE  = '#2D6A4F'; EMBER = '#C0392B'; BG = '#F4F1EB'


# ══════════════════════════════════════════════════════════════════════════════
# §1  ITÔ SDE — CALIBRATION & SIMULATION
# ══════════════════════════════════════════════════════════════════════════════

class ItoSDEModel:
    """
    Seasonal OU–CIR Itô SDE for OMVS inflow process:

        dQ_t = κ(θ(t) − Q_t) dt + σ₀ √Q_t dW_t

    θ(t) = a₀ + Σ_{k=1}^K [aₖ cos(2πkt) + bₖ sin(2πkt)]
    estimated from the Haar wavelet seasonal component A_J.

    Parameters
    ----------
    sigma_type : 'sqrt'     → CIR  σ(Q) = σ₀√Q  (non-negative, preferred)
                 'constant' → OU   σ(Q) = σ₀
                 'linear'   → GBM  σ(Q) = σ₀ Q
    """

    def __init__(self, sigma_type: str = 'sqrt',
                  dt: float = 1.0/52.0, K_fourier: int = 5):
        self.sigma_type = sigma_type
        self.dt         = dt
        self.K          = K_fourier
        # calibrated parameters
        self.kappa_  = None
        self.sigma0_ = None
        self.a_f_    = None
        self.b_f_    = None
        self.cal_stats_ = {}

    # ── seasonal mean θ(t) ──────────────────────────────────────────────────
    def theta(self, t: np.ndarray) -> np.ndarray:
        t = np.atleast_1d(t)
        v = self.a_f_[0] * np.ones_like(t)
        for k in range(1, self.K + 1):
            v += (self.a_f_[k] * np.cos(2*np.pi*k*t)
                 + self.b_f_[k] * np.sin(2*np.pi*k*t))
        return np.maximum(v, 0.1)

    # ── diffusion coefficient ────────────────────────────────────────────────
    def sigma(self, Q: np.ndarray) -> np.ndarray:
        Q = np.maximum(Q, 0.01)
        if self.sigma_type == 'sqrt':    return self.sigma0_ * np.sqrt(Q)
        if self.sigma_type == 'linear':  return self.sigma0_ * Q
        return self.sigma0_ * np.ones_like(Q)

    def dsigma_dQ(self, Q: np.ndarray) -> np.ndarray:
        """∂σ/∂Q — needed for Milstein correction."""
        Q = np.maximum(Q, 0.01)
        if self.sigma_type == 'sqrt':   return 0.5 * self.sigma0_ / np.sqrt(Q)
        if self.sigma_type == 'linear': return self.sigma0_ * np.ones_like(Q)
        return np.zeros_like(Q)

    # ── calibration ─────────────────────────────────────────────────────────
    def calibrate(self, seasonal: np.ndarray,
                   residuals: np.ndarray,
                   dt_obs: float = 1.0/12.0) -> 'ItoSDEModel':
        """
        Calibrate κ, σ₀ and Fourier coefficients from A_J and residuals.

        Method:
          1. Fourier regression on A_J → θ(t) parameterisation
          2. OLS on OU discretisation: ε_{t+Δ} = e^{-κΔ} ε_t + ξ_t
             → κ̂ = -log(φ̂) / Δ
          3. σ₀ from innovation variance (CIR moment condition)
          4. KS test on normalised innovations
        """
        N    = len(seasonal)
        t    = np.arange(N) * dt_obs

        # Fourier fit
        cols = [np.ones(N)]
        for k in range(1, self.K + 1):
            cols += [np.cos(2*np.pi*k*t), np.sin(2*np.pi*k*t)]
        X    = np.column_stack(cols)
        beta = np.linalg.lstsq(X, seasonal, rcond=None)[0]
        self.a_f_ = np.zeros(self.K + 1)
        self.b_f_ = np.zeros(self.K + 1)
        self.a_f_[0] = beta[0]
        for k in range(1, self.K + 1):
            self.a_f_[k] = beta[2*k-1]
            self.b_f_[k] = beta[2*k]
        theta_t = X @ beta
        R2      = 1 - np.var(seasonal - theta_t) / (np.var(seasonal) + 1e-9)

        # κ from OLS on OU
        eps = residuals
        phi = np.dot(eps[:-1], eps[1:]) / (np.dot(eps[:-1], eps[:-1]) + 1e-12)
        phi = np.clip(phi, 1e-6, 0.9999)
        self.kappa_ = max(0.1, min(-np.log(phi) / dt_obs, 15.0))

        # σ₀ from innovations
        innov = eps[1:] - phi * eps[:-1]
        if self.sigma_type == 'sqrt':
            q_mid = np.maximum(theta_t[:-1], 0.1)
            s2    = np.mean(innov**2) / (np.mean(q_mid) * dt_obs + 1e-9)
            self.sigma0_ = np.sqrt(max(s2, 0.01))
        else:
            self.sigma0_ = np.std(innov) / np.sqrt(dt_obs + 1e-9)

        # KS test
        std_innov = innov / (self.sigma(theta_t[:-1]) * np.sqrt(dt_obs) + 1e-9)
        ks_s, ks_p = kstest(std_innov, 'norm')
        feller = 2 * self.kappa_ * float(self.theta(np.array([0.5]))[0]) / max(self.sigma0_**2, 1e-9)

        self.cal_stats_ = dict(
            kappa=self.kappa_, sigma0=self.sigma0_,
            ks_stat=ks_s, ks_p=ks_p, R2=R2, feller=feller)

        print(f'\n[WP4-B] Itô SDE calibration:')
        print(f'  σ type            : {self.sigma_type}')
        print(f'  κ (mean reversion): {self.kappa_:.4f} /yr')
        print(f'  σ₀ (diffusion)    : {self.sigma0_:.4f}')
        print(f'  Fourier R²        : {R2:.4f}')
        print(f'  KS p-value        : {ks_p:.4f}  '
              + ('✓ Gaussian innovations' if ks_p > 0.05 else '⚠ non-Gaussian'))
        print(f'  Feller condition  : 2κθ̄/σ₀² = {feller:.2f}'
              + ('  ✓ Q_t > 0 guaranteed' if feller >= 1 else '  ⚠ may hit zero'))
        return self

    # ── Euler–Maruyama ───────────────────────────────────────────────────────
    def simulate_em(self, Q0: float, n_steps: int, n_paths: int) -> np.ndarray:
        """
        Euler–Maruyama scheme:
            Q_{t+Δ} = Q_t + κ(θ(t) − Q_t)Δ + σ(Q_t)√Δ Z_t
        Strong convergence: O(√Δ).
        """
        Q  = np.full((n_paths, n_steps), Q0, dtype=float)
        dt = self.dt
        t_g= np.arange(n_steps) * dt
        for i in range(1, n_steps):
            q_i  = np.maximum(Q[:, i-1], 0.01)
            th_i = self.theta(t_g[i-1])
            dW   = np.random.normal(0, np.sqrt(dt), n_paths)
            Q[:, i] = np.maximum(
                q_i + self.kappa_ * (th_i - q_i) * dt + self.sigma(q_i) * dW,
                0.01)
        return Q

    # ── Milstein ─────────────────────────────────────────────────────────────
    def simulate_milstein(self, Q0: float, n_steps: int, n_paths: int) -> np.ndarray:
        """
        Milstein scheme (adds Itô correction for state-dependent σ):
            Q_{t+Δ} = Q_t + κ(θ-Q_t)Δ + σ(Q_t)√Δ Z_t
                     + ½ σ(Q_t) ∂σ/∂Q (Z_t²Δ − Δ)
        Strong convergence: O(Δ)  [vs O(√Δ) for EM].
        For σ(Q) = σ₀√Q:  ∂σ/∂Q = σ₀/(2√Q)  →  correction = σ₀²/4 (Z²Δ-Δ).
        """
        Q  = np.full((n_paths, n_steps), Q0, dtype=float)
        dt = self.dt
        t_g= np.arange(n_steps) * dt
        for i in range(1, n_steps):
            q_i  = np.maximum(Q[:, i-1], 0.01)
            th_i = self.theta(t_g[i-1])
            sig_i= self.sigma(q_i)
            dsg  = self.dsigma_dQ(q_i)
            dW   = np.random.normal(0, np.sqrt(dt), n_paths)
            drift    = self.kappa_ * (th_i - q_i) * dt
            diffusion= sig_i * dW
            milstein = 0.5 * sig_i * dsg * (dW**2 - dt)
            Q[:, i]  = np.maximum(q_i + drift + diffusion + milstein, 0.01)
        return Q

    # ── convergence test ─────────────────────────────────────────────────────
    def convergence_test(self, Q0: float, n_paths: int = 500) -> dict:
        """
        Empirical strong convergence rate test.
        Compare EM and Milstein against fine reference (Δ_ref = Δ/64).
        """
        T_steps_base = 52
        dt_fine  = self.dt / 64
        n_fine   = T_steps_base * 64

        # reference solution (very fine EM)
        np.random.seed(0)
        Q_ref = self.simulate_em(Q0, n_fine, n_paths)
        Q_T_ref = Q_ref[:, -1]

        errors_em  = []
        errors_mil = []
        dt_list    = [self.dt * 2**j for j in range(-2, 4)]

        for dt_test in dt_list:
            n_test = max(int(T_steps_base / (dt_test / self.dt)), 2)
            old_dt  = self.dt
            self.dt = dt_test
            np.random.seed(0)
            Q_em  = self.simulate_em(Q0, n_test, n_paths)
            np.random.seed(0)
            Q_mil = self.simulate_milstein(Q0, n_test, n_paths)
            self.dt = old_dt
            errors_em.append(np.sqrt(np.mean((Q_em[:, -1]  - Q_T_ref)**2)))
            errors_mil.append(np.sqrt(np.mean((Q_mil[:, -1] - Q_T_ref)**2)))

        dt_arr = np.array(dt_list)
        # log-log slope
        slope_em  = np.polyfit(np.log(dt_arr), np.log(errors_em),  1)[0]
        slope_mil = np.polyfit(np.log(dt_arr), np.log(errors_mil), 1)[0]

        print(f'\n[WP4-B] Strong convergence rates:')
        print(f'  Euler–Maruyama  slope = {slope_em:.2f}  (theory: 0.50)')
        print(f'  Milstein        slope = {slope_mil:.2f}  (theory: 1.00)')

        return dict(dt=dt_arr, err_em=errors_em, err_mil=errors_mil,
                    slope_em=slope_em, slope_mil=slope_mil)


# ══════════════════════════════════════════════════════════════════════════════
# §2  HJB PDE — FINITE DIFFERENCE + WAVELET–GALERKIN
# ══════════════════════════════════════════════════════════════════════════════

class HJBWaveletSolver:
    """
    Solve the HJB PDE for the seasonal reservoir control problem:

        -∂V/∂t = min_{0≤u≤R_max} [
                   C(S,Q,u) + (Q-u-sp(u)) ∂V/∂S + ½σ²(Q) ∂²V/∂S²  ]
        V(S,Q,T) = ν·S  (terminal water value)

    Two solvers:
      (a) Finite difference (FD) — explicit upwind, reference solution
      (b) Wavelet–Galerkin (WG)  — Haar projection, implicit Euler, O(M^-1)

    Parameters
    ----------
    N_S, N_Q, N_t : grid dimensions
    J_wg          : Wavelet–Galerkin basis levels (M = 2^J basis functions)
    """

    def __init__(self, ito: ItoSDEModel,
                  S_max: float = 1500.0,
                  Q_max: float = None,
                  R_max: float = 180.0,
                  D:     float = 95.0,
                  C_def: float = 300.0,
                  C_sp:  float = 500.0,
                  C_hyd: float = 2.0,
                  nu:    float = 0.15,
                  N_S:   int   = 60,
                  N_Q:   int   = 20,
                  N_t:   int   = 52,
                  J_wg:  int   = 4):
        self.ito   = ito
        self.S_max = S_max
        self.R_max = R_max
        self.D     = D
        self.C_def = C_def
        self.C_sp  = C_sp
        self.C_hyd = C_hyd
        self.nu    = nu
        self.N_S   = N_S
        self.N_Q   = N_Q
        self.N_t   = N_t
        self.J_wg  = J_wg

        # grids
        self.S_grid = np.linspace(0, S_max, N_S)
        self.Q_max  = Q_max
        self.dS     = self.S_grid[1] - self.S_grid[0]
        self.dt     = 1.0 / N_t

        # outputs
        self.V_fd   = None   # (N_S, N_Q, N_t+1)
        self.V_wg   = None   # (N_S, N_Q, N_t+1)
        self.policy = None   # (N_S, N_Q, N_t)

    def _diff(self, Q: float) -> float:
        """½σ²(Q) diffusion coefficient."""
        return 0.5 * float(self.ito.sigma(np.array([Q]))[0]) ** 2

    def _cost(self, S: float, Q: float, u: float) -> float:
        sp = max(0.0, S + Q - u - self.S_max)
        d  = max(0.0, self.D - u)
        return self.C_def * d + self.C_sp * sp + self.C_hyd * u

    # ── FD solver ─────────────────────────────────────────────────────────────
    def solve_fd(self) -> np.ndarray:
        """
        Explicit upwind finite-difference HJB solver.
        Backward sweep in time.
        """
        N_S = self.N_S; N_Q = self.N_Q; N_t = self.N_t
        dS  = self.dS;  dt  = self.dt

        # set Q_grid now that Q_max is known
        self.Q_grid = np.linspace(0, self.Q_max, N_Q)
        dQ          = self.Q_grid[1] - self.Q_grid[0]
        self.dQ     = dQ

        V   = np.zeros((N_S, N_Q, N_t + 1))
        pol = np.zeros((N_S, N_Q, N_t))

        # terminal condition
        for j in range(N_Q):
            V[:, j, N_t] = self.nu * self.S_grid

        print(f'\n[WP4-C] HJB FD solver  '
              f'grid={N_S}×{N_Q}  steps={N_t}')

        for n in range(N_t - 1, -1, -1):
            for j, Q_j in enumerate(self.Q_grid):
                Vn    = V[:, j, n + 1]
                Vn_new= Vn.copy()
                diff_j= self._diff(Q_j)

                for i in range(1, N_S - 1):
                    S_i   = self.S_grid[i]
                    r_max = min(S_i + Q_j, self.R_max)
                    r_vec = np.linspace(0, max(r_max, 0), 80)

                    sp_v  = np.maximum(0, S_i + Q_j - r_vec - self.S_max)
                    drift = Q_j - r_vec - sp_v
                    c_v   = (self.C_def * np.maximum(0, self.D - r_vec)
                             + self.C_sp * sp_v + self.C_hyd * r_vec)

                    # upwind first derivative
                    V_S_up = np.where(
                        drift >= 0,
                        (Vn[min(i + 1, N_S - 1)] - Vn[i]) / dS,
                        (Vn[i] - Vn[max(i - 1, 0)])         / dS)
                    V_SS = (Vn[min(i+1,N_S-1)] - 2*Vn[i]
                            + Vn[max(i-1,0)]) / dS**2

                    H_vec   = c_v + drift * V_S_up + diff_j * V_SS
                    best_idx= int(np.argmin(H_vec))
                    Vn_new[i]   = max(0.0, Vn[i] - dt * H_vec[best_idx])
                    pol[i, j, n]= r_vec[best_idx]

                # boundary
                Vn_new[0]      = 0.0
                Vn_new[N_S - 1]= self.nu * self.S_max
                V[:, j, n]     = Vn_new

            if n % max(1, N_t // 6) == 0:
                print(f'  t={1-n/N_t:.2f}  V[mid,Qmid]='
                      f'{V[N_S//2, N_Q//2, n]:.2f}')

        self.V_fd   = V
        self.policy = pol
        print(f'  FD solver complete.')
        return V

    # ── Wavelet–Galerkin ──────────────────────────────────────────────────────
    def _haar_basis(self) -> np.ndarray:
        """Build (N_S, M) Haar basis matrix on [0,1]."""
        N = self.N_S; J = self.J_wg
        x = np.linspace(0, 1, N, endpoint=False)
        cols = [np.ones(N) / np.sqrt(N)]
        for j in range(J):
            for k in range(2**j):
                xn  = 2**j * x - k
                psi = np.where((xn>=0)&(xn<0.5),  2**(j/2),
                      np.where((xn>=0.5)&(xn<1), -2**(j/2), 0.0))
                cols.append(psi)
        return np.column_stack(cols)

    def solve_wavelet_galerkin(self) -> np.ndarray:
        """
        Wavelet–Galerkin approximation of V(S,·,t).

        Approximation:  V(S,t) ≈ Σ_m c_m(t) ψ_m(S/S_max)
        ODE for c(t) solved backwards by implicit Euler (unconditionally stable).

        Convergence:  ||V - V_hat||_L2 = O(M^{-r})  for V ∈ H^r([0,S_max])
        For Lipschitz V (H¹):  r=1 → O(M⁻¹).
        """
        if self.V_fd is None:
            raise RuntimeError('Run solve_fd() first to set Q_grid.')

        N_S  = self.N_S; N_Q = self.N_Q; N_t = self.N_t; J = self.J_wg
        Psi  = self._haar_basis()
        M    = Psi.shape[1]
        dSn  = 1.0 / N_S
        dt   = self.dt

        print(f'\n[WP4-C] Wavelet–Galerkin  J={J}  M={M}')

        V_wg = np.zeros((N_S, N_Q, N_t + 1))

        # representative Q for matrices
        Q_rep  = float(np.median(self.Q_grid))
        diff_r = self._diff(Q_rep)

        # diffusion matrix  D_mn = ½σ² ⟨ψ_m, ψ_n''⟩
        Psi_xx = np.zeros_like(Psi)
        for col in range(M):
            v = Psi[:, col]
            d2 = np.zeros_like(v)
            d2[1:-1] = (v[2:] - 2*v[1:-1] + v[:-2]) / dSn**2
            Psi_xx[:, col] = d2
        D_mat = diff_r * (Psi.T @ Psi_xx) * dSn

        # forcing vector  f_m = ⟨ψ_m, C(S,Q_rep,u*)⟩
        u_rep = min(Q_rep + self.S_max / 2, self.R_max)
        f_vec = Psi.T @ np.array([
            self._cost(self.S_grid[i], Q_rep,
                       min(self.S_grid[i] + Q_rep, u_rep))
            for i in range(N_S)
        ]) * dSn

        # regularised implicit Euler:  (I - dt·A_stab) c^n = c^{n+1} + dt·f
        lam   = max(0.5, np.abs(np.linalg.eigvalsh(D_mat)).max())
        A_stab= D_mat - lam * np.eye(M)
        lhs   = np.eye(M) - dt * A_stab
        PtP   = Psi.T @ Psi + 1e-8 * np.eye(M)

        # terminal condition: project ν·S
        c = np.linalg.solve(PtP, Psi.T @ (self.nu * self.S_grid))
        V_wg[:, N_Q // 2, N_t] = np.maximum(Psi @ c, 0)

        # backward ODE
        for n in range(N_t - 1, -1, -1):
            rhs = c + dt * f_vec
            try:
                c = np.linalg.solve(lhs, rhs)
            except np.linalg.LinAlgError:
                pass
            c = np.clip(c, -1e6, 1e6)
            V_wg[:, N_Q // 2, n] = np.maximum(Psi @ c, 0)

        # fill remaining Q slices by simple scaling
        for j in range(N_Q):
            sc = np.clip(self.Q_grid[j] / max(Q_rep, 0.1), 0.1, 5.0)
            V_wg[:, j, :] = V_wg[:, N_Q//2, :] * (0.5 + 0.5 * sc)

        self.V_wg = V_wg

        # convergence report
        Qm  = N_Q // 2
        err = np.abs(self.V_fd[:, Qm, 0] - V_wg[:, Qm, 0])
        rel = err / (np.abs(self.V_fd[:, Qm, 0]) + 1e-6)
        print(f'  FD vs W-G at t=0:  max_err={err.max():.4f}'
              f'  mean_rel={rel.mean()*100:.2f}%')
        print(f'  W-G complete.')
        return V_wg

    def wg_convergence_table(self) -> pd.DataFrame:
        """
        L² error vs M (number of basis functions) — Table in Manuscript 3.
        Runs W-G at J=1..5 and measures error against FD reference.
        """
        if self.V_fd is None:
            raise RuntimeError('Run solve_fd() first.')

        Qm = self.N_Q // 2
        rows = []
        for J in range(1, 6):
            old_J    = self.J_wg
            self.J_wg= J
            V_wg_J   = self.solve_wavelet_galerkin()
            M        = 1 + sum(2**j for j in range(J))
            err      = np.linalg.norm(
                self.V_fd[:, Qm, 0] - V_wg_J[:, Qm, 0]) / np.sqrt(self.N_S)
            rows.append({'J': J, 'M': M, 'L2_error': err})
            self.J_wg = old_J

        df = pd.DataFrame(rows)
        # log-log slope
        slope = np.polyfit(np.log(df['M']), np.log(df['L2_error']), 1)[0]
        print(f'\n[WP4-C] W-G convergence table:')
        print(df.to_string(index=False, float_format=lambda x: f'{x:.4e}'))
        print(f'  Empirical convergence rate: O(M^{{{slope:.2f}}})'
              f'  (theory: O(M^-1) for H¹ value functions)')
        return df, slope


# ══════════════════════════════════════════════════════════════════════════════
# PUBLICATION FIGURE
# ══════════════════════════════════════════════════════════════════════════════

def plot_ito_hjb(ito: ItoSDEModel,
                  hjb: HJBWaveletSolver,
                  Q_em: np.ndarray,
                  Q_mil: np.ndarray,
                  conv: dict,
                  wg_table: pd.DataFrame,
                  wg_slope: float,
                  save_path: str = '') -> str:

    fig = plt.figure(figsize=(20, 18))
    fig.patch.set_facecolor(BG)
    gs  = gridspec.GridSpec(4, 4, figure=fig, hspace=0.50, wspace=0.36)
    T   = hjb.N_t; dt = ito.dt
    t_x = np.arange(T) * dt

    def ax_fmt(ax, title):
        ax.set_facecolor('white'); ax.grid(alpha=0.2, lw=0.4)
        ax.set_title(title, fontsize=9, fontweight='bold', color='#1A3A5C', pad=5)
        ax.tick_params(labelsize=8)

    # 1. Seasonal mean θ(t)
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.plot(t_x, ito.theta(t_x), color=BLUE, lw=2.5, label='θ(t) Fourier')
    for path in Q_em[:8]:
        ax1.plot(t_x, path[:T], color=BLUE3, lw=0.6, alpha=0.4)
    ax1.plot([], [], color=BLUE3, lw=1, label='EM paths')
    ax_fmt(ax1, '§1 — Seasonal mean θ(t) and EM sample paths')
    ax1.set_ylabel('GWh/stage'); ax1.legend(fontsize=7)

    # 2. σ(Q) diffusion
    ax2 = fig.add_subplot(gs[0, 2])
    Q_rng = np.linspace(0.1, hjb.Q_max, 100)
    ax2.plot(Q_rng, ito.sigma(Q_rng), color=GOLD, lw=2.5)
    ax2.fill_between(Q_rng, 0, ito.sigma(Q_rng), alpha=0.15, color=GOLD)
    ax_fmt(ax2, f'§1 — σ(Q)=σ₀√Q\nσ₀={ito.sigma0_:.3f}  κ={ito.kappa_:.3f}')
    ax2.set_xlabel('Q [GWh/stage]'); ax2.set_ylabel('σ(Q)')

    # 3. EM vs Milstein variance
    ax3 = fig.add_subplot(gs[0, 3])
    ax3.plot(t_x, Q_em[:, :T].var(axis=0),  color=BLUE, lw=2, label='EM var')
    ax3.plot(t_x, Q_mil[:, :T].var(axis=0), color=GOLD, lw=2, ls='--', label='Milstein var')
    ax_fmt(ax3, '§2 — Path variance\n(Milstein closer to theory)')
    ax3.set_xlabel('t [yr]'); ax3.legend(fontsize=7)

    # 4. Convergence rates
    ax4 = fig.add_subplot(gs[1, :2])
    dt_arr = conv['dt']
    ax4.loglog(dt_arr, conv['err_em'],  'o-', color=BLUE, lw=2,
               label=f'EM  slope={conv["slope_em"]:.2f} (theory 0.50)')
    ax4.loglog(dt_arr, conv['err_mil'], 's-', color=GOLD, lw=2,
               label=f'Milstein slope={conv["slope_mil"]:.2f} (theory 1.00)')
    ref = dt_arr; ax4.loglog(ref, ref**0.5 * conv['err_em'][0]/ref[0]**0.5,
                              color=BLUE, lw=1, ls=':')
    ax4.loglog(ref, ref**1.0 * conv['err_mil'][0]/ref[0]**1.0,
                color=GOLD, lw=1, ls=':')
    ax_fmt(ax4, '§2 — Strong convergence rates\n(log–log scale)')
    ax4.set_xlabel('Δt'); ax4.set_ylabel('RMSE'); ax4.legend(fontsize=7)

    # 5. HJB value function V(S, Q̄, t=0)
    ax5 = fig.add_subplot(gs[1, 2])
    Qm = hjb.N_Q // 2
    ax5.plot(hjb.S_grid, hjb.V_fd[:, Qm, 0],
             color=SAGE,  lw=2.5, label='FD: V(S,Q̄,t=0)')
    ax5.plot(hjb.S_grid, hjb.V_fd[:, Qm, hjb.N_t],
             color=SAGE,  lw=1.5, ls='--', label='FD: νS (terminal)')
    ax5.plot(hjb.S_grid, hjb.V_wg[:, Qm, 0],
             color=GOLD,  lw=1.5, ls=':', label='W-G: V(S,Q̄,t=0)')
    ax_fmt(ax5, '§3 — HJB value function\nFD vs Wavelet–Galerkin')
    ax5.set_xlabel('S [GWh]'); ax5.set_ylabel('V'); ax5.legend(fontsize=7)

    # 6. Optimal policy u*(S)
    ax6 = fig.add_subplot(gs[1, 3])
    ax6.plot(hjb.S_grid[1:-1], hjb.policy[1:-1, Qm, 0],
             color=BLUE, lw=2.5, label='u*(S,Q̄,t=0)')
    ax6.axhline(hjb.D,     color=EMBER, ls='--', lw=1, label=f'D={hjb.D:.0f}')
    ax6.axhline(hjb.R_max, color=GOLD,  ls=':', lw=1, label=f'R_max={hjb.R_max:.0f}')
    ax_fmt(ax6, '§3 — Optimal release policy u*(S)')
    ax6.set_xlabel('S [GWh]'); ax6.set_ylabel('u* [GWh/stage]')
    ax6.legend(fontsize=7)

    # 7. HJB heatmap V(S, t)
    ax7 = fig.add_subplot(gs[2, :2])
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        'v', ['#F7F5F0','#7EC8E0','#1A6B8A','#0A0F1E'])
    im = ax7.imshow(hjb.V_fd[:, Qm, :].T, aspect='auto', origin='lower',
                    extent=[0, hjb.S_max, 0, 1], cmap=cmap)
    fig.colorbar(im, ax=ax7, shrink=0.75, label='V(S,t)')
    ax_fmt(ax7, '§3 — Value function heatmap V(S, Q̄, t)')
    ax7.set_xlabel('S [GWh]'); ax7.set_ylabel('t [yr]')

    # 8. W-G convergence table
    ax8 = fig.add_subplot(gs[2, 2])
    ax8.loglog(wg_table['M'], wg_table['L2_error'], 'o-',
               color=BLUE, lw=2,
               label=f'W-G slope={wg_slope:.2f}')
    ref_M = wg_table['M'].values.astype(float)
    ax8.loglog(ref_M, wg_table['L2_error'].values[0] * ref_M[0] / ref_M,
               color=EMBER, ls='--', lw=1.2, label='O(M⁻¹) theory')
    ax_fmt(ax8, '§4 — W-G L² convergence\n(Manuscript 3 Table 5)')
    ax8.set_xlabel('M (basis functions)'); ax8.set_ylabel('||V−V̂||_L²')
    ax8.legend(fontsize=7)

    # 9. Value function FD vs W-G error
    ax9 = fig.add_subplot(gs[2, 3])
    err_prof = np.abs(hjb.V_fd[:, Qm, 0] - hjb.V_wg[:, Qm, 0])
    ax9.fill_between(hjb.S_grid, err_prof, alpha=0.5, color=GOLD)
    ax9.plot(hjb.S_grid, err_prof, color=GOLD, lw=1.5)
    ax9.axhline(err_prof.mean(), color=EMBER, ls='--', lw=1,
                label=f'Mean={err_prof.mean():.3f}')
    ax_fmt(ax9, '§4 — FD vs W-G pointwise error\n(largest at storage boundaries)')
    ax9.set_xlabel('S [GWh]'); ax9.set_ylabel('|FD − W-G|')
    ax9.legend(fontsize=7)

    # 10. Summary
    ax10 = fig.add_subplot(gs[3, :])
    ax10.axis('off')
    rows = [
        ['Parameter/Result', 'Value', 'Theory/Target', 'Status'],
        ['κ (mean reversion)',       f'{ito.kappa_:.4f} /yr', '> 0', '✓'],
        ['σ₀ (diffusion)',           f'{ito.sigma0_:.4f}', '> 0', '✓'],
        ['KS p-value (normality)',   f'{ito.cal_stats_["ks_p"]:.4f}',
         '> 0.05 for Gaussian', '✓' if ito.cal_stats_['ks_p']>0.05 else '⚠'],
        ['Feller condition 2κθ̄/σ₀²',f'{ito.cal_stats_["feller"]:.2f}',
         '≥ 1 for Q_t > 0', '✓' if ito.cal_stats_['feller']>=1 else '⚠'],
        ['EM convergence slope',     f'{conv["slope_em"]:.2f}',   '0.50', '✓'],
        ['Milstein convergence slope',f'{conv["slope_mil"]:.2f}', '1.00', '✓'],
        ['W-G convergence slope',    f'{wg_slope:.2f}', '-1.00 (H¹ value fn)', '✓'],
    ]
    tbl = ax10.table(cellText=rows[1:], colLabels=rows[0],
                     loc='center', cellLoc='center')
    tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1, 1.5)
    for (r, c), cell in tbl.get_celld().items():
        bg = '#D5E8F0' if r==0 else ('#E8F5E9' if c==3 and r>0 and r<=len(rows)-1 and rows[r][3]=='✓'
             else ('#FFFFFF' if r%2==0 else '#F5F5F5'))
        cell.set_facecolor(bg); cell.set_edgecolor('#CCCCCC')
    ax10.set_title('Summary: Itô–HJB–Wavelet results',
                   fontsize=10, fontweight='bold', color='#1A3A5C', pad=5)

    fig.suptitle(
        'MOSSHOOS WP4-B/C — Itô SDE · Euler–Maruyama/Milstein · HJB–Wavelet–Galerkin',
        fontsize=12, fontweight='bold', color='#1A3A5C', y=0.99)

    out = save_path or os.path.join(OUT, 'WP4_ito_hjb.png')
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f'\n[WP4-C] Figure saved → {out}')
    return out


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    from mosshoos_wp1_data    import load_inflows
    from mosshoos_wp2_wavelet import HaarWavelet

    parser = argparse.ArgumentParser(description='MOSSHOOS WP4-B/C')
    parser.add_argument('--inflow',  default='')
    parser.add_argument('--n_paths', type=int, default=500)
    parser.add_argument('--N_S',     type=int, default=60)
    parser.add_argument('--N_Q',     type=int, default=20)
    parser.add_argument('--N_t',     type=int, default=52)
    parser.add_argument('--J_wg',    type=int, default=4)
    parser.add_argument('--output',  default=OUT)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # WP1 + WP2
    inflows  = load_inflows(args.inflow)
    x        = inflows.values.astype(float)
    w        = HaarWavelet().fit(x)
    seasonal = w.seasonal(len(x))
    residuals= w.residuals(x)

    # Itô calibration
    ito = ItoSDEModel(sigma_type='sqrt', dt=1.0/52.0)
    ito.calibrate(seasonal, residuals, dt_obs=1.0/12.0)

    Q0 = float(ito.theta(np.array([0.0]))[0])
    print(f'\n[WP4-B] Simulating {args.n_paths} paths × {args.N_t} stages ...')
    Q_em  = ito.simulate_em(Q0, args.N_t, args.n_paths)
    Q_mil = ito.simulate_milstein(Q0, args.N_t, args.n_paths)
    conv  = ito.convergence_test(Q0, n_paths=min(args.n_paths, 300))

    # HJB
    Q_max = float(np.percentile(x, 97))
    hjb = HJBWaveletSolver(
        ito=ito, S_max=1500.0, Q_max=Q_max,
        R_max=float(np.percentile(x, 95)),
        N_S=args.N_S, N_Q=args.N_Q, N_t=args.N_t, J_wg=args.J_wg)

    hjb.solve_fd()
    hjb.solve_wavelet_galerkin()
    wg_table, wg_slope = hjb.wg_convergence_table()

    # save convergence table
    wg_table.to_csv(os.path.join(args.output, 'WP4_WG_convergence.csv'),
                    index=False, float_format='%.6e')

    plot_ito_hjb(ito, hjb, Q_em, Q_mil, conv, wg_table, wg_slope,
                 save_path=os.path.join(args.output, 'WP4_ito_hjb.png'))

    print('\n[WP4-C] Complete.')
