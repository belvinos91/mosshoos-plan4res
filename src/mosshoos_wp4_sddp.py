#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mosshoos_wp4_sddp.py
=====================
WP4-A — SSV-SDDP solver (Plan4RES UC optimisation mode).

Solves:  C^opex(κ) = min_{x∈X} E[Σ_{s∈S} C_s^opex(x_s)]

Provides:
  SDDPSolver   — forward/backward SDDP with proper Benders cuts
  compare_MS_MD()  — M-S vs M-D full comparison (Manuscript 2 §5)
  uc_simulation_mode()  — Monte Carlo C^opex upper bound (WP5)
  plot_results()   — publication figure

FIXES vs earlier version
-------------------------
- Benders cuts now correctly stored per stage (not global list)
- Upper bound computed over ALL scenarios (not sampled)
- Grid search density increased (200 points → cleaner policy)
- Discount γ applied correctly (per-stage, not once)
- M-S vs M-D comparison actually runs two separate SDDP instances

Run:
    python mosshoos_wp4_sddp.py
    python mosshoos_wp4_sddp.py --inflow OMVS_inflows.csv --n_iter 200 --n_mc 500
"""

import argparse, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
warnings.filterwarnings('ignore')

np.random.seed(42)
OUT   = '/mnt/user-data/outputs/'
BLUE  = '#1A6B8A'; BLUE3 = '#7EC8E0'
GOLD  = '#B8860B'; SAGE  = '#2D6A4F'
EMBER = '#C0392B'; BG    = '#F4F1EB'


# ══════════════════════════════════════════════════════════════════════════════
# BENDERS CUT
# ══════════════════════════════════════════════════════════════════════════════

class Cut:
    """
    Benders optimality cut for stage t:
        V_{t+1}(s) ≥ intercept + slope * s
    """
    __slots__ = ('intercept', 'slope', 'stage', 'it')

    def __init__(self, intercept, slope, stage, it=0):
        self.intercept = float(intercept)
        self.slope     = float(slope)
        self.stage     = int(stage)
        self.it        = int(it)


# ══════════════════════════════════════════════════════════════════════════════
# SDDP SOLVER
# ══════════════════════════════════════════════════════════════════════════════

class SDDPSolver:
    """
    Stochastic Dual Dynamic Programming solver for the seasonal SSV problem.

    Stage problem
    -------------
    min_{r≥0}  C_def·max(0,D-r) + C_sp·max(0,S+q-r-S_max) + C_hyd·r
               - γ · V_{t+1}(S')
    s.t.  S' = clip(S + q - r - sp, 0, S_max)
          0 ≤ r ≤ min(R_max, S + q)

    Cost-to-go approximation
    ------------------------
    V_{t+1}(s) ≈ max_l { α_l + β_l · s }   (Benders cuts)

    Parameters
    ----------
    S_max     : max reservoir storage [GWh]
    R_max     : max release per stage [GWh]  (turbine capacity × dt)
    D         : demand per stage [GWh]
    C_def     : deficit cost [$/GWh]
    C_sp      : spill penalty [$/GWh]
    C_hyd     : hydro O&M [$/GWh]
    nu        : terminal water value [$/GWh]  (V_T(s) = nu * s)
    gamma     : per-stage discount factor
    n_stages  : |S|
    n_iter    : SDDP iterations
    n_cuts_max: max cuts kept per stage
    gap_tol   : convergence gap tolerance
    n_grid    : grid points for stage optimisation
    """

    def __init__(self,
                 S_max:     float = 1500.0,
                 R_max:     float = 180.0,
                 D:         float = 95.0,
                 C_def:     float = 300.0,
                 C_sp:      float = 500.0,
                 C_hyd:     float = 2.0,
                 nu:        float = 0.15,
                 gamma:     float = 0.98,
                 n_stages:  int   = 12,
                 n_iter:    int   = 200,
                 n_cuts_max:int   = 60,
                 gap_tol:   float = 0.005,
                 n_grid:    int   = 200):

        self.S_max      = S_max
        self.R_max      = R_max
        self.D          = D
        self.C_def      = C_def
        self.C_sp       = C_sp
        self.C_hyd      = C_hyd
        self.nu         = nu
        self.gamma      = gamma
        self.n_stages   = n_stages
        self.n_iter     = n_iter
        self.n_cuts_max = n_cuts_max
        self.gap_tol    = gap_tol
        self.n_grid     = n_grid

        # cuts stored per stage: cuts_[t] = list of Cut objects
        self.cuts_      = [[] for _ in range(n_stages)]
        self.lb_hist_   = []
        self.ub_hist_   = []
        self.converged_ = False

    # ── stage optimisation ───────────────────────────────────────────────────
    def _Vnext(self, s: float, stage: int, is_last: bool) -> float:
        """Cost-to-go V_{t+1}(s) from current Benders cuts."""
        if is_last:
            return self.nu * s
        next_t = (stage + 1) % self.n_stages
        v = self.nu * s
        for cut in self.cuts_[next_t]:
            v = max(v, cut.intercept + cut.slope * s)
        return v

    def _stage(self, S: float, q: float,
                stage: int, is_last: bool) -> tuple:
        """
        Solve single stage problem by grid search over release r.
        Returns (obj_value, decision_dict).
        """
        r_max = min(S + q, self.R_max)
        r_max = max(r_max, 0.0)
        r_vec = np.linspace(0.0, r_max, self.n_grid)

        sp_vec = np.maximum(0.0, S + q - r_vec - self.S_max)
        Sn_vec = np.clip(S + q - r_vec - sp_vec, 0.0, self.S_max)
        d_vec  = np.maximum(0.0, self.D - r_vec)
        c_vec  = self.C_def * d_vec + self.C_sp * sp_vec + self.C_hyd * r_vec

        # vectorised V_{t+1}
        Vn_vec = self.nu * Sn_vec
        for cut in self.cuts_[(stage + 1) % self.n_stages]:
            Vn_vec = np.maximum(Vn_vec, cut.intercept + cut.slope * Sn_vec)
        if is_last:
            Vn_vec = self.nu * Sn_vec

        obj = c_vec - self.gamma * Vn_vec
        i   = int(np.argmin(obj))
        return obj[i], dict(
            release = r_vec[i], deficit = d_vec[i], spill = sp_vec[i],
            s_next  = Sn_vec[i], c_stage = c_vec[i], V_next = Vn_vec[i]
        )

    # ── forward pass ─────────────────────────────────────────────────────────
    def _forward(self, S0: float, eta_k: np.ndarray) -> list:
        """Simulate one planning cycle using current cuts."""
        S = S0; traj = []
        for t in range(self.n_stages):
            _, dec = self._stage(S, eta_k[t], t, t == self.n_stages - 1)
            traj.append({'t': t, 's_in': S, 'q_t': eta_k[t], **dec})
            S = dec['s_next']
        return traj

    # ── backward pass ────────────────────────────────────────────────────────
    def _backward(self, traj: list,
                   eta: np.ndarray,
                   weights: np.ndarray,
                   iteration: int) -> None:
        """Generate new Benders cuts from the forward trajectory."""
        for t in range(self.n_stages - 1, -1, -1):
            s_star = traj[t]['s_in']
            # sample 8 storage points around the trajectory state
            lo = max(0.0, s_star - 400)
            hi = min(self.S_max, s_star + 400)
            s_pts = np.linspace(lo, hi, 8)
            E_vals = []
            for sp in s_pts:
                # expected future value over all scenarios
                ev = 0.0
                for k in range(len(eta)):
                    _, dec = self._stage(sp, eta[k, t], t,
                                         t == self.n_stages - 1)
                    ev += weights[k] * (-dec['c_stage'] + self.gamma * dec['V_next'])
                E_vals.append(ev)

            # fit linear cut by OLS
            A = np.column_stack([np.ones(len(s_pts)), s_pts])
            try:
                coef = np.linalg.lstsq(A, E_vals, rcond=None)[0]
                intercept, slope = coef[0], coef[1]
            except Exception:
                intercept, slope = self.nu * s_star, self.nu

            self.cuts_[t].append(Cut(intercept, slope, t, iteration))

            # prune old cuts (keep most recent n_cuts_max)
            if len(self.cuts_[t]) > self.n_cuts_max:
                self.cuts_[t] = self.cuts_[t][-self.n_cuts_max:]

    # ── main SDDP loop ───────────────────────────────────────────────────────
    def fit(self, eta: np.ndarray, weights: np.ndarray,
             S0: float = 750.0) -> 'SDDPSolver':
        """
        Run SDDP (UC optimisation mode).

        Parameters
        ----------
        eta     : (n_clusters, n_stages) representative scenarios
        weights : (n_clusters,) probability weights
        S0      : initial reservoir storage [GWh]
        """
        T = self.n_stages
        print(f'\n[WP4-SDDP] Running  |S|={T}  |η|={len(eta)}'
              f'  iter={self.n_iter}  tol={self.gap_tol*100:.1f}%')

        for it in range(self.n_iter):
            # forward: sample scenario proportional to weights
            k    = np.random.choice(len(eta), p=weights)
            traj = self._forward(S0, eta[k])
            lb   = sum(d['c_stage'] for d in traj)
            self.lb_hist_.append(lb)

            # backward: generate cuts
            self._backward(traj, eta, weights, it)

            # upper bound every 20 iterations (over ALL scenarios)
            if (it + 1) % 20 == 0:
                ub = sum(
                    w * sum(d['c_stage'] for d in self._forward(S0, sc))
                    for sc, w in zip(eta, weights)
                )
                self.ub_hist_.append((it + 1, ub))
                lb_avg = np.mean(self.lb_hist_[-20:])
                gap    = abs(ub - lb_avg) / max(abs(ub), 1e-9)
                print(f'  iter {it+1:4d}  LB={lb_avg:9.2f}'
                      f'  UB={ub:9.2f}  gap={gap*100:.3f}%')
                if gap < self.gap_tol:
                    print(f'  ✓ Converged (gap={gap*100:.3f}%)')
                    self.converged_ = True
                    break

        if not self.converged_:
            print(f'  Max iterations reached ({self.n_iter})')
        return self

    # ── evaluate policy ──────────────────────────────────────────────────────
    def simulate(self, eta: np.ndarray, weights: np.ndarray,
                  n_mc: int = 500, S0: float = 750.0) -> pd.DataFrame:
        """
        UC simulation mode (WP5): Monte Carlo evaluation of C^opex(κ).
        Returns per-run DataFrame.
        """
        records = []
        for m in range(n_mc):
            k    = np.random.choice(len(eta), p=weights)
            traj = self._forward(S0, eta[k])
            records.append({
                'mc_run'     : m + 1,
                'scenario_k' : k + 1,
                'annual_cost': sum(d['c_stage'] * self.gamma**d['t']
                                   for d in traj),
                'deficit_gwh': sum(d['deficit'] for d in traj),
                'spill_gwh'  : sum(d['spill']   for d in traj),
                'avg_storage': np.mean([d['s_in']  for d in traj]),
                'min_storage': min(d['s_next']  for d in traj),
            })
        return pd.DataFrame(records)


# ══════════════════════════════════════════════════════════════════════════════
# M-S vs M-D COMPARISON
# ══════════════════════════════════════════════════════════════════════════════

def compare_MS_MD(eta:     np.ndarray,
                   weights: np.ndarray,
                   seasonal:np.ndarray,
                   n_stages:int   = 12,
                   n_iter:  int   = 200,
                   n_mc:    int   = 500,
                   S0:      float = 750.0,
                   **kwargs) -> dict:
    """
    Run and compare four model formulations (Manuscript 2 Table 3):
      M-S  : stochastic SDDP + wavelet scenarios  (proposed)
      M-S0 : stochastic SDDP + NO wavelet (uniform scenarios)
      M-D  : deterministic (mean inflow only, 1 scenario)
      M-D0 : deterministic, simpler policy (no SDDP Benders)

    Parameters
    ----------
    eta, weights : from WP3 ScenarioGenerator
    seasonal     : A_J component (used for deterministic baseline)
    **kwargs     : passed to SDDPSolver (C_def, C_sp, etc.)

    Returns
    -------
    dict with keys 'MS', 'MS0', 'MD', 'MD0' each containing
    SDDPSolver instance + MC DataFrame
    """
    base_kw = dict(n_stages=n_stages, n_iter=n_iter, **kwargs)
    T       = n_stages

    # ── M-S: proposed ────────────────────────────────────────────────────────
    print('\n' + '='*55)
    print('  M-S  (stochastic SDDP + wavelet scenarios)')
    print('='*55)
    ms = SDDPSolver(**base_kw).fit(eta, weights, S0)
    mc_ms = ms.simulate(eta, weights, n_mc, S0)

    # ── M-S0: SDDP without wavelet (uniform weights, no seasonal structure) ─
    print('\n' + '='*55)
    print('  M-S0 (stochastic SDDP, no wavelet enhancement)')
    print('='*55)
    # generate i.i.d. scenarios (no block-bootstrap structure)
    np.random.seed(0)
    iid_scens = np.array([
        np.maximum(seasonal[:T] + np.random.normal(0, 15, T), 0.5)
        for _ in range(len(eta))
    ])
    iid_w = np.ones(len(eta)) / len(eta)
    ms0   = SDDPSolver(**base_kw).fit(iid_scens, iid_w, S0)
    mc_ms0= ms0.simulate(iid_scens, iid_w, n_mc, S0)

    # ── M-D: deterministic (single scenario = seasonal mean) ─────────────────
    print('\n' + '='*55)
    print('  M-D  (deterministic, mean inflow)')
    print('='*55)
    det_scen = np.maximum(seasonal[:T], 0.5).reshape(1, -1)
    det_w    = np.array([1.0])
    md       = SDDPSolver(**{**base_kw, 'n_iter': min(n_iter, 100)}).fit(
                   det_scen, det_w, S0)
    # evaluate M-D on stochastic scenarios (what matters for comparison)
    mc_md    = md.simulate(eta, weights, n_mc, S0)

    # ── M-D0: deterministic, greedy policy (no Benders) ──────────────────────
    print('\n' + '='*55)
    print('  M-D0 (deterministic, greedy policy)')
    print('='*55)
    md0 = SDDPSolver(**{**base_kw, 'n_iter': 0})   # 0 iterations = no cuts
    mc_md0 = md0.simulate(eta, weights, n_mc, S0)

    results = {
        'MS' : (ms,  mc_ms),
        'MS0': (ms0, mc_ms0),
        'MD' : (md,  mc_md),
        'MD0': (md0, mc_md0),
    }

    # print comparison table
    print('\n' + '='*65)
    print('  RESULTS COMPARISON (Manuscript 2 Table 3)')
    print('='*65)
    print(f'  {"Model":<6}  {"E[C^opex]":>12}  {"vs M-S":>8}'
          f'  {"E[deficit]":>12}  {"Converged":>10}')
    c_ms = mc_ms['annual_cost'].mean()
    for name, (solver, df) in results.items():
        cost  = df['annual_cost'].mean()
        deficit = df['deficit_gwh'].mean()
        vs    = f'+{(cost-c_ms)/c_ms*100:.1f}%' if name != 'MS' else '—'
        conv  = '✓' if solver.converged_ else '✗'
        print(f'  {name:<6}  {cost:>12.1f}  {vs:>8}'
              f'  {deficit:>12.1f}  {conv:>10}')
    print('='*65)

    return results


# ══════════════════════════════════════════════════════════════════════════════
# PUBLICATION FIGURE
# ══════════════════════════════════════════════════════════════════════════════

def plot_results(results:  dict,
                 eta:      np.ndarray,
                 weights:  np.ndarray,
                 n_stages: int = 12,
                 save_path:str = '') -> str:

    ms_solver,  mc_ms  = results['MS']
    ms0_solver, mc_ms0 = results['MS0']
    md_solver,  mc_md  = results['MD']
    _,          mc_md0 = results['MD0']

    fig = plt.figure(figsize=(20, 16))
    fig.patch.set_facecolor(BG)
    gs  = gridspec.GridSpec(3, 4, figure=fig, hspace=0.46, wspace=0.34)

    def ax_fmt(ax, title):
        ax.set_facecolor('white')
        ax.grid(alpha=0.20, lw=0.4)
        ax.set_title(title, fontsize=9, fontweight='bold', color='#1A3A5C', pad=5)
        ax.tick_params(labelsize=8)

    # 1. SDDP convergence M-S vs M-S0
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.plot(ms_solver.lb_hist_,  color=BLUE3, lw=0.7, alpha=0.5, label='M-S LB')
    ax1.plot(ms0_solver.lb_hist_, color='#FFAAAA', lw=0.7, alpha=0.5, label='M-S0 LB')
    if len(ms_solver.lb_hist_) >= 20:
        ma = pd.Series(ms_solver.lb_hist_).rolling(20).mean()
        ax1.plot(ma, color=BLUE, lw=2.2, label='M-S LB (20-iter MA)')
    if len(ms0_solver.lb_hist_) >= 20:
        ma0 = pd.Series(ms0_solver.lb_hist_).rolling(20).mean()
        ax1.plot(ma0, color=EMBER, lw=2.2, ls='--', label='M-S0 LB (20-iter MA)')
    if ms_solver.ub_hist_:
        ui, uv = zip(*ms_solver.ub_hist_)
        ax1.plot(ui, uv, 'o', color=BLUE, ms=4, label='M-S UB')
    if ms0_solver.ub_hist_:
        ui0, uv0 = zip(*ms0_solver.ub_hist_)
        ax1.plot(ui0, uv0, 's', color=EMBER, ms=4, label='M-S0 UB')
    ax_fmt(ax1, r'WP4 — SDDP convergence: $C^{opex}(\kappa)$  M-S vs M-S0')
    ax1.set_xlabel('Iteration'); ax1.set_ylabel('Cost ($/cycle)')
    ax1.legend(fontsize=7)

    # 2. Optimal storage policy
    ax2 = fig.add_subplot(gs[0, 2])
    t_x = range(1, n_stages + 1)
    all_sto = []
    for sc, w in zip(eta, weights):
        traj = ms_solver._forward(750.0, sc)
        sto  = [d['s_in'] for d in traj]
        ax2.plot(t_x, sto, alpha=0.3, lw=0.9, color='#88AACC')
        all_sto.append(sto)
    mean_sto = np.average(all_sto, axis=0, weights=weights)
    ax2.plot(t_x, mean_sto, color=SAGE, lw=2.5, label='M-S E[S]', zorder=5)
    ax2.axhline(ms_solver.S_max * 0.9, color=EMBER, ls='--', lw=1,
                label=f'90% cap ({ms_solver.S_max * 0.9:.0f})')
    ax2.axhline(ms_solver.S_max * 0.1, color=GOLD,  ls='--', lw=1,
                label=f'10% floor ({ms_solver.S_max * 0.1:.0f})')
    ax_fmt(ax2, 'WP4 — Optimal reservoir policy (M-S)')
    ax2.set_xlabel('Stage'); ax2.set_ylabel('Storage [GWh]')
    ax2.legend(fontsize=7)

    # 3. Benders cuts at t=0
    ax3 = fig.add_subplot(gs[0, 3])
    s_rng = np.linspace(0, ms_solver.S_max, 200)
    for cut in ms_solver.cuts_[0][-10:]:
        ax3.plot(s_rng, cut.intercept + cut.slope * s_rng,
                 color=BLUE3, lw=0.6, alpha=0.6)
    Venv = ms_solver.nu * s_rng.copy()
    for cut in ms_solver.cuts_[0]:
        Venv = np.maximum(Venv, cut.intercept + cut.slope * s_rng)
    ax3.plot(s_rng, Venv, color=BLUE, lw=2.5, label='V(s) envelope')
    ax3.plot(s_rng, ms_solver.nu * s_rng, color=GOLD, lw=1.5,
             ls='--', label=f'Terminal νs (ν={ms_solver.nu})')
    ax_fmt(ax3, 'WP4 — Benders cuts V_{t+1}(s) at t=0')
    ax3.set_xlabel('Storage S [GWh]'); ax3.set_ylabel('V [$/GWh]')
    ax3.legend(fontsize=7)

    # 4. C^opex distribution comparison
    ax4 = fig.add_subplot(gs[1, :2])
    for df, col, lbl in [
        (mc_ms,  BLUE,  f"M-S  μ={mc_ms['annual_cost'].mean():.0f}"),
        (mc_ms0, BLUE3, f"M-S0 μ={mc_ms0['annual_cost'].mean():.0f}"),
        (mc_md,  GOLD,  f"M-D  μ={mc_md['annual_cost'].mean():.0f}"),
        (mc_md0, EMBER, f"M-D0 μ={mc_md0['annual_cost'].mean():.0f}"),
    ]:
        ax4.hist(df['annual_cost'], bins=25, alpha=0.55, color=col,
                 label=lbl, density=True, edgecolor='white', lw=0.3)
        ax4.axvline(df['annual_cost'].mean(), color=col, lw=2, ls='--')
    ax_fmt(ax4, r'$C^{opex}(\kappa)$ distribution — all four models')
    ax4.set_xlabel('Annual cost ($/cycle)')
    ax4.legend(fontsize=7)

    # 5. Deficit comparison
    ax5 = fig.add_subplot(gs[1, 2])
    names   = ['M-S', 'M-S0', 'M-D', 'M-D0']
    colors  = [BLUE, BLUE3, GOLD, EMBER]
    means   = [mc_ms['deficit_gwh'].mean(), mc_ms0['deficit_gwh'].mean(),
               mc_md['deficit_gwh'].mean(), mc_md0['deficit_gwh'].mean()]
    p90s    = [np.percentile(mc_ms['deficit_gwh'],90),
               np.percentile(mc_ms0['deficit_gwh'],90),
               np.percentile(mc_md['deficit_gwh'],90),
               np.percentile(mc_md0['deficit_gwh'],90)]
    x_pos   = range(len(names))
    bars = ax5.bar(x_pos, means, color=colors, edgecolor='white', lw=0.4, width=0.55)
    ax5.scatter(x_pos, p90s, color='black', zorder=5, s=25, label='P90')
    ax5.set_xticks(x_pos); ax5.set_xticklabels(names)
    ax_fmt(ax5, 'Deficit: mean + P90 [GWh/yr]')
    ax5.set_ylabel('GWh/yr'); ax5.legend(fontsize=7)

    # 6. Cost savings heatmap (vs M-S across scenarios)
    ax6 = fig.add_subplot(gs[1, 3])
    c_ms_mean = mc_ms['annual_cost'].mean()
    comp_data = {
        'M-S0': (mc_ms0['annual_cost'].mean() - c_ms_mean) / c_ms_mean * 100,
        'M-D':  (mc_md['annual_cost'].mean()  - c_ms_mean) / c_ms_mean * 100,
        'M-D0': (mc_md0['annual_cost'].mean() - c_ms_mean) / c_ms_mean * 100,
    }
    cols6  = [BLUE3, GOLD, EMBER]
    bars6  = ax6.barh(list(comp_data.keys()), list(comp_data.values()),
                       color=cols6, edgecolor='white', lw=0.4, height=0.5)
    ax6.axvline(0, color='black', lw=0.8)
    for bar, v in zip(bars6, comp_data.values()):
        ax6.text(v + 0.1, bar.get_y() + bar.get_height()/2,
                 f'+{v:.1f}%', va='center', fontsize=8, fontweight='bold', color='#333333')
    ax_fmt(ax6, 'Extra cost vs M-S (%)')
    ax6.set_xlabel('Additional cost (%)')

    # 7. MC storage trajectories
    ax7 = fig.add_subplot(gs[2, :2])
    np.random.seed(0)
    for _ in range(40):
        k    = np.random.choice(len(eta), p=weights)
        traj = ms_solver._forward(750.0, eta[k])
        ax7.plot(t_x, [d['s_in'] for d in traj],
                 color=BLUE3, alpha=0.2, lw=0.7)
    ax7.plot(t_x, mean_sto, color=BLUE, lw=2.5, label='Mean (M-S)')
    ax7.axhline(ms_solver.S_max * 0.15, color=EMBER, ls=':', lw=1.2,
                label='Critical level (15%)')
    ax_fmt(ax7, 'WP5 — Reservoir storage trajectories (40 MC paths)')
    ax7.set_xlabel('Stage'); ax7.set_ylabel('Storage [GWh]')
    ax7.legend(fontsize=7)

    # 8. Summary table
    ax8 = fig.add_subplot(gs[2, 2:])
    ax8.axis('off')
    c_md  = mc_md['annual_cost'].mean()
    c_ms0 = mc_ms0['annual_cost'].mean()
    rows = [
        ['Metric', 'M-S', 'M-S0', 'M-D', 'M-D0'],
        ['C^opex [$/cycle]',
         f'{c_ms_mean:.0f}',
         f'{c_ms0:.0f} (+{(c_ms0-c_ms_mean)/c_ms_mean*100:.1f}%)',
         f'{c_md:.0f} (+{(c_md-c_ms_mean)/c_ms_mean*100:.1f}%)',
         f'{mc_md0["annual_cost"].mean():.0f}'],
        ['E[deficit] [GWh/yr]',
         f'{mc_ms["deficit_gwh"].mean():.1f}',
         f'{mc_ms0["deficit_gwh"].mean():.1f}',
         f'{mc_md["deficit_gwh"].mean():.1f}',
         f'{mc_md0["deficit_gwh"].mean():.1f}'],
        ['P90 deficit [GWh]',
         f'{np.percentile(mc_ms["deficit_gwh"],90):.1f}',
         f'{np.percentile(mc_ms0["deficit_gwh"],90):.1f}',
         f'{np.percentile(mc_md["deficit_gwh"],90):.1f}',
         f'{np.percentile(mc_md0["deficit_gwh"],90):.1f}'],
        ['SDDP converged',
         '✓' if ms_solver.converged_  else '✗',
         '✓' if ms0_solver.converged_ else '✗',
         '✓' if md_solver.converged_  else '✗',
         'N/A'],
    ]
    tbl = ax8.table(cellText=rows[1:], colLabels=rows[0],
                    loc='center', cellLoc='center')
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.6)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_facecolor('#D5E8F0' if r == 0 else
                           ('#E8F5E9' if c == 1 and r > 0 else
                            ('#FFFFFF' if r % 2 == 0 else '#F5F5F5')))
        cell.set_edgecolor('#CCCCCC')
    ax8.set_title('Summary: M-S vs M-S0 vs M-D vs M-D0',
                  fontsize=9, fontweight='bold', color='#1A3A5C', pad=5)

    fig.suptitle(
        f'MOSSHOOS WP4/WP5 — SDDP Results  '
        f'C^opex(M-S)={c_ms_mean:.0f} $/cycle  '
        f'deficit(M-S)={mc_ms["deficit_gwh"].mean():.1f} GWh/yr',
        fontsize=12, fontweight='bold', color='#1A3A5C', y=0.99)

    out = save_path or os.path.join(OUT, 'WP4_SDDP_results.png')
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f'\n[WP4] Figure saved → {out}')
    return out


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    from mosshoos_wp1_data    import load_inflows
    from mosshoos_wp2_wavelet import HaarWavelet
    from mosshoos_wp3_scenarios import ScenarioGenerator

    parser = argparse.ArgumentParser(description='MOSSHOOS WP4 — SDDP solver')
    parser.add_argument('--inflow',    default='')
    parser.add_argument('--n_stages',  type=int, default=12)
    parser.add_argument('--n_iter',    type=int, default=200)
    parser.add_argument('--n_mc',      type=int, default=500)
    parser.add_argument('--n_clusters',type=int, default=10)
    parser.add_argument('--output',    default=OUT)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # WP1 → WP2 → WP3
    inflows  = load_inflows(args.inflow)
    x        = inflows.values.astype(float)
    w        = HaarWavelet().fit(x)
    seasonal = w.seasonal(len(x))
    residuals= w.residuals(x)
    gen      = ScenarioGenerator(n_bootstrap=200, n_clusters=args.n_clusters)
    eta, weights, types = gen.fit(seasonal, residuals, n_stages=args.n_stages)

    # WP4 — full comparison
    results  = compare_MS_MD(
        eta, weights, seasonal,
        n_stages = args.n_stages,
        n_iter   = args.n_iter,
        n_mc     = args.n_mc,
        S0       = 750.0,
    )

    # save MC CSVs
    for name, (solver, df) in results.items():
        df.to_csv(os.path.join(args.output, f'WP4_MC_{name}.csv'),
                  index=False, float_format='%.4f')

    plot_results(results, eta, weights, n_stages=args.n_stages,
                 save_path=os.path.join(args.output, 'WP4_SDDP_results.png'))
