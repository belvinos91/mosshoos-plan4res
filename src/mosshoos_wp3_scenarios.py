#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mosshoos_wp3_scenarios.py
==========================
WP3 — Stochastic scenario generation for Plan4RES SSV.

Provides:
  ScenarioGenerator        — block-bootstrap + K-means
  validate_scenarios()     — KS test, ACF, energy score (Manuscript 1 §4)
  plot_scenarios()         — publication figure

Run standalone:
    python mosshoos_wp3_scenarios.py
    python mosshoos_wp3_scenarios.py --inflow OMVS_inflows.csv --n_clusters 10

Import:
    from mosshoos_wp3_scenarios import ScenarioGenerator
    gen = ScenarioGenerator(n_clusters=10, block_size=4, n_bootstrap=200)
    eta, weights, types = gen.fit(seasonal, residuals)
"""

import argparse, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.cluster import KMeans
from scipy.stats import ks_2samp, kstest, norm
warnings.filterwarnings('ignore')

np.random.seed(42)
OUT   = '/mnt/user-data/outputs/'
BLUE  = '#1A6B8A'; BLUE2 = '#2E9BBF'; BLUE3 = '#7EC8E0'
GOLD  = '#B8860B'; SAGE  = '#2D6A4F'; EMBER = '#C0392B'; BG = '#F4F1EB'


# ══════════════════════════════════════════════════════════════════════════════
# SCENARIO GENERATOR
# ══════════════════════════════════════════════════════════════════════════════

class ScenarioGenerator:
    """
    Generate the stochastic scenario set η for the Plan4RES SSV model.

    Method
    ------
    1. Block-bootstrap resampling of wavelet residuals (block length b=4)
       → preserves empirical ACF at lags τ < b (Künsch 1989)
    2. Stage-varying Gaussian perturbation calibrated to σ_ε(t)
       → captures seasonal heteroskedasticity of residuals
    3. K-means clustering on bootstrap ensemble
       → reduces N_bootstrap → N_c representative scenarios
    4. Probability weights w_k = |C_k| / N_bootstrap
       → equivalent to minimising L² Wasserstein distance

    Parameters
    ----------
    n_bootstrap : int   — number of block-bootstrap resamples before K-means
    n_clusters  : int   — |η| representative scenarios
    block_size  : int   — block length b (in stages)
    perturb_frac: float — Gaussian perturbation fraction (default 0.25)
    random_state: int
    """

    def __init__(self,
                 n_bootstrap:  int   = 200,
                 n_clusters:   int   = 10,
                 block_size:   int   = 4,
                 perturb_frac: float = 0.25,
                 random_state: int   = 42):
        self.n_bootstrap  = n_bootstrap
        self.n_clusters   = n_clusters
        self.block_size   = block_size
        self.perturb_frac = perturb_frac
        self.random_state = random_state

        # fitted attributes
        self.centers_     = None   # (n_clusters, n_stages)
        self.weights_     = None   # (n_clusters,)
        self.types_       = None   # 'DRY' | 'median' | 'WET'
        self.vol_annual_  = None   # annual volume per scenario
        self.bootstrap_   = None   # full bootstrap ensemble (n_bootstrap, n_stages)

    def fit(self,
            seasonal:  np.ndarray,
            residuals: np.ndarray,
            n_stages:  int = 12) -> tuple:
        """
        Parameters
        ----------
        seasonal  : (N,) seasonal component A_J
        residuals : (N,) stochastic residuals x - A_J
        n_stages  : stages per planning cycle (12=monthly, 52=weekly)

        Returns
        -------
        centers : (n_clusters, n_stages)
        weights : (n_clusters,)
        types   : (n_clusters,) str array
        """
        T      = n_stages
        n_full = len(residuals) // T
        if n_full < 2:
            raise ValueError(
                f'Not enough data: need at least 2 full cycles of {T} stages, '
                f'got {len(residuals)} stages ({n_full} full cycles). '
                f'Reduce n_stages or provide longer inflow series.')

        R         = residuals[:n_full * T].reshape(n_full, T)
        seas_cyc  = seasonal[:T]
        sigma_t   = R.std(axis=0)          # stage-wise std of residuals
        B         = self.block_size

        np.random.seed(self.random_state)
        scenarios = []
        for _ in range(self.n_bootstrap):
            eps = np.zeros(T)
            n_blocks = T // B
            for b in range(n_blocks):
                yr  = np.random.randint(0, n_full)
                st  = np.random.randint(0, max(1, T - B + 1))
                src = R[yr, st:st + B]
                if len(src) < B:
                    src = np.pad(src, (0, B - len(src)), mode='edge')
                eps[b * B:(b + 1) * B] = src
            # Gaussian perturbation (stage-varying)
            noise = np.random.normal(0, sigma_t * self.perturb_frac)
            scenarios.append(np.maximum(seas_cyc + eps + noise, 0.1))

        scenarios = np.array(scenarios)    # (n_bootstrap, T)
        self.bootstrap_ = scenarios

        # K-means reduction
        km      = KMeans(n_clusters=self.n_clusters,
                          random_state=self.random_state,
                          n_init=20, max_iter=500)
        labels  = km.fit_predict(scenarios)
        centers = km.cluster_centers_
        weights = np.array([np.sum(labels == k) / self.n_bootstrap
                             for k in range(self.n_clusters)])

        # classify DRY / median / WET
        vols   = centers.sum(axis=1)
        p10    = np.percentile(vols, 10)
        p90    = np.percentile(vols, 90)
        types  = np.where(vols <= p10, 'DRY',
                 np.where(vols >= p90, 'WET', 'median'))

        self.centers_    = centers
        self.weights_    = weights
        self.types_      = types
        self.vol_annual_ = vols

        self._print_summary(T)
        return centers, weights, types

    def _print_summary(self, T):
        vols = self.vol_annual_
        print(f'\n[WP3] Scenario generation:')
        print(f'  Bootstrap N       : {self.n_bootstrap}')
        print(f'  Block size b      : {self.block_size} stages')
        print(f'  K-means clusters  : {self.n_clusters}')
        print(f'  Weighted mean vol : {np.dot(vols, self.weights_):.1f} GWh/cycle')
        print(f'  P10 / P50 / P90   : '
              f'{np.percentile(vols,10):.0f} / '
              f'{np.percentile(vols,50):.0f} / '
              f'{np.percentile(vols,90):.0f} GWh/cycle')
        print(f'  {"Scen":>4}  {"Vol":>8}  {"w_k":>6}  {"type":>8}')
        for k in range(self.n_clusters):
            print(f'  {k+1:>4}  {vols[k]:>8.0f}  {self.weights_[k]:>6.3f}  '
                  f'{self.types_[k]:>8}')


# ══════════════════════════════════════════════════════════════════════════════
# STATISTICAL VALIDATION  (Manuscript 1 §4)
# ══════════════════════════════════════════════════════════════════════════════

def _acf(x, max_lag):
    x   = x - x.mean()
    std = np.std(x)
    if std == 0:
        return np.zeros(max_lag)
    return np.array([
        np.corrcoef(x[:-l], x[l:])[0, 1] if l > 0 else 1.0
        for l in range(max_lag)
    ])


def _energy_score(scenarios: np.ndarray,
                  weights:   np.ndarray,
                  obs:       np.ndarray) -> float:
    """
    Energy score (proper scoring rule for multivariate forecasts):
        ES(P, x) = E_P||X - x|| - (1/2) E_P||X - X'||
    Approximated via importance sampling from weighted discrete measure.
    """
    N = len(scenarios)
    # E_P||X - obs||
    term1 = np.sum(weights * np.linalg.norm(scenarios - obs, axis=1))
    # (1/2) E_P||X - X'||
    diff = 0.0
    for i in range(N):
        for j in range(N):
            diff += weights[i] * weights[j] * np.linalg.norm(
                scenarios[i] - scenarios[j])
    return term1 - 0.5 * diff


def validate_scenarios(gen: ScenarioGenerator,
                        inflows: np.ndarray,
                        n_stages: int = 12) -> dict:
    """
    Statistical validation of scenario set η against historical record.

    Metrics (Manuscript 1 §4):
      KS statistic    — Kolmogorov–Smirnov on annual volumes
      ACF error       — max |ρ_scen(τ) - ρ_obs(τ)| for τ = 1..6
      Energy score    — proper scoring rule (lower = better)
      P10/P90 error   — percentile reproduction accuracy

    Returns dict of all metrics (also printed).
    """
    centers = gen.centers_
    weights = gen.weights_
    T       = n_stages

    # observed annual volumes
    n_full   = len(inflows) // T
    obs_vols = inflows[:n_full * T].reshape(n_full, T).sum(axis=1)

    # weighted scenario annual volumes
    scen_vols = centers.sum(axis=1)

    # ── KS test (annual volumes) ──────────────────────────────────────────────
    # Build weighted sample from scenario set
    n_ks     = 1000
    k_samp   = np.random.choice(len(centers), p=weights, size=n_ks)
    ks_stat, ks_p = ks_2samp(obs_vols, scen_vols[k_samp])

    # ── ACF error ─────────────────────────────────────────────────────────────
    obs_resid  = inflows - inflows.mean()
    acf_obs    = _acf(obs_resid, 7)[1:]   # lags 1..6
    acf_scen_  = []
    for sc, w in zip(centers, weights):
        acf_scen_.append(_acf(sc - sc.mean(), 7)[1:])
    acf_scen   = np.average(acf_scen_, axis=0, weights=weights)
    acf_err    = np.abs(acf_scen - acf_obs).max()

    # ── Energy score ──────────────────────────────────────────────────────────
    # Use first full cycle as observation
    obs_cycle  = inflows[:T]
    # compare wavelet vs i.i.d. bootstrap baseline
    np.random.seed(42)
    iid_scens  = np.array([
        np.maximum(np.random.choice(inflows, T), 0.1)
        for _ in range(len(centers))
    ])
    iid_weights= np.ones(len(centers)) / len(centers)
    es_wavelet = _energy_score(centers,   weights,     obs_cycle)
    es_iid     = _energy_score(iid_scens, iid_weights, obs_cycle)
    es_improvement = (es_iid - es_wavelet) / es_iid * 100

    # ── Percentile errors ─────────────────────────────────────────────────────
    obs_p10 = np.percentile(obs_vols, 10)
    obs_p90 = np.percentile(obs_vols, 90)
    p10_err = abs(np.percentile(scen_vols, 10) - obs_p10) / obs_p10 * 100
    p90_err = abs(np.percentile(scen_vols, 90) - obs_p90) / obs_p90 * 100

    results = dict(
        ks_stat         = ks_stat,
        ks_p            = ks_p,
        acf_max_err     = acf_err,
        energy_score_wv = es_wavelet,
        energy_score_iid= es_iid,
        es_improvement  = es_improvement,
        p10_error_pct   = p10_err,
        p90_error_pct   = p90_err,
        acf_obs         = acf_obs,
        acf_scen        = acf_scen,
    )

    print('\n[WP3] Statistical validation (Manuscript 1 §4):')
    print(f'  KS statistic D            : {ks_stat:.4f}  (p={ks_p:.3f})'
          + ('  ✓ accept H0' if ks_p > 0.05 else '  ✗ reject H0'))
    print(f'  ACF max error (lags 1-6)  : {acf_err:.4f}')
    print(f'  Energy score (wavelet)    : {es_wavelet:.2f}')
    print(f'  Energy score (i.i.d.)     : {es_iid:.2f}')
    print(f'  ES improvement            : {es_improvement:.1f}%')
    print(f'  P10 volume error          : {p10_err:.1f}%')
    print(f'  P90 volume error          : {p90_err:.1f}%')

    return results


# ══════════════════════════════════════════════════════════════════════════════
# PUBLICATION FIGURE
# ══════════════════════════════════════════════════════════════════════════════

def plot_scenarios(gen:         ScenarioGenerator,
                   inflows:     np.ndarray,
                   validation:  dict,
                   n_stages:    int   = 12,
                   save_path:   str   = '') -> str:

    centers = gen.centers_
    weights = gen.weights_
    types   = gen.types_
    T       = n_stages

    fig = plt.figure(figsize=(18, 12))
    fig.patch.set_facecolor(BG)
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.46, wspace=0.35)

    def ax_fmt(ax, title):
        ax.set_facecolor('white')
        ax.grid(alpha=0.20, lw=0.4)
        ax.set_title(title, fontsize=9, fontweight='bold', color='#1A3A5C', pad=5)
        ax.tick_params(labelsize=8)

    t_x = range(1, T + 1)

    # 1. Scenario set η
    ax1 = fig.add_subplot(gs[0, :2])
    cmap = {'DRY': EMBER, 'WET': BLUE, 'median': '#AAAAAA'}
    for k, (sc, w, tp) in enumerate(zip(centers, weights, types)):
        ax1.plot(t_x, sc, color=cmap[tp],
                 alpha=0.35 + 0.55 * w / weights.max(), lw=1.2)
    mean_sc = np.average(centers, axis=0, weights=weights)
    ax1.plot(t_x, mean_sc, color='#1A3A5C', lw=2.8,
             label=r'Weighted mean $\bar{q}$', zorder=5)
    from matplotlib.lines import Line2D
    leg = [Line2D([0],[0], color=EMBER, lw=2, label='DRY (≤P10)'),
           Line2D([0],[0], color=BLUE,  lw=2, label='WET (≥P90)'),
           Line2D([0],[0], color='#AAAAAA', lw=2, label='Median'),
           Line2D([0],[0], color='#1A3A5C', lw=2.8, label='Weighted mean')]
    ax1.legend(handles=leg, fontsize=7)
    ax_fmt(ax1, r'WP3 — Scenario set $\eta$  (block-bootstrap + K-means)')
    ax1.set_xlabel('Stage'); ax1.set_ylabel('GWh/stage')

    # 2. Scenario weights
    ax2 = fig.add_subplot(gs[0, 2])
    bar_c = [cmap[tp] for tp in types]
    ax2.bar(range(1, len(weights) + 1), weights, color=bar_c,
            edgecolor='white', lw=0.4)
    ax2.axhline(1 / len(weights), color='black', ls='--', lw=1,
                label='Uniform')
    ax_fmt(ax2, r'Scenario weights $w_k$')
    ax2.set_xlabel('Scenario k'); ax2.set_ylabel('Weight')
    ax2.legend(fontsize=7)

    # 3. Annual volume distribution
    ax3 = fig.add_subplot(gs[1, 0])
    n_full   = len(inflows) // T
    obs_vols = inflows[:n_full * T].reshape(n_full, T).sum(axis=1)
    ax3.hist(obs_vols, bins=15, color=BLUE3, alpha=0.7, density=True,
             label='Observed', edgecolor='white', lw=0.3)
    ax3.vlines(centers.sum(axis=1), 0,
               ax3.get_ylim()[1] if ax3.get_ylim()[1] > 0 else 0.01,
               color=bar_c, lw=2, alpha=0.7, label='Scenarios')
    ax_fmt(ax3, 'Annual volume distribution\n(KS test validation)')
    ax3.set_xlabel('Annual inflow [GWh/cycle]')
    ks = validation['ks_stat']; ksp = validation['ks_p']
    ax3.set_title(
        f'Annual volume dist  KS={ks:.3f} p={ksp:.3f}',
        fontsize=9, fontweight='bold', color='#1A3A5C', pad=5)
    ax3.legend(fontsize=7)

    # 4. ACF comparison
    ax4 = fig.add_subplot(gs[1, 1])
    lags     = range(1, 7)
    acf_obs  = validation['acf_obs']
    acf_scen = validation['acf_scen']
    ax4.bar([l - 0.2 for l in lags], acf_obs,  0.35, color=BLUE,  label='Observed ACF')
    ax4.bar([l + 0.2 for l in lags], acf_scen, 0.35, color=GOLD,  label='Scenario ACF', alpha=0.8)
    ax4.axhline(0, color='black', lw=0.5)
    ci = 1.96 / np.sqrt(len(inflows))
    ax4.axhline( ci, color=EMBER, ls='--', lw=1, label='95% CI')
    ax4.axhline(-ci, color=EMBER, ls='--', lw=1)
    ax4.set_title(
        f'ACF preservation  max_err={validation["acf_max_err"]:.4f}',
        fontsize=9, fontweight='bold', color='#1A3A5C', pad=5)
    ax4.set_xlabel('Lag'); ax4.legend(fontsize=7)
    ax4.grid(alpha=0.2, lw=0.4); ax4.set_facecolor('white')
    ax4.tick_params(labelsize=8)

    # 5. Energy score comparison
    ax5 = fig.add_subplot(gs[1, 2])
    methods = ['Wavelet\n(proposed)', 'i.i.d.\nbootstrap']
    scores  = [validation['energy_score_wv'], validation['energy_score_iid']]
    bars    = ax5.bar(methods, scores, color=[SAGE, EMBER],
                      edgecolor='white', lw=0.4, width=0.5)
    ax5.set_title(
        f'Energy score  improvement={validation["es_improvement"]:.1f}%',
        fontsize=9, fontweight='bold', color='#1A3A5C', pad=5)
    ax5.set_ylabel('Energy score (lower=better)')
    ax5.grid(alpha=0.2, lw=0.4, axis='y'); ax5.set_facecolor('white')
    ax5.tick_params(labelsize=8)
    for bar, s in zip(bars, scores):
        ax5.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                 f'{s:.1f}', ha='center', fontsize=8, fontweight='bold')

    fig.suptitle(
        f'MOSSHOOS WP3 — Stochastic Scenario Generation  |η|={len(centers)}',
        fontsize=12, fontweight='bold', color='#1A3A5C', y=0.99)

    out = save_path or os.path.join(OUT, 'WP3_scenarios.png')
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f'[WP3] Figure saved → {out}')
    return out


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    from mosshoos_wp1_data    import load_inflows
    from mosshoos_wp2_wavelet import HaarWavelet, fourier_fit_seasonal

    parser = argparse.ArgumentParser(description='MOSSHOOS WP3 — Scenario generation')
    parser.add_argument('--inflow',      default='')
    parser.add_argument('--n_clusters',  type=int, default=10)
    parser.add_argument('--n_bootstrap', type=int, default=200)
    parser.add_argument('--block_size',  type=int, default=4)
    parser.add_argument('--n_stages',    type=int, default=12,
                        help='Stages per planning cycle (12=monthly, 52=weekly)')
    parser.add_argument('--output',      default=OUT)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # WP1
    inflows = load_inflows(args.inflow)
    x       = inflows.values.astype(float)

    # WP2
    from mosshoos_wp2_wavelet import HaarWavelet
    w         = HaarWavelet().fit(x)
    seasonal  = w.seasonal(len(x))
    residuals = w.residuals(x)

    # WP3
    gen = ScenarioGenerator(
        n_bootstrap  = args.n_bootstrap,
        n_clusters   = args.n_clusters,
        block_size   = args.block_size,
    )
    eta, weights, types = gen.fit(seasonal, residuals, n_stages=args.n_stages)

    # validation
    val = validate_scenarios(gen, x, n_stages=args.n_stages)

    # save scenario CSV
    T      = args.n_stages
    eta_df = pd.DataFrame(
        eta,
        index=[f'scenario_{k+1}' for k in range(len(eta))],
        columns=[f'stage_{t+1}' for t in range(T)])
    eta_df['weight'] = weights
    eta_df['type']   = types
    out_csv = os.path.join(args.output, 'WP3_scenarios_eta.csv')
    eta_df.to_csv(out_csv, float_format='%.4f')
    print(f'[WP3] Scenarios saved → {out_csv}')

    # save validation CSV
    val_out = {k: v for k, v in val.items()
               if not isinstance(v, np.ndarray)}
    pd.Series(val_out).to_csv(
        os.path.join(args.output, 'WP3_validation.csv'), header=False)

    plot_scenarios(gen, x, val, n_stages=args.n_stages,
                   save_path=os.path.join(args.output, 'WP3_scenarios.png'))
