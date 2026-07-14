#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mosshoos_wp2_wavelet.py
========================
WP2 — Multi-level Haar DWT decomposition of OMVS inflow series.

Provides:
  HaarWavelet              — decompose / reconstruct / energy analysis
  fourier_fit_seasonal()   — Fourier parameterisation of A_J for Itô SDE
  plot_decomposition()     — 6-panel publication figure

Run standalone:
    python mosshoos_wp2_wavelet.py
    python mosshoos_wp2_wavelet.py --inflow OMVS_inflows.csv --levels 4

Import:
    from mosshoos_wp2_wavelet import HaarWavelet
    w = HaarWavelet(levels=4).fit(inflows.values)
    seasonal  = w.seasonal(len(inflows))
    residuals = w.residuals(inflows.values)
"""

import argparse, os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr
warnings.filterwarnings('ignore')

np.random.seed(42)
OUT = '/mnt/user-data/outputs/'
BLUE='#1A6B8A'; BLUE2='#2E9BBF'; BLUE3='#7EC8E0'
GOLD='#B8860B'; SAGE='#2D6A4F'; EMBER='#C0392B'; BG='#F4F1EB'


# ══════════════════════════════════════════════════════════════════════════════
# HAAR WAVELET CLASS
# ══════════════════════════════════════════════════════════════════════════════

class HaarWavelet:
    """
    Multi-level orthonormal Haar Discrete Wavelet Transform.

    Mathematical basis
    ------------------
    Single-level QMF pair (Haar filters h, g):
        A_{j+1}(n) = (A_j(2n) + A_j(2n+1)) / sqrt(2)   [approximation]
        D_{j+1}(n) = (A_j(2n) - A_j(2n+1)) / sqrt(2)   [detail]

    Full decomposition:  x_t = A_J(t) + sum_{j=1}^J D_j(t)

    Energy fraction per level:
        eta_j = ||D_j||^2 / ||x||^2
        eta_0 = ||A_J||^2 / ||x||^2
        sum_j eta_j = 1  (Parseval)

    Scale selection criterion:
        J* = argmin_J { J : eta_0(J) >= eta_min }
    """

    def __init__(self, levels: int = 4):
        self.levels   = levels
        self.coeffs_  = []        # [A_J, D_J, ..., D_1]
        self.energy_  = []        # energy fraction per level
        self.labels_  = []

    # ── single-level ─────────────────────────────────────────────────────────
    @staticmethod
    def _fwd(x):
        n = len(x); h = n // 2
        a = (x[:2*h:2] + x[1:2*h:2]) / np.sqrt(2)
        d = (x[:2*h:2] - x[1:2*h:2]) / np.sqrt(2)
        return a, d

    @staticmethod
    def _inv(a, d):
        n  = len(a); rec = np.zeros(2 * n)
        rec[0::2] = (a + d) / np.sqrt(2)
        rec[1::2] = (a - d) / np.sqrt(2)
        return rec

    # ── multi-level ──────────────────────────────────────────────────────────
    def fit(self, signal: np.ndarray) -> 'HaarWavelet':
        """Decompose signal. Returns self for method chaining."""
        x = signal.astype(float).copy()
        cs = []
        for _ in range(self.levels):
            if len(x) < 2:
                break
            a, d = self._fwd(x)
            cs.append(d)
            x = a
        cs.append(x); cs.reverse()          # [A_J, D_J, ..., D_1]
        self.coeffs_ = cs

        tot   = sum(np.sum(c**2) for c in cs)
        self.energy_ = [np.sum(c**2) / tot * 100 for c in cs]
        J = len(cs) - 1
        self.labels_ = (
            ['A_J (seasonal trend)'] +
            [f'D_{J-i} (scale 2^{J-i} stages)' for i in range(J)]
        )
        return self

    def seasonal(self, N: int) -> np.ndarray:
        """Reconstruct A_J seasonal component (zero out all details)."""
        low = [self.coeffs_[0]] + [np.zeros_like(c) for c in self.coeffs_[1:]]
        cur = low[0].copy()
        for d in low[1:]:
            m   = min(len(cur), len(d))
            cur = self._inv(cur[:m], d[:m])
        out = np.zeros(N)
        m2  = min(N, len(cur))
        out[:m2] = cur[:m2]
        if m2 < N:
            out[m2:] = cur[-1]
        return out

    def residuals(self, signal: np.ndarray) -> np.ndarray:
        """Stochastic residuals: x - A_J."""
        return signal - self.seasonal(len(signal))

    def reconstruct(self) -> np.ndarray:
        """Full reconstruction from all coefficient levels."""
        cur = self.coeffs_[0].copy()
        for d in self.coeffs_[1:]:
            m   = min(len(cur), len(d))
            cur = self._inv(cur[:m], d[:m])
        return cur

    def print_summary(self):
        print(f'\n[WP2] Haar DWT  levels={self.levels}')
        for lab, pct in zip(self.labels_, self.energy_):
            bar = '\u2588' * max(1, int(pct / 3))
            print(f'  {lab:<36} {pct:5.1f}%  {bar}')
        print(f'  Seasonal energy A_J : {self.energy_[0]:.1f}%')

    @classmethod
    def optimal_levels(cls, signal: np.ndarray,
                        eta_min: float = 75.0,
                        max_levels: int = 7) -> int:
        """Return smallest J such that A_J captures >= eta_min % energy."""
        for J in range(1, max_levels + 1):
            w = cls(J).fit(signal)
            if w.energy_[0] >= eta_min:
                return J
        return max_levels


# ══════════════════════════════════════════════════════════════════════════════
# FOURIER FIT OF SEASONAL COMPONENT  (feeds Itô SDE in WP4-B)
# ══════════════════════════════════════════════════════════════════════════════

def fourier_fit_seasonal(seasonal: np.ndarray,
                          dt: float = 1.0 / 12.0,
                          K: int = 5) -> dict:
    """
    Fit truncated Fourier series to A_J seasonal component:
        theta(t) = a_0 + sum_{k=1}^K [a_k cos(2pi k t) + b_k sin(2pi k t)]

    Parameters
    ----------
    seasonal : A_J reconstructed series
    dt       : time step in years (1/12 for monthly, 1/52 for weekly)
    K        : number of harmonics to retain

    Returns
    -------
    dict with keys: a (array K+1), b (array K+1), R2, rmse
         and callable theta(t_array) -> values
    """
    N    = len(seasonal)
    t    = np.arange(N) * dt

    cols = [np.ones(N)]
    for k in range(1, K + 1):
        cols += [np.cos(2 * np.pi * k * t),
                 np.sin(2 * np.pi * k * t)]
    X    = np.column_stack(cols)
    beta = np.linalg.lstsq(X, seasonal, rcond=None)[0]

    a_f = np.zeros(K + 1)
    b_f = np.zeros(K + 1)
    a_f[0] = beta[0]
    for k in range(1, K + 1):
        a_f[k] = beta[2 * k - 1]
        b_f[k] = beta[2 * k]

    fitted = X @ beta
    ss_res = np.sum((seasonal - fitted) ** 2)
    ss_tot = np.sum((seasonal - seasonal.mean()) ** 2)
    R2     = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    rmse   = np.sqrt(ss_res / N)

    def theta(t_arr):
        t_arr = np.atleast_1d(t_arr)
        val   = a_f[0] * np.ones_like(t_arr)
        for k in range(1, K + 1):
            val += (a_f[k] * np.cos(2 * np.pi * k * t_arr)
                   + b_f[k] * np.sin(2 * np.pi * k * t_arr))
        return np.maximum(val, 0.1)

    result = dict(a=a_f, b=b_f, K=K, R2=R2, rmse=rmse, theta=theta)
    print(f'[WP2] Fourier fit  K={K}  R²={R2:.4f}  RMSE={rmse:.2f} GWh/stage')
    return result


# ══════════════════════════════════════════════════════════════════════════════
# ACF HELPER
# ══════════════════════════════════════════════════════════════════════════════

def acf(x: np.ndarray, max_lag: int = 24) -> np.ndarray:
    """Sample autocorrelation function up to max_lag."""
    x   = x - x.mean()
    std = np.std(x)
    if std == 0:
        return np.zeros(max_lag)
    return np.array([
        np.corrcoef(x[:-l], x[l:])[0, 1] if l > 0 else 1.0
        for l in range(max_lag)
    ])


# ══════════════════════════════════════════════════════════════════════════════
# PUBLICATION FIGURE
# ══════════════════════════════════════════════════════════════════════════════

def plot_decomposition(inflows: pd.Series,
                        wavelet: HaarWavelet,
                        seasonal: np.ndarray,
                        residuals: np.ndarray,
                        fourier: dict,
                        save_path: str = '') -> str:
    """
    6-panel WP2 figure:
      1. Inflow + A_J seasonal + Fourier fit
      2. Detail levels D_1..D_J
      3. Wavelet energy spectrum
      4. Residuals
      5. Residual ACF  (shows non-trivial autocorrelation)
      6. Partial wavelet reconstruction test
    """
    fig = plt.figure(figsize=(18, 14))
    fig.patch.set_facecolor(BG)
    gs  = gridspec.GridSpec(3, 3, figure=fig, hspace=0.48, wspace=0.36)
    t   = inflows.index

    def ax_fmt(ax, title):
        ax.set_facecolor('white')
        ax.grid(alpha=0.20, lw=0.4)
        ax.set_title(title, fontsize=9, fontweight='bold', color='#1A3A5C', pad=5)
        ax.tick_params(labelsize=7)

    # 1. Full series + seasonal + Fourier
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.plot(t, inflows.values, color=BLUE3, lw=0.7, alpha=0.7, label='Raw')
    ax1.plot(t, seasonal,       color=BLUE,  lw=2.0, label='A_J seasonal (IDWT)')
    t_yr = np.arange(len(inflows)) / 12.0
    ax1.plot(t, fourier['theta'](t_yr), color=GOLD, lw=1.5,
             ls='--', label=f'Fourier fit (K={fourier["K"]}, R²={fourier["R2"]:.3f})')
    ax_fmt(ax1, 'WP2 \u2014 OMVS inflow & wavelet seasonal component A_J')
    ax1.set_ylabel('GWh/month')
    ax1.legend(fontsize=7)

    # 2. Detail levels
    ax2 = fig.add_subplot(gs[0, 2])
    J   = wavelet.levels
    colours_d = [EMBER, GOLD, BLUE2, SAGE]
    for idx in range(1, min(J + 1, 5)):
        d = wavelet.coeffs_[idx]
        # stretch detail to original length for plotting
        d_full = np.repeat(d, 2 ** (idx - 0) * 1)[:len(inflows) // 2]
        lab = wavelet.labels_[idx][:12]
        if idx - 1 < len(colours_d):
            ax2.plot(d_full, color=colours_d[idx - 1], lw=0.9,
                     alpha=0.8, label=lab)
    ax_fmt(ax2, 'Detail levels D_j\n(stochastic residual components)')
    ax2.set_xlabel('Stage index'); ax2.legend(fontsize=6)

    # 3. Energy spectrum
    ax3 = fig.add_subplot(gs[1, 0])
    short_lbl = [l.split(' ')[0] for l in wavelet.labels_]
    bar_cols  = [BLUE] + [BLUE2, BLUE3, GOLD, '#CCCCCC', '#DDDDDD'][:J]
    ax3.bar(short_lbl, wavelet.energy_, color=bar_cols[:len(short_lbl)],
            edgecolor='white', lw=0.5)
    ax3.axhline(75, color=EMBER, ls='--', lw=1, label='75% threshold')
    ax_fmt(ax3, 'Wavelet energy spectrum (%)')
    ax3.set_ylabel('% signal energy')
    ax3.legend(fontsize=7)
    ax3.tick_params(axis='x', rotation=30)

    # 4. Residuals
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.fill_between(t, residuals, alpha=0.4, color=BLUE3)
    ax4.plot(t, residuals, color=BLUE, lw=0.6)
    ax4.axhline(0, color='black', lw=0.8, ls=':')
    ax_fmt(ax4, f'Stochastic residuals\n(x \u2212 A_J)  \u03c3={residuals.std():.1f} GWh/month')
    ax4.set_ylabel('GWh/month')

    # 5. ACF of residuals
    ax5 = fig.add_subplot(gs[1, 2])
    lag_max = min(25, len(residuals) // 4)
    acf_r   = acf(residuals, lag_max)
    ci      = 1.96 / np.sqrt(len(residuals))
    ax5.bar(range(lag_max), acf_r, color=BLUE, alpha=0.75, width=0.7)
    ax5.axhline( ci, color=EMBER, ls='--', lw=1, label='95% CI')
    ax5.axhline(-ci, color=EMBER, ls='--', lw=1)
    ax5.axhline(0,   color='black', lw=0.5)
    ax_fmt(ax5, 'Residual ACF\n(non-trivial \u2192 justifies block-bootstrap)')
    ax5.set_xlabel('Lag (months)')
    ax5.legend(fontsize=7)

    # 6. Reconstruction test
    ax6 = fig.add_subplot(gs[2, :])
    rec = wavelet.reconstruct()
    n   = min(len(inflows), len(rec))
    err = np.abs(inflows.values[:n] - rec[:n])
    ax6.plot(t[:n], inflows.values[:n], color=BLUE3, lw=1.2,
             label='Original')
    ax6.plot(t[:n], rec[:n], color=EMBER, lw=1.0, ls='--',
             label=f'Reconstructed  max_err={err.max():.3f}')
    ax_fmt(ax6, 'Parseval reconstruction test (original vs IDWT)')
    ax6.set_ylabel('GWh/month')
    ax6.legend(fontsize=8)

    fig.suptitle(
        f'MOSSHOOS WP2 \u2014 Haar DWT Decomposition  J={wavelet.levels}  '
        f'A_J energy={wavelet.energy_[0]:.1f}%',
        fontsize=12, fontweight='bold', color='#1A3A5C', y=0.99)

    out = save_path or os.path.join(OUT, 'WP2_wavelet.png')
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f'[WP2] Figure saved \u2192 {out}')
    return out


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    from mosshoos_wp1_data import load_inflows, validate_inflows

    parser = argparse.ArgumentParser(description='MOSSHOOS WP2 \u2014 Wavelet decomposition')
    parser.add_argument('--inflow', default='',
                        help='Path to OMVS monthly inflow CSV')
    parser.add_argument('--levels', type=int, default=0,
                        help='DWT levels (0 = auto-select by 75%% energy criterion)')
    parser.add_argument('--output', default=OUT)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # WP1
    inflows = load_inflows(args.inflow)
    _       = validate_inflows(inflows)
    x       = inflows.values.astype(float)

    # auto-select levels
    J = args.levels if args.levels > 0 else HaarWavelet.optimal_levels(x)
    print(f'[WP2] Using J={J} levels')

    # WP2
    w         = HaarWavelet(J).fit(x)
    w.print_summary()
    seasonal  = w.seasonal(len(x))
    residuals = w.residuals(x)
    fourier   = fourier_fit_seasonal(seasonal, dt=1.0/12.0, K=5)

    # save outputs
    out_df = pd.DataFrame({
        'date'      : inflows.index,
        'inflow'    : x,
        'seasonal'  : seasonal,
        'residual'  : residuals,
        'fourier'   : fourier['theta'](np.arange(len(x)) / 12.0),
    }).set_index('date')
    out_csv = os.path.join(args.output, 'WP2_decomposition.csv')
    out_df.to_csv(out_csv, float_format='%.4f')
    print(f'[WP2] Decomposition saved \u2192 {out_csv}')

    plot_decomposition(inflows, w, seasonal, residuals, fourier,
                       save_path=os.path.join(args.output, 'WP2_wavelet.png'))

    print('\n[WP2] Summary for manuscripts:')
    print(f'  J* (optimal levels) : {J}')
    print(f'  A_J energy fraction : {w.energy_[0]:.1f}%')
    print(f'  Residual std        : {residuals.std():.2f} GWh/month')
    print(f'  ACF lag-1           : {acf(residuals, 2)[1]:.3f}')
    print(f'  Fourier R²          : {fourier["R2"]:.4f}')
