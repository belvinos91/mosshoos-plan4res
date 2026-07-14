#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mosshoos_wp1_data.py
====================
WP1 — Data preparation for MOSSHOOS.

Provides:
  load_inflows()         — load real CSV or generate synthetic West African inflows
  validate_inflows()     — basic quality checks
  plot_inflows()         — exploratory data figure

HOW TO USE WITH REAL DATA
--------------------------
Your CSV must have two columns:
    date,inflow_GWh
    1994-01-01,8.2
    1994-02-01,5.1
    ...

Monthly frequency. At least 10 years. Source: SENELEC / OMVS Secretariat.

Then run:
    python mosshoos_wp1_data.py --inflow OMVS_inflows.csv --plot

Or import:
    from mosshoos_wp1_data import load_inflows
    inflows = load_inflows('OMVS_inflows.csv')

HOW TO USE WITHOUT REAL DATA (synthetic fallback)
--------------------------------------------------
    python mosshoos_wp1_data.py --plot
    # Uses realistic synthetic Sahelian hydrology
"""

import argparse
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
warnings.filterwarnings('ignore')

np.random.seed(42)

OUT = '/mnt/user-data/outputs/'

# ── colours ──────────────────────────────────────────────────────────────────
BLUE  = '#1A6B8A'
BLUE3 = '#7EC8E0'
GOLD  = '#B8860B'
SAGE  = '#2D6A4F'
EMBER = '#C0392B'
BG    = '#F4F1EB'


# ══════════════════════════════════════════════════════════════════════════════
# REAL DATA LOADER
# ══════════════════════════════════════════════════════════════════════════════

def load_inflows(csv_path: str = '',
                 n_years_synthetic: int = 30,
                 freq: str = 'MS') -> pd.Series:
    """
    Load OMVS monthly inflow series.

    Parameters
    ----------
    csv_path : str
        Path to CSV with columns [date, inflow_GWh].
        If empty or file not found, synthetic data is generated.
    n_years_synthetic : int
        Years of synthetic data to generate if no CSV.
    freq : str
        Pandas frequency string. 'MS' = month start (default).

    Returns
    -------
    pd.Series  index=DatetimeIndex, values=GWh/month, name='inflow_GWh'
    """
    if csv_path and os.path.isfile(csv_path):
        print(f'[WP1] Loading real inflows from {csv_path}')
        df = pd.read_csv(csv_path, index_col=0, parse_dates=True)

        # accept any column name for inflow
        col = df.columns[0]
        series = df[col].dropna().astype(float)

        # enforce monotone DatetimeIndex
        series.index = pd.to_datetime(series.index)
        series       = series.sort_index()
        series.name  = 'inflow_GWh'

        # resample to monthly if finer
        if series.index.inferred_freq not in ('MS', 'M', 'BMS'):
            series = series.resample('MS').mean()

        n = len(series)
        print(f'[WP1] Loaded {n} monthly records  '
              f'({series.index[0].date()} \u2192 {series.index[-1].date()})')
        print(f'[WP1] Mean={series.mean():.1f}  Std={series.std():.1f}  '
              f'Min={series.min():.1f}  Max={series.max():.1f} GWh/month')
        return series

    print('[WP1] No real CSV provided \u2014 generating synthetic West African inflows')
    return _generate_synthetic(n_years_synthetic)


# ══════════════════════════════════════════════════════════════════════════════
# SYNTHETIC GENERATOR
# ══════════════════════════════════════════════════════════════════════════════

def _generate_synthetic(n_years: int = 30) -> pd.Series:
    """
    Realistic synthetic monthly inflow for OMVS / Senegal River.

    Physical basis:
      - Sharp Sahelian rainy season: peak August (month 7), onset July,
        recession September\u2013October. Dry Jan\u2013May \u2248 0.
      - ENSO-like interannual cycle:  4.3-year period, amplitude 8 GWh/month.
      - Atlantic multi-decadal:       7.2-year period, amplitude 4 GWh/month.
      - AR(1) stochastic residual:    \u03c1 = 0.42, \u03c3 = 12 GWh/month.
      - Long-term trend:              +0.5 %/yr (slight positive from
                                      OMVS infrastructure improvements).
    """
    T  = n_years * 12
    t  = np.arange(T)
    m  = t % 12          # month index 0=Jan, 7=Aug

    # seasonal envelope
    seas = (8.0
            + 60.0 * np.exp(-0.5 * ((m - 7.0) / 1.8) ** 2)   # Aug peak
            + 18.0 * np.exp(-0.5 * ((m - 8.0) / 1.3) ** 2)   # Sep tail
            + 8.0  * np.exp(-0.5 * ((m - 6.0) / 1.0) ** 2))  # Jul onset

    # interannual cycles (ENSO + AMO)
    yr   = t / 12.0
    enso = 8.0  * np.sin(2 * np.pi * yr / 4.3 + 1.1)
    amo  = 4.0  * np.sin(2 * np.pi * yr / 7.2 + 0.3)

    # long-term trend
    trend = 1.0 + 0.005 * yr

    # AR(1) residuals
    ar1 = np.zeros(T)
    for i in range(1, T):
        ar1[i] = 0.42 * ar1[i - 1] + np.random.normal(0, 12.0)

    inflows = np.maximum(seas * trend + enso + amo + ar1, 0.5)

    dates  = pd.date_range('1994-01', periods=T, freq='MS')
    series = pd.Series(inflows, index=dates, name='inflow_GWh')

    print(f'[WP1] Synthetic:  {T} months  '
          f'Mean={series.mean():.1f}  Std={series.std():.1f}  '
          f'Min={series.min():.1f}  Max={series.max():.1f} GWh/month')
    return series


# ══════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ══════════════════════════════════════════════════════════════════════════════

def validate_inflows(series: pd.Series) -> dict:
    """
    Basic quality checks on the inflow series.

    Returns dict of statistics; prints warnings.
    """
    stats = {
        'n'          : len(series),
        'n_years'    : len(series) / 12,
        'mean'       : series.mean(),
        'std'        : series.std(),
        'cv'         : series.std() / series.mean(),
        'min'        : series.min(),
        'max'        : series.max(),
        'p10'        : series.quantile(0.10),
        'p50'        : series.quantile(0.50),
        'p90'        : series.quantile(0.90),
        'n_zeros'    : (series <= 0).sum(),
        'n_missing'  : series.isna().sum(),
        # annual total
        'ann_mean'   : series.resample('Y').sum().mean(),
        'ann_std'    : series.resample('Y').sum().std(),
    }

    print('\n[WP1] Validation report:')
    print(f'  Records         : {stats["n"]} months ({stats["n_years"]:.1f} years)')
    print(f'  Mean/Std        : {stats["mean"]:.1f} / {stats["std"]:.1f} GWh/month')
    print(f'  CV              : {stats["cv"]:.3f}')
    print(f'  P10/P50/P90     : {stats["p10"]:.1f} / {stats["p50"]:.1f} / {stats["p90"]:.1f}')
    print(f'  Annual mean     : {stats["ann_mean"]:.1f} GWh/yr')
    print(f'  Annual std      : {stats["ann_std"]:.1f} GWh/yr')

    if stats['n_zeros'] > 0:
        print(f'  WARNING: {stats["n_zeros"]} zero/negative values found')
    if stats['n_missing'] > 0:
        print(f'  WARNING: {stats["n_missing"]} missing values found')
    if stats['n_years'] < 10:
        print(f'  WARNING: less than 10 years of data \u2014 wavelet decomposition may be unreliable')

    return stats


# ══════════════════════════════════════════════════════════════════════════════
# EXPLORATORY PLOT
# ══════════════════════════════════════════════════════════════════════════════

def plot_inflows(series: pd.Series,
                 stats:  dict,
                 save_path: str = '') -> str:
    """
    4-panel exploratory figure:
      1. Full time series
      2. Annual totals bar
      3. Monthly climatology (mean ± 1σ)
      4. Autocorrelation function
    """
    fig = plt.figure(figsize=(16, 10))
    fig.patch.set_facecolor(BG)
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.32)

    def ax_fmt(ax, title):
        ax.set_facecolor('white')
        ax.grid(alpha=0.22, lw=0.5)
        ax.set_title(title, fontsize=10, fontweight='bold', color='#1A3A5C', pad=5)
        ax.tick_params(labelsize=8)

    # 1. Time series
    ax1 = fig.add_subplot(gs[0, :])
    ax1.fill_between(series.index, series.values, alpha=0.3, color=BLUE3)
    ax1.plot(series.index, series.values, color=BLUE, lw=0.8, label='Monthly inflow')
    ann = series.resample('Y').sum() / 12
    ax1.plot(ann.index, ann.values, color=GOLD, lw=2.5, label='Annual mean', zorder=5)
    ax_fmt(ax1, 'WP1 \u2014 OMVS Monthly Inflows [GWh/month]')
    ax1.set_ylabel('GWh/month')
    ax1.legend(fontsize=8)

    # 2. Annual totals
    ax2 = fig.add_subplot(gs[1, 0])
    annual = series.resample('Y').sum()
    colors = [EMBER if v < annual.quantile(0.25)
              else (SAGE if v > annual.quantile(0.75) else BLUE)
              for v in annual.values]
    ax2.bar(annual.index.year, annual.values, color=colors,
            edgecolor='white', lw=0.4, width=0.7)
    ax2.axhline(annual.mean(), color=GOLD, lw=2, ls='--',
                label=f'Mean={annual.mean():.0f} GWh/yr')
    ax_fmt(ax2, 'Annual inflow totals [GWh/yr]')
    ax2.set_ylabel('GWh/yr')
    ax2.legend(fontsize=7)
    ax2.tick_params(axis='x', rotation=45)

    # 3. Monthly climatology
    ax3 = fig.add_subplot(gs[1, 1])
    monthly_grp = series.groupby(series.index.month)
    m_mean = monthly_grp.mean()
    m_std  = monthly_grp.std()
    months = range(1, 13)
    mlbl   = ['J','F','M','A','M','J','J','A','S','O','N','D']
    ax3.fill_between(months,
                     m_mean - m_std, m_mean + m_std,
                     alpha=0.25, color=BLUE3, label='\u00b11\u03c3')
    ax3.plot(months, m_mean.values, color=BLUE, lw=2.5, marker='o',
             ms=4, label='Monthly mean')
    ax3.set_xticks(months); ax3.set_xticklabels(mlbl, fontsize=8)
    ax_fmt(ax3, 'Monthly climatology (mean \u00b1 1\u03c3) [GWh/month]')
    ax3.set_ylabel('GWh/month')
    ax3.legend(fontsize=7)

    fig.suptitle('MOSSHOOS WP1 \u2014 OMVS Inflow Data Analysis',
                 fontsize=13, fontweight='bold', color='#1A3A5C', y=0.99)

    out = save_path or os.path.join(OUT, 'WP1_inflows.png')
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=BG)
    plt.close()
    print(f'[WP1] Figure saved \u2192 {out}')
    return out


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='MOSSHOOS WP1 \u2014 Data preparation')
    parser.add_argument('--inflow',  default='',
                        help='Path to OMVS monthly inflow CSV (date, inflow_GWh)')
    parser.add_argument('--plot',    action='store_true',
                        help='Produce exploratory figure')
    parser.add_argument('--output',  default=OUT)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    series = load_inflows(args.inflow)
    stats  = validate_inflows(series)

    # save validated series
    out_csv = os.path.join(args.output, 'WP1_inflows_validated.csv')
    series.to_csv(out_csv, header=True)
    print(f'[WP1] Validated series saved \u2192 {out_csv}')

    if args.plot:
        plot_inflows(series, stats, save_path=os.path.join(args.output, 'WP1_inflows.png'))
