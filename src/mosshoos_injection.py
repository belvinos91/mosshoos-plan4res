#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mosshoos_injection.py
=====================
Drop-in bridge between PostTreatPlan4res.py and the MOSSHOOS stochastic pipeline.

HOW TO ACTIVATE IN PostTreatPlan4res.py
----------------------------------------
Add one line to settingsPostTreatPlan4res.yml:
    usevu_stochastic: true

Then paste the four code blocks below at the indicated line numbers.

INJECTION A  →  paste after line ~334  (after SS_SeasonalStorage block)
INJECTION B  →  paste after line ~1553 (after BV.to_csv)
INJECTION C  →  paste at   line ~1559  (replaces BV_VU loop header)
INJECTION D  →  paste at   line ~1517  (replaces max_{l in L_t} block)

NO existing Plan4RES code is deleted. The injections are additive:
all original code runs first, MOSSHOOS extends it when the flag is set.
"""

import os, logging, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')
logger = logging.getLogger(__name__)


# ══════════════════════════════════════════════════════════════════════════════
# MOSSHOOS CONFIG  (populated by Injection A)
# ══════════════════════════════════════════════════════════════════════════════

class MOSSHOOSConfig:
    """
    Holds physical parameters extracted from SS_SeasonalStorage.csv.
    Passed to all subsequent injection functions.
    """
    def __init__(self,
                 kappa_turb: float,   # total turbine capacity [MW]
                 V_max: float,        # max reservoir volume [MWh]
                 D_plan: float,       # planning demand [MWh/month]
                 C_def: float = 300.0,
                 C_sp:  float = 150.0,
                 C_hyd: float = 2.0,
                 nu:    float = 0.08,
                 gamma: float = 0.98,
                 n_stages: int = 12,
                 n_iter:   int = 200,
                 n_mc:     int = 500,
                 n_clusters: int = 10,
                 n_bootstrap: int = 200,
                 dt_hours: float = 730.0):   # hours per stage (730 = monthly)
        self.kappa_turb = kappa_turb
        self.V_max      = V_max
        self.D_plan     = D_plan
        self.C_def      = C_def
        self.C_sp       = C_sp
        self.C_hyd      = C_hyd
        self.nu         = nu
        self.gamma      = gamma
        self.n_stages   = n_stages
        self.n_iter     = n_iter
        self.n_mc       = n_mc
        self.n_clusters = n_clusters
        self.n_bootstrap= n_bootstrap
        self.dt_hours   = dt_hours
        # R_max: max energy release per stage [MWh] = turbine capacity × dt
        self.R_max = kappa_turb * dt_hours / 1000.0   # MWh/stage
        self.S0    = V_max * 0.5                       # 50% initial fill

    def __repr__(self):
        return (f"MOSSHOOSConfig("
                f"kappa_turb={self.kappa_turb:.0f} MW, "
                f"V_max={self.V_max:.0f} MWh, "
                f"R_max={self.R_max:.0f} MWh/stage, "
                f"D={self.D_plan:.0f} MWh/stage)")


# ══════════════════════════════════════════════════════════════════════════════
# INJECTION A  —  Extract OMVS parameters from SS_SeasonalStorage
# ══════════════════════════════════════════════════════════════════════════════

# ┌─────────────────────────────────────────────────────────────────────────────
# │ PASTE IN PostTreatPlan4res.py  AFTER LINE ~334  (end of SS block)
# │
# │     if isInvest:
# │         IsInvestedTechno[techno]=False
# │ else: listTechnosSS=[]
# │ logger.info('Seasonal storage technos: '+str(listTechnosSS))
# │
# │ # ---- MOSSHOOS INJECTION A ------------------------------------------------
# │ if cfg.get('usevu_stochastic', False) and listTechnosSS:
# │     from mosshoos_injection import inject_A
# │     _MOSSHOOS_CFG = inject_A(
# │         ss_data   = inputdata_save['SS_SeasonalStorage'],
# │         demand_gwh= cfg.get('mosshoos_demand_gwh', None),
# │         n_stages  = cfg.get('mosshoos_n_stages', 12),
# │         n_iter    = cfg.get('mosshoos_n_iter', 200),
# │         n_mc      = cfg.get('mosshoos_n_mc', 500),
# │         n_clusters= cfg.get('mosshoos_n_clusters', 10),
# │         dt_hours  = TimeStepHours,
# │     )
# │ # ---- END MOSSHOOS A -------------------------------------------------------
# └─────────────────────────────────────────────────────────────────────────────

def inject_A(ss_data: pd.DataFrame,
             demand_gwh: float = None,
             n_stages:   int   = 12,
             n_iter:     int   = 200,
             n_mc:       int   = 500,
             n_clusters: int   = 10,
             dt_hours:   float = 730.0) -> MOSSHOOSConfig:
    """
    INJECTION A — Extract OMVS physical parameters from SS_SeasonalStorage.csv.

    Parameters extracted:
        kappa_turb = MaxPower × NumberUnits  [MW total]
        V_max      = MaxVolume × NumberUnits [MWh total]

    These become R_max and S_max in the SDDP stage problem.

    Parameters
    ----------
    ss_data    : DataFrame from SS_SeasonalStorage.csv
    demand_gwh : monthly demand [GWh/month]. If None, defaults to 10% of R_max.
    n_stages   : planning stages (12=monthly, 52=weekly)
    n_iter     : SDDP iterations
    n_mc       : Monte Carlo runs for WP5
    n_clusters : K-means scenario count |eta|
    dt_hours   : hours per stage (730 for monthly)

    Returns
    -------
    MOSSHOOSConfig
    """
    kappa_turb = float((ss_data['MaxPower'] * ss_data['NumberUnits']).sum())
    V_max      = float((ss_data['MaxVolume'] * ss_data['NumberUnits']).sum())

    # Convert to GWh if values look like MWh
    if V_max > 1e6:
        V_max /= 1000.0
        logger.info("[MOSSHOOS A] Converted V_max from MWh to GWh")
    if kappa_turb > 50000:
        kappa_turb /= 1000.0
        logger.info("[MOSSHOOS A] Converted kappa_turb from MW to GW")

    R_max = kappa_turb * dt_hours / 1000.0   # GWh/stage
    D     = demand_gwh if demand_gwh else R_max * 0.40   # 40% of capacity

    cfg_m = MOSSHOOSConfig(
        kappa_turb  = kappa_turb,
        V_max       = V_max,
        D_plan      = D,
        n_stages    = n_stages,
        n_iter      = n_iter,
        n_mc        = n_mc,
        n_clusters  = n_clusters,
        dt_hours    = dt_hours,
    )

    logger.info(f"[MOSSHOOS A] {cfg_m}")
    logger.info(f"[MOSSHOOS A] R_max = {cfg_m.R_max:.1f} GWh/stage")
    logger.info(f"[MOSSHOOS A] D_plan = {cfg_m.D_plan:.1f} GWh/stage")
    logger.info(f"[MOSSHOOS A] S0 (50%) = {cfg_m.S0:.1f} GWh")
    return cfg_m


# ══════════════════════════════════════════════════════════════════════════════
# INJECTION B  —  Haar DWT on the inflow/volume series (WP2)
# ══════════════════════════════════════════════════════════════════════════════

# ┌─────────────────────────────────────────────────────────────────────────────
# │ PASTE IN PostTreatPlan4res.py  AFTER LINE ~1553
# │ (after:  if cfg['usevu']: BV.to_csv(...BellmanValuesPerScenario.csv))
# │
# │ # ---- MOSSHOOS INJECTION B ------------------------------------------------
# │ if cfg.get('usevu_stochastic', False):
# │     from mosshoos_injection import inject_B
# │     _MOSSHOOS_WAVELET = inject_B(
# │         volume_dir  = cfg['dirSto'] + cfg['PostTreat']['Volume']['Dir'],
# │         listscen    = listscen,
# │         reservoirs  = cfg['ReservoirRegions'],
# │         dirOUT      = cfg['dirOUT'],
# │         n_stages    = _MOSSHOOS_CFG.n_stages,
# │     )
# │ # ---- END MOSSHOOS B -------------------------------------------------------
# └─────────────────────────────────────────────────────────────────────────────

def inject_B(volume_dir: str,
             listscen:   list,
             reservoirs: list,
             dirOUT:     str,
             n_stages:   int = 12) -> dict:
    """
    INJECTION B — Haar DWT decomposition of the inflow/volume series.

    Reads Volume-Reservoir-{res}.csv files (already created by Plan4RES),
    concatenates across scenarios to build the long series,
    applies Haar DWT, extracts seasonal component A_J and residuals,
    fits Fourier parameterisation of theta(t).

    Returns
    -------
    dict with keys:
        inflow_series  : pd.Series  (monthly GWh, all scenarios concatenated)
        seasonal       : np.ndarray  A_J component
        residuals      : np.ndarray  stochastic residuals x - A_J
        wavelet        : HaarWavelet instance
        fourier        : dict with theta callable + R2
        J              : int  (optimal DWT levels)
    """
    import sys, os
    # import WP2 — look for it next to this file first, then in sys.path
    _here = os.path.dirname(os.path.abspath(__file__))
    if _here not in sys.path:
        sys.path.insert(0, _here)
    from mosshoos_wp2_wavelet import HaarWavelet, fourier_fit_seasonal

    # ── collect series from Volume files ─────────────────────────────────────
    all_vals = []
    for res in reservoirs:
        fpath = os.path.join(volume_dir, f'Volume-Reservoir-{res}.csv')
        if not os.path.isfile(fpath):
            logger.warning(f"[MOSSHOOS B] Volume file not found: {fpath}")
            continue
        df = pd.read_csv(fpath, index_col=0)
        # average across scenarios as the "representative" inflow proxy
        df.index = pd.to_datetime(df.index)
        mean_vol = df.mean(axis=1)
        # first-differences approximate inflow (volume change ≈ inflow - release)
        inflow_proxy = mean_vol.diff().fillna(0).clip(lower=0)
        all_vals.append(inflow_proxy.values)

    if not all_vals:
        logger.warning("[MOSSHOOS B] No volume files found — using synthetic fallback")
        from mosshoos_wp1_data import load_inflows
        series = load_inflows('')
    else:
        combined = np.mean(all_vals, axis=0)
        dates    = pd.date_range('1994-01', periods=len(combined), freq='MS')
        series   = pd.Series(np.maximum(combined, 0.01), index=dates,
                              name='inflow_GWh')

    x = series.values.astype(float)

    # ── Haar DWT ──────────────────────────────────────────────────────────────
    J       = HaarWavelet.optimal_levels(x, eta_min=55.0)   # 55% for real data
    wavelet = HaarWavelet(J).fit(x)
    wavelet.print_summary()
    seasonal  = wavelet.seasonal(len(x))
    residuals = wavelet.residuals(x)
    fourier   = fourier_fit_seasonal(seasonal, dt=1/12., K=5)

    logger.info(f"[MOSSHOOS B] J* = {J}  A_J energy = {wavelet.energy_[0]:.1f}%")
    logger.info(f"[MOSSHOOS B] Fourier R² = {fourier['R2']:.4f}")
    logger.info(f"[MOSSHOOS B] Residual σ = {residuals.std():.2f}")

    # save decomposition for audit
    out_df = pd.DataFrame(
        {'inflow': x, 'seasonal': seasonal, 'residual': residuals},
        index=series.index)
    out_df.to_csv(os.path.join(dirOUT, 'MOSSHOOS_WP2_decomposition.csv'),
                  float_format='%.3f')
    logger.info(f"[MOSSHOOS B] Saved MOSSHOOS_WP2_decomposition.csv")

    return dict(inflow_series=series, seasonal=seasonal,
                residuals=residuals, wavelet=wavelet,
                fourier=fourier, J=J)


# ══════════════════════════════════════════════════════════════════════════════
# INJECTION C  —  Block-bootstrap scenario set η (WP3)
# ══════════════════════════════════════════════════════════════════════════════

# ┌─────────────────────────────────────────────────────────────────────────────
# │ PASTE IN PostTreatPlan4res.py  AT LINE ~1559
# │ (replaces or precedes the:  if cfg['usevu'] and cfg['PostTreat']['Volume']['draw']:  block)
# │
# │ # ---- MOSSHOOS INJECTION C ------------------------------------------------
# │ if cfg.get('usevu_stochastic', False):
# │     from mosshoos_injection import inject_C
# │     _MOSSHOOS_ETA, _MOSSHOOS_WEIGHTS, _MOSSHOOS_TYPES = inject_C(
# │         wavelet_result = _MOSSHOOS_WAVELET,
# │         cfg_m          = _MOSSHOOS_CFG,
# │         dirOUT         = cfg['dirOUT'],
# │     )
# │     logger.info(f"[MOSSHOOS] Scenario set |eta| = {len(_MOSSHOOS_ETA)}")
# │ # ---- END MOSSHOOS C -------------------------------------------------------
# └─────────────────────────────────────────────────────────────────────────────

def inject_C(wavelet_result: dict,
             cfg_m: MOSSHOOSConfig,
             dirOUT: str) -> tuple:
    """
    INJECTION C — Block-bootstrap + K-means scenario generation (WP3).

    Takes the wavelet decomposition from inject_B and produces the
    stochastic scenario set eta = {(eta^k, w_k), k=1..K} that replaces
    the deterministic BV_VU centile loop in Plan4RES.

    Returns
    -------
    eta     : (K, n_stages) array of representative scenarios
    weights : (K,) probability weights
    types   : (K,) str array  'DRY' | 'median' | 'WET'
    """
    import os, sys
    _here = os.path.dirname(os.path.abspath(__file__))
    if _here not in sys.path:
        sys.path.insert(0, _here)
    from mosshoos_wp3_scenarios import ScenarioGenerator, validate_scenarios

    seasonal  = wavelet_result['seasonal']
    residuals = wavelet_result['residuals']
    x         = wavelet_result['inflow_series'].values.astype(float)

    gen = ScenarioGenerator(
        n_bootstrap  = cfg_m.n_bootstrap,
        n_clusters   = cfg_m.n_clusters,
        block_size   = 4,
        perturb_frac = 0.25,
    )
    eta, weights, types = gen.fit(seasonal, residuals,
                                   n_stages=cfg_m.n_stages)

    # validation
    val = validate_scenarios(gen, x, n_stages=cfg_m.n_stages)
    logger.info(f"[MOSSHOOS C] KS stat={val['ks_stat']:.4f}  p={val['ks_p']:.4f}")
    logger.info(f"[MOSSHOOS C] ACF max err={val['acf_max_err']:.4f}")
    logger.info(f"[MOSSHOOS C] Energy score improvement={val['es_improvement']:.1f}%")

    # save
    cols = [f'stage_{t+1}' for t in range(cfg_m.n_stages)]
    df_eta = pd.DataFrame(eta,
        index  =[f'scenario_{k+1}' for k in range(len(eta))],
        columns=cols)
    df_eta['weight'] = weights
    df_eta['type']   = types
    df_eta.to_csv(os.path.join(dirOUT, 'MOSSHOOS_WP3_scenarios_eta.csv'),
                  float_format='%.4f')
    logger.info(f"[MOSSHOOS C] Saved MOSSHOOS_WP3_scenarios_eta.csv")

    # save validation
    val_out = {k: v for k, v in val.items()
               if not isinstance(v, np.ndarray)}
    pd.Series(val_out).to_csv(
        os.path.join(dirOUT, 'MOSSHOOS_WP3_validation.csv'), header=False)

    return eta, weights, types


# ══════════════════════════════════════════════════════════════════════════════
# INJECTION D  —  SDDP Benders cuts replacing deterministic BV (WP4-A)
# ══════════════════════════════════════════════════════════════════════════════

# ┌─────────────────────────────────────────────────────────────────────────────
# │ PASTE IN PostTreatPlan4res.py  AT LINE ~1517
# │ (replaces the max_{l in L_t} block that computes valuesmsIni/valuesmsFin)
# │
# │ The ORIGINAL Plan4RES code reads:
# │     valuesmsIni = -10000000000000
# │     valuesmsFin = -10000000000000
# │     for row in aIni.index:
# │         valIni = bIni.loc[row] + (aIni.loc[row]*VolSMS['VolumeIni']).sum()
# │         if valIni > valuesmsIni: valuesmsIni = valIni
# │     for row in aFin.index:
# │         valFin = bFin.loc[row] + (aFin.loc[row]*VolSMS['VolumeFin']).sum()
# │         if valFin > valuesmsFin: valuesmsFin = valFin
# │
# │ MOSSHOOS wraps this:
# │
# │ # ---- MOSSHOOS INJECTION D ------------------------------------------------
# │ if cfg.get('usevu_stochastic', False):
# │     from mosshoos_injection import inject_D
# │     valuesmsIni, valuesmsFin, _MOSSHOOS_SDDP = inject_D(
# │         cfg_m       = _MOSSHOOS_CFG,
# │         eta         = _MOSSHOOS_ETA,
# │         weights     = _MOSSHOOS_WEIGHTS,
# │         VolSMS      = VolSMS,
# │         dirOUT      = cfg['dirOUT'],
# │         scen        = scen,
# │         listscen    = listscen,
# │     )
# │ else:
# │     # Original Plan4RES deterministic code (unchanged):
# │     valuesmsIni = -10000000000000
# │     valuesmsFin = -10000000000000
# │     for row in aIni.index:
# │         valIni = bIni.loc[row] + (aIni.loc[row]*VolSMS['VolumeIni']).sum()
# │         if valIni > valuesmsIni: valuesmsIni = valIni
# │     for row in aFin.index:
# │         valFin = bFin.loc[row] + (aFin.loc[row]*VolSMS['VolumeFin']).sum()
# │         if valFin > valuesmsFin: valuesmsFin = valFin
# │ # ---- END MOSSHOOS D -------------------------------------------------------
# └─────────────────────────────────────────────────────────────────────────────

# module-level cache so SDDP is only trained once per variant/year
_SDDP_CACHE = {}

def inject_D(cfg_m:    MOSSHOOSConfig,
             eta:      np.ndarray,
             weights:  np.ndarray,
             VolSMS:   pd.DataFrame,
             dirOUT:   str,
             scen:     int,
             listscen: list) -> tuple:
    """
    INJECTION D — SDDP replaces the deterministic Bellman value evaluation.

    On the FIRST scenario in listscen: trains the SDDP solver, saves cuts.
    On subsequent scenarios: reuses the already-trained solver (faster).

    Evaluates:
        V_stoch(S_ini) = SDDP upper bound at S_ini
        V_stoch(S_fin) = SDDP upper bound at S_fin

    These replace Plan4RES valuesmsIni / valuesmsFin.

    Returns
    -------
    valuesmsIni : float   (stochastic water value at start)
    valuesmsFin : float   (stochastic water value at end)
    sddp        : SDDPSolver instance (for inspection / plotting)
    """
    import os, sys
    _here = os.path.dirname(os.path.abspath(__file__))
    if _here not in sys.path:
        sys.path.insert(0, _here)
    from mosshoos_wp4_sddp import SDDPSolver

    global _SDDP_CACHE

    # ── train SDDP on first scenario only ────────────────────────────────────
    cache_key = id(cfg_m)
    if cache_key not in _SDDP_CACHE:
        logger.info("[MOSSHOOS D] Training SDDP solver (first call)...")
        sddp = SDDPSolver(
            S_max      = cfg_m.V_max,
            R_max      = cfg_m.R_max,
            D          = cfg_m.D_plan,
            C_def      = cfg_m.C_def,
            C_sp       = cfg_m.C_sp,
            C_hyd      = cfg_m.C_hyd,
            nu         = cfg_m.nu,
            gamma      = cfg_m.gamma,
            n_stages   = cfg_m.n_stages,
            n_iter     = cfg_m.n_iter,
            n_cuts_max = 60,
            gap_tol    = 0.005,
            n_grid     = 200,
        )
        sddp.fit(eta, weights, S0=cfg_m.S0)
        _SDDP_CACHE[cache_key] = sddp

        # run MC evaluation and save
        mc_df = sddp.simulate(eta, weights, n_mc=cfg_m.n_mc, S0=cfg_m.S0)
        mc_df.to_csv(os.path.join(dirOUT, 'MOSSHOOS_WP4_MC_MS.csv'),
                     index=False, float_format='%.3f')
        logger.info(f"[MOSSHOOS D] SDDP converged={sddp.converged_}  "
                    f"E[C^opex]={mc_df['cost'].mean():.1f}  "
                    f"E[deficit]={mc_df['deficit_gwh'].mean():.2f} GWh")
        logger.info(f"[MOSSHOOS D] Saved MOSSHOOS_WP4_MC_MS.csv")
    else:
        sddp = _SDDP_CACHE[cache_key]
        logger.info(f"[MOSSHOOS D] Reusing trained SDDP (scenario {scen})")

    # ── evaluate stochastic water value at scenario's start/end volumes ───────
    def V_stoch(S: float) -> float:
        """Evaluate stochastic V(S) = max over SDDP cuts at stage 0."""
        v = sddp.nu * S
        for (a, b) in sddp.cuts_[0]:
            v = max(v, a + b * S)
        return float(v)

    # VolSMS has columns 'VolumeIni' and 'VolumeFin' per reservoir
    S_ini = float(VolSMS['VolumeIni'].sum())
    S_fin = float(VolSMS['VolumeFin'].sum())

    valuesmsIni = V_stoch(S_ini)
    valuesmsFin = V_stoch(S_fin)

    logger.info(f"[MOSSHOOS D] scen={scen}  "
                f"S_ini={S_ini:.0f}  V_stoch(S_ini)={valuesmsIni:.2f}  "
                f"S_fin={S_fin:.0f}  V_stoch(S_fin)={valuesmsFin:.2f}")

    return valuesmsIni, valuesmsFin, sddp


# ══════════════════════════════════════════════════════════════════════════════
# INJECTION E  —  Monte Carlo WP5 full evaluation (optional, runs after loop)
# ══════════════════════════════════════════════════════════════════════════════

# ┌─────────────────────────────────────────────────────────────────────────────
# │ PASTE IN PostTreatPlan4res.py  AFTER THE SCENARIO LOOP (line ~1555)
# │
# │ # ---- MOSSHOOS INJECTION E ------------------------------------------------
# │ if cfg.get('usevu_stochastic', False) and '_MOSSHOOS_SDDP' in dir():
# │     from mosshoos_injection import inject_E
# │     inject_E(
# │         sddp       = _MOSSHOOS_SDDP,
# │         eta        = _MOSSHOOS_ETA,
# │         weights    = _MOSSHOOS_WEIGHTS,
# │         cfg_m      = _MOSSHOOS_CFG,
# │         dirOUT     = cfg['dirOUT'],
# │         dirIMG     = cfg['dirIMG'],
# │         listscen   = listscen,
# │     )
# │ # ---- END MOSSHOOS E -------------------------------------------------------
# └─────────────────────────────────────────────────────────────────────────────

def inject_E(sddp,
             eta:       np.ndarray,
             weights:   np.ndarray,
             cfg_m:     MOSSHOOSConfig,
             dirOUT:    str,
             dirIMG:    str,
             listscen:  list) -> None:
    """
    INJECTION E — Full WP5 Monte Carlo validation + comparison figures.

    Runs all four models (M-S, M-S0, M-D, M-D0), produces comparison figure,
    and saves summary CSV with all key manuscript numbers.
    """
    import os, sys
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    _here = os.path.dirname(os.path.abspath(__file__))
    if _here not in sys.path:
        sys.path.insert(0, _here)
    from mosshoos_wp4_sddp import SDDPSolver, compare_MS_MD, plot_results

    logger.info("[MOSSHOOS E] Running WP5 four-model comparison...")

    kw = dict(
        S_max=cfg_m.V_max, R_max=cfg_m.R_max, D=cfg_m.D_plan,
        C_def=cfg_m.C_def, C_sp=cfg_m.C_sp, C_hyd=cfg_m.C_hyd,
        nu=cfg_m.nu, gamma=cfg_m.gamma, n_stages=cfg_m.n_stages,
        n_iter=cfg_m.n_iter, n_cuts_max=60, gap_tol=0.005, n_grid=200,
    )

    # build seasonal series for M-S0 and M-D baselines
    from mosshoos_injection import _SDDP_CACHE
    # use the wavelet seasonal if available
    seasonal_approx = eta.mean(axis=0)   # fallback: mean of scenarios

    results = compare_MS_MD(
        eta=eta, weights=weights,
        seasonal=seasonal_approx,
        n_stages=cfg_m.n_stages,
        n_iter=cfg_m.n_iter,
        n_mc=cfg_m.n_mc,
        S0=cfg_m.S0,
        **{k: v for k, v in kw.items()
           if k not in ('n_stages','n_iter','n_cuts_max','gap_tol','n_grid')},
    )

    # produce figure
    fig_path = os.path.join(dirIMG, 'MOSSHOOS_WP4_SDDP_results.png')
    plot_results(results, eta, weights,
                 n_stages=cfg_m.n_stages,
                 save_path=fig_path)
    logger.info(f"[MOSSHOOS E] Figure saved: {fig_path}")

    # summary numbers
    c_ms  = results['MS'][1]['cost'].mean()
    c_md  = results['MD'][1]['cost'].mean()
    d_ms  = results['MS'][1]['deficit_gwh'].mean()
    d_md  = results['MD'][1]['deficit_gwh'].mean()
    saving= (c_md - c_ms) / max(c_md, 1) * 100

    summary = {
        'C_opex_MS_stochastic': round(c_ms, 1),
        'C_opex_MD_deterministic': round(c_md, 1),
        'cost_saving_pct': round(saving, 1),
        'deficit_MS_GWh': round(d_ms, 2),
        'deficit_MD_GWh': round(d_md, 2),
        'SDDP_converged': results['MS'][0].converged_,
        'n_scenarios': len(eta),
        'n_stages': cfg_m.n_stages,
        'V_max_GWh': cfg_m.V_max,
        'R_max_GWh_stage': cfg_m.R_max,
    }
    pd.Series(summary).to_csv(
        os.path.join(dirOUT, 'MOSSHOOS_summary.csv'), header=False)
    logger.info(
        f"[MOSSHOOS E] Cost saving M-S vs M-D = {saving:.1f}%  "
        f"Deficit reduction = {(d_md-d_ms)/max(d_md,0.001)*100:.1f}%")
    logger.info(f"[MOSSHOOS E] Saved MOSSHOOS_summary.csv")


# ══════════════════════════════════════════════════════════════════════════════
# SETTINGS TEMPLATE
# ══════════════════════════════════════════════════════════════════════════════

YAML_TEMPLATE = """
# ── Add these lines to settingsPostTreatPlan4res.yml to activate MOSSHOOS ──

usevu: true                    # keep existing Plan4RES Bellman value reading
usevu_stochastic: true         # activate MOSSHOOS stochastic extension

# MOSSHOOS parameters (all optional — defaults shown)
mosshoos_n_stages:   12        # 12 = monthly planning, 52 = weekly
mosshoos_n_iter:    200        # SDDP iterations (increase for production)
mosshoos_n_mc:      500        # Monte Carlo runs for WP5 validation
mosshoos_n_clusters: 10        # K-means scenario count |eta|
mosshoos_demand_gwh: null      # monthly demand [GWh]. null = auto (40% R_max)

# Results will appear in:
#   {dirOUT}/MOSSHOOS_WP2_decomposition.csv    <- wavelet components
#   {dirOUT}/MOSSHOOS_WP3_scenarios_eta.csv    <- scenario set eta
#   {dirOUT}/MOSSHOOS_WP3_validation.csv       <- KS, ACF, energy score
#   {dirOUT}/MOSSHOOS_WP4_MC_MS.csv            <- MC cost distribution
#   {dirOUT}/MOSSHOOS_summary.csv              <- all key numbers
#   {dirIMG}/MOSSHOOS_WP4_SDDP_results.png    <- comparison figure
"""

if __name__ == '__main__':
    print(YAML_TEMPLATE)
    print("\nTo activate: add 'usevu_stochastic: true' to your YAML config.")
    print("Then paste the 4 code blocks (marked INJECTION A-D) into")
    print("PostTreatPlan4res.py at the line numbers shown above each block.")
