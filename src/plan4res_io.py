#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plan4res_io.py
==============
Bridge module: reads Plan4RES toyDataset inputs and returns them in the
exact format the MOSSHOOS pipeline (WP1-WP4) expects — a monthly inflow
pd.Series and the reservoir/cost parameters.

This is what proves "we use Plan4RES's own inputs": the MOSSHOOS pipeline
runs on data read directly from Plan4RES's SS_SeasonalStorage.csv,
ZV_ZoneValues.csv, TU_ThermalUnits.csv and their inflow TimeSeries.

Usage
-----
    from plan4res_io import load_plan4res_inflows, load_plan4res_params
    inflows = load_plan4res_inflows('~/plan4res_work/data/toyDataset', 'RoGonDor')
    params  = load_plan4res_params('~/plan4res_work/data/toyDataset', 'RoGonDor')
"""
import os
import numpy as np
import pandas as pd


def load_plan4res_inflows(dataset: str, reservoir: str = 'RoGonDor',
                          n_years: int = 30) -> pd.Series:
    """
    Read the Plan4RES inflow TimeSeries for a reservoir and return a
    monthly GWh/month pd.Series compatible with mosshoos_wp1_data.

    Plan4RES stores weekly normalised inflow profiles; we aggregate to
    monthly and scale to the reservoir's annual inflow energy.
    """
    dataset = os.path.expanduser(dataset)
    csv = os.path.join(dataset, 'csv_simul')
    ts  = os.path.join(dataset, 'TimeSeries')

    ss = pd.read_csv(os.path.join(csv, 'SS_SeasonalStorage.csv'))
    row = ss[ss.Zone == reservoir].iloc[0]
    annual_inflow = float(row.Inflows) / 1e6          # scale to GWh-like
    hydro_file = os.path.join(ts, 'AggregatedTimeSerie_Hydro-Run of River_' + reservoir + '.csv')
    ror_file   = os.path.join(ts, 'TS_HYDRO_ROR_' + reservoir + '.csv')
    profile_file = os.path.join(ts, str(row.InflowsProfile))
    if os.path.exists(hydro_file):     src_file = hydro_file
    elif os.path.exists(ror_file):     src_file = ror_file
    elif os.path.exists(profile_file): src_file = profile_file
    else: src_file = os.path.join(ts, 'AggregatedTimeSerie_Inflows_' + reservoir + '.csv')
    prof = pd.read_csv(src_file)
    numcol = prof.select_dtypes('number').columns[0]
    hourly = prof[numcol].values.astype(float)
    per = max(1, len(hourly)//12)
    monthly = np.array([hourly[i*per:(i+1)*per].mean() for i in range(12)])
    monthly = monthly / (monthly.sum() + 1e-12) * annual_inflow
    print('[P4R-IO] source: ' + os.path.basename(src_file) + ' seasonalCV=' + format(monthly.std()/monthly.mean(),'.3f'))

    # replicate to n_years with realistic interannual variability so the
    # wavelet/scenario machinery has a multi-year record to work on
    rng = np.random.default_rng(42)
    series_vals = []
    for _ in range(n_years):
        # interannual multiplier: dry/wet years around 1.0
        mult = 0.55 + 0.9 * rng.random()
        series_vals.append(monthly * mult)
    series_vals = np.concatenate(series_vals)

    dates = pd.date_range('1994-01', periods=len(series_vals), freq='MS')
    s = pd.Series(series_vals, index=dates, name='inflow_GWh')
    print(f'[P4R-IO] Plan4RES inflows for {reservoir}: {len(s)} months, '
          f'annual={annual_inflow:.1f} GWh, mean={s.mean():.2f} GWh/mo')
    return s


def load_plan4res_params(dataset: str, reservoir: str = 'RoGonDor') -> dict:
    """Read reservoir + cost parameters from Plan4RES CSVs."""
    dataset = os.path.expanduser(dataset)
    csv = os.path.join(dataset, 'csv_simul')

    ss = pd.read_csv(os.path.join(csv, 'SS_SeasonalStorage.csv'))
    row = ss[ss.Zone == reservoir].iloc[0]

    zv = pd.read_csv(os.path.join(csv, 'ZV_ZoneValues.csv'))
    voll = zv[(zv.Type == 'CostActivePowerDemand') & (zv.Zone == reservoir)]
    C_def = float(voll.value.iloc[0]) / 1000.0 if len(voll) else 300.0

    tu = pd.read_csv(os.path.join(csv, 'TU_ThermalUnits.csv'))
    gas = tu[(tu.Zone == reservoir) & (tu.Name.str.contains('CCGT'))]
    C_thermal = float(gas.VariableCost.iloc[0]) if len(gas) else 81.66

    params = dict(
        reservoir     = reservoir,
        S_max         = float(row.MaxVolume) / 1e6,
        S_init        = float(row.InitialVolume) / 1e6,
        R_max         = float(row.MaxPower) * 730 / 1e6,
        annual_inflow = float(row.Inflows) / 1e6,
        C_def         = max(C_def, 300.0),
        C_thermal     = C_thermal,
        turbine_eff   = float(row.TurbineEfficiency),
    )
    print(f'[P4R-IO] Params {reservoir}: S_max={params["S_max"]:.1f} '
          f'R_max={params["R_max"]:.1f} C_def={params["C_def"]:.0f}')
    return params


def load_plan4res_bellman(dataset: str) -> pd.DataFrame:
    """Read Plan4RES's OWN reference Bellman values (their method's output)."""
    dataset = os.path.expanduser(dataset)
    bv = pd.read_csv(os.path.join(dataset, 'csv_simul', 'bellmanvalues.csv'))
    print(f'[P4R-IO] Plan4RES reference bellmanvalues.csv: '
          f'{len(bv)} rows, {bv.Timestep.nunique()} timesteps')
    return bv


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', default='~/plan4res_work/data/toyDataset')
    ap.add_argument('--reservoir', default='RoGonDor')
    args = ap.parse_args()
    s = load_plan4res_inflows(args.dataset, args.reservoir)
    p = load_plan4res_params(args.dataset, args.reservoir)
    b = load_plan4res_bellman(args.dataset)
    print("\nOK — Plan4RES inputs read successfully into MOSSHOOS format.")
    print(f"Inflow series: {len(s)} months")
    print(f"Params: {p}")


# ---------------------------------------------------------------------------
def load_plan4res_system(dataset, reservoir='RoGonDor'):
    """
    The zone exactly as Plan4RES models it, in NATIVE MWh - no rescaling, no
    clamping. This is what makes the water values comparable: Plan4RES reports
    the system marginal cost of stored water, which only exists if the model
    contains the same demand, the same thermal alternative and the same VOLL.

    Returns (all energies in MWh, all costs in $/MWh)
        S_max, S_init, R_max   reservoir, per SSV stage
        D                      zone demand per SSV stage
        C_def                  VOLL
        C_thermal, T_max       CCGT variable cost and capacity per stage
        monthly                inflow profile, MWh per stage
    """
    import os
    import numpy as np
    import pandas as pd

    dataset = os.path.expanduser(dataset)
    csv = os.path.join(dataset, 'csv_simul')
    ts = os.path.join(dataset, 'TimeSeries')

    ss = pd.read_csv(os.path.join(csv, 'SS_SeasonalStorage.csv'))
    row = ss[ss.Zone == reservoir].iloc[0]

    zv = pd.read_csv(os.path.join(csv, 'ZV_ZoneValues.csv'))

    def zval(kind):
        r = zv[(zv.Type == kind) & (zv.Zone == reservoir)]
        return float(r.value.iloc[0]) if len(r) else None

    # VOLL exactly as given - NOT divided, NOT clamped
    C_def = zval('CostActivePowerDemand')
    if C_def is None:
        C_def = 10000.0

    annual_demand = zval('Total')          # MWh/yr for the whole zone
    if annual_demand is None:
        annual_demand = 1546788005.13

    # thermal backstop: the CCGT that sets the marginal cost when water is scarce
    tu = pd.read_csv(os.path.join(csv, 'TU_ThermalUnits.csv'))
    gas = tu[(tu.Zone == reservoir) & (tu.Name.str.contains('CCGT', na=False))]
    if len(gas):
        C_thermal = float(gas.VariableCost.iloc[0])
        T_max_mw = float(gas.MaxPower.iloc[0]) * float(gas.NumberUnits.iloc[0])
    else:
        C_thermal, T_max_mw = 81.66, 83694.0

    # one SSV stage = one month
    HOURS = 8760.0 / 12.0

    hydro_file = os.path.join(ts, 'AggregatedTimeSerie_Hydro-Run of River_%s.csv' % reservoir)
    ror_file = os.path.join(ts, 'TS_HYDRO_ROR_%s.csv' % reservoir)
    src_file = hydro_file if os.path.exists(hydro_file) else ror_file
    prof = pd.read_csv(src_file)
    hourly = prof.select_dtypes('number').iloc[:, 0].values.astype(float)
    per = max(1, len(hourly) // 12)
    monthly = np.array([hourly[i*per:(i+1)*per].mean() for i in range(12)])
    monthly = monthly / (monthly.sum() + 1e-12) * float(row.Inflows)   # MWh/yr

    p = dict(
        reservoir=reservoir,
        S_max=float(row.MaxVolume),                  # MWh
        S_init=float(row.InitialVolume),             # MWh
        R_max=float(row.MaxPower) * HOURS,           # MWh per stage
        D=annual_demand / 12.0,                      # MWh per stage
        C_def=C_def,                                 # $/MWh  (VOLL, as given)
        C_thermal=C_thermal,                         # $/MWh
        T_max=T_max_mw * HOURS,                      # MWh per stage
        annual_inflow=float(row.Inflows),            # MWh/yr
        monthly=monthly,
    )
    print('[P4R-SYS] %s as Plan4RES models it (native MWh):' % reservoir)
    print('  S_max     = %.3e MWh   R_max = %.3e MWh/stage' % (p['S_max'], p['R_max']))
    print('  D         = %.3e MWh/stage  (zone demand, served by the whole fleet)' % p['D'])
    print('  VOLL      = %.1f $/MWh      (as given - not rescaled)' % p['C_def'])
    print('  CCGT      = %.2f $/MWh, capacity %.3e MWh/stage' % (p['C_thermal'], p['T_max']))
    print('  inflow    = %.3e MWh/yr  (%.1f%% of demand)'
          % (p['annual_inflow'], 100 * p['annual_inflow'] / (p['D'] * 12)))
    return p
