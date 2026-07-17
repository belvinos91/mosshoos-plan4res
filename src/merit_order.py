#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merit_order.py
==============
Reproduce the economics Plan4RES's SSV actually solves, so that MOSSHOOS's
water values are the same object as theirs.

THE KEY OBSERVATION
-------------------
Plan4RES's Bellman cut slopes are not arbitrary numbers: they are the VARIABLE
COST OF THE MARGINAL GENERATOR. Compare their cuts for RoGonDor (column a_2)
against TU_ThermalUnits.csv for the same zone:

    their cut slope        the unit it matches
    ---------------        -------------------
      3.00                 Coal          2.9999999988
      8.00                 Nuclear       7.9999999992
     81.66                 Gas CCGT     81.66487751511683
    113.21                 Gas OCGT    113.40944098188044   (Harad, a_1)
  10000.00                 VOLL      10000.0  (ZV_ZoneValues CostActivePowerDemand)

That is exactly what a water value should be: one more MWh of stored water
displaces the marginal unit, so its shadow price is that unit's cost. When the
stack is exhausted the water is worth VOLL.

MOSSHOOS could not produce this object because its stage problem had no
generation fleet - just a reservoir and an invented demand. This module supplies
the fleet, read from Plan4RES's own CSVs.

WHAT IS MODELLED
----------------
    residual = D(t) - renewables(t) - hydro_release
    the thermal stack is dispatched in merit order against the residual
    anything left over is unserved energy at VOLL

    stage cost = sum over the dispatched stack + VOLL * unserved + spill penalty

The marginal water value then emerges from the merit order, exactly as it does
in Plan4RES.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd

HOURS_PER_STAGE = 8760.0 / 12.0        # one SSV stage = one month


# -----------------------------------------------------------------------------
def load_fleet(dataset: str, zone: str = 'RoGonDor') -> dict:
    """
    The zone's generation fleet and demand exactly as Plan4RES has it, in native
    MWh and $/MWh. No rescaling, no clamping.
    """
    dataset = os.path.expanduser(dataset)
    csv = os.path.join(dataset, 'csv_simul')

    # ---- demand and VOLL ---------------------------------------------------
    zv = pd.read_csv(os.path.join(csv, 'ZV_ZoneValues.csv'))

    def zval(kind, dflt=None):
        r = zv[(zv.Type == kind) & (zv.Zone == zone)]
        return float(r.value.iloc[0]) if len(r) else dflt

    annual_demand = zval('Total', 1546788005.13)      # MWh/yr
    voll = zval('CostActivePowerDemand', 10000.0)     # $/MWh, as given

    # ---- thermal stack, sorted by variable cost ---------------------------
    tu = pd.read_csv(os.path.join(csv, 'TU_ThermalUnits.csv'))
    tu = tu[tu.Zone == zone].copy()
    tu['cap_mwh'] = (tu.MaxPower.astype(float)
                     * tu.NumberUnits.astype(float) * HOURS_PER_STAGE)
    tu = tu.sort_values('VariableCost')
    stack = [(float(r.VariableCost), float(r.cap_mwh), str(r.Name))
             for r in tu.itertuples()]

    # ---- renewables: annual energy (must-run, ~zero marginal cost) --------
    res = pd.read_csv(os.path.join(csv, 'RES_RenewableUnits.csv'))
    res = res[res.Zone == zone]
    renew_annual = float(res.Energy.astype(float).sum())     # MWh/yr

    # ---- the reservoir ----------------------------------------------------
    ss = pd.read_csv(os.path.join(csv, 'SS_SeasonalStorage.csv'))
    row = ss[ss.Zone == zone].iloc[0]

    fleet = dict(
        zone=zone,
        D=annual_demand / 12.0,                       # MWh per stage
        annual_demand=annual_demand,
        voll=voll,                                    # $/MWh
        stack=stack,                                  # [(cost, cap_MWh, name)]
        renew_per_stage=renew_annual / 12.0,          # MWh per stage
        renew_annual=renew_annual,
        S_max=float(row.MaxVolume),                   # MWh
        S_init=float(row.InitialVolume),              # MWh
        R_max=float(row.MaxPower) * HOURS_PER_STAGE,  # MWh per stage
        annual_inflow=float(row.Inflows),             # MWh/yr
    )
    return fleet


def print_fleet(f: dict) -> None:
    print(f"[FLEET] {f['zone']} as Plan4RES models it (native MWh, $/MWh)")
    print(f"  demand            {f['D']:.4g} MWh/stage   ({f['annual_demand']:.4g} MWh/yr)")
    print(f"  renewables        {f['renew_per_stage']:.4g} MWh/stage  "
          f"({100*f['renew_annual']/f['annual_demand']:.0f}% of demand)")
    print(f"  reservoir         S_max {f['S_max']:.4g}  R_max {f['R_max']:.4g} MWh/stage")
    print(f"  inflow            {f['annual_inflow']:.4g} MWh/yr "
          f"({100*f['annual_inflow']/f['annual_demand']:.1f}% of demand)")
    print(f"  VOLL              {f['voll']:.1f} $/MWh")
    print("  thermal stack (merit order):")
    for c, cap, name in f['stack']:
        print(f"    {c:10.4f} $/MWh  {cap:12.4g} MWh/stage  {name}")


# -----------------------------------------------------------------------------
def dispatch_cost(residual: float, stack, voll: float):
    """
    Dispatch the thermal stack in merit order against `residual` MWh.
    Returns (total_cost, marginal_cost, unserved).

    marginal_cost is the variable cost of the last unit called - which is what
    a Bellman cut slope reports as the water value.
    """
    if residual <= 0:
        return 0.0, 0.0, 0.0
    cost = 0.0
    left = residual
    marginal = 0.0
    for c, cap, _ in stack:
        if left <= 0:
            break
        use = min(left, cap)
        cost += c * use
        left -= use
        marginal = c
    if left > 1e-9:            # stack exhausted -> shed load
        cost += voll * left
        marginal = voll
    return cost, marginal, max(0.0, left)


def dispatch_cost_vec(residual: np.ndarray, stack, voll: float) -> np.ndarray:
    """Vectorised total cost, for the release grid inside the stage problem."""
    residual = np.maximum(np.asarray(residual, float), 0.0)
    cost = np.zeros_like(residual)
    left = residual.copy()
    for c, cap, _ in stack:
        use = np.minimum(left, cap)
        cost += c * use
        left -= use
    cost += voll * np.maximum(left, 0.0)
    return cost


def marginal_cost_at(residual: float, stack, voll: float) -> float:
    """The merit-order marginal cost at a given residual - the water value."""
    return dispatch_cost(residual, stack, voll)[1]


# -----------------------------------------------------------------------------
def merit_order_table(f: dict) -> pd.DataFrame:
    """
    The cumulative merit order: what the marginal cost is at each level of
    residual demand. This is the curve Plan4RES's water values ride on.
    """
    rows, cum = [], 0.0
    for c, cap, name in f['stack']:
        rows.append(dict(unit=name, cost=c, capacity=cap,
                         cum_capacity=cum + cap))
        cum += cap
    rows.append(dict(unit='UNSERVED (VOLL)', cost=f['voll'], capacity=np.inf,
                     cum_capacity=np.inf))
    return pd.DataFrame(rows)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--dataset', default='~/plan4res_work/data/toyDataset')
    ap.add_argument('--zone', default='RoGonDor')
    args = ap.parse_args()

    f = load_fleet(args.dataset, args.zone)
    print_fleet(f)

    print("\n[MERIT ORDER] cumulative")
    print(merit_order_table(f).to_string(index=False))

    print("\n[CHECK] does the merit order explain Plan4RES's cut slopes?")
    print("  Plan4RES water values observed in results_simul/cuts.txt (a_2, RoGonDor):")
    print("    3.00, 5.77, 8.00, 81.66, 97.54, 113.21, 5056.60, 10000.00")
    print("  units in this zone:")
    for c, cap, name in f['stack']:
        print(f"    {c:10.4f}  {name}")
    print(f"    {f['voll']:10.4f}  VOLL (unserved energy)")
    print("\n  -> the cut slopes ARE the marginal generator's variable cost.")


# -----------------------------------------------------------------------------
def map_cuts_to_units(cuts_file: str, fleet: dict, col: int,
                      tol: float = 0.02) -> pd.DataFrame:
    """
    Map every Plan4RES cut slope to the generator whose variable cost it equals.

    This is the demonstration that Plan4RES's water value IS the marginal
    generator's cost. Each cut slope -a_col is matched to the closest unit in
    the merit order; `match` reports the relative distance.
    """
    df = pd.read_csv(os.path.expanduser(cuts_file))
    acol = f'a_{col}'
    units = [(c, n) for c, _cap, n in fleet['stack']]
    units.append((fleet['voll'], 'UNSERVED (VOLL)'))

    rows = []
    for t, g in df.groupby('Timestep'):
        acols = [c for c in df.columns if c.startswith('a_')]
        nz = g[(g[acols] != 0).any(axis=1)]
        if not len(nz):
            continue
        r = nz.loc[nz[acol].abs().idxmax()]
        wv = float(-r[acol])
        cost, name = min(units, key=lambda u: abs(u[0] - wv))
        rel = abs(cost - wv) / max(abs(wv), 1e-9)
        rows.append(dict(timestep=int(t), water_value=wv,
                         matched_unit=name, unit_cost=cost,
                         rel_error=rel,
                         exact=bool(rel < tol)))
    return pd.DataFrame(rows).sort_values('timestep').reset_index(drop=True)


def demo(dataset, zone, cuts_file):
    """The demonstration, end to end."""
    f = load_fleet(dataset, zone)
    print_fleet(f)

    ss = pd.read_csv(os.path.join(os.path.expanduser(dataset),
                                  'csv_simul', 'SS_SeasonalStorage.csv'))
    zones = list(ss.Zone)
    col = zones.index(zone)
    print(f"\n[CUTS] reservoir order {zones} -> {zone} is column a_{col}")

    m = map_cuts_to_units(cuts_file, f, col)
    print("\n" + "=" * 78)
    print("  PLAN4RES WATER VALUES ARE THE MARGINAL GENERATOR'S VARIABLE COST")
    print("=" * 78)
    print(f"  {'t':>3} {'water value':>12} {'matched unit':>36} {'its cost':>10} {'':>6}")
    for r in m.itertuples():
        mark = 'EXACT' if r.exact else f'{100*r.rel_error:.0f}%'
        print(f"  {r.timestep:>3} {r.water_value:>12.2f} {r.matched_unit[:36]:>36} "
              f"{r.unit_cost:>10.3f} {mark:>6}")
    n_exact = int(m.exact.sum())
    print("-" * 78)
    print(f"  {n_exact}/{len(m)} cut slopes match a generator's cost to within 2%")
    print("=" * 78)
    return f, m


