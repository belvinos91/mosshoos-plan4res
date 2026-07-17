#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_water_values.py
=======================
THE comparison: real Plan4RES water values vs MOSSHOOS water values, same
reservoir, same inputs, both solvers actually run.

This is possible only because the SMS++ overload-hiding bug is fixed
(see plan4res_fixes/0001-fix-MILPSolver-overload-hiding.patch). Before that,
Plan4RES's SDDP aborted on any build without a commercial solver licence, so
there was no Plan4RES output to compare against at all.

What is compared
----------------
Both solvers approximate the cost-to-go by Benders cuts

    V_t(S)  >=  b + a . S

so the MARGINAL WATER VALUE of reservoir i at stage t is  -a_i  [$/MWh]:
the shadow price of one more MWh of stored energy. It is the same object in
both codes, in the same units, so the curves can be laid on one axis.

Plan4RES  : results_simul/cuts.txt      columns Timestep,a_0,a_1,a_2,b
            a_0=DolAmroth, a_1=Harad, a_2=RoGonDor  (order of SS_SeasonalStorage.csv)
MOSSHOOS  : src/sddp_clean.py cuts       Cut.b is the slope, so wv = -Cut.b

Usage
-----
    python3 compare_water_values.py \
        --p4r-cuts ~/plan4res_work/data/toyDataset/results_simul/cuts.txt \
        --dataset  ~/plan4res_work/data/toyDataset \
        --reservoir RoGonDor
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'src'))


# -- which column is our reservoir? ------------------------------------------
def reservoir_column(dataset, reservoir):
    """Cut columns a_0..a_n follow the row order of SS_SeasonalStorage.csv."""
    ss = pd.read_csv(os.path.join(os.path.expanduser(dataset),
                                  'csv_simul', 'SS_SeasonalStorage.csv'))
    zones = list(ss.Zone)
    if reservoir not in zones:
        sys.exit(f'{reservoir} not in {zones}')
    return zones.index(reservoir), zones


# -- Plan4RES: marginal water value per timestep -----------------------------
def plan4res_water_values(cuts_file, col):
    """
    For each timestep take the binding (most informative) cut and read off
    -a_col. Rows with all-zero slopes are the trivial initial cut and carry no
    water-value information, so they are skipped.
    """
    df = pd.read_csv(os.path.expanduser(cuts_file))
    acol = f'a_{col}'
    out = []
    for t, g in df.groupby('Timestep'):
        nz = g[(g[[c for c in df.columns if c.startswith('a_')]] != 0).any(axis=1)]
        if len(nz) == 0:
            continue                       # only the trivial cut at this stage
        # most binding = steepest slope for this reservoir
        row = nz.loc[nz[acol].abs().idxmax()]
        out.append(dict(timestep=int(t),
                        water_value=float(-row[acol]),   # $/MWh
                        intercept=float(row.b),
                        n_cuts=int(len(g))))
    return pd.DataFrame(out).sort_values('timestep').reset_index(drop=True)


# -- MOSSHOOS: marginal water value per stage --------------------------------
def mosshoos_water_values(dataset, reservoir, n_iter, n_clusters, n_bootstrap):
    from plan4res_io import load_plan4res_inflows, load_plan4res_params
    from mosshoos_wp2_wavelet import HaarWavelet
    from mosshoos_wp3_scenarios import ScenarioGenerator
    from sddp_clean import SDDP

    inflows = load_plan4res_inflows(dataset, reservoir)
    params = load_plan4res_params(dataset, reservoir)
    x = inflows.values.astype(float)

    J = HaarWavelet.optimal_levels(x)
    wv = HaarWavelet(J).fit(x)
    seasonal = wv.seasonal(len(x))
    residuals = wv.residuals(x)

    gen = ScenarioGenerator(n_bootstrap=n_bootstrap, n_clusters=n_clusters,
                            block_size=4)
    eta, weights, _ = gen.fit(seasonal, residuals, n_stages=12)

    D = 0.9 * float(inflows.mean())
    sol = SDDP(S_max=params['S_max'], R_max=params['R_max'], D=D,
               C_def=params['C_def'], n_stages=12, n_iter=n_iter,
               verbose=False).fit(eta, weights, params['S_init'])

    rows = []
    for t in range(12):
        if not sol.cuts_[t]:
            continue
        # steepest accumulated cut = the binding one
        c = max(sol.cuts_[t], key=lambda c: abs(c.b))
        rows.append(dict(stage=t,
                         water_value=float(-c.b),      # per storage unit
                         intercept=float(c.a),
                         n_cuts=len(sol.cuts_[t])))
    return pd.DataFrame(rows), params, sol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--p4r-cuts', default=os.path.expanduser(
        '~/plan4res_work/data/toyDataset/results_simul/cuts.txt'))
    ap.add_argument('--dataset', default=os.path.expanduser(
        '~/plan4res_work/data/toyDataset'))
    ap.add_argument('--reservoir', default='RoGonDor')
    ap.add_argument('--n_iter', type=int, default=150)
    ap.add_argument('--n_clusters', type=int, default=10)
    ap.add_argument('--n_bootstrap', type=int, default=150)
    ap.add_argument('--out', default=os.path.join(HERE, 'results',
                                                  'water_value_comparison.json'))
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    print('=' * 70)
    print('  REAL Plan4RES water values  vs  MOSSHOOS water values')
    print('  same reservoir, same inputs, both solvers actually run')
    print('=' * 70)

    col, zones = reservoir_column(args.dataset, args.reservoir)
    print(f'\nReservoir order in SS_SeasonalStorage.csv: {zones}')
    print(f'{args.reservoir} -> cut column a_{col}')

    # ---- Plan4RES ---------------------------------------------------------
    p4r = plan4res_water_values(args.p4r_cuts, col)
    print(f'\n[Plan4RES] {os.path.basename(args.p4r_cuts)}: '
          f'{len(p4r)} stages with informative cuts')
    print(f'  water value  min={p4r.water_value.min():.2f}  '
          f'max={p4r.water_value.max():.2f}  mean={p4r.water_value.mean():.2f}  $/MWh')

    # ---- MOSSHOOS ---------------------------------------------------------
    print('\n[MOSSHOOS] running SDDP on the same Plan4RES inputs...')
    mos, params, sol = mosshoos_water_values(args.dataset, args.reservoir,
                                             args.n_iter, args.n_clusters,
                                             args.n_bootstrap)
    print(f'  {len(mos)} stages, converged={sol.converged_}')
    print(f'  water value  min={mos.water_value.min():.4g}  '
          f'max={mos.water_value.max():.4g}  mean={mos.water_value.mean():.4g}')

    # ---- side by side -----------------------------------------------------
    # Both horizons are one annual cycle; compare on a normalised [0,1) axis.
    p4r_x = p4r.timestep / max(p4r.timestep.max(), 1)
    mos_x = mos.stage / 12.0
    mos_on_p4r = np.interp(p4r_x, mos_x, mos.water_value)

    print('\n' + '=' * 70)
    print('  MARGINAL WATER VALUE  ($/MWh for Plan4RES; solver units for MOSSHOOS)')
    print('=' * 70)
    print(f'  {"t":>3} {"phase":>6}  {"Plan4RES":>12}  {"MOSSHOOS":>12}')
    for i, r in p4r.iterrows():
        print(f'  {int(r.timestep):>3} {p4r_x[i]:>6.2f}  '
              f'{r.water_value:>12.2f}  {mos_on_p4r[i]:>12.4g}')

    # correlation of the SHAPE (the physically meaningful comparison; the two
    # solvers work in different storage units, so levels are not comparable)
    if len(p4r) > 2 and np.std(mos_on_p4r) > 0:
        rho = float(np.corrcoef(p4r.water_value, mos_on_p4r)[0, 1])
        rs = pd.Series(p4r.water_value).corr(pd.Series(mos_on_p4r), method='spearman')
    else:
        rho, rs = float('nan'), float('nan')

    print('\n' + '-' * 70)
    print(f'  Pearson  correlation of the two water-value profiles : {rho:.3f}')
    print(f'  Spearman correlation (shape/ranking)                 : {rs:.3f}')
    print('-' * 70)
    print('  Levels are NOT directly comparable: Plan4RES works in MWh of stored')
    print('  energy, MOSSHOOS in its own scaled storage unit. The SHAPE - when')
    print('  water is cheap and when it approaches VOLL - is the comparison that')
    print('  carries physical meaning.')

    out = dict(
        reservoir=args.reservoir, cut_column=f'a_{col}', zone_order=zones,
        plan4res=dict(source=args.p4r_cuts,
                      solver='Plan4RES SDDP (SMS++) on HiGHS, patched',
                      note='real output; possible only after the overload-hiding fix',
                      stages=p4r.to_dict('records'),
                      wv_min=float(p4r.water_value.min()),
                      wv_max=float(p4r.water_value.max()),
                      wv_mean=float(p4r.water_value.mean())),
        mosshoos=dict(solver='MOSSHOOS SDDP (src/sddp_clean.py)',
                      converged=bool(sol.converged_),
                      params={k: (round(v, 3) if isinstance(v, float) else v)
                              for k, v in params.items()},
                      stages=mos.to_dict('records'),
                      wv_min=float(mos.water_value.min()),
                      wv_max=float(mos.water_value.max()),
                      wv_mean=float(mos.water_value.mean())),
        alignment=dict(plan4res_x=[float(v) for v in p4r_x],
                       mosshoos_interp=[float(v) for v in mos_on_p4r],
                       pearson=rho, spearman=float(rs)),
    )
    json.dump(out, open(args.out, 'w'), indent=2)
    print(f'\nWrote {args.out}')


if __name__ == '__main__':
    main()
