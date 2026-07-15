# MOSSHOOS x Plan4RES - Validated Results

Generated from `results/comparison_results.json` and `results/storage_sweep.json`.
Every number below is produced by the pipeline; none is hand-entered.

## How to verify

```bash
python3 src/sddp_clean.py         # SDDP: monotone LB, gap closes, VSS>=0, valid cuts
python3 src/sde_convergence.py    # EM ~0.5, Milstein ~1.0 (path-coupled)
python3 src/hjb_clean.py          # HJB: stable, U-shaped V, monotone policy
python3 src/wavelet_galerkin.py   # Haar basis: O(M^-1) rate

./run_all.sh --omvs data/OMVS_inflows.csv
python3 storage_sweep.py --inflow data/OMVS_inflows.csv
python3 dashboard_server.py --port 8050      # live dashboard
```

## Data

- **Real OMVS Manantali**: Real OMVS data: data/OMVS_inflows.csv, 660 months (1961-2015), CV = 1.36,
  mean 102.6 GWh/month, converted from daily discharge with the documented
  factor 13.77324 MWh/day per m3/s.
- **Plan4RES toyDataset**: RoGonDor reservoir, read directly from their
  `SS_SeasonalStorage.csv`, `ZV_ZoneValues.csv`, hourly hydro series and
  `bellmanvalues.csv`.

## Headline result (OMVS-Manantali, S_max = 308 GWh, D = 92 GWh/mo)

| Model | Mean cost | P90 cost | Deficit (GWh) |
|-------|----------:|---------:|--------------:|
| M-S  MOSSHOOS stochastic | 338,723 | 416,282 | 178.4 |
| M-S0 stochastic, no wavelet | 338,764 | 416,282 | 178.4 |
| M-D  det.-equiv. deterministic | 338,749 | 416,282 | 178.4 |
| M-D0 greedy, no water value | 365,237 | 443,283 | 178.8 |

- **M-S vs greedy (M-D0): 7.3%** - the value of water-value-based operation
- **M-S vs deterministic-equivalent (M-D): 0.0%** - the additional value of modelling uncertainty

## Storage-sensitivity (the defensible contribution)

VSS as a function of usable storage, real 55-year record, corrected solver:

| S / mean_mo | S_max (GWh) | vs greedy % | vs determ. % |
|------------:|------------:|------------:|-------------:|
|  0.10 |      10.3 |      1.53 |      0.00 |
|  0.25 |      25.7 |      2.26 |      0.00 |
|  0.50 |      51.3 |      3.20 |      0.00 |
|  1.00 |     102.6 |      4.51 |     -0.11 |
|  1.50 |     153.9 |      6.09 |     -0.00 |
|  2.00 |     205.3 |      8.63 |      0.00 |
|  3.00 |     307.9 |     13.30 |      0.00 |
|  4.00 |     410.5 |     17.11 |     -0.00 |
|  6.00 |     615.8 |     17.01 |      9.58 |

- Against greedy the gain rises monotonically with storage, peaking at
  **17.1%** at 4.0x mean monthly inflow.
- Against the deterministic water-value policy the gain is ~0 across the
  operational range and reaches **9.6%** only at 6.0x.

**Interpretation.** When the reservoir is frequently binding the release is forced by
physical constraints, so the water value has little discretion to exploit and stochastic
and deterministic policies coincide. Stochastic value emerges only once storage is large
enough for genuine multi-period hedging. The value of stochastic seasonal storage
valuation in the OMVS system is therefore governed by the storage-to-inflow regime and
the cascade topology - not by a single headline number.

## Plan4RES toyDataset

VSS = 0. This is **correct**: their dataset is a single synthetic year (LotR zone names)
with a seasonal cycle (CV 0.10-0.25) but no interannual variability - no droughts, no
floods. With no uncertainty to hedge, stochastic and deterministic policies coincide.
This validates that the method does not manufacture value where none exists.

## Continuous-time strand

**Ito SDE (CIR) calibration.** kappa = 15.0 /yr, sigma0 = 21.769,
Fourier R2 = 0.8467, Feller ratio = 13.847
(Q_t > 0 guaranteed).
KS p = 0.0 - innovations are **heavy-tailed**.
Real Sahelian inflows are not Gaussian; droughts and floods produce fat tails. This is a finding, and it is what motivates the CIR treatment over a Gaussian assumption.

**Strong convergence** (path-coupled Brownian increments):
EM slope 0.929 (theory 0.5), Milstein slope 1.174 (theory 1.0).
On the calibrated parameters kappa saturates, so the SDE is drift-dominated and both
schemes approach first order; the controlled unit test (moderate kappa, Feller-safe) is
the proper demonstration of the 0.5 / 1.0 split and gives EM 0.657, Milstein 1.031.

**HJB PDE.** stable = True, U-shaped = True,
V minimum at S = 127.4 GWh
(V empty = 9382.5, V min = 1339.7, V full = 12105.8).
The U-shape is the physically correct signature: too little water causes deficits, too
much causes spill.

**Wavelet-Galerkin.** Empirical rate O(M^-1.029) against the corrected FD
reference; theory O(M^-1) for value functions in H1. V has kinks in S where the release
regime switches, so V is in H1 but not H2 - a smoother basis cannot beat O(M^-1), which
is exactly why the Haar basis is the right choice.

## Defects found and corrected

Six real defects were found and fixed. Each is documented because each would have
invalidated results under review.

1. **SDDP bound and cuts.** The lower bound was the cost of a single randomly sampled
   forward pass, so it oscillated wildly under high inflow variance; cuts were fitted by
   OLS through sampled points and then pruned, so the value function could get worse
   between iterations. Replaced with accumulating supporting hyperplanes (never pruned)
   and a deterministic first-stage lower bound. `src/sddp_clean.py`.

2. **HJB sign error.** The backward HJB is -dV/dt = H, so the backward step is
   V(t) = V(t+dt) + dt*H. The code used V - dt*H, flipping the diffusion sign into a
   backward heat equation (anti-diffusion) and diverging to |V| ~ 1e116.

3. **HJB diffusion on the wrong state.** Randomness lives in the inflow Q, not in storage
   S (S evolves deterministically given u and Q). The code applied 0.5*sigma^2*V_SS;
   corrected to kappa*(theta-Q)*V_Q + 0.5*sigma^2(Q)*V_QQ with S carrying only the drift
   (Q-u-spill)*V_S. `src/hjb_clean.py`.

4. **Convergence test path coupling.** The RNG was re-seeded before each discretization,
   so every dt level was driven by a *different* Brownian path. Strong convergence is
   pathwise, so the measured error was dominated by path differences: both slopes capped
   near 0.3 and Milstein's first order was hidden. Fine increments are now generated once
   and summed into coarser levels. The Milstein implementation was correct all along.
   `src/sde_convergence.py`.

5. **Reservoir parameters never reached the solver.** `compare_plan4res.py` called
   `compare_MS_MD` without S_max, R_max or D, so every run used the solver defaults
   (S_max=1500, R_max=180, D=95) regardless of the reservoir being studied. On RoGonDor
   this meant demand was 14x the inflow and every policy was forced to release everything,
   making all models identical. Fixed; the demand is now sized to the reservoir's
   deliverable energy on both branches (the zone total in ZV_ZoneValues.csv is met by the
   whole fleet, not by the reservoir, and must not be used as D).

6. **Bound scale mismatch.** The lower bound was discounted and included the terminal
   water value (both inherited through the cut recursion); the upper bound was an
   undiscounted sum of stage costs with no terminal term. The gap therefore could never
   close (a ~10.5% floor at gamma=0.98 over 12 stages) and could even go negative.
   Both bounds are now on the same scale. A hardcoded sigma=15 in the i.i.d. ablation
   was also scaled to the data.

## Remaining work

- **Cascade.** The model is single-reservoir. The OMVS cascade
  (Manantali -> Gouina -> Felou, downstream run-of-river plants fed by upstream releases)
  is not represented. The station parameter file assigns Manantali's head and capacity to
  the BafingMakana *station* with zero active storage, which deletes the reservoir from the
  model; this must be corrected before cascade results are reported.
- **Plan4RES cross-validation.** Pending resolution of the upstream
  `oneStepBackward: no linearization is available` issue on their stock toyDataset
  (see `plan4res_fixes/FIXES_APPLIED.md`).
- **Cyclic terminal value.** nu is an approximation; a fixed-point treatment would close
  the residual ~4% LB-UB gap.
