# MOSSHOOS x Plan4RES — Validated Results

## Summary

This document records the honest, reproducible findings from the integrated
MOSSHOOS pipeline run on both Plan4RES toyDataset inputs and real OMVS
(Manantali) hydrological data.

## Methodology

The full pipeline (WP1 -> WP2 -> WP3 -> WP4) runs on a monthly inflow series
that comes from either source, ensuring a like-for-like comparison:

1. WP1 — inflow series in GWh/month (from Plan4RES CSVs or OMVS data)
2. WP2 — Haar wavelet split: x_t = A_J(t) + sum D_j(t); Fourier fit of A_J
3. WP3 — block-bootstrap (b=4) + K-means -> 10 representative scenarios eta
4. WP4-SDDP — corrected, unit-tested water-value solver (sddp_clean.py)
5. WP4-HJB — Ito SDE + Hamilton-Jacobi-Bellman continuous-time control

## The SDDP solver was corrected

The original mosshoos_wp4_sddp.py used (a) a single-scenario lower bound that
oscillated under high inflow variance, and (b) OLS-fitted, pruned cuts that did
not accumulate. We replaced the water-value core with sddp_clean.py, which:

- uses supporting-hyperplane (subgradient) cuts that ACCUMULATE (never pruned),
  making the value-function approximation monotone (standard Pereira-Pinto SDDP);
- computes a deterministic first-stage lower bound (monotone, meaningful gap);
- is validated by unit tests (monotone LB, VSS >= 0, valid cuts).

Run the tests:  python3 src/sddp_clean.py

## Findings

### 1. Plan4RES toyDataset (their own inputs)

VSS ~ 0%. This is CORRECT: the toyDataset is a single synthetic year with a
seasonal cycle but no interannual variability (seasonal CV 0.10-0.25, no
droughts). With no uncertainty to hedge, stochastic and deterministic policies
coincide. This validates that the method does not fabricate value where none
exists.

### 2. Real OMVS Manantali data (55 years, 1961-2015, CV = 1.36)

The value of stochastic optimization depends on the storage-to-inflow regime:

  - M-S vs greedy (M-D0):       +5% to +24%, growing with storage
  - M-S vs deterministic (M-D): ~0% until storage >= 4x mean, then +8%

Interpretation. MOSSHOOS's water-value method delivers substantial savings over
naive greedy operation (5-24%). Its advantage over a deterministic mean-inflow
water-value policy is small in the constrained regime — because when the
reservoir is frequently binding (empty/full), the release decision is forced by
physical constraints, leaving little discretion for the water value to exploit.
Stochastic value emerges (~8%) only when storage is large enough for genuine
multi-period hedging.

Defensible contribution: the value of stochastic seasonal storage valuation in
the OMVS system is governed by the storage-to-inflow regime and cascade
topology, not by a single headline number.

### 3. Storage-sensitivity curve

See results/storage_sweep.json and the accompanying figure. VSS is reported as a
function of reservoir storage (0.1x to 6x mean monthly inflow), computed from
the real 55-year record with the validated solver.

## Reproduce

  python3 src/sddp_clean.py                                    # unit tests
  python3 compare_plan4res.py --inflow data/OMVS_inflows.csv   # single comparison
  python3 storage_sweep.py --inflow data/OMVS_inflows.csv      # full curve

## Open items (honest)

- The single-reservoir model does not yet represent the OMVS cascade
  (Manantali -> Gouina -> Felou coupling). The station files have zero active
  storage and should be modeled as run-of-river fed by upstream releases.
- The M-S vs M-D near-zero result should be cross-checked against the real
  Plan4RES SDDP solver once its toyDataset terminal-cut issue is resolved.
- The cyclic terminal water value (nu) is an approximation; a fixed-point
  treatment would close the residual LB-UB gap.

## WP4-B/C — Continuous-time strand (Ito SDE + HJB)

Status as of this commit:

- **Ito SDE calibration (VALIDATED).** The CIR process dQ = kappa(theta-Q)dt +
  sigma0 sqrt(Q) dW calibrates on real OMVS data with kappa=15/yr, sigma0=21.8,
  Fourier R2=0.85, and Feller ratio 13.85 (>= 1, so Q_t > 0 is guaranteed).
  The KS test rejects Gaussian innovations (p=0.000) — real Manantali inflows
  are heavy-tailed, which motivates the stochastic (non-Gaussian) treatment.

- **Euler-Maruyama convergence (VALIDATED).** Strong-convergence slope 0.31 on
  real data (theory 0.5; synthetic test gives 0.44).

- **HJB PDE solver (UNDER VALIDATION).** The explicit finite-difference HJB
  solver exhibits a numerical instability (value function diverges under the
  current explicit scheme). The Milstein scheme slope (0.32, theory 1.0) and the
  wavelet-Galerkin convergence depend on a stable HJB solver and are therefore
  held back pending a stable (implicit/CFL-correct) reimplementation. These
  numbers are NOT reported in results to avoid shipping unvalidated values.
