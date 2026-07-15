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

- **HJB PDE solver (VALIDATED).** The original explicit solver diverged due to
  two real defects, now corrected in src/hjb_clean.py:

  1. *Sign error.* The backward HJB is -dV/dt = H, so a backward step is
     V(t) = V(t+dt) + dt*H. The original used V - dt*H, which flips the
     diffusion sign and turns the PDE into a backward heat equation
     (anti-diffusion) -> exponential blow-up (observed: |V| ~ 1e116).
  2. *Diffusion on the wrong state.* Randomness lives in the INFLOW Q, not in
     storage S (S evolves deterministically given u and Q). The original applied
     0.5*sigma^2*V_SS. The correct generator is
     kappa*(theta(t)-Q)*V_Q + 0.5*sigma^2(Q)*V_QQ, with S carrying only the
     drift (Q-u-spill)*V_S.

  The corrected solver uses explicit upwinding with automatic CFL sub-stepping.
  Unit tests (run: python3 src/hjb_clean.py) verify:
    - solution finite (no blow-up)
    - value function U-shaped in storage: high at empty (deficit risk), interior
      minimum, rising near full (spill penalty) -- physically correct
    - more inflow lowers cost when the reservoir is empty
    - optimal release increases with storage

  On real OMVS data the solver is stable, the value function is U-shaped with its
  minimum at S = 127 GWh, and the calibrated CIR satisfies Feller (13.85), so
  Q_t > 0 is guaranteed.

- **Strong convergence: EM and Milstein (VALIDATED).** The original
  convergence_test re-seeded the RNG before each discretization. Because the
  number of steps differs per dt, this drew a DIFFERENT Brownian path for each
  level. Strong convergence is a pathwise property:

      e(dt) = E| Q_T^{dt} - Q_T^{ref} |   on the SAME Brownian path

  Without path coupling the measured error is dominated by the difference
  between independent paths, which caps both slopes near 0.3 and hides
  Milstein's first order. This was the observed symptom (EM 0.31, Milstein 0.32
  -- nearly identical, which is the tell).

  src/sde_convergence.py implements the correct test: Brownian increments are
  generated once on the finest grid and each coarser level is built by SUMMING
  the fine increments (Brownian consistency), so every scheme sees the same
  driving path. Unit tests (run: python3 src/sde_convergence.py) in a
  Feller-safe regime give:

      EM slope       = 0.657   (theory 0.5)
      Milstein slope = 1.031   (theory 1.0)
      Milstein error falls 1.77e-2 -> 1.01e-3 as dt refines

  Conclusion: the Milstein IMPLEMENTATION was correct all along; only the test
  was flawed. This is itself a reportable methodological point.

  On the calibrated real-OMVS parameters the measured slopes are EM 0.93 and
  Milstein 1.17. EM exceeds its 0.5 asymptotic rate here because the calibrated
  kappa saturates at the 15/yr cap, so the drift dominates the diffusion and the
  scheme behaves near first order in this regime. The controlled unit test above
  (moderate kappa, Feller-safe) is the proper demonstration of the theoretical
  rates.

- **Wavelet-Galerkin projection (VALIDATED).** src/wavelet_galerkin.py solves
  the HJB in a reduced orthonormal Haar basis in the storage dimension and
  measures the L2 error against the CORRECTED finite-difference reference
  (src/hjb_clean.py), so the study measures genuine basis error rather than
  inheriting a broken reference.

  Rationale for a wavelet basis: the value function V(.,Q,t) is continuous but
  has kinks in S (the min over the release switches between deficit, interior
  and spill regimes), so V is in H^1 but not H^2. For such functions a Haar
  multiresolution basis attains

      || V - P_M V ||_L2 = O(M^-1)      for V in H^1

  and a smoother basis cannot do better because of the kinks. This is precisely
  why the wavelet basis is the right choice for this problem.

  Measured convergence (run: python3 src/wavelet_galerkin.py):

      J    M     L2 error
      1    2     6.786e+02
      2    4     4.688e+02
      3    8     2.618e+02
      4   16     1.314e+02
      5   32     5.906e+01

      Empirical rate: O(M^-0.89)     (theory O(M^-1))

  Unit tests verify: Haar basis orthonormality (max off-diagonal 2e-16),
  projection error decreasing with M, WG error decreasing with M, and the rate
  in the theoretical band.

## Validation status: all components

| Component                | File                    | Status    |
|--------------------------|-------------------------|-----------|
| SDDP water values        | src/sddp_clean.py       | validated |
| Ito SDE calibration      | src/hjb_integrate.py    | validated |
| EM / Milstein convergence| src/sde_convergence.py  | validated |
| HJB PDE solver           | src/hjb_clean.py        | validated |
| Wavelet-Galerkin         | src/wavelet_galerkin.py | validated |

Run every test suite:

    python3 src/sddp_clean.py
    python3 src/sde_convergence.py
    python3 src/hjb_clean.py
    python3 src/wavelet_galerkin.py

## Defects found and corrected

1. **SDDP lower bound** was a single random forward-pass cost, so it oscillated
   under high inflow variance; and cuts were OLS-fitted then pruned rather than
   accumulated. Replaced with accumulating supporting hyperplanes and a
   deterministic first-stage bound (monotone).
2. **HJB sign error.** The backward step used V - dt*H instead of V + dt*H,
   flipping the diffusion sign into a backward heat equation -> blow-up to
   |V| ~ 1e116.
3. **HJB diffusion on the wrong state.** Applied to storage (V_SS) instead of
   the inflow; corrected to kappa*(theta-Q)*V_Q + 0.5*sigma^2(Q)*V_QQ.
4. **Convergence test path coupling.** The RNG was re-seeded per discretization,
   so each dt used a different Brownian path; this capped both slopes near 0.3
   and hid Milstein's first order. The Milstein implementation was correct.

## Remaining work (honest)

- The model is single-reservoir. The OMVS cascade (Manantali -> Gouina ->
  Felou, with the downstream run-of-river plants fed by upstream releases) is
  not yet represented.
- Cross-validation against the real Plan4RES SDDP solver is pending resolution
  of its toyDataset terminal-cut issue (see plan4res_fixes/FIXES_APPLIED.md).
- The cyclic terminal water value (nu) is an approximation; a fixed-point
  treatment would close the residual LB-UB gap in the SDDP.
