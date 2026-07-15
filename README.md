# MOSSHOOS × Plan4RES

**Stochastic seasonal storage valuation for the Senegal River / WAPP system.**

A stochastic extension of the [Plan4RES](https://github.com/plan4res) seasonal-storage
framework, combining Haar wavelet decomposition, block-bootstrap scenario generation,
stochastic dual dynamic programming, and Hamilton–Jacobi–Bellman control — validated on
55 years of real OMVS inflow data.

S. B. Affognon · B. M. Ndiaye · P. Mendy · C. M. F. Kebe
UCAD Dakar · University of Nairobi · CTSS

---

## Install and run - one paste

**macOS / Linux** (VS Code terminal):
```bash
git clone https://github.com/belvinos91/mosshoos-plan4res.git && cd mosshoos-plan4res && bash install.sh
```

**Windows** (VS Code PowerShell terminal):
```powershell
git clone https://github.com/belvinos91/mosshoos-plan4res.git; cd mosshoos-plan4res; .\install.ps1
```

That is everything: it finds a usable Python, builds an isolated `.venv` (nothing on your
system is touched), installs the dependencies, runs all four unit-test suites, runs the
pipeline on the real OMVS data, builds `dashboard.html` and opens it.

Nothing to configure. No solver licence. No Plan4RES installation required.

> The repository is private, so `git clone` will ask you to sign in. Ask the author to add
> you as a collaborator first, or use the `dashboard.html` file directly - it is standalone
> and needs nothing installed at all.

## Quick start (manual)

```bash
git clone https://github.com/belvinos91/mosshoos-plan4res.git
cd mosshoos-plan4res
./setup.sh                                    # detect, install, verify
./run_all.sh --omvs data/OMVS_inflows.csv     # run the pipeline
python3 dashboard_server.py --port 8050       # open http://localhost:8050
```

**Plan4RES is optional.** MOSSHOOS needs no external solver — the SDDP and HJB are
self-contained, so everything here runs whether or not Plan4RES is installed and with or
without a commercial solver licence.

**A direct comparison against Plan4RES output is not yet possible.** Their seasonal
storage module (SSV) does not currently run on open-source-solver builds — see
[`plan4res_fixes/PLAN4RES_SDDP_BUG.md`](plan4res_fixes/PLAN4RES_SDDP_BUG.md). Until that
is resolved, the benchmarks here are internal (M-S vs M-S0 vs M-D vs M-D0) and **no claim
of superiority over Plan4RES is made.**

---

## Verify it yourself

Every component ships with unit tests that assert the mathematical properties that make
the results trustworthy:

```bash
python3 src/sddp_clean.py        # monotone lower bound, gap closes, VSS>=0, cuts valid
python3 src/sde_convergence.py   # Euler-Maruyama ~0.5, Milstein ~1.0 (path-coupled)
python3 src/hjb_clean.py         # solution finite, V U-shaped, policy monotone
python3 src/wavelet_galerkin.py  # Haar basis orthonormal, O(M^-1) convergence
```

Nothing in the results is hand-entered: `RESULTS.md` is **generated** from the computed
JSON by `write_results.py`, so it cannot drift from the numbers.

---

## What is compared

| Label | Method |
|-------|--------|
| **M-S** | MOSSHOOS — stochastic: wavelet + block-bootstrap + SDDP |
| **M-S0** | ablation — stochastic SDDP without the wavelet step |
| **M-D** | benchmark — our SDDP on a single mean-inflow scenario (deterministic-equivalent) |
| **M-D0** | baseline — greedy, no water value |

VSS = C(M-D) − C(M-S)

> **Important.** M-D is **not** Plan4RES output. Plan4RES's seasonal storage module (SSV)
> is itself SDDP - the binary is `sddp_solver` and its `bellmanvalues.csv` holds Benders
> cuts over scenarios. It is a *stochastic* method, not a deterministic one. M-D here is
> **our own SDDP fitted on a single mean-inflow scenario**: a deterministic-equivalent
> benchmark that isolates the value of representing uncertainty at all.
>
> A true MOSSHOOS-vs-Plan4RES comparison requires running their solver, which is currently
> blocked by a bug in their SDDP on open-source-solver builds
> (see [`plan4res_fixes/PLAN4RES_SDDP_BUG.md`](plan4res_fixes/PLAN4RES_SDDP_BUG.md)).
> Until that is resolved, **no claim of superiority over Plan4RES is made here.**
>
> The comparison that does test the MOSSHOOS contribution is **M-S vs M-S0** - our
> wavelet + block-bootstrap scenarios against plain scenarios, same solver. On the current
> runs that difference is ~0%, which is an open question we are investigating, not a result
> we are claiming.
. We report the comparison against **both** M-D and the greedy
baseline, because they answer different questions: M-D0 measures the value of water-value
operation at all; M-D measures the additional value of modelling uncertainty.

---

## Headline finding

On the real 55-year Manantali record (CV = 1.36), at the reservoir's actual scale:

- **M-S vs greedy: +7.3%** — the value of water-value-based operation
- **M-S vs deterministic: 0.0%** — no additional value from uncertainty *at this storage*

The full picture is a **curve**, not a number (`storage_sweep.py`):

| Storage (× mean monthly inflow) | vs greedy | vs deterministic |
|--------------------------------:|----------:|-----------------:|
| 0.1× | +1.5% | 0% |
| 1× | +4.5% | 0% |
| 3× (Manantali) | +13.3% | 0% |
| 4× | +17.1% | 0% |
| 6× | +17.0% | **+9.6%** |

**Interpretation.** When the reservoir binds frequently, the release is forced by physical
constraints and the water value has little discretion to exploit — so stochastic and
deterministic policies coincide. Stochastic value emerges only once storage is large
enough for genuine multi-period hedging. The value of stochastic seasonal storage
valuation is governed by the storage-to-inflow regime and the cascade topology, not by a
single headline number.

On the **Plan4RES toyDataset**, VSS = 0. That is correct: their dataset is a single
synthetic year with a seasonal cycle but no interannual variability. With no uncertainty
to hedge, the two policies coincide — which validates that the method does not manufacture
value where none exists.

See [RESULTS.md](RESULTS.md) for the full validated numbers and the defects found and
corrected along the way.

---

## Dashboard

Two ways, same content:

```bash
python3 dashboard_server.py --port 8050   # live: change parameters, re-run, download
python3 make_static.py                    # standalone dashboard.html, offline, no server
```

The live server exposes a small API (`/api/results`, `/api/run`, `/api/download/<file>`)
and lets you move sliders — iterations, scenarios K, Monte-Carlo paths, bootstrap cycles —
press **Run**, watch the solver log stream in, and download the new output. The static
export bakes the current results into a single file you can email.

Tabs: Overview · Plan4RES · MOSSHOOS method (equations) · Comparison · Scenarios ·
Wavelet · Continuous-time · Run models · Data & downloads.

---

## Method

**WP1** — monthly inflow series, from Plan4RES CSVs *or* real OMVS data.
Both sources funnel into one format, which is what makes the comparison like-for-like.

**WP2** — Haar wavelet split `x_t = A_J(t) + Σ_j D_j(t)`; Fourier fit of `A_J` gives the
seasonal drift `θ(t)` used by the SDE.

**WP3** — block-bootstrap (block b = 4 months, preserving autocorrelation) + K-means
reduction to K representative scenarios; validated by KS, ACF and energy score.

**WP4-A** — SDDP water values: `Q_t(S) = min_u { C(S,q,u) + γ E_q[Q_{t+1}(S')] }`, with
`Q_{t+1}(S') ≥ max_ℓ {a_ℓ + b_ℓ S'}`. Supporting hyperplanes accumulate, so the
approximation is a valid lower bound and improves monotonically.

**WP4-B** — Itô SDE for the inflow: `dQ = κ(θ(t) − Q)dt + σ₀√Q dW` (CIR keeps `Q > 0`
when Feller holds).

**WP4-C** — HJB control:

```
-∂V/∂t = min_u { C(S,Q,u) + (Q−u−sp)·∂V/∂S + κ(θ(t)−Q)·∂V/∂Q + ½σ²(Q)·∂²V/∂Q² }
```

Randomness lives in the inflow `Q`; storage evolves deterministically given `(u,Q)`. `V`
has kinks in `S` where the release regime switches, so `V ∈ H¹` — which is exactly why a
Haar wavelet basis attains the optimal `O(M⁻¹)` rate and a smoother basis cannot.

---

## Layout

```
├── setup.sh                  detect environment, install, verify
├── run_all.sh                one command: pipeline + dashboard
├── compare_plan4res.py       runs WP1-WP4 on Plan4RES inputs or OMVS data
├── storage_sweep.py          VSS vs storage-to-inflow ratio
├── write_results.py          generates RESULTS.md from the computed JSON
├── dashboard_server.py       live dashboard + re-run API (stdlib only)
├── dashboard_app.html        the dashboard
├── make_static.py            bake it into a standalone offline file
├── src/
│   ├── sddp_clean.py         validated SDDP water values  (+ unit tests)
│   ├── hjb_clean.py          validated HJB PDE solver     (+ unit tests)
│   ├── sde_convergence.py    path-coupled convergence     (+ unit tests)
│   ├── wavelet_galerkin.py   Haar reduced-basis HJB       (+ unit tests)
│   ├── hjb_integrate.py      runs the continuous-time strand from pipeline inputs
│   ├── plan4res_io.py        reads Plan4RES's own inputs into MOSSHOOS format
│   ├── mosshoos_wp1_data.py  … wp2_wavelet … wp3_scenarios … wp4_sddp … wp4_ito_hjb
│   └── PostTreatPlan4res_MOSSHOOS.py
├── data/OMVS_inflows.csv     real Manantali monthly inflow energy, 1961-2015
├── results/                  generated: comparison_results.json/.csv, storage_sweep.json
└── plan4res_fixes/
    ├── FIXES_APPLIED.md      MPI init, USE_MPI flag, ssv.sh paths
    └── PLAN4RES_SDDP_BUG.md  full diagnosis of the HiGHS "no linearization" failure
```

---

## Data

**Real OMVS Manantali** — 660 months (1961–2015) of daily discharge from the OMVS record,
converted to monthly energy with the documented factor 13.77324 MWh/day per m³/s.
Mean 102.6 GWh/month, CV = 1.36, including the Sahel droughts.

**Plan4RES toyDataset** — read directly from their `SS_SeasonalStorage.csv`,
`ZV_ZoneValues.csv`, hourly hydro series and `bellmanvalues.csv`, so the "same inputs"
claim is literal.

---

## On running Plan4RES itself

Plan4RES's own `sddp_solver` currently **fails on HiGHS-only builds** with
`SDDPOptimizer::oneStepBackward: no linearization is available`. We traced it: the inner
unit-commitment subproblem is solved as a MIP despite `intRelaxIntVars 1`; a MIP has no
duals; no duals means no linearization; SDDP cannot build a cut. It is reproducible on
their own toyDataset and affects every user without a commercial licence.

Full trace, including three further defects we fixed (MPI initialization, a missing
`USE_MPI` compile flag, and six bare-filename paths in `ssv.sh`), is in
[`plan4res_fixes/PLAN4RES_SDDP_BUG.md`](plan4res_fixes/PLAN4RES_SDDP_BUG.md).

This is **why M-D is reimplemented here** rather than invoked: it is not that we did not
run their solver, it is that their solver does not currently run without Gurobi. Our M-D
follows their documented method — our SDDP on a single mean-inflow scenario — on their own
input files.

---

## Requirements

Python ≥ 3.10 with numpy, pandas, scipy, scikit-learn. Nothing else.
`./setup.sh` checks and installs.

On srv-ortad use `/usr/bin/python3` (3.12); the `/usr/local` 3.10 build lacks `_ctypes`,
which breaks pandas.
