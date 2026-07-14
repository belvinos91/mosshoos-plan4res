# MOSSHOOS × Plan4RES — Stochastic Seasonal Storage Valuation

Stochastic extension of the [Plan4RES](https://github.com/plan4res) seasonal-storage
framework for the Senegal / WAPP power system, combining **Haar wavelet decomposition**,
**block-bootstrap scenario generation**, **stochastic dual dynamic programming (SDDP)**,
and **Hamilton–Jacobi–Bellman (HJB) control**.

**Authors:** S. B. Affognon, B. M. Ndiaye, P. Mendy, C. M. F. Kebe
(UCAD Dakar · University of Nairobi · CTSS)

---

## What this repository proves

The donor/supervisor requirement is: **use the same Plan4RES inputs to produce
comparable outputs.** This repository does exactly that — the full MOSSHOOS
pipeline (WP1–WP4) runs on data read **directly from Plan4RES's own toyDataset
CSVs**, and produces a side-by-side comparison of:

| Label | Method | |
|-------|--------|---|
| **M-S**  | MOSSHOOS | stochastic: wavelet + block-bootstrap + SDDP |
| **M-S0** | ablation | stochastic SDDP **without** the wavelet step |
| **M-D**  | Plan4RES | deterministic mean-inflow Benders cuts |
| **M-D0** | baseline | greedy (no water value) |

The **Value of Stochastic Solution** is VSS = C(M-D) − C(M-S).

---

## Quick start — verify it yourself

```bash
# 1. install dependencies
pip install -r requirements.txt          # or: pip install --break-system-packages -r requirements.txt

# 2. run on Plan4RES toyDataset inputs (default path ~/plan4res_work/data/toyDataset)
./run_all.sh

# 3. or point to an explicit dataset / reservoir
./run_all.sh /path/to/toyDataset RoGonDor

# 4. or run on real OMVS data
./run_all.sh --omvs data/OMVS_inflows.csv

# 5. open the result
xdg-open dashboard.html        # or just open dashboard.html in a browser
```

`run_all.sh` does two things:
1. `compare_plan4res.py` — runs WP1–WP4, writes `results/comparison_results.json`
2. `build_dashboard.py` — regenerates `dashboard.html` **from that JSON**

**Every number in the dashboard comes from the JSON the pipeline just computed.
Nothing is hardcoded** — re-run and everything updates. This is what makes the
results verifiable.

---

## Repository layout

```
mosshoos-plan4res/
├── run_all.sh               # one-command entry point
├── compare_plan4res.py      # runs MOSSHOOS on Plan4RES inputs, writes results JSON
├── build_dashboard.py       # generates dashboard.html from results JSON (no hardcoding)
├── dashboard.html           # generated comparison dashboard (open in browser)
├── requirements.txt
├── results/                 # comparison_results.json + .csv (generated)
├── src/                     # the MOSSHOOS scientific pipeline
│   ├── mosshoos_wp1_data.py       # WP1  data loading & validation
│   ├── mosshoos_wp2_wavelet.py    # WP2  Haar DWT + Fourier seasonal fit
│   ├── mosshoos_wp3_scenarios.py  # WP3  block-bootstrap + K-means scenarios
│   ├── mosshoos_wp4_sddp.py       # WP4-A SDDP (M-S vs M-D vs M-D0)
│   ├── mosshoos_wp4_ito_hjb.py    # WP4-B/C Itô SDE + HJB + wavelet-Galerkin
│   ├── mosshoos_injection.py      # Plan4RES PostTreat injection points
│   ├── mosshoos_run_all.py        # standalone full MOSSHOOS pipeline
│   ├── PostTreatPlan4res_MOSSHOOS.py  # patched Plan4RES post-treatment
│   └── plan4res_io.py             # reads Plan4RES toyDataset -> MOSSHOOS format
└── plan4res_fixes/          # fixes to make Plan4RES's own solver build & run
    └── FIXES_APPLIED.md
```

## Method (equations)

**WP2 — Haar wavelet:**  x_t = A_J(t) + Σ_j D_j(t)

**WP2 — Fourier seasonal:**  θ(t) = a₀ + Σ_k [a_k cos(2πkt) + b_k sin(2πkt)]

**WP3 — scenarios:**  block-bootstrap (block b=4) + K-means → {(η^k, w_k)}

**WP4-A — SDDP Bellman cuts:**  V_t(S) = max_ℓ { b_ℓ + a_ℓ·S }

**WP4-B — Itô SDE (OU / CIR):**
  dQ = κ(μ−Q)dt + σ dW   /   dQ = κ(θ−Q)dt + σ√Q dW

**WP4-C — HJB:**
  ρJ = sup_u { R + (Q−u)J_V + κ(μ−Q)J_Q + ½σ² J_QQ }
  u*(V,Q) = Π[0,u_max]( (b − J_V)/c )

## Key finding

On Plan4RES's **low-variability toy data**, deterministic planning is nearly
adequate (small VSS), but the greedy baseline is clearly worse — confirming the
water-value mechanism. On **real OMVS data** (Senegal River, CV ≈ 1.02, severe
droughts), stochastic optimisation delivers a substantial cost reduction.
**Stochastic seasonal storage valuation matters precisely in high-variability
systems like the Senegal River basin.**

## Requirements

Python ≥ 3.10 with numpy, pandas, scipy, scikit-learn, matplotlib.
On the srv-ortad server use `/usr/bin/python3` (3.12) — the /usr/local python3.10
lacks the `_ctypes` module.
