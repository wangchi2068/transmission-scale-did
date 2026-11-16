# Experiment design — implementation of `theory/simulation_spec.md` (ADEMP) and the application

This file records how the code implements the pre-specified design and every place where it departs from
it, with the reason. Nothing here is a result; all numbers quoted are design quantities, timings or pilot
diagnostics (pilot outputs are in `logs/` and are labelled "not paper numbers").

## 1. Aims, DGM, estimands, methods, performance (as specified)

* **Aims / DGM**: 25 scenarios exactly as in the spec table (`../code/scenarios.py::build_scenarios`): core
  block 1–16 (adoption feedback/random × spillover yes/no × N ∈ {500, 2000} × T ∈ {12, 26}), null block
  17–20 (τ = 0, k ∈ {10, ∞}), misspecification block 21–25 (X1 wrong generation interval, X2 reporting drop
  after adoption, X3 heterogeneous τ_i, X4 susceptible depletion, X5 dynamic effects). Generator:
  `../code/dgp.py` (multiplicative renewal with network pressure, S = 2 lags, negative-binomial innovations
  with size k (Poisson for k = ∞); geometry fixed per N: uniform points, symmetrised 6-nearest-neighbour
  graph, G = 50 strip "states"; adoption calibration ḡ in `calibration.json` from `../code/calibrate.py`).
* **Estimands**: τ, κ, φ (generating values); event-time effects for X5; for X3 the per-replication switcher
  mean of τ_i and the first-order pseudo-true value of Proposition 3 (`aggregate.x3_targets`); ψ† per
  replication (true θ₀ in the ν recursion); population ψ from 5,000 common-random-number pairs per scenario.
* **Methods**: M0, M0-ET (E = 6), M0-OI (two-step GMM + Hansen J), M0′, M1, M1-BC, M2, M3, M4, Oracle,
  all hand-written in NumPy/SciPy (`../code/qdgmm.py`, `../code/comparators.py`). Failures are recorded per method
  (`ok`, `reason`) and never replaced by another method; convergence rates are reported with denominators.
* **Performance** (`../code/aggregate.py`): bias, empSE, model SE, relative SE error, RMSE, coverage of
  t_{G−1} intervals, CI width, rejection rate of H₀: θ = 0, convergence rate, each with MCSE (Morris, White
  & Crowther 2019); results JSON in the 03-results-protocol schema plus `metrics.csv`.
* **Random numbers**: seed = 20261112 + 10⁵·scenario + rep; two generator streams `default_rng([seed, 0])`
  (structure) and `default_rng([seed, 1])` (innovations); ψ-truth pairs use seed + 9·10⁸. All methods use the same data set within a replication.

## 2. Deviations from the spec (all decided before any full run)

| # | item | spec | implementation | reason / evidence |
|---|---|---|---|---|
| 1 | unit-effect spread | a_g, b_i ~ N(0, 0.2²), 0.3·π_g | multiplied by `ALPHA_SCALE = 0.25` | with the spec values and no depletion, persistent growth heterogeneity produced county-week counts up to ~10¹² at T = 26; pilot table `logs/pilot_alpha_tails.txt` |
| 2 | additional estimator **M0-G** | — | M0 with global root selection: profile the φ-equation on a grid φ ∈ [0, 3], refine every sign change, add the boundary candidate φ = 0, choose the smallest over-identified criterion Q (all moments incl. the M0-OI instruments); inference conditional on the chosen root, boundary choices flagged | the just-identified equations can have several roots when the true φ is at or near 0; M0 (local solver) can then settle on a spurious interior root (pilot, 30 reps: `logs/pilot/pilot_root_selection.txt`). M0 stays the pre-specified primary estimator; M0-G is reported alongside it |
| 3 | R | 1000 | `run_plan.json` (benchmark-based, ≥ 1000; waves 200 → 1000 → R) | step instruction: choose R for ≈ 4 h on 16 cores; see §4 |
| 4 | M1 offset | ω per unit | ω may be (T, N); `poisson_fe` subsets it in half panels (M1-BC) | bug fix: M1-BC previously failed on half panels |
| 5 | `poisson_fe` options | — | `with_kappa`, `phi_fixed` | needed to reproduce the N3 pseudo-true values in the tests; the simulation uses the default (κ, φ free) |
| 6 | M0-H, Callaway–Sant'Anna | optional | not run | optional in the spec; not needed for any pre-specified criterion |

### 2b. Post-hoc deviation (decided in run-sim-application, *after* interim results were seen)

| # | item | spec | implementation | reason / evidence |
|---|---|---|---|---|
| 7 | additional estimator **M0-GI** (post hoc) | — | M0-G with a different selection criterion: same candidate roots (grid φ ∈ [0, 3], sign-change refinement, boundary φ = 0), but Q is built from the M0 moments plus the φ-instrument evaluated at other reference values, N/Λ(0.05) and N/Λ(2.0) (`qdgmm.fit_global(selection="inv")`, `qdgmm.spec_selection_inv`). Run by `../code/simulate_supp.py` on the *same seeds and data sets* as the main run (rows `sim/rows_m0gi.jsonl`, merged by (sid, rep) in `../code/aggregate.py`, data fingerprint checked); measures use the replications it reached (`scenarios[sid].supp_R`) | interim R = 200 results showed M0 and the pre-specified M0-G under-covering τ in every no-spillover core scenario (M0 cov. 0–38 %, M0-G 46–78 %). Diagnosis (a local diagnostic script, not distributed): no code defect; with φ0 = 0 the M0 equations have a spurious root that tracks φ̄ (global identification at the boundary, flagged as unproven in theory). Instruments at other reference values are valid at the true root for any reference value but not at a φ̄-dependent spurious root. A 40-rep diagnostic (local, not distributed) motivated the choice. **M0 remains the pre-specified primary estimator and M0-G the pre-specified remedy; both are reported unchanged. M0-GI is reported separately and labelled post hoc everywhere (Table 8, figures 1–2, application table).** |

## 3. Application (`../code/application.py`, data in `data_contract.md`)

Main specification: M0 with Ferretti generation interval, adjacency network, ω = pop/10⁵, φ̄ = 0.5, all
comparators (M0-G, M0-OI with J, M0′, M0-ET, M1, M1-BC, M2, M3, M4), φ-profile, ψ for M0 and M0-G.
Pre-specified sensitivities (each with M0 and M0-G): Nishiura serial interval; w₁ grid {0.6, …, 1.0};
commuting network; policy lagged one week; ω × 0.25 and × 4; φ̄ ∈ {0.25, 1}; excluding Kansas and Texas
(county opt-outs not coded); lead-indicator regression.

Implementation choices fixed here (not results):
* **Generation interval**: daily distribution discretised to weeks by averaging over a uniform infection
  day within the week; lags 1–2 kept and renormalised; dropped same-week and > 2-week mass reported.
* **ψ by horizon**: besides the 26-week totals, ψ is reported cumulatively over the first h ∈ {4, 8, 13}
  adoption weeks, because the ν recursion compounds exp(−κ D̄) through renewal and contagion and its
  total can be dominated by the last weeks (seen in the smoke run; interpretation is left to the results
  step).
* **Lead placebo**: under adoption feedback the lead coefficient is not zero even when the model is
  correct (Proposition 2(c) and numeric check N1), so the lead regression is reported for comparison only;
  the specification test of the paper is Hansen's J from M0-OI.
* **Figures**: Times New Roman is obtained by `fetch_fonts.sh` (Microsoft core-fonts installer, verified by
  sha256, unpacked with 7-Zip, never executed) because the style kit refuses substitutes.

## 4. Compute plan (`../code/benchmark.py` → `run_plan.json`)

Benchmark on the 4-core build pod (5 timed reps + 20 ψ pairs per scenario, 4 concurrent workers):
≈ 17 core-seconds per replication of all 25 scenarios, ≈ 1,700 core-seconds for the 25 × 5,000 ψ pairs,
peak worker RSS ≈ 0.2 GB. R is the largest value on 1000, 1500, …, 8000 whose projected wall time on
16 cores (85 % parallel efficiency) is ≤ 75 % of 4 h; the grid is capped at 8000 because aggregate.py holds
all rows in memory (~18 kB per row). Plan: see `run_plan.json` (R, projected hours, memory). Coverage MCSE
at p = 0.95 is √(0.95·0.05/R).

Priority: wave 1 brings every scenario (core feedback+spillover → core feedback no-spill → core random →
null → misspecification, cheaper N first) to 200 reps and its ψ truth; wave 2 to 1000; wave 3 to R. A run
stopped early therefore leaves a balanced design. Resume: rerun the same command.

## 5. Verification done in this step

* `python -I ../theory/numeric_checks.py`: reproduces the Web Appendix B checks (S1, S3, S4, S5, S7, S8)
  that are distributed with the paper. The development-time pytest suite that cross-checked these against
  the implementation was not part of the distribution; the checks retained here are self-contained and
  deterministic.
* Smoke run (`smoke/`, not paper numbers): kill −9 of the whole process group during the simulation,
  truncated line, resume, aggregate/tables/figures with partial scenarios; resumed vs uninterrupted runs
  agree to 3e-15 relative.
