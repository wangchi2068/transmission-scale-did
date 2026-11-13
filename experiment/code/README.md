# experiment/ — code and data for the QD-GMM paper

Run everything with `bash run_all.sh` from `experiment/` (stages, environment variables, priority and resume
are documented in its header). Design and deviations: `experiment_design.md`; data: `data_contract.md`;
compute plan: `run_plan.json`. Python ≥ 3.10, `requirements.txt`.

## Layout

| path | role |
|---|---|
| `src/qdgmm.py` | **proposed estimator**: QD-GMM (M0), M0′, M0-ET, M0-OI (two-step + Hansen J), M0-G (global root selection); cluster sandwich, t_{G−1} intervals |
| `src/psi.py` | counterfactual recursion ν(θ), cases averted ψ (treated / network) with delta-method SE |
| `src/comparators.py` | M1 Poisson FE (conditional likelihood), M1-BC split-panel jackknife, M2 Cori log R DiD, M3 log(Y+1) TWFE, M4 pooled QMLE, Oracle |
| `src/dgp.py`, `src/scenarios.py` | data-generating mechanism and the 25-scenario table; seeds; worker count |
| `src/simulate.py` | resumable multiprocessing Monte Carlo runner → `sim/rows.jsonl` |
| `src/aggregate.py` | performance measures + MCSE from whatever rows exist → `results.json`, `sim/metrics.csv` |
| `src/make_tables.py`, `src/make_figures.py`, `src/style_kit.py` | LaTeX tables and PDF/PNG figures from `results.json` (+ manifests) |
| `src/data_prep.py` | raw files → `data/processed/panel.npz` (run with `python3 -I`) |
| `src/application.py` | application: main analysis, comparators, sensitivities → `application/application_results.json` |
| `src/benchmark.py` | timing/memory per scenario → `run_plan.json` |
| `src/calibrate.py`, `src/pilot_alpha_tails.py` | pilot calibration (`calibration.json`) and the ALPHA_SCALE evidence (`logs/`) |
| `tests/` | pytest suite (`python3 -m pytest`; ~7 s) |
| `fetch_data.sh`, `fetch_fonts.sh` | data via the platform data library with sha256 checks; Times New Roman for figures |
| `smoke/` | tiny end-to-end run — **not paper numbers** |

## Quick reference

```bash
python3 -m pytest                                   # all tests
SMOKE=1 QDG_WORKERS=4 bash run_all.sh               # tiny end-to-end check into smoke/
QDG_WORKERS=16 bash run_all.sh                      # full study (R from run_plan.json), resumable
STAGES="agg tables figures" bash run_all.sh         # refresh outputs from the rows written so far
```

Worker count: `QDG_WORKERS`, otherwise the usable cores (CPU affinity / cgroup quota). Containers whose CPU
limit is invisible to the process (the 4-core build pod reports 36 CPUs) must set `QDG_WORKERS` explicitly.
BLAS/OpenMP threads are pinned to 1.

## Results files

* `results.json`: `task`, `dataset`, `metrics`, `seeds`, `statistics{std_ddof}`, `provenance{mode, label,
  run_complete, …}`, `summary{method:{metric:{cell:{mean, std, n_seeds}}}}` plus full per-scenario
  performance tables and the application block.
* `tables/manifest.json` lists each table with sha256; `figures/manifest.json` lists basenames only.
