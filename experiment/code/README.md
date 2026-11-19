# experiment/ — code and data for the QD-GMM paper

Python 3.11; dependencies are pinned in `requirements.txt`. Data provenance, SHA-256 checksums and
German county varicella / school-holiday coding are in `../data/data_contract.md`.

## Layout

| path | role |
|---|---|
| `qdgmm.py` | **proposed estimator**: QD-GMM (M0), M0′, M0-ET, M0-OI (two-step + Hansen J), M0-G (global root selection), M0-GI (invariance root selection) |
| `psi.py` | counterfactual recursion ν(θ), cases averted ψ (treated / network) with delta-method SE |
| `comparators.py` | M1 Poisson FE, M1-BC split-panel jackknife, M2 log R DiD, M3 log(Y+1) TWFE, M4 pooled QMLE, Oracle |
| `dgp.py`, `scenarios.py` | simulation data-generating mechanism and the 25-scenario table |
| `simulate.py`, `aggregate.py` | resumable Monte Carlo runner and metrics aggregator → `../analysis/results.json` |
| `make_tables.py`, `make_figures.py`, `style_kit.py` | simulation LaTeX tables and publication figures |
| `data_prep.py` | verifies SHA-256 checksums of raw German files (`../data/raw/germany/`) → `../data/processed/panel.npz` |
| `application.py` | German varicella application (M0, M0-G, M0-GI, comparators, bootstraps, jackknife, Stock–Wright sets, event-time, sensitivities, NUTS-2, ascertainment) → `../analysis/germany_application_results.json` |
| `calibrated_sim_germany.py` | calibrated finite-cluster simulation on the 400-county × 261-week German panel ($R=500 \times 2$) → `../analysis/germany_calibrated_sim.json` |
| `make_figure_germany.py`, `make_tables_germany.py` | generates Figure 3 (`figures/fig_03_application.png/.pdf`) and all 8 application LaTeX tables |
| `../theory/numeric_checks.py` | numerical checks of the derivations in Web Appendix B |

## Reproducing the German Varicella Application

```bash
python experiment/data/fetch_germany.py --check   # or run without --check to re-download
python experiment/code/data_prep.py
python experiment/code/application.py
python experiment/code/calibrated_sim_germany.py
python experiment/code/make_figure_germany.py
python experiment/code/make_tables_germany.py
```
