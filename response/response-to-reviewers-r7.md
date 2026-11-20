# Response to Seventh-Round Review (`re-review` on `49b1102`)

We sincerely thank the referee for the exceptionally thorough code-and-decomposition audit on `49b1102`. Every issue raised in Major Comments 1--6, Minor Comments 1--6, and the Plan B finite-cluster / multi-start / jackknife unification recommendations in Section 4.5(2) and Section 5 has now been implemented, verified in full re-runs, and synchronised across `main.pdf` (38 pages), `supporting_information.pdf` (50 pages), `paper_zh/main.pdf` (26 pages), `bibliography.bib`, and the replication suite.

---

## Response to Major Comments

### Major Comment 1: Known-Truth Comparator Claims & Mechanism Decomposition (`M2` Prior Pseudo-Counts)
- **Changes made (`sections/abstract.tex` p.1, `sections/application.tex` §4.2 p.19, `sections/discussion.tex` §5 p.22, `sections/appendix.tex` Web E p.33, `experiment/code/application.py`, `experiment/code/calibrated_sim_germany.py`, `tables/table_12_calibrated_sim.tex` Web p.41):**
  - **Abstract & §5:** Removed `"confirmed"` from `sections/abstract.tex` (now `"evaluated against design-calibrated simulation truth"`) and replaced `"flips sign"` in `§5` (`main.pdf` p.22) with `"remains negative across county and NUTS-2 scales where two-step EpiEstim shifts from $+0.064$ to $-0.102$"`.
  - **Empirical `+0.064` vs. simulated distribution (`§4.2`, `main.pdf` p.19):** Explicitly reported that simulated county `M2` has mean `-0.015` (SD `0.019`), placing the observed county `M2` contrast (`+0.064`) **`+4.13` simulated SDs** above the calibrated M2 mean (`+2.68` SDs in high noise), so calibrated working dynamics do not alone explain the observed sign reversal without additional mechanisms such as holiday reporting variation.
  - **Prior-pseudocount decomposition & sensitivity (`calibrated_sim_germany.py`, `application.py`, `Web Table 11(c)` Web p.41):** Added deterministic noise-free evaluation ($Y=\mu$ with realized $L$) and weaker Gamma priors directly into `calibrated_sim_germany.py` ($R=500 \times 2$) and `application.py`:
    - Under default EpiEstim $(a=1, b=5)$, noise-free `M2` already retains **`+0.212`** bias (`+0.213` high-noise), confirming that sampling noise and zero realizations add only **`+0.012`** (`+0.014` high-noise).
    - Weakening the prior to $(a=0.1, b=5)$ halves simulated county `M2` bias from `+0.224` to **`+0.099`** (`+0.114` high-noise), reduces NUTS-2 `M2` bias to **`-0.001`** (`+0.002` high-noise), and shifts the observed county `M2` estimate from `+0.064` to **`-0.093`** (SE `0.017`); ultra-weak $(a=0.01, b=100)$ reduces simulated bias further to **`+0.016`** (`+0.045` high-noise).
  - **M3 estimand clarification (`§4.2`, `Web Table 11(c)`):** Reported county `M3` bias explicitly (`+0.159` calibrated, `+0.164` high-noise) and noted that `M3` targets the $\log(Y+1)$ TWFE slope rather than $\tau+\kappa$.

### Major Comment 2: NUTS-2 `M0` Convergence Counts, Non-Oracle Multi-Starts & Global Root Check
- **Changes made (`calibrated_sim_germany.py`, `sections/application.tex` §4.2–4.3 p.19–20, `sections/appendix.tex` Web E–F p.33–34, `tables/table_12_calibrated_sim.tex` Web p.41):**
  - Replaced `beta_true` starting values on NUTS-2 with strictly non-oracle multi-starts (`NUTS_STARTS`) and added a full 9-point non-negative global root check (`X2.global_roots`, $\phi \ge 0$).
  - Explicitly reported in `§4.2` and `Web Table 11(c)` that unconstrained non-oracle NUTS-2 `M0` converges on **`449/500`** calibrated (`tau+kappa` bias `+0.013`, $\hat\phi$ mean `0.113`, bias `-0.212`) and **`446/500`** high-noise (`tau+kappa` bias `-0.008`, $\hat\phi$ mean `0.170`, bias `-0.156`) panels, with remaining unconstrained iterates failing by crossing $\Lambda(\phi)\le 0$.
  - On the non-negative global grid ($\phi\ge 0$), all **`500/500`** NUTS-2 panels resolve in both arms, with **`33.8%`** (calibrated) and **`25.0%`** (high-noise) selecting the $\phi=0$ boundary ($\tau+\kappa$ bias `+0.014` / `-0.008`, $\hat\phi$ mean `0.120` / `0.160`).
  - Restricted the `"100% unique interior convergence"` statement in `§4.3` and Web Appendix F.2 explicitly to **county-level `M0`**.

### Major Comment 3 & Section 4.5(2): Finite-Cluster Premise, Unified Jackknife, Table 3 Contrast, and Nested Pairs Simulation (Plan B)
- **Changes made (`calibrated_sim_germany.py`, `application.py`, `make_tables_germany.py`, `tables/table_envelope.tex` p.37, `tables/table_12_calibrated_sim.tex` Web p.41, `sections/application.tex` p.19–20, `sections/discussion.tex` p.22, `sections/appendix.tex` Web p.33):**
  - **Calibration premise:** Explicitly stated in `§4.3`, `Table 3` footnote, `Web Table 11` footnote, and Web Appendix E that all calibrated rates assume the maintained working model (Assumption A7 and parallel trends hold at $\tau=-0.300, \kappa=0.060, \phi=0.326$).
  - **Unified delete-one-state jackknife:** Unified `calibrated_sim_germany.py` and `application.py` so both drop any isolated county (`deg == 0`, i.e.\ Berlin `11000` when Brandenburg is omitted), use identical non-oracle fixed starts (`[[-0.3, 0.0, 0.3], [0.0, 0.0, 0.5]]`), and check fold convergence explicitly (`16/16` observed; **`500/500`** calibrated and high-noise). Unified jackknife coverage is **`0.970`** calibrated (SE ratio `1.142`) and **`0.960`** high-noise (SE ratio `1.126`), explicitly noted as slightly conservative ($\approx +13$--$14\%$ SE overstatement).
  - **Nested `pairs cluster bootstrap-t` simulation (Plan B):** Added a fast vectorized `_pairs_bootstrap_fast` routine inside `calibrated_sim_germany.py` (verified to `1.24e-13` against `WB.pairs_moment_bootstrap`) and ran $R=200$, $B=199$ panels in both arms (`MCSE \approx 0.020`, `198.9/199` mean inner convergence). Empirical coverage is **`0.915`** calibrated and **`0.945`** high-noise (`Web Table 11(a)` Web p.41, `§4.3` p.20).
  - **Discussion & Table 3 alignment:** Reconciled `§5` (`main.pdf` p.21–22) and `Table 3` (`main.pdf` p.37) to highlight `WCR-score` (`0.934`--`0.946`) and unified `jackknife` (`0.960`--`0.970`) as the primary finite-cluster checks alongside `pairs cluster bootstrap-t` (`0.915`--`0.945`) as a calibrated supplementary check. Added `vs.\ M0 0.94--0.97` directly to Table 3's `M0-GI (post hoc)` weak-contagion row.

### Major Comment 4: Lag-1 Working Coding & Reporting Lag 2/3 Estimates in Main Text
- **Changes made (`sections/application.tex` §4.1 p.17–18, `paper_zh/main.tex`):**
  - Replaced the biological justification for lag 1 with an explicit statement that lag-1 coding is a post hoc baseline working specification, noting varicella's 10--21 day (mean $\approx$14--16 day) incubation window plus reporting delay implies biological transmission effects could align with lags 2--3.
  - Reported lag 2 (`-0.331`, `95%` CI `[-0.714, 0.052]`), lag 3 (`-0.541`, `[-0.822, -0.260]`), BW/BY-excluded lag 2 (`+0.033`), and all-holidays lag 2 (`+0.060`) directly in `§4.1` (`main.pdf` p.18), explaining how this sensitivity reinforces the strictly illustrative scope.

### Major Comment 5: Citation Accuracy, `Adda (2016)` Section IV, `Rieck et al. (2020)` Vaccination Data & `Elenev et al. (2025)`
- **Changes made (`sections/application.tex` p.17 & p.21, `sections/discussion.tex` p.22, `sections/appendix.tex` Web F.7, `bibliography.bib`, `paper_zh/bibliography.bib`):**
  - Separated general school-holiday transmission studies (`cauchemez2008school`, `jackson2014varicella`) from staggered German summer break studies (`isphording2021schools`, `vonbismarckosten2022schools`).
  - Rewrote the discussion of `\citet{adda2016economic}` (`§4.1`, `main.pdf` p.17) to accurately highlight Section IV's dynamic spatial transmission model (3-week varicella lag, cross-regional spillover, Table VI Panel A) alongside Section III's outcome regressions, contrasting QD-GMM's forward quasi-differencing and predetermined feedback instruments.
  - Replaced `suzuki2023varicella` in the vaccination sentence (`§4.4`, `main.pdf` p.21) with **`rieck2020impfquoten`** (2015--2019 contemporaneous RKI-KV and school-entry uptake surveillance, DOI `10.25646/7027`), retaining `siedler2010varicella` and `reuss2010varicella` specifically for early programme rollout and East--West uptake history, and noting whether state-level uptake trends interacting with season are absorbed by common week effects; kept `suzuki2023varicella` strictly in `§5` for temperature/holiday seasonality.
  - Rewrote `elenev2025staggered` in `§5` (`main.pdf` p.22) as paralleling cross-unit staggered policy spillovers rather than separating force of infection, and explicitly justified the log-$R$ scale relative to `roth2023parallel`.

### Major Comment 6: Canonical Re-Download Checksum & OpenHolidays Licence
- **Changes made (`experiment/data/fetch_germany.py`, `.gitignore`, `experiment/data/data_contract.md`, `sections/backmatter.tex` p.23, `sections/appendix.tex` Web D p.26–27):**
  - Removed volatile `extracted_date`/`archive_baseline_date` keys from `rki_varicella_2014_2019.json` (writing run metadata to `fetch_metadata.json`, gitignored), sorted `population_fallback` keys canonically, and matched `ensure_ascii=True, indent=1` in `fetch_holidays`.
  - Fixed the checksum-mismatch warning message to avoid falsely blaming SurvStat retroactive corrections.
  - **Live verification:** `python experiment/data/fetch_germany.py --verify-download` re-downloaded all three raw files from RKI SurvStat, OpenHolidays, and GitHub and matched the archived LF SHA-256 checksums **`3/3 OK`** (`2,475,214` / `261,944` / `1,289,615` bytes).
  - Added the OpenHolidays MIT / open data licence explicitly to Data Availability (`main.pdf` p.23), Web Appendix D (`Web p.26–27`), and `paper_zh/main.tex`.

---

## Response to Minor Comments

1. **Discussion phrasing (`§5` p.22):** Qualified WCR-score coverage gains with the maintained-working-model premise and adjusted `"preserves a negative working contrast"` to `"remains negative"`.
2. **Dynamic generation of Table 3 (`experiment/code/make_tables_germany.py`):** Added `write_table_envelope()` reading dynamically from `experiment/analysis/results.json`, `verification/root_selection/full_summary_tsummary.json`, `verification/same_week/results.json`, and `germany_calibrated_sim.json`.
3. **§4.1 root-selection phrasing (`§4.1` p.18):** Clarified `"unique interior root"` with an immediate qualification noting `M0-G` selects the $\phi=0$ boundary root (`9.113` vs. `9.120`).
4. **Chinese companion draft (`paper_zh/main.pdf`, 26 pages):** Synchronised all Major Comments 1--6 and Plan B simulation figures cleanly, fixed `vonbismarckosten2022schools` placement in F.7, and removed `成对 bootstrap` `"essential"` inconsistency.
5. **Chinese overfull box:** Replaced `\texttt{[ReportingDate]...}` with `\path{...}` in `paper_zh/main.tex`; `paper_zh/main.log` now has **`0` overfull boxes**.
6. **BibTeX DOIs:** Added verified DOIs across `bibliography.bib` and `paper_zh/bibliography.bib` (`adda2016economic`, `vonbismarckosten2022schools`, `isphording2021schools`, `siedler2010varicella`, `reuss2010varicella`, `rieck2020impfquoten`, `suzuki2023varicella`, `roth2023parallel`, `elenev2025staggered`).

---

## Verified Repository & Build Checks

- `python experiment/data/fetch_germany.py --check`: **`3/3 OK`**.
- `python experiment/data/fetch_germany.py --verify-download`: **`3/3 OK`** live upstream match against archived SHA-256 (`0735f021...`, `4f51fd93...`, `8a0d0bf2...`).
- `python -u experiment/code/application.py`: `16/16` unified jackknife convergence, `999/999` pairs convergence, plus `M2_prior_sensitivity` (`372.2s`).
- `python -u experiment/code/calibrated_sim_germany.py`: `500/500` county `M0`, `500/500` unified jackknife (`0.970` / `0.960`), `449/500` & `446/500` unconstrained NUTS-2 `M0` (`500/500` global $\phi\ge 0$), `200/200` nested pairs (`R=200, B=199`: `0.915` / `0.945`), plus `M2` weak/noise-free decomposition (`2575.51s`, 16 workers).
- `python experiment/code/make_tables_germany.py`: regenerated all **10** tables dynamically (`table_envelope.tex` + `table_7..18`).
- `python scripts/scan_control_chars.py`: **`total 0`**; `qa_checks.py`: abstract **`205 <= 225`** words, **`0` unresolved `??`** across `main.pdf` (38 pp; `30 / 25` counted), `supporting_information.pdf` (50 pp), and `paper_zh/main.pdf` (26 pp); `UNCITED: []`, `MISSING: []`; **`0` undefined references/citations** and **`0` overfull boxes** across all three `.log` files.
