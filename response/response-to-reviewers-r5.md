# Response to Fifth-Round Review (`re-review-ab6d07e.md`)

We thank the referee for the thorough fifth-round verification of commit `ab6d07e`, including exact line checks across `main.pdf`, `supporting_information.pdf`, `paper_zh/main.pdf`, and full re-execution of the replication scripts. Every factual point, methodological concern, empirical diagnostic, and reproducibility gap identified in Round 5 has been addressed directly below and in the revised manuscript (`main.pdf`, 38 pages total, 31 Biometrics-counted pages), Supporting Information (`supporting_information.pdf`, 45 pages, Web Appendices A--F and Web Tables 1--16), and Chinese companion (`paper_zh/main.pdf`, 26 pages).

---

## Overarching Disclosure & Strategic Repositioning

1. **Replacement of US COVID-19 Mask-Mandate Application (Honest Disclosure).**
   In Round 4, the empirical application used US county COVID-19 case counts and state mask mandates, where short generation times ($w_0 \gtrsim 0.3$) and immediate ascertainment/behavioural simultaneity violated the weekly domain boundary established in Section 3.2.1. Between Rounds 4 and 5, we replaced that COVID-19 dataset completely with a 2015--2019 German county-level weekly varicella (\emph{Windpocken}) panel paired with federal-state summer school holidays ($N=400$ counties, $G=16$ states, $T=261$ weeks, 115,129 cases), where varicella's 14-day mean serial interval yields negligible within-week transmission ($w_0 = 2.6 \times 10^{-6}$). We apologize for omitting an explicit disclosure of this empirical switch in `response-to-reviewers-r4.md`; that historical file is retained as an archive, and this Round 5 letter documents the complete provenance, protocol status, and post hoc formulation of the German application.

2. **Repositioning Section 4 as a Diagnostic Stress Test & Relaxing Forced Compression.**
   Rather than stretching the German holiday estimates into an unconfounded "21.2% causal reduction", we follow the referee's diagnostic findings directly: **Sections 2 and 3 remain the core mature contribution** (identification, forward quasi-differencing, martingale difference scores, and 25-scenario pre-specified Monte Carlo evidence), while **Section 4 is repositioned as a deliberate stress test** of both the estimator and its built-in failure-mode diagnostics (`sections/abstract.tex`, `sections/introduction.tex`, `sections/application.tex`, `sections/discussion.tex`, and `paper_zh/main.tex`). Instead of trimming safety caveats or compressing into an artificial 25-page envelope, the manuscript runs 31 counted pages (38 total PDF pages) with explicit scope boundaries, full protocol departures disclosure, and all sensitivity and placebo evidence intact.

---

## Point-by-Point Response to Major Comments

### Major Comment 1: Total-Scale Event-Time Onset Effect ($\hat\tau_0 + c_0\hat\kappa = -0.208$)
- **Issue.** On the identified policy-level total scale $\hat\tau_e + c_e\hat\kappa$, holiday onset week $e_0$ is $-0.208$ (SE $0.097$, $p=0.049$, 95\% CI $[-0.416,-0.001]$), not $-0.008$.
- **Action.** `experiment/code/application.py` now computes both direct $\hat\tau_e$ and total-scale $\hat\tau_e + c_e\hat\kappa$ (with within-bin average $\bar D_{it}$ factor $c_e$) plus the joint reopening test for $e_6, e_7$ ($p_F=0.172$). Web Table 14 (`tables/table_14_event_weights.tex`, SI p.39), Figure 3(b) (`figures/fig_03_application.png`), Section 4.3 (`main.pdf` p.20--21), and Web Appendix F.5 (SI p.38--39) now report $\hat\tau_0 + c_0\hat\kappa = -0.208$ verbatim and explicitly acknowledge that an immediate $-0.208$ drop in *Meldewoche* onset conflicts with pure post-incubation 14-day biological transmission, pointing to holiday clinic-visitation or notification compression.
- **Ascertainment Check (Web Table 15, SI p.40).** In addition to lagged-$D$ tipping points ($\rho^*_{\mathrm{CI}}=9.4\%$, $\rho^*_0=17.9\%$), we added contemporaneous-$D^{(0)}$ onset-moment matching: eliminating the contemporaneous onset drop requires $\rho^*_{\mathrm{on}}=21.5\%$, which drives the baseline lag-1 M0 estimate to $+0.007$ (SE $0.063$) and leaves lag 2 at $-0.213$ (SE $0.162$, covering zero).

### Major Comment 2: Non-Rotating Southern States, Seasonality, and Vaccination Uptake
- **Issue.** Bavaria (BY) and Baden-W\"urttemberg (BW) never rotate summer dates, regional seasonality and vaccination heterogeneity were absent from text, and 4-region $\times$ week FE covered zero.
- **Action.** `experiment/code/application.py` and Web Table 10 (`tables/table_11_regional.tex`, SI p.36) now report simultaneous exclusion of BW and BY ($N=260, G=14$) at lag 1 ($-0.211$, $[-0.431,0.008]$) and lag 2 ($+0.033$, $[-0.158,0.224]$), 4-region $\times$ week FE at lag 1 ($-0.169$, $[-0.376,0.038]$) and lag 2 ($-0.072$, $[-0.348,0.204]$), state-specific linear time trends ($-0.315$ lag 1; $-0.373$ lag 2), and state-specific 1-harmonic seasonal curves ($-0.224$ lag 1; $-0.318$, $[-0.694,0.057]$ lag 2). Section 4.1 (`main.pdf` p.18) and Section 5 (`main.pdf` p.22) explicitly discuss both late-summer seasonal troughs and southern differences in childhood varicella vaccination uptake (`siedler2010varicella`, `reuss2010varicella`, `suzuki2023varicella`), concluding that rotating holiday variation alone does not isolate an unconfounded effect.
- **Placebos & Randomisation (Web Table 16, SI p.41).** New placebo runs cover Easter, Autumn, Winter, and Pentecost breaks at lags 0--2, calendar shifts ($-4$ to $+4$ weeks), and 200 cross-state schedule permutations. All four non-summer breaks exhibit immediate negative lag-0 dips ($-0.198$ to $-0.497$, all $p<0.05$) followed by positive lag-2 rebounds (Easter $+0.241$ $[0.119,0.363]$; Autumn $+0.150$ $[0.054,0.247]$), directly illustrating reporting compression and post-break catch-up.

### Major Comment 3: Contradiction Between Causal Claims, Rejected $J_{LM,\mathrm{oi}}$, and Stock--Wright Sets
- **Issue.** $J_{LM,\mathrm{oi}}$ rejected ($p_{\mathrm{cal}}=0.004$ under 0-imputed populations; lag-1 growth $6.74$, $p_F=0.020$), `M0-G` sits on the $\phi=0$ boundary under $S_{3,\mathrm{oi}}$, $S_{3,\mathrm{oi}}(\chi^2)$ was empty for all $\bar\phi\in\{0.25,0.5,1.0\}$, and $F$-calibrated $S_3(F)$ spanned `[0.00, 3.00]`.
- **Action.** Removed all unqualified causal claims from the Abstract (`225` words), Introduction, Section 4, Discussion, and Chinese text. Web Table 13 (`tables/table_16_robust_sets.tex`, SI p.38), Section 4.2 (`main.pdf` p.19--20), and Web Appendix F.4 (SI p.38) now disclose verbatim that $S_{3,\mathrm{oi}}$ selects $\phi=0$ ($9.113$ vs $9.120$) so `M0-G` hits the boundary, $S_{3,\mathrm{oi}}(\chi^2)$ is empty across all three $\bar\phi$, and $S_{3,\mathrm{oi}}(F)$, $S_{3,\mathrm{inv}}(F)$, and 97.5\% $S_1(F)$ cover the entire $[0.00,3.00]$ grid.

### Major Comment 4: Selective Summary & All-Holiday Lag-1 Disclosure
- **Issue.** Selective reporting of all-holiday lag 0 while omitting all-holiday lag 1, single-year null/borderline $J_{LM}$ flags, and sensitivity of $\hat\psi^\dagger$ to $\omega_i$.
- **Action.** `main.pdf` Section 4 (`main.pdf` p.18--21), Web Table 9 (`tables/table_7b_application_sensitivity.tex`, SI p.35), and Web Table 15 (SI p.40) now explicitly report all-holiday lag 1 ($-0.081$, $[-0.183,0.021]$) and lag 2 ($+0.060$, $[-0.022,0.142]$), lag-2 summer ($-0.331$, $[-0.714,0.052]$), single-year 2015 ($+0.046$) and 2017 ($J_{LM,\mathrm{inv}}=8.01, p_F=0.050$), and the $>4\times$ sensitivity of counterfactual counts averted $\hat\psi^\dagger$ to $\omega_i$ scaling ($5{,}355$ treated / $15{,}615$ network-wide at $\omega_i\times 0.25$ vs.\ $1{,}585$ / $1{,}633$ at $\omega_i\times 4.0$).

### Major Comment 5: Finite-Cluster Undercoverage, WCR-Score Calibration, and Remark 4 Non-PSD Covariance
- **Issue.** Baseline $t_{15}$ and unrestricted wild score intervals undercovered $\tau+\kappa$, and Remark 4 border-pair adjustment produced an indefinite covariance matrix (`NaN` SEs on $\tau,\kappa$).
- **Action.** Implemented null-imposed restricted score wild cluster bootstrap (`WCR-score`, Webb6 and Rademacher) in both `application.py` ($B=999$) and `calibrated_sim_germany.py` ($R=500 \times 2$, $B=399$). On the calibrated DGP ($k=4.96, \sigma=0.2025$), baseline $t_{15}$ coverage on $\tau+\kappa$ is $0.902$ (`se_ratio` $=0.855$; high-noise $0.904$, `se_ratio` $=0.830$) and unrestricted Mammen/Webb6 cover $0.872/0.894$, whereas **WCR-score restores near-nominal coverage ($0.936$ Webb6, $0.934$ Rademacher; high-noise $0.946$ / $0.946$)**. Observed WCR-score 95\% CIs are $[-0.419,-0.139]$ (Webb6, $p=0.002$) and $[-0.429,-0.119]$ (Rademacher, $p=0.003$). In `remark4_adjacent_dependence`, non-PSD is explicitly flagged (`min_eig = -0.020`, negative diagonal variances on $\tau,\kappa$, formal quadratic SE $0.0482$), and Section 4.2, Remark 4 (SI p.17), Web Appendix E (SI p.32), and Section 5 remove all "slightly conservative" wording.

### Major Comment 6: Data Availability, Runnable Download Script, Checksums, and *Meldewoche* Contract
- **Issue.** Missing runnable download script, Windows CRLF checksum mismatch, median population imputation for 3 counties, and reporting-week vs symptom-onset clarification.
- **Action.** Created runnable `experiment/data/fetch_germany.py` (fetching RKI SurvStat@RKI 2.0 SOAP 1.2, OpenHolidays API, and AGS county GeoJSON; verified `--check` 3/3 OK) and verified via `GetAllDimensions` that SurvStat@RKI 2.0 exposes only `ReportingDate` (*Meldewoche*, no symptom onset date) and that all-case totals (`22,260` in 2017) match unfiltered `[ReferenzDefinition]`. Locked `experiment/data/raw/** binary` in `.gitattributes`, normalized raw JSONs to LF, recorded exact SHA-256 checksums in `experiment/data/checksums.tsv`, and added official neighbouring-year RKI `"population_fallback"` values for the three zero-2017 counties (`09188`, `09763`, `16069`), yielding **`n_imputed_pop_counties = 0`** across all 400 counties (`data_prep.py`, Web Appendix D SI p.26--31).

### Major Comment 7: Page Limit & Transparency
- **Issue.** Biometrics 25-page envelope vs full caveat disclosure.
- **Action.** Kept full caveats, boundaries, Protocol & Departures, and diagnostics un-truncated (`31 / 25` Biometrics-counted pages, `38` total PDF pages in `main.pdf`), explicitly reported in `README.md` and above.

### Major Comment 8: Restoring Protocol vs Post Hoc Disclosure & Boundary Heuristics
- **Issue.** Deletion of "Departures from the protocol", relabelling `M0-G` as pre-specified, and overstating boundary-root rules.
- **Action.** Restored explicit distinction in `sections/methods.tex` (§2.6), `sections/simulation.tex` (`Departures from the protocol`, clarifying that `M0-G` / `M0-GI` were formulated **post hoc** on observing boundary multiplicity in the zero-contagion scenarios, and that all German application coding and sensitivities were formulated post hoc during revision), and `paper_zh/main.tex`, using explicit heuristic/diagnostic wording rather than universal selection guarantees.

---

## Point-by-Point Response to Minor Comments

1. **Minor 1 (`S_1(0)` and single-year $J_{LM,\mathrm{inv}}$):** Corrected to `5.78` (`5.777` unrounded) in text/Table and corrected 2017 `J_LM,inv = 8.01, p_F = 0.050` across Section 4.2 and Web Appendix F.
2. **Minor 2 (Jackson et al. 2014 external benchmark):** Removed the mismatched "aligns with European varicella benchmarks (22%--31%)" claim from Abstract, Introduction, Section 4, Discussion, and Chinese text; cited `adda2016economic`, `isphording2021schools`, and `vonbismarckosten2022schools` accurately on school-closure literature.
3. **Minor 3 (Serial vs generation interval):** Clarified throughout Section 4 and Web Appendix D (`SI p.31`) that Vink et al. (2014) fit a normal serial-interval summary (mean 14.0, SD 2.5 days), which we discretise via a log-normal surrogate as our working generation-interval weights.
4. **Minor 4 ($\omega_i$ & zero-2017 population imputation):** Eliminated median imputation completely (`n_imputed_pop_counties = 0`) using documented official neighbouring-year RKI incidence/count population figures (`09188`: 136,083 [2018]; `09763`: 66,934 [2015]; `16069`: 63,694 [2018]).
5. **Minor 5 (Remark 4 `NaN`):** Synchronised Remark 4 in `sections/appendix_moved.tex` (`SI p.17`), `sections/appendix.tex` (Web Appendix E, `SI p.33`), Section 4.2 (`main.pdf` p.19), and `paper_zh/main.tex` to flag `non_psd = True` (`min eig = -0.020`) and negative diagonal variances on $\tau,\kappa$.
6. **Minor 6 ($\ge 4$ days threshold justification & $\ge 3/\ge 5$ results):** Clarified the majority-of-week ($\ge 4$ of 7 days) baseline rule in Section 4 (`main.pdf` p.17) and reported $\ge 3$-day ($-0.223$, $[-0.332,-0.113]$) and $\ge 5$-day ($-0.190$, $[-0.282,-0.098]$) results in Section 4 and Web Table 9 (`SI p.35`).
7. **Minor 7 (Author/funding placeholders):** Retained `Apevon Science` and `[Funding statement to be completed before submission.]` as requested by the author for pre-submission completion, noted in `README.md`.
8. **Minor 8 (Figure 3(b) caption):** Rewrote Figure 3 (`main.pdf` p.35) to plot and label the total-policy scale $\hat\tau_e + c_e\hat\kappa$ with $\hat\tau_0 + c_0\hat\kappa = -0.208$ annotated cleanly without overlap.
9. **Minor 9 (Outdated citations & literature):** Updated `goldstein2024wastewater` to *Biometrics* 80(3):ujae074 (2024) and added `adda2016economic`, `vonbismarckosten2022schools`, `isphording2021schools`, `siedler2010varicella`, `reuss2010varicella`, and `suzuki2023varicella`.
10. **Minor 10 (Developer comment in `application.py`):** Removed `# Wait: check exact convention used in prior application` from `experiment/code/application.py`.
11. **Minor 11 (Rounding of $w$):** Synchronised discrete weights as $w = (0.141, 0.720, 0.138)$ (`0.1414, 0.7205, 0.1381` in Web Appendix D) across all files.

---

## Replication & Verification Summary (Verified After Final Edits)

- `python experiment/data/fetch_germany.py --check`: **3/3 OK**
  - `rki_varicella_2014_2019.json`: 2,475,214 bytes, SHA-256 `0735f0211c4d51389a20cb93dbcd4d5b48f8f7790d1c68176c8841c12d0b2489`
  - `germany_school_holidays_2014_2019.json`: 261,944 bytes, SHA-256 `4f51fd93e04f0f7ae37a5d78aea0966d12f74acf15b7293a795f155826fb695e`
  - `counties_ags.geojson`: 1,289,615 bytes, SHA-256 `8a0d0bf29ff31b6e536622613ce6e2f6198a6999bd725a631b4b1707780a6abd`
- `python experiment/code/data_prep.py`: **passed** (`n_imputed_pop_counties = 0`, $N=400$, $T=261$, 115,129 cases).
- `python -u experiment/code/application.py`: **passed** in 627.67s (`experiment/analysis/germany_application_results.json`).
- `python -u experiment/code/calibrated_sim_germany.py`: **passed** in 140.95s (`R=500 \times 2`, 1,000/1,000 converged, including `WCR-score` Webb6 and Rademacher coverage).
- `python experiment/code/make_figure_germany.py` & `make_tables_germany.py`: **passed**, regenerated `fig_03_application.png/.pdf` and all 9 LaTeX tables (`tables/table_7_application.tex`, `table_7b_application_sensitivity.tex`, `table_11_regional.tex`--`table_18_placebos.tex`).
- Build & QA after final `.tex` edits:
  - `python scripts/scan_control_chars.py`: `total 0`
  - `python scripts/qa_checks.py`: running head `'DiD on the transmission scale'` (29 chars $\le 50$), abstract `225` words ($\le 225$), `main.pdf` 38 pages (`[]` unresolved `??`), `supporting_information.pdf` 45 pages (`[]` unresolved `??`), `paper_zh/main.pdf` 26 pages (`[]` unresolved `??`).
  - `python scripts/section_pages.py`: `Biometrics`-counted pages `31 / 25`, zero undefined references or citations across all `.log` files.
