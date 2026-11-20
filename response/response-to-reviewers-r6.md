# Response to Sixth-Round Review (`re-review` on `1dfc773`)

We thank the referee for the extraordinarily thorough verification across `main.pdf`, `supporting_information.pdf`, `paper_zh/main.pdf`, `bibliography.bib`, `experiment/data/fetch_germany.py`, and the replication suite. Every technical bug, reproducibility defect, main-text/appendix discrepancy, citation omission, and numerical typo raised in Round 6 has been fixed below and verified directly on the freshly compiled outputs (`main.pdf` 36 pages, `supporting_information.pdf` 49 pages, `paper_zh/main.pdf` 25 pages).

---

## Response to Major Comments

### Major Comment 1: Reconciling Main-Text $\phi$ Identification & Diagnostic Findings
- **Fixes made (`main.pdf` p.18 & p.20, §4.1 & §4.4; Table 4 p.36; Web Appendix F.4--F.7 Web p.37--39 & p.41--44):**
  - **§4.1 (`main.pdf` p.18):** Replaced unqualified `"excluding zero contagion coupling"` with an explicit qualification noting that while local M0 converges to $\hat\phi=0.326$ where exact-moment $S_1(F)=[0.21,2.40]$ excludes 0, **the pre-specified global rule `M0-G` selects the $\phi=0$ boundary ($\hat\tau+\hat\kappa=-0.215$, Table 4) because $S_{3,\mathrm{oi}}$ marginally favours $\phi=0$ ($9.113$ vs.\ $9.120$)**, and `M0-GI` retains the interior root only under $S_{3,\mathrm{inv}}$.
  - **§4.4 (`main.pdf` p.20):** Added a direct diagnostic summary citing **Web Appendices E--F, Web Tables 13 & 16 (Web p.41 & p.44), and Web Figure 5 (Web p.38)**: explicitly stating that (i) lagged working moments are rejected by $J_{LM,\mathrm{oi}}=9.12$ ($p_{\mathrm{cal}}=0.004$), $S_{3,\mathrm{oi}}(\chi^2)$ is empty across all $\bar\phi$, and $F$-calibrated joint/invariance sets span $[0,3]$; and (ii) non-summer holiday placebos show lag-0 drops and lag-2 rebounds alongside summer onset $-0.208$ ($p=0.049$), placing the panel in Table 3's **Assumption A7-violating regime** (where simulation X2 also exhibited $+0.424$ upward bias on $\hat\phi$, meaning $\hat\phi=0.326$ may itself absorb reporting variation).

### Major Comment 2: Known-Truth Comparator Benchmark in Design-Calibrated Simulations
- **Fixes made (`experiment/code/calibrated_sim_germany.py`, Web Table 11(c) Web p.40, §4.2 `main.pdf` p.18, §5 `main.pdf` p.21):**
  - Expanded `calibrated_sim_germany.py` ($R=500 \times 2$ arms, 1,000/1,000 converged) to evaluate **county-level `M2` (EpiEstim log-$R$), county-level `M3` ($\log(Y+1)$ TWFE), and aggregated 38-region `NUTS-2` (`M0` and `M2`) directly against the known generating truth ($\tau+\kappa=-0.239$, $\phi=0.326$)**.
  - **Results (`Web Table 11(c)`, Web p.40):** Under simulated county zero shares (`54.5%` calibrated, `57.8%` high-noise), county `M0` has bias `-0.003` (`-0.016` high-noise), whereas county `M2` exhibits severe positive bias **`+0.224`** (`+0.227` high-noise) and flips strictly positive in **`22.8%`** (`32.6%` high-noise) of draws. Aggregating identical simulated data to NUTS-2 (`2.3%--3.6%` zeros) reduces `M2` bias by two-thirds to **`+0.076`** (`+0.083` high-noise), while `NUTS-2 M0` bias remains small (`+0.013` calibrated, `-0.008` high-noise).
  - Updated §4.2 (`main.pdf` p.18) and §5 (`main.pdf` p.21) to report the observed panel as *"consistent with Theorem 3(a)'s small-count distortion $B_J$, which is directly confirmed against known truth in Web Table 11(c)"*, replacing over-strong claims (`isolating`, `clear comparative advantages`).

### Major Comment 3: Operating Envelope (Table 3, `main.pdf` p.35) & §4.3 Few-Cluster Claims (`main.pdf` p.19)
- **Fixes made (`tables/table_envelope.tex`, `calibrated_sim_germany.py`, `make_tables_germany.py`, §4.3 `main.pdf` p.19, Web Appendix E Web p.33):**
  - **Boundary row:** Labelled `M0-GI (post hoc)`, removed `+ S1(phi)` from the coverage claim, and explicitly reported the weak-contagion coverage cost (`0.88--0.95` at $\phi_0=0$; `0.90--0.93` at $\phi_0\in[0.05,0.10]$).
  - **Finite-cluster row & delete-one-state jackknife:** Added full **16-fold delete-one-state jackknife ($t_{15}$)** inside `calibrated_sim_germany.py` across all $R=500 \times 2$ panels (**`500/500` converged in both arms**): calibrated jackknife $t_{15}$ coverage is **`0.966`** (`high-noise`: **`0.964`**), directly backing the Table 3 recommendation alongside `WCR-score` (`0.936` Webb6 / `0.934` Rademacher; `high-noise` `0.946` / `0.946`). Table 3 footnote and §4.3 now explicitly report `MCSE ~ 0.010` at $R=500$, plus the $R=2000$ validation designs (`Web Appendix C`). In §4.3 (`main.pdf` p.19), replaced `"agreeing closely"` with exact interval widths (`"similar in width to ...; coverage of pairs bootstrap not simulated"`).
  - **Assumption A7 row:** Added **`phi (+0.42)` bias in X2** explicitly (`biases \tau+\kappa (\approx -0.35) & \phi (+0.42) in X2`).

### Major Comment 4: Holiday Placebos, Calendar-Shift Lead Alignment, and Permutation Caveats
- **Fixes made (`make_tables_germany.py`, Web Table 16 Web p.44, Web Appendix F.5 & F.7 Web p.37--39, Web Figure 5 Web p.38, §4.4 `main.pdf` p.20):**
  - **Calendar-shift lead labels (`make_tables_germany.py` line 592; Web Table 16(b) Web p.44):** Corrected `1 + sh >= 0 -> lag (1+sh)` and `1 + sh < 0 -> lead -(1+sh)`, so shifts `-4, -3, -2, -1` now read `-4 wk (lead 3)`, `-3 wk (lead 2)`, `-2 wk (lead 1)`, and `-1 wk (lag 0)` cleanly.
  - **Pre-trend vs. calendar-shift reconciliation (`Web Appendix F.5 & F.7`, Web Figure 5 caption Web p.38):** Explained explicitly why overlapping shifted schedule indicators at leads 3--1 produce positive pre-window slopes (`+0.296, +0.248, +0.328`), whereas event-time bins $e_{-2}, e_{-1}$ (`p_F = 0.765`) condition directly on relative holiday onset; both are treated jointly as flagging calendar/reporting confounding rather than clearing causal validity.
  - **Permutation caveat & baseline lag-1 rationale:** Clarified in Web Appendix F.7 and Web Table 16 footnote that state permutation ($p=0.010$ at $B=200$, minimum attainable two-sided $p$-value) lacks strict exchangeability because BW and BY never rotate and cannot separate transmission from calendar-aligned reporting. Explained in §4.1 (`main.pdf` p.17) that lag-1 baseline coding tracks the ~14-day serial interval in weekly *Meldewoche* counts while noting lag-2 sensitivity (where BW/BY exclusion gives `+0.033` and all-holidays gives `+0.060`).

### Major Comment 5: `fetch_germany.py` Live Parsing Bug, Extraction Dates, and Data Availability
- **Fixes made (`experiment/data/fetch_germany.py`, `data_contract.md`, `sections/backmatter.tex` p.22, Web Appendix D Web p.26--32, `paper_zh/main.tex`, `README.md`):**
  - Fixed `_parse_olap` in `fetch_germany.py` to read `b:ColumnName` and `b:RowName` (instead of `b:Id`), added strict empty-parse guards (`RuntimeError` if any year returns `<400` rows or `0` cases), corrected UTF-8 state names (`Baden-Württemberg`, `Thüringen`) and 1-space holiday indentation, and added post-download SHA-256 verification.
  - **Live test verified:** `_parse_olap(_soap_olap(2017, ...))` fetched **412 rows $\times$ 52 weeks, total = 22,260 cases**, matching `rki_varicella_2014_2019.json` to the single case.
  - Recorded exact archive extraction date **`2026-11-19`**, GitHub persistence URL (`https://github.com/wangchi2068/transmission-scale-did`), and open data licences (RKI SurvStat@RKI 2.0 public data, OpenHolidays MIT/open schedules, `geofeatures-ags-germany` BKG/GeoBasis-DE Data licence Germany) in `data_contract.md`, `backmatter.tex`, Web Appendix D, and `paper_zh/main.tex`.

### Major Comment 6: 25-Page Biometrics Limit
- **Author decision retained:** Per explicit author instruction (*"质量为主，不强制压缩 25 页"*), the repository holds the un-compressed working manuscript (`36` pages total, `28 / 25` counted pages, `4` Tables + `3` Figures) so full diagnostics, Operating Envelope guidance, and protocol disclosures remain intact. `README.md` explicitly states that this is the un-compressed full working manuscript and a page-limited journal submission version will be formatted separately before ScholarOne submission.

### Major Comment 7: Literature Citations & Unused `.bib` Entries
- **Fixes made (`sections/application.tex` p.17 & p.20, `sections/discussion.tex` p.21, `sections/appendix.tex` Web p.38--39, `bibliography.bib`, `paper_zh/bibliography.bib`):**
  - **Zero uncited `.bib` entries (`UNCITED: []`, `MISSING: []`):** Now cited `adda2016economic`, `isphording2021schools`, and `vonbismarckosten2022schools` in §4.1 and Web Appendix F.7; cited `siedler2010varicella` and `reuss2010varicella` in §4.4 and Web Appendix F.7 (highlighting German routine vaccination since 2004 and North/South state uptake heterogeneity); cited `suzuki2023varicella` in §4.4 and §5; cited `roth2023parallel` and `elenev2025staggered` in §5; and pruned the 10 obsolete US COVID-19 mask entries (`cdc2021maskorders`, `dong2020dashboard`, `ferretti2020quantifying`, `ginther2021kansas`, `guy2021masks`, `lyu2020masks`, `nishiura2020serial`, `vandyke2020kansas`, `conley2011inference`, `mackinnon2018wild`).

---

## Response to Minor Comments

1. **Five Web Appendix numbers + Remark 4 `16%`:** All fixed to match `germany_application_results.json` and `germany_calibrated_sim.json` (`0.431` -> `0.451` Web p.34; state trends `[-0.471,-0.160]` -> `[-0.550,-0.081]` Web p.34; 1-harmonic `[-0.371,-0.077]` -> `[-0.352,-0.096]` Web p.34; `make_tables_germany.py` NUTS-2 footnote now reads JSON dynamically: `7.21 (0.064)` and `4.99 (0.041)` Web p.35; F.5 `e_1--e_5` `-0.374 to -0.298` Web p.38; Remark 4 `16%` -> `14.5%` / SE ratio `0.855` Web p.17 and `paper_zh/main.tex`).
2. **Serial interval vs. generation interval:** Distinguished throughout §4, §3.2.1 (`serial or generation intervals`), Web Appendix D, Web Table 12 (`table_13_gi_timing.tex`), Table 4 footnote, and Web Appendix F.3 (`log-normal working generation-interval surrogate fitted to Vink et al. 2014 serial-interval moments`).
3. **Lag-1 rationale:** Documented in §4.1 (`main.pdf` p.17).
4. **Author & funding placeholders:** Retained intentionally per author instruction until formal submission.
5. **`suzuki2023varicella` & unused `.bib`:** Fixed `Ayaka` -> `Ayako` in both `.bib` files; pruned unused legacy entries (`UNCITED: []`).
6. **Vaccination uptake heterogeneity:** Integrated with `siedler2010varicella` and `reuss2010varicella` in §4.4 (`main.pdf` p.20) and Web Appendix F.7 (`Web p.39`).
7. **Chinese draft (`paper_zh/main.pdf`, 25 pages):** Fixed `0.431` -> `0.451` and `约 16%` -> `14.5%`; added `M0-G` boundary root (`\hat\tau+\hat\kappa=-0.215, \phi=0`), calendar-shift positive leads, `jackknife 0.966/0.964`, and known-truth `M2` bias (`+0.224 -> +0.076`); toned down all `"显著"` / `"稳健"` wording.
8. **Web Figure 5 reference:** Cited explicitly in §4.4 (`main.pdf` p.20) and Web Appendix F.5 & F.7 (`Web p.37--39`); Web Figure 5 (`Web p.38`) caption updated to reconcile `p_F=0.765` bins vs. positive overlapping pre-window calendar shifts.
9. **README font & runtime notes:** Added Times New Roman `style_kit.py` requirement and 4-worker vs. single-thread runtime note in `README.md`.
10. **Table 3 `M0-GI` label:** Marked `M0-GI (post hoc)` explicitly in `tables/table_envelope.tex` (`main.pdf` p.35).

---

## Response to Section 10 Suggestions (Scope & Feasibility for Beyond 7/10)

We appreciate the referee's constructive suggestions on how the German empirical implementation could be extended further toward a structural causal application and an exhaustive finite-cluster benchmark. After careful consideration, we retain our current scope for the following reasons:

1. **Symptom-onset data and joint reporting--transmission cascade modelling:** Publicly available `SurvStat@RKI 2.0` exports are aggregated strictly by notification week (*Meldewoche*) rather than symptom-onset date, and routine county-level varicella surveillance lacks an auxiliary ascertainment anchor (such as sentinel testing or serological panels) needed to separately identify time-varying reporting delays from biological transmission without strong parametric delay assumptions. As the referee notes in Section 1, reframing §4 as an *illustrative sparse-count implementation* rather than a causal holiday-effect study resolves the prior causal overclaim cleanly: the application now fulfills its methodological purpose by validating Theorem 3(a)'s small-count `M2` distortion against known truth in design-calibrated simulations (`Web Table 11(c)`, Web p.40) and demonstrating how onset-week and non-summer placebo diagnostics expose Assumption A7 violations in practice (§4.4 `main.pdf` p.20, Web Appendix F.7 Web p.38--39). Joint nowcasting and ascertainment-delay cascades remain a natural standalone extension already highlighted in §5 (`main.pdf` p.21).
2. **Nested Monte Carlo simulation of nonlinear pairs cluster bootstrap-$t$:** In commit `3252fb6`, we implemented full 16-fold delete-one-state jackknife ($t_{15}$) inference inside `experiment/code/calibrated_sim_germany.py` across all $R=500 \times 2$ design-calibrated panels (`500/500` converged in both arms), establishing empirical coverage of `0.966` (`0.964` high-noise). Together with `WCR-score` (`0.934`--`0.946`), **every finite-cluster inference method recommended in Table 3 (`main.pdf` p.35) now has direct Monte Carlo coverage backing**. By contrast, `pairs cluster bootstrap-t` is *not* recommended in Table 3 and appears in §4.3 (`main.pdf` p.19) solely as an observed-sample sensitivity check explicitly flagged as `"coverage not simulated"`. Nesting $B$ nonlinear GMM profile optimizations inside 1,000 outer Monte Carlo panels (where the current simulation including 16-fold jackknife already takes `488.6s`) would multiply computation by orders of magnitude while offering limited methodological return, especially because `pairs` resampling across $G=16$ strongly unbalanced states can yield degenerate cluster weights.

---

## Verified Final Repository & Build Checks

- `python experiment/data/fetch_germany.py --check`: **3/3 OK** (`2,475,214` / `261,944` / `1,289,615` bytes).
- Live RKI SOAP test (`2017`): `412 rows, 52 cols, total = 22,260` (exact match).
- `python -u experiment/code/calibrated_sim_germany.py`: **500/500** converged in both arms (`488.6s`), adding simulated delete-one-state jackknife (`0.966` / `0.964`) and known-truth M2/M3/NUTS-2 (`M2` county bias `+0.224`, `M2` NUTS-2 `+0.076`).
- `python experiment/code/make_tables_germany.py`: regenerated all 9 tables dynamically from JSON.
- `python scripts/scan_control_chars.py`: `total 0`.
- `python scripts/qa_checks.py`: running head `29 <= 50`, abstract `202 <= 225`, `main.pdf` 36 pages (`[]`), `supporting_information.pdf` 49 pages (`[]`), `paper_zh/main.pdf` 25 pages (`[]`), `UNCITED: []`, `MISSING: []`, zero undefined references or citations across all logs.
