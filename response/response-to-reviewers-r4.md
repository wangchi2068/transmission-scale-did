# Response to the Fourth-Round Referee

**Manuscript:** "Difference-in-differences on the transmission scale: sequential moment estimation of policy effects on the reproduction number under feedback and spatial contagion"

**Journal:** Biometrics

**Authors:** [author names to be supplied before submission]

We thank the referee for the fourth-round report. All five major points (§4.1), the actionable minor items (§4.2), and the documentation/errata items (§4.3) have been addressed. Rather than re-compressing prose to force an artificial 25-page squeeze (which caused the regressions noted in §4.1 items 4–5), we restored full qualifiers and the standalone `Data Availability Statement` heading; the resulting manuscript is 27 Biometrics-counted pages (PDF pp. 2–28, printed pp. 1–27; abstract 222 words ≤ 225; zero overfull boxes and zero unresolved cross-references). In addition, a pre-submission self-audit identified and resolved several technical, numerical, and bibliographic issues (listed in Section 4 below).

**Page conventions.** "PDF p." refers to the page of `main.pdf` (equals printed page + 1); "Web p." refers to the page of `supporting_information.pdf`.

---

## 1. Major comments (§4.1)

### 1.1 Scoping the 0.783 same-week coverage to §3.2.1 (Discussion, PDF p. 21 / printed p. 20)
**Agreed and rewritten.** The Discussion now explicitly scopes `0.783` to the §3.2.1 design and contrasts it with the application's precision:
> "That calibration does not include same-week transmission: in the same-week design of Section 3.2.1, where the bias (−0.053) exceeds the standard error (0.046), the interval covers 0.783 and J has no power; at the application's precision (SE 0.31), a comparable bias would barely move coverage, so the concern is the target—a pseudo-true value—rather than interval width."

Synced to `paper_zh/main.tex`. Erratum for `response-to-reviewers-r3.md` (lines 21, 88): `0.783` comes from Section 3.2.1 (`verification/same_week/results.json`), not Web Table 7 (which reports the root-selection validation).

### 1.2 Rewriting the Introduction same-week sentence and unifying 0.30–0.39 (§1, PDF pp. 3–4; §5, PDF p. 21)
**Agreed and rewritten.** The Introduction sentence now reads:
> "The generation interval places 0.30–0.39 of transmission in the same week, a regime in which the weekly model is biased and J has no power (Section 3.2.1); the rejection by J and Q in Section 4 is therefore a separate finding."

The Discussion (§5, PDF p. 21) now also uses `0.30–0.39` instead of `about a third`, and both Introduction and Discussion use `division × week effects` to match the Abstract (§4.2-6).

### 1.3 Qualifying the empty robust set in the Abstract (Abstract, PDF p. 1)
**Agreed and updated.** The Abstract now reads: `"and the robust set for the coupling parameter is empty with over-identifying moments."` With a small wording trim earlier in the Abstract (`"Excluding same-week cross-cluster transmission"`), the abstract remains 222 words (≤ 225).

### 1.4 Restoring the Proposition 3 qualifiers and separating Remark 2 (§2.5, PDF p. 11)
**Agreed and restored.** Section 2.5 now restores the first-round `m5` qualifiers and separates Remark 2 from the non-negative weights:
> "a first-order expansion with (κ, φ) held at their true values giving non-negative weights under time-invariant unit heterogeneity (Proposition 3; Goodman-Bacon, 2021; de Chaisemartin and D'Haultfœuille, 2020), whereas under dynamic effects pooled weights can be negative, motivating the event-time specification (Remark 2)"

It also restores the distinction between the Jensen gap and Nickell-type feedback term and the classical fixed-T inconsistency of Poisson fixed effects.

### 1.5 Restoring standalone `Data Availability Statement` and correcting deletion statements (Back Matter, PDF p. 21; README)
**Agreed and restored.** `sections/backmatter.tex` now restores standalone unnumbered `\section*{Data Availability Statement}` (pointing to Web Appendix D and `experiment/data`), `\section*{Supporting Information}`, `\section*{Conflict of Interest}`, and `\section*{Funding}`. `README.md` has been corrected to drop the inaccurate claim that nothing was removed in the prior round and to document the actual edits.

---

## 2. Minor comments (§4.2)

1. **§2.6 step (4) and Web Appendix F.4 alignment (PDF p. 12; Web p. 36):** §2.6 now recommends `"report the profiled S set for φ and, if reported, the projection for τ+κ with its strong-identification condition"`. Web Appendix F.4 and Web Table 13's footnote now state the strong-identification condition explicitly—and also correct an error in the previous draft (see §4.1 below: the *profiled* $S$ test itself relies on strong identification of $(\delta,\tau,\kappa)$ at the tested $\phi$).
2. **`negative thereafter` vs. $e_5 = +0.009$ (§4, PDF p. 18):** Changed to `"mostly negative thereafter (+0.009 at e5, down to −0.471 for e6+)"`, and aligned the linearisation-check contrast with its matching evaluation point (see §4.2 below).
3. **`0.074` vs. `0.075` (§4, PDF p. 18):** Replaced the grid-point value `0.074` with the exact M0 root `0.075` (`0.07470`), matching Web Table 9 (`tab:app-all`), and cited Web Table 9 explicitly (erratum for `response-to-reviewers-r3.md` M5, which cited Web Table 13).
4. **Unequal cluster size (§4, PDF p. 16):** Replaced `"so the largest state carries much of the weight in every moment"` with `"from 1 to 254 counties (median 64)"` (verified directly on `panel.npz`: Texas is 8.2% of counties and 10.1% of cases, close to California 9.8% and Florida 8.1%).
5. **Reporting sentence wording (§4, PDF p. 19):** Restored `"a share 1 − ρ of what they would report untreated"` and added `"respectively"`.
6. **`regional` vs. `division × week` (§1, PDF p. 3; §5, PDF p. 21):** Unified to `division × week effects` across Abstract, Introduction, and Discussion.
7. **Chen & Roth (2024):** Remains cited in Web Appendix A.9.
8. **Author and funding placeholders:** Retained as placeholders to be completed by the authors directly before ScholarOne upload.

---

## 3. Response letter and README errata (§4.3)

- **README.md:** Fixed printed page breakdown (main text + Back Matter printed pp. 1–20, References printed pp. 21–27, endfloat tables/figures printed pp. 28–33; counted 27 pages), unified jackknife excess to `14%–30%` (`13.6%–30.1%`), fixed Web Table numbering (`Web Tables 9–15` in F, `Web Table 9` for `split_app_table.py`), removed the erroneous `6/10` round label, and updated version notes.
- **Auxiliary docs:** Updated `experiment/data/data_contract.md` and `experiment/code/README.md` to remove stale references to `src/data_prep.py`, `fetch_data.sh`, and `tests/test_application_data.py`, and clarified that `run_plan.json` is written by `benchmark.py`.
- **R3 response letter errata:** Noted above (`0.783` is from §3.2.1, not Web Table 7; `0.075` is in Web Table 9, not Web Table 13; `response-to-reviewers-r3.md` addressed the third-round re-review of `b57ce4e`, whose cluster-size edit appeared on printed p. 15 / PDF p. 16).

---

## 4. Additional pre-submission fixes from self-audit

1. **Profiled Stock–Wright $S$ set condition (§2.6, PDF p. 12; §4, PDF p. 18; Web Appendix F.4, Web p. 36; Web Table 13):** Removed the incorrect statement that the profiled $S$ set for $\phi$ does not need strong identification of $(\tau,\kappa)$. Because `robust.py` profiles out $(\delta,\tau,\kappa)$ via $\psi_g = g_{t,g} - \hat\Gamma_{t,n}\hat\Gamma_{n,n}^{-1} g_{n,g}$, both the profiled $S$ test (Stock & Wright, 2000, Theorem 3) and the Bonferroni projection rely on strong identification of $(\delta,\tau,\kappa)$ at the tested $\phi$, which fails near $\phi=0$. Clarified the degrees of freedom in Web Table 13's footnote ($d=1$ for $S_1$, $d=3$ for $S_3$).
2. **Event-time contrast evaluation-point alignment (§4, PDF p. 18):** The previous draft quoted M0-evaluated weights (`2.109` to `−0.471`) followed by the M0-ET2-evaluated contrast (`0.079`, SE `0.240`, `p = 0.742`). The text now reports the matching M0-evaluated contrast (`−0.066`, SE `0.268`, `p = 0.806`) alongside the M0-ET2-evaluated value (`0.079`, SE `0.240`, `p = 0.742`).
3. **Calibrated simulation DGP disclosure (§4, PDF p. 16; §5, PDF p. 20; Web Appendix E & F.2):** Clarified that the calibrated simulation satisfies the weekly conditional mean A1 with fixed (exogenous) adoption dates (no feedback) and spatially correlated state × week shocks (mean adjacent correlation 0.21, a mild departure from A4).
4. **Methodological and proof polish:** Added the staggered-adoption requirement of the default instrument $s_{it}D_{it}$ in §2.4 (PDF p. 9), noted the 1–2 week case-reporting delay in §4 (PDF p. 19), clarified Figure 3(b)'s caption that M0-G compares $Q$ only over root candidates and $\phi=0$, defined $M_t$ explicitly in Theorem 1(b)'s proof and tightened the $\eta$-derivative bounds in Lemma 5 (Web Appendix A), moved Web Table 9 (`tab:app-all`) out of simulation appendix C.2 into the start of Web Appendix F (Web p. 33, preserving Web Table numbers 1–15), corrected `"five runs"` to `"four sensitivity runs (five including the main fit)"` in §4, `"within ±0.006"` to specify treated units (`±0.007` network-wide) in §3.2, and `"53.6%"` to `"53.5%"` in Web Appendix C (matching Web Table 6).
5. **Bibliography (`bibliography.bib` and `paper_zh/bibliography.bib`):** Added missing coauthor Peter Teunis to `wallinga2004different` (now correctly rendered as `Wallinga and Teunis, 2004`), fixed surname parsing for `Cameron, A. Colin` and `Santos Silva, J. M. C.`, added middle initial `C. C.` to `Heyde`, and added missing end page numbers to six entries (`arellano1991tests`, `hansen1982gmm`, `hausman1984count`, `newey1991uniform`, `nickell1981biases`, `rothenberg1971identification`).
