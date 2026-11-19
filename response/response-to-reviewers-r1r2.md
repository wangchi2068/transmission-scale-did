# Response to the Referee (first and second rounds)

**Manuscript:** "Difference-in-differences on the transmission scale: sequential moment estimation of policy effects on the reproduction number under feedback and spatial contagion"

**Journal:** Biometrics

**Authors:** [author names to be supplied by the authors]

This file answers the first-round report (major revision, 5/10). Its page references have been re-checked against the current build, commit `2244b15`, in which the main text occupies printed pages 1–19 and the Supporting Information 43 pages. The second-round report (minor revision, 7/10) is answered in `response-to-reviewers-r3.md`.

We thank the referee for a close and constructive report. The report reran parts of our archive, and several of its own computations (the Q reference values, the adjacent-cluster score check and the graph diameter) became analyses in the revision. Below we answer every point in the order of the report: weaknesses W1–W7, the numbered requests under Comments 1–7, minor comments m1–m15, presentation items t1–t6 and the questions of Section 4.4.

**Conventions.**
- "p." is the printed page of the main manuscript (PDF page minus one; figures and tables follow the references, Biometrics style).
- "SI p." is the printed page of the Supporting Information (43 pages).
- Section, table and figure numbers refer to the revised files.
- Every number quoted below is taken from the result files of the revision archive. Each new Web Table names its result file in its note.
- "Post hoc" marks analyses that were added in revision after seeing the first-round results. All revision analyses are post hoc, and the manuscript labels them so (§3.1, p. 11; §4, p. 14; Web Appendix F).

## Summary of the main changes

1. **The application is now presented as a diagnosis of the weekly working model.** The model is still rejected after adding region × week or division × week effects and region-specific coupling: J_LM is 16.9–17.8, against calibrated 95th percentiles of 6.38–6.66 (Web Table 10). The Abstract, §1, §4 and §5 therefore report τ+κ = −0.466 (95% CI −1.091 to 0.159) as the pseudo-true value of the rejected weekly model. It is not presented as a causal effect of mask orders.
2. **New application diagnostics (Web Appendix F, Web Tables 10–15).** These are:
   - regional week effects;
   - an application-calibrated simulation, which gives reference distributions for J, Q_G and Q_GI and evaluates the three interval methods;
   - generation-interval reassignment of same-week mass;
   - Stock–Wright identification-robust sets;
   - the implicit event-time weights, with a joint-sandwich contrast;
   - an ascertainment tipping point.
3. **Theory scope is stated where it applies.**
   - Theorem 2 covers interior, globally identified cases. φ0 = 0 violates A6(i), and post-selection inference is outside the theorem (§2.5, p. 8).
   - Lemma 3 is a large-G device that does not describe the 49-state graph (§2.5, p. 8).
   - M0-G and M0-GI are described as selection heuristics without a stated guarantee (§2.6, p. 10; §3.2, p. 12).
4. **Positioning.** The contribution is stated as the integration of log-R parallel trends, the renewal mean with imported pressure and quasi-differenced moments, together with its diagnostics (§1, pp. 1–2). Kosfeld et al., Sun and Xiao, and Zhang and Crawford are cited and contrasted. Theorem 3, Proposition 3 and Theorem 4 are labelled as adaptations of known results.
5. **Space.** Comparator and misspecification detail moved to Web Appendix C. The freed space in §4 carries the diagnostics. The main text stays within the 25-page Biometrics count.

---

## Major weaknesses

### W1. The application rejects the working model

We agree, and the application has been reframed.

- **Rejection persists.** Rejection survives every regional specification we fitted (Comment 1.1).
- **Q is outside its calibrated range.**
  - In an application-calibrated simulation (R = 600), the application's Q_G = 40.64 exceeds all draws (maximum 35.13).
  - J = 17.87 also exceeds all draws (95th percentile 6.66, maximum 12.45).
  - Q_GI = 27.62 is exceeded in 8 of 600 draws.
  - §4 states this directly (p. 16; Web Table 11, SI p. 35).
- **Dependence on φ̄ is reported.**
  - The interior component of the S set moves with φ̄: [0.53, 0.85] at φ̄ = 0.25, [0.58, 0.97] at φ̄ = 0.5 and [0.64, 1.15] at φ̄ = 1.
  - The text states that the interior root and the φ estimate depend on φ̄ (§4, p. 16; Web Table 13).
- **Reframed reading.**
  - The Abstract now says the diagnostics reject the weekly model.
  - The total τ+κ = −0.466 is described as a summary of the weekly fit, not a causal effect (unnumbered title page; §1, p. 2; §4, p. 15; §5, p. 18).

### W2. Same-week stress regime; the cited sub-weekly extension does not exist

Partly addressed.

- **Wording corrected.** We removed the reference to a sub-weekly extension.
- **Disclosure in the text.** The text now states, in §3.2.1 (p. 13) and Web Appendix A.12 (SI p. 18):
  - the application lies in the w0 = 0.3 stress regime (same-week mass about 0.30 under Ferretti, 0.39 under Nishiura);
  - in this regime τ+κ under-covers and J has no power against the misspecification;
  - a sub-weekly mean is needed, and this paper does not provide one.
- **New analysis (Comment 2.2).** Reassigning same-week mass moves τ+κ only between −0.476 and −0.444 (Web Table 12).
- **Not done: the sub-weekly estimator.** We did not implement a sub-weekly estimator. Comment 2.1 gives the reason.

### W3. Boundary and root-selection inference lie outside Theorem 2

Addressed.

- **Scope statement.** §2.5 (p. 9) states:
  - Theorem 2 covers interior, globally identified cases;
  - φ0 = 0 is a boundary that violates A6(i);
  - in the zero-contagion designs the population moments have a second root that moves with φ̄;
  - post-selection inference is not covered.
- **Population second root.** It appears in all 40 large-N draws of the four zero-contagion designs (Web Table 8, SI p. 31).
- **Discussion sentence replaced** (§5, pp. 18-19). See Comment 3.2.
- **Identification-robust S sets added** (§2.6, p. 10; §4, p. 16; Web Table 13).

### W4. The dynamic-effects narrative is stronger than the evidence

Addressed.

- **Abstract.** It no longer makes a dynamic-effects claim.
- **Derivation.** Web Appendix A.10 (SI p. 15) now derives the implicit weights of the pooled estimand.
- **§4 (p. 16).** It now reports:
  - the joint F test (F = 2.089 on 7 and 30 df, p = 0.076);
  - the pooled–dynamic contrasts with joint-sandwich standard errors;
  - J on the event-time sample.
- **Wording.** The divergence is described as consistent with dynamic weighting but not distinguishable from sampling error.

### W5. Finite-cluster inference is not calibrated

Addressed with an application-calibrated simulation (Comment 5.1).

- **Coverage.** For τ+κ, t48 sandwich coverage is 0.935, wild 0.925 and pairs 0.933 (R = 600). In a stress design matched to the application's standard error, the figures are 0.950, 0.930 and 0.923 (R = 300).
- **Primary interval.** The t48 sandwich interval is declared primary. The text says this designation was made in revision and is post hoc (§4, p. 14).
- **Few-cluster literature.** Ferman and Pinto (2019) and Canay, Santos and Shaikh (2021) are discussed (§5, p. 18; Web Appendix E, SI p. 33).
- **Lemma 3.** Its large-G scope is stated (§2.5, p. 8).

### W6. Novelty and positioning

Addressed.

- **Contribution.** The Introduction states the contribution as an integration of known components together with its diagnostics (§1, pp. 1–2).
- **Related work.** It cites and contrasts:
  - Kosfeld et al. (2021): spatial DiD on growth rates with regional spillovers;
  - Sun and Xiao (2025): network DiD for US county mask mandates;
  - Zhang and Crawford (2025): time-varying confounding in epidemic interventions.
- **Adaptations.** Theorem 3, Proposition 3 and Theorem 4 are labelled as adaptations (§2.5, p. 8; Web Appendices A.9–A.11).

### W7. Robustness to ascertainment is narrower than implied

Addressed.

- **Scope of Proposition 1.** A paragraph after Proposition 1 (§2.2, p. 5) states that, on the connected US adjacency graph with φ0 > 0, the proposition requires a nationally homogeneous unit reporting factor. It covers testing-time trends, but not county differences in reporting or policy-linked reporting.
- **Tipping point.** The new analysis is in §4 (p. 16) and Web Table 15: an 18.7% treatment-linked drop in reporting would offset the point estimate.
- **Testing data.** Not used (Comment 7.2 gives the reason).

---

## Comment 1: the application after J rejects

**1.1 Region × week effects or region-specific φ.**
- **What we fitted.** We extended the estimator to grouped week effects and grouped coupling. The code was validated against the paper's M0 fit. Specifications:
  - region × week;
  - division × week;
  - region-specific φ;
  - region-specific φ with division × week.
- **Result.** J_LM is 16.93–17.77 in every specification, against calibrated 95th percentiles of 6.38–6.66 (all p < 0.001).
- **Other quantities.**
  - τ+κ ranges from −0.698 to −0.466, and every interval covers zero.
  - Equal regional φ is not rejected (p = 0.505 and 0.630).
  - Regional φ estimates range from 0.27 to 1.08.
- **Q.** In the paper's form Q is undefined when grouped week moments outnumber the 49 clusters, because the estimated Ω is singular. It is reported where it is defined; J_LM is used otherwise.
- **Location.** §4, p. 15; Web Table 10, SI p. 34; Web Appendix F.1.

**1.2 Reframe §4 and the Abstract if rejection persists.**
- Rejection persists, and we reframed accordingly.
- §4 is now titled "Diagnosing a Weekly Transmission Model for State Mask Orders".
- The Abstract ends by stating that −0.47 (95% CI −1.09 to 0.16) summarises the weekly fit and is not a causal effect.
- §4 reads −0.466 as the pseudo-true value of a rejected model.
- **Location.** Unnumbered title page; §1, p. 2; §4, p. 15; §5, p. 18.

**1.3 Report Q against its simulated reference distribution.**
- **Main text (§4, p. 15).** In the application-calibrated design (R = 600):
  - Q_G = 40.64 exceeds every draw (95th percentile 24.93, maximum 35.13);
  - Q_GI = 27.62 is exceeded in 8 of 600 draws;
  - J = 17.87 exceeds every draw.
- **Stress design (R = 300).** Q_G and J exceed every draw; Q_GI is exceeded in 5 of 300.
- **Location.** Web Table 11, SI p. 35, which also gives the scenario 1 and 12 references the referee computed.

**Remark 4 adjacent-cluster check (point in our favour).**
- We reproduced the check with our code. There are 112 adjacent state pairs, with mean score correlations 0.082 (τ), 0.088 (κ) and −0.092 (φ).
- Adding the adjacent cross-products changes the τ+κ standard error from 0.3108 to 0.2915.
- §4 states that the rejected component is the mean model, not the clustered variance.
- **Location.** §4, p. 15; Web Appendix A.11, Remark 4, SI p. 17.

## Comment 2: same-week transmission

**2.1 Implement a sub-weekly renewal mean, or explain why the weekly mean is adequate.**

Not done in this revision.

- **Why it is not done here.**
  - The archive contains daily cumulative JHU counts.
  - Daily county counts carry day-of-week reporting effects and reporting delays. These break the "deterministic reporting without delays" condition of Proposition 1.
  - A daily renewal QD-GMM therefore needs new moment conditions and its own simulation validation. That is a separate methods contribution.
- **Weekly mean adequacy.** We do not claim the weekly mean is adequate for these data.
- **What the paper states instead.**
  - The application lies in the stress regime, where the weekly estimator is biased by about one-sixth of τ+κ and J has no power.
  - A sub-weekly mean is needed, and the paper does not provide one.
- **Location.** §3.2.1, p. 13; §5, p. 18; Web Appendix A.12, SI p. 18. Comment 2.2 is the substitute analysis.

**2.2 Alternative generation-interval discretisations.**
- **Variants.** Ferretti and Nishiura generation intervals, each with:
  - the base discretisation;
  - same-week mass folded into lag 1;
  - infector-time shifts of 1.75, 3.5 and 7 days.
- **Results.**
  - Across the ten variants τ+κ ranges from −0.476 to −0.444, and the φ estimate from 0.19 to 0.82.
  - J rejects in 9 of 10 variants. The exception is the Nishiura 3.5-day shift (J = 0.67, p = 0.715), a single post hoc variant.
  - M0-G departs from M0 to a near-zero root under the 1.75- and 7-day shifts.
  - The existing w1 sensitivity runs (0.6–1.0) give τ+κ between −0.476 and −0.366. They are now connected to the same-week argument.
- **Location.** §4, p. 15; Web Table 12, SI p. 36; Web Appendix F.3.

**2.3 Correct the wording that cites a non-existent extension.** Done. The sentence now says that a sub-weekly mean is needed and is not provided here (§3.2.1, p. 13; §4, p. 15; §5, p. 18; Web Appendix A.12, SI p. 18).

## Comment 3: root selection and boundary inference

**3.1 Scope of Theorem 2.** Done.
- §2.5 (p. 9) states:
  - interior, globally identified cases only;
  - φ0 = 0 violates A6(i);
  - the population second root moves with φ̄;
  - post-selection inference is not covered.
- Web Appendix A.1 (SI p. 2) repeats the A6 condition.
- The population-root check now covers all four zero-contagion designs: the second root appears in 40 of 40 large-N draws (Web Table 8, SI p. 31).

**3.2 Replace the Discussion sentence on global search.** Done (§5, pp. 18-19). The sentence now reads as follows:
- pre-specified M0-G under-covers at φ0 = 0 (52.9%–71.0%);
- post hoc M0-GI repairs this (87.5%–95.2%) but over-selects the boundary when contagion is small and positive. At φ0 = 0.05 coverage drops from 0.960–0.970 to 0.921–0.935;
- both are heuristics, and the selected root should be reported with the profile and the robust sets.

**3.3 Identification-robust confidence sets.** Provided.
- **Method.** We invert profiled Stock–Wright S statistics on a 304-point φ grid. Week effects, τ and κ are profiled out, and the result is projected to τ+κ by Bonferroni.
- **Results at φ̄ = 0.5.**
  - The just-identified set for φ is [0, 0.002] ∪ [0.58, 0.97].
  - With the two over-identifying moments the set is empty, in line with J.
  - The τ+κ projection is [−3.55, 6.29], or [−1.19, 0.26] on the interior component.
- **Location.** §2.6, p. 10, where the set is step (4) of the recommended workflow; §4, p. 16; §5, p. 18; Web Table 13, SI p. 36; Web Appendix F.4.

**3.4 If robust sets are out of scope, restrict claims and present M0-G/M0-GI as heuristics.** Robust sets are provided (3.3). M0-G and M0-GI are also described as selection heuristics without a stated guarantee (§2.6, p. 10; §3.2, p. 12; §5, p. 18).

## Comment 4: event-time results

**4.1 Derive Remark 2; report the implicit weights.** Done (Web Appendix A.10, SI p. 15; Web Table 14, SI p. 37).
- **Derivation.**
  - A first-order expansion of the pooled moment gives pooled τ ≈ Σ_b ω_b τ*_b.
  - ω is the first row of A^(−1) H_b, with A = [H_b 1, H_κ, H_φ].
  - The weights sum to one and can be negative. The pooled κ is also contaminated when bin effects differ.
- **Estimated weights.**
  - At the M0-ET2 estimate they are 1.79, −0.407, −0.177, −0.124, 0.173, −0.101 and −0.154 for bins e0 to e6+.
  - At the pooled estimate they are 2.109, −0.294, −0.333, −0.001, −0.019, 0.009 and −0.471.
  - The pooled estimand is therefore not a convex average of the bin effects.

**4.2 Pooled–dynamic contrast with a joint-sandwich standard error.** Done. Standard errors come from stacked influence functions with 37 clusters and t36 critical values. Results:
- Pooled τ estimate: 0.056 (SE 0.181).
- Pooled minus ω-weighted bins: 0.079 (SE 0.240, p = 0.742).
- Pooled minus equal-weighted bins: 0.733 (SE 0.457, p = 0.117).
- Pooled minus cell-share-weighted bins: 0.899 (SE 0.592, p = 0.138).
- **Location.** §4, p. 16; Web Table 14, SI p. 37.

**4.3 J and the φ-profile for the M0-ET2 fit.** Done on the identical 2,656-county, 37-cluster sample.
- J (df 2) is 14.37 for M0 and 13.85 for M0-ET2, so the dynamic specification does not explain the rejection.
- The φ-profile is reported in Web Table 14.
- **Location.** §4, p. 16; Web Table 14, SI p. 37.

**4.4 Rephrase the Abstract and p. 17.** Done.
- The Abstract no longer mentions dynamic weighting.
- §4 reports that the bins are all negative but imprecise (joint F p = 0.076), and that the divergence is "consistent with dynamic weighting but not distinguishable from sampling error".
- **Location.** §4, p. 16; Table 3 note, p. 31; §5, p. 18.

## Comment 5: finite-cluster inference

**5.1 Application-calibrated simulation of the three interval methods.** Done.
- **Design.** The generator keeps the observed 49-state panel, county counts, adoption dates, network and generation weights. The truth is the M0 estimate. It adds NB2 noise (k = 15.92) and spatially correlated state × week shocks (σ = 0.23).
- **Calibrated design (R = 600).**
  - τ+κ: bias 0.007, SD 0.156, mean SE 0.142.
  - Coverage: t48 0.935, score wild 0.925, pairs 0.933.
- **Stress design (σ = 0.50; R = 300).** This design matches the application's SE of 0.311.
  - Bias 0.015, SD 0.366, mean SE 0.348.
  - Coverage: 0.950, 0.930 and 0.923.
- **τ alone.** All three methods over-cover (0.980–0.993).
- **Assumption.** Both designs assume the working model holds, so they evaluate the interval methods, not the model.
- **Location.** §4, p. 14; Web Appendix E, SI p. 32; Web Table 11, SI p. 35; Web Appendix F.2.

**5.2 Declare the primary interval.**
- The cluster sandwich with t48 critical values is primary, because it is the best calibrated of the three in both designs.
- The same sentence states that the designation was made in revision, after the calibrated simulation, and is post hoc.
- **Location.** §4, p. 14.

**5.3 Few-treated-cluster and heteroskedastic-cluster inference.**
- Ferman and Pinto (2019), Conley and Taber (2011), MacKinnon and Webb (2018) and Canay, Santos and Shaikh (2021) are cited.
- The text notes that the wild bootstrap's validity depends on conditions on the treated clusters that may not hold here. Jackknife and pairs intervals are therefore reported as sensitivity checks.
- **Location.** §5, p. 18; Web Appendix E, SI p. 33.

**5.4 Lemma 3 is a large-G device.**
- Under A4′ the dependence radius of Lemma 3 is at least 28 cluster steps. This exceeds the diameter of 11 of the 49-state graph.
- (L) is therefore an assumption in the application.
- **Location.** §2.5, p. 8; Web Appendix A.2, SI p. 4.

## Comment 6: novelty and related work

**6.1 State the contribution as the integration with its diagnostics.** Done (§1, p. 2).

**6.2 Cite and contrast Kosfeld et al., Sun and Xiao, and Zhang and Crawford.** Done (§1, p. 1; §4, p. 14).
- Both Kosfeld et al. and Sun and Xiao work on the outcome scale, without renewal dynamics or feedback adoption.
- Zhang and Crawford is cited with Imai and Kim for the feedback problem.
- **No numerical comparison.** The estimands differ: growth-rate spatial DiD and network DiD on case outcomes, versus log-R effects. The comparison is therefore qualitative.

**6.3 Mark Theorem 3, Proposition 3 and Theorem 4 as adaptations.** Done.
- **Theorem 3.** Adapts known results on predetermined regressors in count panels (Hausman, Hall and Griliches 1984; Wooldridge 1997, 1999; Blundell, Griffith and Windmeijer 2002). It is stated to provide a common footing for the comparators, not new inconsistency results.
- **Proposition 3.** Adapts Goodman-Bacon (2021) and de Chaisemartin and D'Haultfœuille (2020).
- **Remark 2.** Adapts Sun and Abraham (2021) and Borusyak, Jaravel and Spiess (2024).
- **Theorem 4(b).** Labelled a standard delta-method consequence of Theorem 2. None of the four references the report suggests is its source, so we did not attach them there.
- **Location.** §2.5, p. 8; Web Appendices A.9–A.11, SI pp. 13–16.

**6.4 Move comparator-bias material to the supplement.** Done.
- §3.2 keeps one sentence of comparator biases.
- M1-BC, M3 and Oracle results and the misspecification detail are in Web Appendix C.
- **Location.** §3.2, p. 12; Web Appendix C, SI p. 20.

## Comment 7: case ascertainment

**7.1 Application-level sensitivity to treatment-linked reporting.** Done.
- **Method.** Under Y_obs = c_i c_t (1−ρ)^D X, we corrected counts, rebuilt both pressures and refitted M0 over ρ ∈ [−0.6, 0.6].
- **Results.**
  - τ+κ reaches zero at ρ = 0.187, an 18.7% treatment-linked drop in reporting.
  - The upper 95% limit reaches zero at ρ = −0.076, and the lower limit at ρ = 0.410.
  - J rejects over the whole grid.
- **Location.** §4, p. 16; §5, p. 18; Web Table 15, SI p. 38; Web Appendix F.6.

**7.2 Test positivity or testing volume as a check.** Not done.
- **Reasons.**
  - The archive holds no testing data.
  - County-level testing volume for 2020 is not consistently available.
  - State-level series would need a separate measurement model linking testing to the reporting factor c_it.
- **Substitute.** The tipping point of 7.1 shows how large a reporting change would have to be.
- **Future data.** Web Appendix F.6 notes wastewater surveillance as a reporting-independent signal for future applications.

**7.3 State the national-homogeneity implication of Proposition 1 in the main text.** Done (§2.2, p. 5).

---

## Minor comments

**m1.** The jackknife range is corrected to 14%–30%. The archived excesses are 30.1%, 13.6% and 16.5% (§4, p. 14).

**m2.** All five M0-G departures are now listed: policy lag (+0.421), ω×4 (+0.395), and w1 = 0.9, 0.7 and 0.6 (+0.013, −0.260 and −0.261). The text adds that the policy-lag run is the only one whose M0 interval excludes zero (−0.764; −1.346 to −0.182), while M0-G gives +0.421 there (§4, p. 16; Web Table 9 note, SI p. 32).

**m3.** Null scenarios 17–20 are added to the list of designs flagged for φ estimation. Their φ bias is at most 0.0027 (§3.2, p. 12).

**m4.** Population-ψ coverage (72.1%–92.5% under spillover) is reported next to ψ† coverage (85.7%–95.4%), with the reason: the intervals omit process noise. The application's four-week estimate is labelled ψ† = 341,000 (SE 397,000) (§3.2, p. 12; §4, p. 16).

**m5.** The Proposition 3 summary now states that (κ, φ) are held at their true values and that the weights are first-order (§2.5, p. 8).

**m6.** Theorem 3(b) is presented as an illustration of the classical fixed-T inconsistency. The text notes that the bias is −0.004 in the randomised design (§2.5, p. 8; Web Appendix A.9, SI p. 14).

**m7.** "Crosses zero cleanly" is removed.
- The text now reports that the profile is also near zero at the boundary (about −0.0009, where the τ and κ estimates are −2.00 and 3.37).
- Q is given at both candidates: Q_G 40.64 versus 41.34, and Q_GI 27.62 versus 40.32.
- Figure 3 was not regenerated. Its caption states the boundary candidate and Q at both candidates.
- **Location.** §4, p. 14; Figure 3 caption, p. 28.

**m8.** §2.6 now states that the solver requires only Λ > 0. Fits with a negative φ estimate are reported as such; they occurred in about 0.1% of replications in one zero-contagion design (§2.6, p. 10).

**m9.** The score wild bootstrap perturbs the cluster influence terms with the estimated Γ held fixed and does not impose the null. The pairs bootstrap re-solves the estimator. Both are now evaluated in the calibrated simulation, where wild coverage of τ+κ is 0.925 and 0.930 (§4, p. 14; Web Appendix E; Web Table 11).

**m10.** Both statements were correct. The code uses a 0-based window index 2 (ENGINE_FROM = 2), which is the third window week. The text and the verification README now say so (§3.2.1, p. 13; Web Appendix A.12).

**m11.** Pre-specified and post hoc elements are kept visibly distinct.
- One list in §3.1 (p. 12) gives the five protocol departures and names the revision analyses as post hoc.
- §4 (p. 14) labels the primary-interval designation and the new application analyses as post hoc.
- Table 2 (p. 30) and §5 (p. 18) use the same labels.

**m12.** A2(b) is linked to sequential ignorability in the time-series cross-section literature (Blackwell and Glynn 2018) (§2.2, p. 5).

**m13. Partly done.**
- Chen and Roth (2024) is cited where the log(Y+1) TWFE comparator is analysed (Web Appendix A.9, SI p. 14).
- Adding it at the M3 description in §3.1 pushed the references to a 26th counted page, beyond the 25-page limit, so the main-text citation was not kept.

**m14.** §4 (p. 14) now states that state-level orders make 61.9% of counties interior, so τ and κ are separately identified only from 1,177 border counties. τ+κ is stable across sensitivity runs (−0.764 to −0.366), while the τ estimate is positive in the main run and four others. The total is therefore the quantity the design identifies.

**m15.** Holtz et al. (2020) is cited at A2(b). Coordination effects from interdependent state responses are part of the estimand through κ (§2.2, p. 5).

## Typographical and presentation issues

**t1.** Corrected to "Web Table 5" (§3.2, p. 12).

**t2.** M1-BC failures are now given per scenario: 85.3% (scenario 2), 76.1% (scenario 4), 3.6% (scenario 6) and 0.1% (scenario 1) of 8000 runs, with convergence in all other scenarios. Coverage of 25.7%–82.7% is over converged runs (Web Appendix C, SI p. 21).

**t3.** The Table 2 footnote is rewritten in the table file and in the table generator (Table 2, p. 30).

**t4. Partly done.** Author names, affiliations and the funding statement must come from the authors. The placeholders remain in the title block and the Funding section (Back matter, p. 19) and will be completed before submission.

**t5.** The duplicate key callaway2023pandemic is removed. The previously uncited entries are now cited where they belong:
- stock2000weak and andrews2019weak at the robust sets;
- holtz2020interdependence at A2(b);
- chen2024logs in Web Appendix A.9;
- xu2023did and haber2021impact in the Introduction.

All 98 bibliography entries are cited.

**t6.** Same change as Comment 6.4. The freed space carries the application diagnostics (§3.2, p. 12; §4, p. 15; Web Appendix C, SI p. 20).

---

## Questions in Section 4.4

1. **Regional effects; causal statement.**
   - Regional week effects and region-specific φ do not bring J back into its calibrated range (J_LM 16.9–17.8 against 6.38–6.66).
   - We therefore make no causal statement about mask orders.
   - −0.466 is reported as the pseudo-true value of the rejected weekly model (Comments 1.1–1.3; §4, p. 15; Web Table 10).
2. **Sub-weekly mean; reassignment.**
   - The sub-weekly mean is not implemented here (Comment 2.1 gives the reason).
   - Reassigning same-week mass moves τ+κ only between −0.476 and −0.444 (Comment 2.2; §4, p. 15; Web Table 12).
3. **Identification-robust sets.** Provided for φ and τ+κ (Comment 3.3; §4, p. 16; Web Table 13).
4. **Implicit weights; contrast SE.** Derived and estimated. The contrasts are 0.079 (SE 0.240) and 0.733 (SE 0.457) (Comments 4.1–4.2; §2.4, p. 8; §4, p. 16; Web Table 14).
5. **Interval methods; primary interval.**
   - Coverage of τ+κ is 0.935 (t48), 0.925 (wild) and 0.933 (pairs) in the calibrated design.
   - The t48 sandwich is primary, a designation made post hoc (Comments 5.1–5.2; §4, p. 14; Web Table 11).
6. **Comparison with Sun and Xiao and with Kosfeld et al.** Qualitative, because the estimands differ (Comment 6.2; §1, p. 1; §4, p. 14).
7. **Reporting change that explains the estimate.** An 18.7% treatment-linked drop in reporting (Comment 7.1; §4, p. 16; Web Table 15).
8. **Requested corrections.** m1–m4, t1–t3 and t5 are done. t4 awaits the authors' information (see above).

---

## Items to be completed by the authors before submission

- **Author names and affiliations.** The title block currently carries a placeholder.
- **Funding statement.** The Funding section currently carries a placeholder.
- **Corresponding author and contact details** for this letter.

[Corresponding author name, affiliation and contact details]
