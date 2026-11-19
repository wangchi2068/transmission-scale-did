# Response to the Second-Round Referee

**Manuscript:** "Difference-in-differences on the transmission scale: sequential moment estimation of policy effects on the reproduction number under feedback and spatial contagion"

**Journal:** Biometrics

**Authors:** [author names to be supplied by the authors]

We thank the referee for the second report. Its conclusion was that the remaining problems could be answered by changes of wording rather than new estimation, and we have followed that in every case. No estimate, standard error, coverage figure or test statistic changed in this round; all changes are to how the existing results are described, plus two rendering faults and one double-rounding error. Below we answer the six questions of §4.4 and the five minor items of §4.2 in order.

**Conventions.** "p." is the printed page of the main manuscript (25 of 25 Biometrics-counted pages; figures and tables follow the references). "Web p." is a page of the Supporting Information (43 pages). The abstract sits on the unnumbered title page.

## Main issues

### M1. The discussion overstated what the calibrated simulation supports (§3.1 M1; p. 18)

**Agreed, and corrected.** The sentence now carries the condition and the same-week figure:

> "with 49 states, of which 24 switch, an application-calibrated simulation **in which the weekly model holds** supports $t_{G-1}$ sandwich intervals, with coverage 0.935 for τ+κ in the base design and 0.950 once the state-by-week shocks are scaled up …, whereas inference with few treated clusters can fail more generally. **That calibration does not cover the same-week regime this application occupies, where the same interval covers 0.783** and J has no power (Section 3.2.1)."

The 0.783 is the τ+κ coverage reported in §3.2.1 for the w0 = 0.3 regime (Web Table 7). We did not add a same-week arm to the application-calibrated simulation. Doing so would require a daily generator for the observed 49-state panel and its adoption dates, which is a different simulation study from the one the referee's own C5.1 asked for, and we preferred a correctly conditioned sentence to a new design run in the revision window.

### M2. Two sentences in the introduction exceeded the diagnosis (§3.1 M2; p. 3)

**Agreed, both rewritten.**

- Same-week regime: the sentence now reads "The generation-interval weights put about 0.30 of transmission in the same week as the infector, where the weekly model cannot work and J cannot detect the aggregation (Section 3.2.1); **that arithmetic, not a diagnostic test, rules the weekly specification out**." The claim about what identifies the regime is now attached to the generation-interval calculation, which is where it comes from.
- Robust sets: "identification-robust sets avoid root selection" became "identification-robust sets **for φ** do not depend on the selected root", and the application paragraph now states that the set has a boundary and an interior component and is emptied by the over-identifying moments, "so here the sets confirm the rejection rather than replace root selection with a usable interval."

### M3. The validity of the τ+κ projection needed stating (§3.1 M3; p. 17; Web pp. 34, 36)

**Agreed, stated in the text.** The paper now says the projection combines a 97.5% set for φ with a 97.5% Wald interval at each φ, so its second step presumes (τ, κ) strongly identified there, and that on the boundary component this fails: at φ = 0 the fit gives τ̂ = −2.00, κ̂ = 3.37, and the conditional standard error of τ+κ is 2.13 (source: `robust_sets/results.json`, `S1_at_0.tpk_se_phi_fixed` = 2.1276). The projection is described as "a bound computed under that condition, not an interval with guaranteed coverage", and the φ S set, which does not need the condition, is recommended in its place.

We did not compute a fixed-τ+κ S-type statistic. The referee's own assessment was that the projection "is not therefore invalid, because it is wide anyway"; a conditional statement plus a pointer to the set that does not need the condition answers the question without a new inferential procedure.

### M4. Event-time: the ω-weighted contrast is near zero by construction (§3.1 M4; p. 17)

**Agreed, both changes made.** The text now gives the non-convex weights in the main manuscript — "2.109 on the adoption bin and negative thereafter, down to −0.471 for e6+" (source: `event_time/results.json`, `implicit_weights.at_pooled_estimate.omega`) — and labels the ω-weighted contrast as a check: "a check of the linearisation rather than a comparison of estimands, since the two sides agree to first order by construction". The equally weighted contrast is described as the informative one.

### M5. Two disclosures were incomplete in the main text (§3.1 M5; pp. 14, 17)

**Agreed, both added.**

- Unequal cluster size is now in §4's data paragraph: "cluster size is unequal, from 1 to 254 counties, so the largest state carries much of the weight in every moment." This sits beside the 24-of-49 switching-cluster statement, where the reader meets the few-cluster issue.
- The movement of τ̂ with φ̄ is now in the same paragraph as the movement of φ̂: "τ̂ moves with it, from 0.074 at φ̄ = 0.25 to 0.018 at φ̄ = 1, while τ̂+κ̂ stays between −0.476 and −0.466." The three values are the minimum-S1 grid points in `robust_sets/results.json` for the three reference values and match Web Table 13.

### Abstract scope: "also with regional week effects"

**Changed to "division × week effects".** R4 (region × week) has no calibrated reference in Web Table 10 because the calibration design does not carry regional week effects; A1, A3 and A2 do, and their maximum 95th percentile is 6.66. The abstract now names the specification that has a reference. The conclusion is unchanged: R4's J_LM = 17.77 exceeds every calibrated maximum in the table.

### J_LM's inverted version at the 5% edge

**Noted in Web Table 10's note.** The note now says the two versions carry different evidence, that the invariance version sits at the 5% edge (0.046 to 0.074) and is exceeded in 31 of 600 calibrated draws, and that the rejection rests on the over-identifying version.

## Minor items (§4.2)

1. **Double rounding of 0.292.** Corrected to 0.291 (source 0.291479), in the main text and the Chinese copy.
2. **Table 3's note.** Now reads "mean of e0–e5 −0.677", which is the quantity the row reports. The seven-bin mean is −0.712 and is not quoted.
3. **"Web Appendix C.1" not printed.** The two subsubsection headings in Web Appendix C now carry their numbers, "C.1 Validation of the root-selection rules on independent designs" and "C.2 Population roots in the zero-contagion designs".
4. **Code variable in Web Appendix A.12.** Replaced by prose: "from the third window week, that is, from index 2 counting the first window week as 0."
5. **Overfull line on Web p. 33.** The two long `\texttt` paths were replaced by "the qdx module of the replication package"; the Supporting Information now compiles with no overfull box.
6. **Chen & Roth cited only in the appendix.** Unchanged, for the page-limit reason the referee accepted. Status remains partial in the referee's table.
7. **Author and funding placeholders.** Unchanged; these need the authors.

## Typographic and presentation

1. The response letter's "Abstract, p. 1" now reads "unnumbered title page", here and in the earlier round's letter.
2. The new sentences and the cuts that make room for them leave a net addition of about 170 words across Sections 1–5, which still fits inside 25 of 25 counted pages. No number, condition or disclosure was removed. The acknowledgements and data-availability statements now share a paragraph, and the conflict-of-interest and funding statements share a heading.

## What changed, in one line each

| Location | Change |
|---|---|
| Abstract | "regional week effects" → "division × week effects"; the empty φ set is now named |
| §1 | same-week regime attributed to the generation-interval arithmetic; robust-set sentence scoped to φ and to what it delivers here |
| §4 data | unequal cluster sizes (1 to 254 counties) added |
| §4 over-identification | 0.292 → 0.291 |
| §4 robust sets | conditional-identification statement for the projection, with the φ = 0 numbers; τ̂ movement across φ̄ added |
| §4 event time | non-convex implicit weights in the text; ω-weighted contrast labelled a linearisation check |
| §5 | calibration sentence conditioned on the working model, with the same-week coverage 0.783 |
| Web Table 10 note | two J_LM versions carry different evidence |
| Web Table 11 | "stress design" → "high-noise design" |
| Web Appendix A.12 | code variable replaced by prose |
| Web Appendix C | C.1 and C.2 numbers printed |
| Web Appendix F | long paths removed; overfull box gone |
| Chinese copy | all of the above, plus one typo |

**Sources.** Every figure quoted here is from the revision archive under `verification/revision/` (`robust_sets`, `event_time`, `regional_fe`, `calibrated_sim`). The same-week coverage is from `verification/same_week/results.json` and is reported in §3.2.1 of the main text and Web Table 7.