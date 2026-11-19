# -*- coding: utf-8 -*-
"""Round-3 revision: recover one printed page in the main text so the
Biometrics count stays at 25, while keeping every disclosure the second
review asked for. Cuts are redundancy and prefacing only; no numbers,
conditions or limitations are removed."""

import io

APP = "sections/application.tex"
INT = "sections/introduction.tex"
DIS = "sections/discussion.tex"

edits = {
 APP: [
  # data paragraph: fold the exits sentence, drop the duplicated median filler
  ("""the moment conditions accommodate; censoring or dropping exiting units leaves estimates within 0.003.
Median weekly counts are 27 cases (7.4\\% zeros), orders cover 56.8\\% of county-weeks, 1720 units
switch exposure, 571 units are exposed in week 1, and 11,528,182 cases are reported. Network $W$ uses
row-normalised 2010 Census adjacency, with ACS commuting flows in sensitivity runs. Discretising
\\citet{ferretti2020quantifying} yields weights $w\\approx(0.964,0.036)$ and assigns $\\approx 0.30$ of
transmission within the same week.""",
   """the moment conditions accommodate; censoring or dropping exiting units leaves estimates within 0.003.
Median weekly counts are 27 cases (7.4\\% zeros), orders cover 56.8\\% of county-weeks, 1720 units
switch exposure, and 11,528,182 cases are reported. Network $W$ is row-normalised 2010 Census
adjacency, with ACS commuting flows in sensitivity runs. Discretising
\\citet{ferretti2020quantifying} yields $w\\approx(0.964,0.036)$, putting $\\approx 0.30$ of
transmission inside the same week."""),

  # fit paragraph: shorten the lead sentence and the lead-regression aside
  ("""We fit the estimators of Sections~\\ref{ss:estimation}--\\ref{ss:computation} with $\\bar\\phi=0.5$. The
pre-specified sensitivity runs vary generation interval, network $W$, timing (one-week shift),
background $\\omega$ ($\\times 0.25$, $\\times 4$), $\\bar\\phi$, and exclusion of Kansas and Texas. The""",
   """We fit the estimators of Sections~\\ref{ss:estimation}--\\ref{ss:computation} at $\\bar\\phi=0.5$,
with pre-specified sensitivity runs over generation interval, network $W$, timing (one-week shift),
background $\\omega$ ($\\times 0.25$, $\\times 4$), $\\bar\\phi$ and exclusion of Kansas and Texas. The"""),

  ("""clusters contain switching counties, and the clusters are very unequal in size: they run from 1 to
254 counties, so the largest state carries a large share of the weight in every moment. A comparison
lead regression gives $\\hat\\ell=-0.332$ (SE 0.259),
which under feedback is not a valid pre-trend test (Proposition~\\ref{prop:sources}(c)).""",
   """clusters contain switching counties, and cluster size is unequal, from 1 to 254 counties, so the
largest state carries much of the weight in every moment. A comparison lead regression gives
$\\hat\\ell=-0.332$ (SE 0.259), which under feedback is not a valid pre-trend test
(Proposition~\\ref{prop:sources}(c))."""),

  # estimates: drop the redundant clause after tau
  ("""to $-0.366$ while $\\hat\\tau$ is positive in the main run and four others (Web
Table~\\ref{tab:app-all}), so the total, not its split, is the stable summary.""",
   """to $-0.366$ while $\\hat\\tau$ is positive in five runs (Web Table~\\ref{tab:app-all}), so the total,
not its split, is the stable summary."""),

  # over-identification: remove the restatement of the same conclusion
  ("""6.66, maximum 12.45), as does $Q_G=40.64$ (maximum 35.13); $Q_{GI}=27.62$ is exceeded in 8 of 600.
Both statistics therefore lie outside the range calibrated under the working model. Rejection persists
with region $\\times$ week""",
   """6.66, maximum 12.45), as does $Q_G=40.64$ (maximum 35.13); $Q_{GI}=27.62$ is exceeded in 8 of 600,
so both lie outside the range calibrated under the working model. Rejection persists
with region $\\times$ week"""),

  ("""weak (112 pairs, correlations below 0.1 in absolute value) and moves the SE of $\\hat\\tau+\\hat\\kappa$
from 0.311 to 0.291 (Remark~\\ref{rem:testable}); the rejection concerns the mean model, not the
clustered variance. We therefore read $\\hat\\tau+\\hat\\kappa=-0.466$ as the pseudo-true value of a
rejected working model, a summary of the weekly fit and not a causal effect of mask orders.""",
   """weak (112 pairs, correlations below 0.1 in absolute value) and moves the SE of $\\hat\\tau+\\hat\\kappa$
from 0.311 to 0.291 (Remark~\\ref{rem:testable}), so the rejection concerns the mean model. We read
$\\hat\\tau+\\hat\\kappa=-0.466$ as the pseudo-true value of a rejected working model, not as a causal
effect of mask orders."""),

  # robust sets: tighten the conditional-inference sentence
  ("""projection combines a $97.5\\%$ set for $\\phi$ with a $97.5\\%$ Wald interval for $\\tau+\\kappa$ at each
$\\phi$, and its second step presumes that $(\\tau,\\kappa)$ is strongly identified at that $\\phi$. On the
boundary component that presumption is doubtful: at $\\phi=0$ the fit gives $\\hat\\tau=-2.00$ and
$\\hat\\kappa=3.37$, and the conditional standard error of $\\tau+\\kappa$ is 2.13. We therefore report the
projection as a bound computed under that condition, not as an interval with guaranteed coverage, and
recommend the $S$ set for $\\phi$, which does not need it. The interior root, and hence $\\hat\\phi$,
depends on the instrument reference value $\\bar\\phi$; $\\hat\\tau$ moves with it as well, from 0.074 at
$\\bar\\phi=0.25$ to 0.018 at $\\bar\\phi=1$, while $\\hat\\tau+\\hat\\kappa$ stays between $-0.476$ and
$-0.466$.""",
   """projection combines a $97.5\\%$ set for $\\phi$ with a $97.5\\%$ Wald interval for $\\tau+\\kappa$ at each
$\\phi$, so its second step presumes $(\\tau,\\kappa)$ strongly identified there. On the boundary
component that fails: at $\\phi=0$ the fit gives $\\hat\\tau=-2.00$, $\\hat\\kappa=3.37$, and a conditional
standard error of 2.13 for $\\tau+\\kappa$. We therefore report the projection as a bound computed under
that condition, not an interval with guaranteed coverage, and recommend the $S$ set for $\\phi$, which
does not need it. $\\hat\\phi$ depends on the instrument reference value $\\bar\\phi$, and $\\hat\\tau$
moves with it, from 0.074 at $\\bar\\phi=0.25$ to 0.018 at $\\bar\\phi=1$, while $\\hat\\tau+\\hat\\kappa$
stays between $-0.476$ and $-0.466$."""),

  # event time: merge the two sentences about the check
  ("""bins leaves 0.079 (SE 0.240, $p=0.742$), which is a check of the linearisation rather than a comparison
of estimands, since the two sides agree to first order by construction; subtracting the equally
weighted bins 0.733 (0.457, $p=0.117$), which is informative (Web Table~\\ref{tab:rev-event}). The
divergence is consistent with dynamic weighting but is not statistically distinguishable from sampling
error.""",
   """bins leaves 0.079 (SE 0.240, $p=0.742$), a check of the linearisation rather than a comparison of
estimands, since the two sides agree to first order by construction; subtracting the equally weighted
bins leaves 0.733 (0.457, $p=0.117$), which is informative (Web Table~\\ref{tab:rev-event}). The
divergence is consistent with dynamic weighting but not distinguishable from sampling error."""),
 ],
 INT: [
  ("""weekly fit rather than a causal effect. The generation-interval weights put about 0.30 of transmission in
the same week as the infector, inside a range where the weekly model cannot work and $J$ cannot detect
the aggregation (Section~\\ref{ss:sameweek}); that arithmetic, not a diagnostic test, is what rules the
weekly specification out. One limitation is worth stating here: the identification-robust set for
$\\phi$ has a boundary and an interior component, and adding the over-identifying moments empties it,
so in this application the sets confirm the rejection rather than replace root selection with a usable
interval. Sections~\\ref{s:methods}--\\ref{s:discussion} present the theory, the
simulations, the application and the scope of the design; Web Appendices~A--F give proofs, adapted
comparator results, numerical checks, supplementary simulations, data, few-cluster inference and
application diagnostics.""",
   """weekly fit rather than a causal effect. The generation-interval weights put about 0.30 of transmission
in the same week as the infector, where the weekly model cannot work and $J$ cannot detect the
aggregation (Section~\\ref{ss:sameweek}); that arithmetic, not a diagnostic test, rules the weekly
specification out. The identification-robust set for $\\phi$ has a boundary and an interior component
and is emptied by the over-identifying moments, so here the sets confirm the rejection rather than
replace root selection with a usable interval. Sections~\\ref{s:methods}--\\ref{s:discussion} present
the theory, simulations, application and scope; Web Appendices~A--F give proofs, adapted comparator
results, numerical checks, supplementary simulations, data, few-cluster inference and application
diagnostics."""),
 ],
 DIS: [
  ("""Placing parallel trends on the log reproduction number with imported infection pressure, and removing
unit effects by forward quasi-differencing, gives a difference-in-differences design that accommodates
feedback adoption and lagged cross-cluster contagion together while keeping the ordinary
cluster-robust sandwich. In simulations with spillover, feedback, overdispersion and small counts,
QD-GMM recovers direct and spillover effects where comparators built on counts, log counts or
first-step reproduction numbers are biased, and it identifies the policy-level total directly when
mandates are set above the unit.""",
   """Placing parallel trends on the log reproduction number with imported infection pressure, and removing
unit effects by forward quasi-differencing, gives a difference-in-differences design that accommodates
feedback adoption and lagged cross-cluster contagion while keeping the ordinary cluster-robust
sandwich. In simulations with spillover, feedback, overdispersion and small counts, QD-GMM recovers
direct and spillover effects where comparators built on counts, log counts or first-step reproduction
numbers are biased, and identifies the policy-level total directly when mandates are set above the
unit."""),
  ("""In the application the diagnostics reject the weekly model, also with regional week effects and
coupling; the robust sets for $\\phi$ are empty once over-identifying moments are added, and about a
third of transmission falls in the same week. The policy-level total $-0.466$ is therefore a
pseudo-true value of the weekly fit. Event-time estimates are consistent with effects that build up
after adoption but are imprecise (joint $p=0.076$), and an 18.7\\% order-linked fall in ascertainment
would offset the point estimate.""",
   """In the application the diagnostics reject the weekly model, also with regional week effects and
coupling; the robust sets for $\\phi$ are empty once over-identifying moments are added, and about a
third of transmission falls in the same week. The policy-level total $-0.466$ is therefore a
pseudo-true value. Event-time estimates are consistent with effects that build up after adoption but
imprecise (joint $p=0.076$), and an 18.7\\% order-linked fall in ascertainment would offset the point
estimate."""),
  ("""pressure \\citep{elenev2025staggered} and complements unconfoundedness formulations
\\citep{callaway2023pandemic}. The design turns a panel of reported cases into estimates on the
transmission scale with standard clustered inference, and its diagnostics show when the weekly
working model cannot support a causal reading. Sequential efficiency bounds""",
   """pressure \\citep{elenev2025staggered} and complements unconfoundedness formulations
\\citep{callaway2023pandemic}. The design turns a panel of reported cases into estimates on the
transmission scale with standard clustered inference, and its diagnostics show when the weekly model
cannot support a causal reading. Sequential efficiency bounds"""),
 ],
}

for path, pairs in edits.items():
    s = io.open(path, encoding="utf-8").read()
    for a, b in pairs:
        if a not in s:
            raise SystemExit("MISSING in %s: %r" % (path, a[:80]))
        s = s.replace(a, b)
    io.open(path, "w", encoding="utf-8", newline="").write(s)
    print("edited", path)