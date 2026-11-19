# -*- coding: utf-8 -*-
"""Round-3 revision: recover the remaining page overflow. Only redundant
phrasing and repeated definitions are cut; every number, condition and
disclosure stays."""

import io

edits = {
 "sections/backmatter.tex": [
  ("""All data are publicly available; sources, checksums and download instructions are in Web Appendix D
and \\texttt{experiment/data}. Reproduction code accompanies this article.""",
   """All data are publicly available; sources, checksums and download instructions are in Web Appendix D.
Reproduction code accompanies this article."""),
  ("""Web Appendices A--F (regularity conditions, proofs and adapted comparator results, numerical checks,
supplementary simulations, data sources and cleaning, few-cluster inference, application diagnostics), referenced in Sections~\\ref{s:intro}--\\ref{s:discussion},
are available as Supporting Information accompanying this article.""",
   """Web Appendices A--F (regularity conditions, proofs, adapted comparator results, numerical checks,
supplementary simulations, data sources, few-cluster inference and application diagnostics), referenced
in Sections~\\ref{s:intro}--\\ref{s:discussion}, are available as Supporting Information."""),
 ],
 "sections/simulation.tex": [
  ("""The design follows \\citet{morris2019ademp}; aims, generators, estimands, comparators and success
criteria were pre-specified before the full runs.""",
   """The design follows \\citet{morris2019ademp}; aims, generators, estimands, comparators and success
criteria were pre-specified."""),
  ("""Counties sit uniformly on the unit square, linked to their six nearest neighbours via a symmetrised,
row-normalised adjacency matrix $W$ and grouped into 50 contiguous state clusters. Unit effects
combine state, county and adoption-propensity components, correlating adoption with baseline
transmission. Outcomes follow \\eqref{eq:mean} with $S=2$, generation weights $(0.85,0.15)$,
population-proportional background input $\\omega$, negative binomial ($k=10$) or Poisson ($k=\\infty$)
innovations, and week effects calibrated for median county-week counts between 5 and 50. Under
random adoption, states draw adoption weeks uniformly on $\\{2,\\dots,T+2\\}$; under feedback adoption,
states and counties adopt sequentially once lagged incidence exceeds propensity-dependent thresholds.
Adoption is absorbing throughout.""",
   """Counties sit uniformly on the unit square, linked to their six nearest neighbours by a symmetrised,
row-normalised adjacency $W$ and grouped into 50 contiguous state clusters. Unit effects combine
state, county and adoption-propensity components, correlating adoption with baseline transmission.
Outcomes follow \\eqref{eq:mean} with $S=2$, generation weights $(0.85,0.15)$, population-proportional
background $\\omega$, negative binomial ($k=10$) or Poisson ($k=\\infty$) innovations, and week effects
calibrated for median county-week counts of 5 to 50. Under random adoption, states draw adoption weeks
uniformly on $\\{2,\\dots,T+2\\}$; under feedback, states and counties adopt sequentially once lagged
incidence exceeds propensity-dependent thresholds. Adoption is absorbing."""),
  ("""Estimands are generating $(\\tau,\\kappa,\\phi)$, the switcher mean of $\\tau_i$ in X3, event-time
effects in X5, and $\\psi^\\dagger$ of Theorem~\\ref{thm:psi} (with population $\\psi$ approximated from
5000 paired paths). On identical $R=8000$ data sets without replacement of failed runs, we evaluate""",
   """Estimands are generating $(\\tau,\\kappa,\\phi)$, the switcher mean of $\\tau_i$ in X3, event-time
effects in X5, and $\\psi^\\dagger$ of Theorem~\\ref{thm:psi} (population $\\psi$ from 5000 paired paths).
On identical $R=8000$ data sets, without replacing failed runs, we evaluate"""),
  ("""true $\\alpha_i$. Intervals use $t_{G-1}$ critical values clustered by state. Success required bias
within 2 MCSE of zero for $(\\tau,\\kappa,\\phi)$ and coverage in $[93\\%,97\\%]$ under feedback.""",
   """true $\\alpha_i$. Intervals use $t_{G-1}$ critical values clustered by state, and success required
bias within 2 MCSE of zero for $(\\tau,\\kappa,\\phi)$ and coverage in $[93\\%,97\\%]$ under feedback."""),
  ("""sets. The validation designs of Web Appendix~\\ref{app:rsel}, the population root checks, the
application-calibrated simulation and the application analyses added in revision are also post hoc.""",
   """sets. The validation designs of Web Appendix~\\ref{app:rsel}, the population root checks, the
application-calibrated simulation and the revision's application analyses are also post hoc."""),
 ],
 "sections/methods.tex": [],
}

for path, pairs in edits.items():
    if not pairs:
        continue
    s = io.open(path, encoding="utf-8").read()
    for a, b in pairs:
        if a not in s:
            raise SystemExit("MISSING in %s: %r" % (path, a[:80]))
        s = s.replace(a, b)
    io.open(path, "w", encoding="utf-8", newline="").write(s)
    print("edited", path)