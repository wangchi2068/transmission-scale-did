# Plan for boundary-robust inference (reviewer M1)

The simulation shows two separate failures at phi0 = 0: the local solver may select a
non-consistent root, and the Wald interval is unreliable even at the correct root (the oracle
covers 80.2% to 94.5% at T = 26 with random adoption). Root selection addresses the first only.

Proposed procedure, implemented on top of the existing profile:

1. For each phi on the implementation grid, profile (gamma, tau, kappa) out and record the GMM
   objective Q(phi) using the M0-OI instruments.
2. The profile-based confidence set for (tau, kappa) is the projection
   C = {(tau, kappa) : min_phi Q(phi; tau, kappa) <= chi2_{q-p, 1-alpha}}.
   At kappa fixed this is an Anderson-Rubin-style test that does not assume the boundary
   asymptotics of the Wald interval.
3. Report C for phi0 = 0 and phi0 = 0.05 scenarios, alongside the Wald interval, in the Web
   Appendix. If C is unbounded in tau, report that as the outcome rather than a number.

Status: the procedure is specified and the profile objects it needs already exist in
`experiment/code/qdgmm.py` (`fit_global`, `phi_profile`). The Monte Carlo evaluation has not been
run; the manuscript therefore states that the boundary case needs inference robust to weak
identification instead of claiming a fix.
