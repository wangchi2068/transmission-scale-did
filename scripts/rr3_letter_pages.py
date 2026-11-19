# -*- coding: utf-8 -*-
"""Align the round-3 response letter's page references with the built PDF."""

import io

P = "D:/codeg_work/论文相关/peer-review-of-transmission-scale-did-f4dda8b8/response-to-reviewers-r3.md"
s = io.open(P, encoding="utf-8").read()

pairs = [
 ("### M1. The discussion overstated what the calibrated simulation supports (§3.1 M1; p. 19)",
  "### M1. The discussion overstated what the calibrated simulation supports (§3.1 M1; p. 18)"),
 ("### M3. The validity of the τ+κ projection needed stating (§3.1 M3; p. 18; Web p. 34)",
  "### M3. The validity of the τ+κ projection needed stating (§3.1 M3; p. 17; Web pp. 34, 36)"),
 ("### M4. Event-time: the ω-weighted contrast is near zero by construction (§3.1 M4; p. 18)",
  "### M4. Event-time: the ω-weighted contrast is near zero by construction (§3.1 M4; p. 17)"),
 ("### M5. Two disclosures were incomplete in the main text (§3.1 M5; p. 4, p. 18)",
  "### M5. Two disclosures were incomplete in the main text (§3.1 M5; pp. 14, 17)"),
]
for a, b in pairs:
    if a not in s:
        raise SystemExit("MISSING: %r" % a[:70])
    s = s.replace(a, b)

io.open(P, "w", encoding="utf-8", newline="").write(s)
print("page refs aligned")