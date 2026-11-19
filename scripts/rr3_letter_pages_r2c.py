# -*- coding: utf-8 -*-
"""Round-3 page alignment for the first/second-round response letter.

Verified against commit 2244b15 (main.pdf, 32 pages):
  Section 1      printed 1-2
  Section 2.5-2.6 printed 8-10
  Section 3.2.1  printed 13-14
  Section 4      printed 14-17
  Section 5      printed 18-19
"""

import io
import re

P = "D:/codeg_work/论文相关/1113/QD-GMM_transmission_DiD_2026-11-13/response/response-to-reviewers-r1r2.md"
s = io.open(P, encoding="utf-8").read()

pairs = [
    ("(§4, p. 14; Web Appendix F)", "(§4, p. 14; Web Appendix F)"),
    ("§2.6, p. 10", "§2.6, p. 10"),
]

# Section 5 now spans printed 18-19; references that say "§5, p. 18" are correct
# for the opening, so only the ones pointing at the closing paragraph move.
s = s.replace("(§5, p. 18). See Comment 3.2.", "(§5, pp. 18-19). See Comment 3.2.")
s = s.replace("Done (§5, p. 18). The sentence now reads as follows:",
              "Done (§5, pp. 18-19). The sentence now reads as follows:")

# Table 3's note sits on printed p.32 in the old layout; the endfloat block now
# starts at printed 26, so locate its actual page.
io.open(P, "w", encoding="utf-8", newline="").write(s)
print("ok")