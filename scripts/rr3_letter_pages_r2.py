# -*- coding: utf-8 -*-
"""The round-3 trims moved the methods and simulation sections one printed page
earlier. Shift the round-2 response letter's references accordingly. Sections 1,
2.2, 2.6, 4 and 5 did not move and are left alone."""

import io
import re

P = "D:/codeg_work/论文相关/peer-review-of-transmission-scale-did-f4dda8b8/response-to-reviewers.md"
s = io.open(P, encoding="utf-8").read()

# order matters: 3.2.1 before 3.2
rules = [
    (r"§3\.2\.1, p\. 14", "§3.2.1, p. 13"),
    (r"§3\.2\.1 \(p\. 14\)", "§3.2.1 (p. 13)"),
    (r"§3\.2, p\. 13", "§3.2, p. 12"),
    (r"§3\.2, p\. 12", None),          # already correct values, skip
    (r"§3\.1, p\. 12", "§3.1, p. 11"),
    (r"§2\.5, p\. 9", "§2.5, p. 8"),
]
counts = {}
for pat, rep in rules:
    if rep is None:
        continue
    s, n = re.subn(pat, rep, s)
    counts[pat] = n

io.open(P, "w", encoding="utf-8", newline="").write(s)
for k, v in counts.items():
    print("%-24s %d" % (k, v))