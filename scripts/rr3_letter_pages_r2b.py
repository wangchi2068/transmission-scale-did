# -*- coding: utf-8 -*-
"""Shift the round-2 response letter's Section 4 references one printed page
earlier. The round-3 trims pulled Methods and Simulation up by one page, so the
application now occupies printed pp. 14-17 instead of 15-17."""

import io
import re

P = "D:/codeg_work/论文相关/peer-review-of-transmission-scale-did-f4dda8b8/response-to-reviewers.md"
s = io.open(P, encoding="utf-8").read()

# only rewrite §4-references; leave §3.2.1 / §5 / SI pages untouched
before = s
mapping = [("§4, p. 17", "§4, p. 16"),
           ("§4 (p. 17)", "§4 (p. 16)"),
           ("§4, p. 16", "§4, p. 15"),
           ("§4 (p. 16)", "§4 (p. 15)"),
           ("§4, p. 15", "§4, p. 14"),
           ("§4 (p. 15)", "§4 (p. 14)")]

# two passes with placeholders so a chain does not cascade
s = s.replace("§4, p. 17", "@@A@@").replace("§4 (p. 17)", "@@B@@")
s = s.replace("§4, p. 16", "@@C@@").replace("§4 (p. 16)", "@@D@@")
s = s.replace("§4, p. 15", "@@E@@").replace("§4 (p. 15)", "@@F@@")
s = s.replace("@@A@@", "§4, p. 16").replace("@@B@@", "§4 (p. 16)")
s = s.replace("@@C@@", "§4, p. 15").replace("@@D@@", "§4 (p. 15)")
s = s.replace("@@E@@", "§4, p. 14").replace("@@F@@", "§4 (p. 14)")

# the one reference that already pointed at printed p.14 (§4 opening) stays
io.open(P, "w", encoding="utf-8", newline="").write(s)
print("changed" if s != before else "no change")
for m in re.finditer(r"§4[, ]*\(?p\. 1[0-9]\)?", s):
    pass
print(len(re.findall(r"§4[, ]*\(?p\. 1[0-9]\)?", s)), "§4 refs")