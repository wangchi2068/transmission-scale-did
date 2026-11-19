# -*- coding: utf-8 -*-
"""Locate manuscript anchors for the response letters: given a list of short
distinctive phrases, report the PDF page and printed page (PDF - 1)."""

import io
import re
import sys
import unicodedata

import pypdf

PDF = "main.pdf"
PHRASES = [
    ("discussion calibration sentence", "application-calibrated simulation in which the weekly model holds"),
    ("intro same-week arithmetic", "rules the weekly specification out"),
    ("application cluster sizes", "cluster size is unequal"),
    ("application robust-set condition", "strongly identified there"),
    ("application event-time weights", "2.109 on the adoption bin"),
    ("application data paragraph", "3091 contiguous counties"),
    ("application estimates", "policy-level total"),
    ("protocol departures", "M0-GI was added after interim boundary failures"),
    ("remark testable", "adjacent states"),
    ("table 3 note anchor", "mean of"),
]


def norm(s: str) -> str:
    s = unicodedata.normalize("NFKD", s)
    s = s.replace("\u2013", "-").replace("\u2014", "-").replace("\u2019", "'")
    return re.sub(r"\s+", "", s).lower()


reader = pypdf.PdfReader(PDF)
pages = [norm(p.extract_text() or "") for p in reader.pages]

for label, phrase in PHRASES:
    target = norm(phrase)
    hits = [i + 1 for i, t in enumerate(pages) if target in t]
    # tolerate interleaved whitespace-free extraction differences
    if not hits:
        # fall back: match on first 40 normalised chars
        short = target[:40]
        hits = [i + 1 for i, t in enumerate(pages) if short in t]
    print("%-34s %s" % (label, hits if hits else "NOT FOUND"))

print("\nsection starts:")
for name, pat in [("Intro", "1.Introduction"), ("Methods", "2.Methods"),
                  ("Simulation", "3.SimulationStudy"), ("Application", "4.Application"),
                  ("Discussion", "5.Discussion"), ("Backmatter", "Acknowledgements"),
                  ("References", "References")]:
    for i, t in enumerate(pages):
        if pat.lower().replace(" ", "") in t:
            print("  %-12s PDF %d  printed %d" % (name, i + 1, i))
            break