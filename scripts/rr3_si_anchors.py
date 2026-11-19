# -*- coding: utf-8 -*-
"""Locate supporting-information anchors for the response letter."""

import re
import unicodedata

import pypdf


def norm(s):
    s = unicodedata.normalize("NFKD", s)
    s = s.replace("\u2013", "-").replace("\u2014", "-").replace("\u2019", "'")
    return re.sub(r"\s+", "", s).lower()


reader = pypdf.PdfReader("supporting_information.pdf")
pages = [norm(p.extract_text() or "") for p in reader.pages]

PH = [
    ("robust sets table", "identi-fication-robustconfidencesets"),
    ("robust sets text", "zero-transfertsandwichintervals"),
    ("C.1 heading", "C.1Validationoftheroot-selection"),
    ("A.12 window week", "thirdwindowweek"),
    ("F intro", "Allanalysesinthisappendixwereaddedinrevision"),
    ("event weights table", "implicitweightsoftheevent-tim"),
    ("calibrated sim table", "Application-calibratedsimulation"),
]

for label, pat in PH:
    t = norm(pat)
    hits = [i + 1 for i, p in enumerate(pages) if t in p]
    if not hits:
        hits = [i + 1 for i, p in enumerate(pages) if t[:30] in p]
    print("%-26s %s" % (label, hits or "NOT FOUND"))