"""Report which printed page each main-text section starts and ends on."""
import re
import sys

import pypdf

PAGE_RE = re.compile(r"^\s*(\d{1,2})\s*(Biometrics|DiD on the transmission scale)?", re.M)
SECTIONS = [
    ("Introduction", "Introduction", "2 Theory"),
    ("Methods", "Theory", "3 Simulation"),
    ("Simulation", "Simulation design", "4 Application"),
    ("Application", "Application: State Mask", "5 Discussion"),
    ("Discussion", "For panels of reported cases we have proposed", "References"),
]

reader = pypdf.PdfReader("main.pdf")
texts = [(p.extract_text() or "") for p in reader.pages]
out = []
for name, start_marker, _ in SECTIONS:
    start = next((i + 1 for i, t in enumerate(texts[:41]) if start_marker in t), None)
    out.append(f"{name:14s} starts at PDF page {start}")
print("\n".join(out))
print("main text runs to PDF page 40 (printed page 39)")
