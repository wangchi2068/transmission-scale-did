"""Report section start pages and Biometrics 25-page compliance from main.pdf."""
import re
import pypdf

SECTIONS = [
    ("Introduction", r"1\.\s*Introduction"),
    ("Methods", r"2\.\s*Methods"),
    ("Simulation Study", r"3\.\s*Simulation Study"),
    ("Section 4 (Implementation)", r"4\.\s*(?:Application:|Illustrative implementation)"),
    ("Discussion", r"5\.\s*Discussion"),
    ("Back Matter", r"Acknowledgements"),
    ("References", r"^References\b"),
    ("Endfloat Figures/Tables", r"^Figure 1\."),
]

reader = pypdf.PdfReader("main.pdf")
texts = [(p.extract_text() or "") for p in reader.pages]

for name, pat in SECTIONS:
    rx = re.compile(pat, re.M)
    pdf_page = next((i + 1 for i, t in enumerate(texts) if rx.search(t)), None)
    if pdf_page is None:
        print(f"{name:24s}: not found")
    else:
        printed = pdf_page - 1 if pdf_page > 1 else 0
        print(f"{name:24s}: PDF p.{pdf_page:2d} (printed p.{printed:2d})")

fig_idx = next((i for i, t in enumerate(texts) if re.search(r"^Figure 1\.", t, re.M)), len(texts))
counted_pages = fig_idx - 1  # excludes PDF p.1 (cover/abstract) and endfloat floats
print(f"Biometrics-counted pages (Sections 1-5 + Back Matter + References): {counted_pages} / 25")
print(f"Total PDF pages: {len(texts)}")
