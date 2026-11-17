"""QA checks for Biometrics compliance and cross-document PDF integrity."""
import re
from pathlib import Path
import pypdf

ROOT = Path(__file__).resolve().parents[1]

main_tex = (ROOT / "main.tex").read_text(encoding="utf-8")
m = re.search(r"\\title\[([^\]]*)\]", main_tex)
head = m.group(1) if m else ""
print(f"running head ({len(head)} chars, <=50): {head!r}")
assert len(head) <= 50, "Running head exceeds 50 characters"

abstract = (ROOT / "sections/abstract.tex").read_text(encoding="utf-8")
body = abstract.split(r"\begin{abstract}")[1].split(r"\end{abstract}")[0]
n_words = len(re.findall(r"[A-Za-z0-9][A-Za-z0-9\-\.]*", body))
print(f"abstract words ({n_words}, <=225)")
assert n_words <= 225, "Abstract exceeds 225 words"

for pdf_name in ("main.pdf", "supporting_information.pdf", "paper_zh/main.pdf"):
    p = ROOT / pdf_name
    if not p.exists():
        print(f"{pdf_name}: MISSING")
        continue
    reader = pypdf.PdfReader(str(p))
    q_pages = [i + 1 for i, pg in enumerate(reader.pages) if "??" in (pg.extract_text() or "")]
    print(f"{pdf_name}: {len(reader.pages)} pages, unresolved '??' pages={q_pages}")
    assert not q_pages, f"Unresolved references in {pdf_name}: {q_pages}"
