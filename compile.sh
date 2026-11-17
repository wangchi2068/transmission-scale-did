#!/bin/bash
set -euo pipefail

# Compile Biometrics manuscript (main.pdf) and Supporting Information (supporting_information.pdf)
# with bidirectional xr cross-references, plus optional source/PDF QA checks.

if command -v python3 >/dev/null 2>&1; then
    PY=python3
elif command -v python >/dev/null 2>&1; then
    PY=python
else
    PY=""
fi

if [ -n "$PY" ]; then
    "$PY" scripts/scan_control_chars.py
fi

echo "[1/5] Initial pass on main.tex..."
xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null

echo "[2/5] Compiling supporting_information.tex (BibTeX + 2 passes)..."
xelatex -interaction=nonstopmode -halt-on-error supporting_information.tex >/dev/null
bibtex supporting_information >/dev/null
xelatex -interaction=nonstopmode -halt-on-error supporting_information.tex >/dev/null
xelatex -interaction=nonstopmode -halt-on-error supporting_information.tex >/dev/null

echo "[3/5] Settling main.tex (BibTeX + 2 passes)..."
xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
bibtex main >/dev/null
xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null

echo "[4/5] Final xr sync pass on supporting_information.tex..."
xelatex -interaction=nonstopmode -halt-on-error supporting_information.tex >/dev/null

if [ -f paper_zh/main.tex ]; then
    echo "[5/5] Compiling Chinese companion paper_zh/main.tex..."
    (cd paper_zh && xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null && bibtex main >/dev/null && xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null && xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null)
fi

if [ -n "$PY" ]; then
    "$PY" scripts/qa_checks.py
    "$PY" scripts/section_pages.py
fi

if [ "${1:-}" = "--clean" ]; then
    rm -f main.aux main.bbl main.blg main.fff main.log main.out main.ttt \
          supporting_information.aux supporting_information.bbl supporting_information.blg \
          supporting_information.log supporting_information.out \
          paper_zh/main.aux paper_zh/main.bbl paper_zh/main.blg paper_zh/main.log paper_zh/main.out
fi

echo "Done: main.pdf and supporting_information.pdf are up to date."
