"""Scan the LaTeX sources for control characters that replaced a backslash.

A shell heredoc can turn '\\t' into a tab and '\\tau' into a tab plus 'au'. This reports every
remaining tab or other stray control character with its line, so they can be fixed before a build.
"""
from pathlib import Path

FILES = [Path("paper_zh/main.tex"), Path("main.tex")] + sorted(Path("sections").glob("*.tex"))
BAD = {9: "TAB", 7: "BEL", 8: "BS", 11: "VT", 12: "FF"}

found = 0
for f in FILES:
    for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
        for code, name in BAD.items():
            if chr(code) in line:
                found += 1
                print(f"{f}:{i}: {name} -> {line.strip()[:120]}")
print("total", found)
