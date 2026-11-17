"""Small QA helpers for the manuscript sources (no arguments needed).

Checks: running-head length, abstract word count, presence of the data-contract path in the
application table note, and whether the README still claims that internal file names were removed
from table notes.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

main = (ROOT / "main.tex").read_text(encoding="utf-8")
m = re.search(r"\\title\[([^\]]*)\]", main)
head = m.group(1)
print("running head:", len(head), repr(head))

abstract = (ROOT / "sections/abstract.tex").read_text(encoding="utf-8")
body = abstract.split(r"\begin{abstract}")[1].split(r"\end{abstract}")[0]
print("abstract words:", len(re.findall(r"[A-Za-z0-9][A-Za-z0-9\-\.]*", body)))

table = (ROOT / "tables/table_7_application.tex").read_text(encoding="utf-8")
print("table note cites data contract:", "data_contract" in table)

readme = (ROOT / "README.md").read_text(encoding="utf-8")
print("README claims internal names removed:", "内部文件名已从表注中移除" in readme)
