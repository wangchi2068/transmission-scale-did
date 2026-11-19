"""Round-3 revision: symbolic-number and typographic fixes flagged in the
second review (section 4.2, items 1-5)."""

import io
import re

def edit(path, pairs, required=True):
    s = io.open(path, encoding="utf-8").read()
    for a, b in pairs:
        if a not in s:
            if required:
                raise SystemExit("MISSING in %s: %r" % (path, a[:70]))
            continue
        s = s.replace(a, b)
    io.open(path, "w", encoding="utf-8", newline="").write(s)
    print("edited", path)

# item 1: double rounding of the adjacent-score SE (source 0.291479)
edit("sections/application.tex", [("from 0.311 to 0.292", "from 0.311 to 0.291")])

# item 2: the mean in Table 3's note is over e0-e5, not over all seven bins
edit("tables/table_7_application.tex", [
    ("across all seven bins ($-0.343$ to $-0.924$, mean $-0.677$,",
     "across all seven bins ($-0.343$ to $-0.924$, mean of $e_0$--$e_5$ $-0.677$,"),
])

# item 3: print the C.1 / C.2 numbers in the subsubsection headings
edit("sections/appendix.tex", [
    ("\\subsubsection{Validation of the root-selection rules on independent designs}",
     "\\subsubsection{C.1 Validation of the root-selection rules on independent designs}"),
    ("\\subsubsection{Population roots in the zero-contagion designs}",
     "\\subsubsection{C.2 Population roots in the zero-contagion designs}"),
])
print("done")