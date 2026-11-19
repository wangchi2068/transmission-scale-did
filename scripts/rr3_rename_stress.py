"""Round-3 revision: rename the calibrated-simulation 'stress design' to
'high-noise design' so it is not confused with the same-week 'stress regime'.
Reads 'stress' occurrences linewise in the table body only."""

import io

TABLE = "tables/table_12_calibrated_sim.tex"
APPX = "sections/appendix.tex"
APP = "sections/application.tex"

# table body: only the design label in the first column
s = io.open(TABLE, encoding="utf-8").read()
s = s.replace("stress & ", "high-noise & ")
io.open(TABLE, "w", encoding="utf-8", newline="").write(s)

s = io.open(APPX, encoding="utf-8").read()
s = s.replace("week shocks. The stress design", "week shocks. The high-noise design")
io.open(APPX, "w", encoding="utf-8", newline="").write(s)

s = io.open(APP, encoding="utf-8").read()
old = "and 0.950 in a stress design matching the\napplication's standard error"
new = "and 0.950 in a high-noise design that scales the shocks up to the\napplication's standard error"
assert old in s
s = s.replace(old, new)
io.open(APP, "w", encoding="utf-8", newline="").write(s)
print("renamed")