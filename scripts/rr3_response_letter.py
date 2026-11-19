# -*- coding: utf-8 -*-
"""Round-3 revision: update the response letter for the second round."""

import io

P = "D:/codeg_work/论文相关/peer-review-of-transmission-scale-did-f4dda8b8/response-to-reviewers.md"
s = io.open(P, encoding="utf-8").read()

old = "not a causal effect (Abstract, p. 1; §1, p. 2; §4, p. 16; §5, p. 18)."
new = "not a causal effect (unnumbered title page; §1, p. 2; §4, p. 16; §5, p. 18)."
assert old in s
s = s.replace(old, new)

old2 = "- **Location.** Abstract, p. 1; §1, p. 2; §4, p. 16; §5, p. 18."
new2 = "- **Location.** Unnumbered title page; §1, p. 2; §4, p. 16; §5, p. 18."
assert old2 in s
s = s.replace(old2, new2)

io.open(P, "w", encoding="utf-8", newline="").write(s)
print("updated")