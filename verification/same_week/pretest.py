"""Pre-test: calibrate the within-week re-entry g of same_week_dgp.py on the core scenario.

g controls how fast the current week's infections re-enter the same week (rate c = w0 g / 7 per day).  g = 1
is the plain within-week process (c = w0 / 7, within-week renewal A = w0); with the adoption feedback, a
lower peak incidence on the days before the lagged feedback arrives reduces the number of treated units, so
the calibration matches two moments of the weekly model on weeks 3-12 of the window (same seeds):
  * mean outcome ratio engine / weekly model,
  * mean treated-share ratio engine / weekly model.
Both are near-linear in g; g is solved by linear interpolation from the pre-test grid and written to
calib.json (read by same_week_dgp).  Run:  python -I pretest.py     (REPS env var, default 40)
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, "..", "..", "experiment", "code"))
sys.path.insert(0, HERE)
sys.path.insert(0, CODE)

import numpy as np  # noqa: E402

import dgp  # noqa: E402
import same_week_dgp as sw  # noqa: E402
import scenarios as sc  # noqa: E402

GRID = (1.0, 1.5, 2.0, 3.0)


def moments(g, reps):
    sw.CAL["g"] = g
    p = sc.build_scenarios()[1].params
    geo = dgp.make_geography(p.N)
    rat, tre = [], []
    for r in range(reps):
        seed = sc.seed_for(1, r)
        A = dgp.simulate(p, geo, [seed, 0], [seed, 1])
        B = sw.simulate_same_week(p, geo, [seed, 0], [seed, 1], 0.3)
        rat.append(B["Yall"][4:].mean() / A["Yall"][4:].mean())
        tre.append(B["Dall"][4:].mean() / A["Dall"][4:].mean())
    return (float(np.mean(rat)), float(np.std(rat, ddof=1) / np.sqrt(reps)),
            float(np.mean(tre)), float(np.std(tre, ddof=1) / np.sqrt(reps)))


if __name__ == "__main__":
    reps = int(os.environ.get("REPS", "40"))
    rows = []
    for g in GRID:
        y, yse, t, tse = moments(g, reps)
        rows.append(dict(g=g, ratio=y, ratio_mcse=yse, treated_ratio=t, treated_ratio_mcse=tse))
        print(rows[-1], flush=True)
    gs = np.array([r["g"] for r in rows])
    # mean of the two normalised gaps (ratio - 1 and treated_ratio - 1), both should be 0
    y = np.array([0.5 * (r["ratio"] - 1 + r["treated_ratio"] - 1) for r in rows])
    # linear fit, then root
    b, a = np.polyfit(gs, y, 1)
    g_hat = float(-a / b) if b != 0 else float(gs[np.argmin(np.abs(y))])
    sel = rows[int(np.argmin(np.abs(gs - 1))) : int(np.argmin(np.abs(gs - 1))) + 2]
    if len(sel) >= 2:
        y2 = [0.5 * (s["ratio"] - 1 + s["treated_ratio"] - 1) for s in sel]
        g_hat = sel[0]["g"] + (sel[1]["g"] - sel[0]["g"]) * (-y2[0]) / (y2[1] - y2[0])
    print("calibrated g = %.3f" % g_hat)
    cal = dict(g=round(float(g_hat), 3), pretest=rows, n_reps=reps,
               note="g from the panel check: mean of (mean ratio - 1) and (treated-share ratio - 1) = 0")
    with open(os.path.join(HERE, "calib.json"), "w") as fh:
        json.dump(cal, fh, indent=1)
