"""Pilot calibration of the DGP constants that simulation_spec.md leaves to the pilot.

1. gbar (constant in gamma_t) per (T, adoption): median county-week count in the estimation window = 20
   (spec: must lie in [5, 50]); evaluated on N = 500, (phi, kappa) = (0.4, -0.1), tau0 = log 0.8, k = 10,
   R_CAL replications with fixed seeds (common random numbers across bisection steps).
2. depletion_pop for X4: P_i = depletion_pop * scale_i such that the median (over replications) of the
   maximum county cumulative attack rate (pre-sample + window infections / P_i) is 0.20 (spec: "max ~20%").
Writes calibration.json next to src/. Also records count quantiles for every calibrated cell.
Run:  python src/calibrate.py   (about 1-3 min on 4 cores)
"""
from __future__ import annotations

import json
import os
import sys
import time
from dataclasses import replace

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dgp  # noqa: E402
import scenarios as sc  # noqa: E402

R_CAL = 24
TARGET_MEDIAN = 20.0
TARGET_ATTACK = 0.20


def _params(T, adopt, gbar):
    return sc._base(500, T, adopt, 0.4, -0.1, sc.TAU0, 10.0, gbar)


def median_count(p, geo, R=R_CAL, seed=777):
    Ys = []
    for r in range(R):
        s = dgp.simulate(p, geo, [seed, r, 0], [seed, r, 1])
        Ys.append(s["Yall"][dgp.S_LAGS:].ravel())
    Y = np.concatenate(Ys)
    return float(np.median(Y)), Y


def calibrate_gbar(T, adopt, geo):
    lo, hi = -0.8, 0.8
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        m, _ = median_count(_params(T, adopt, mid), geo)
        if m > TARGET_MEDIAN:
            hi = mid
        else:
            lo = mid
    g = round(0.5 * (lo + hi), 4)
    m, Y = median_count(_params(T, adopt, g), geo, seed=778)       # check on fresh seeds
    return g, dict(median=m, q90=float(np.quantile(Y, 0.9)), q99=float(np.quantile(Y, 0.99)),
                   q999=float(np.quantile(Y, 0.999)), max=float(Y.max()), frac_zero=float((Y == 0).mean()))


def max_attack(p, geo, R=R_CAL, seed=779):
    vals = []
    for r in range(R):
        s = dgp.simulate(p, geo, [seed, r, 0], [seed, r, 1])
        vals.append(float(np.max(s["Iall"].sum(0) / s["pop"])))
    return float(np.median(vals))


def calibrate_depletion(gbar, geo2000):
    base = replace(sc._base(2000, 26, "feedback", 0.4, -0.1, sc.TAU0, 10.0, gbar))
    lo, hi = np.log(1e2), np.log(1e7)
    for _ in range(20):
        mid = 0.5 * (lo + hi)
        a = max_attack(replace(base, depletion_pop=float(np.exp(mid))), geo2000, R=24)
        if a > TARGET_ATTACK:
            lo = mid
        else:
            hi = mid
    c = float(np.round(np.exp(0.5 * (lo + hi)), 1))
    return c, max_attack(replace(base, depletion_pop=c), geo2000, R=24, seed=780)


def main(out_path=sc.CALIB_PATH):
    t0 = time.time()
    geo = dgp.make_geography(500)
    res = {"gbar": {}, "gbar_check": {}, "target_median": TARGET_MEDIAN, "R_cal": R_CAL,
           "alpha_scale": sc.ALPHA_SCALE, "design": "N=500, phi=0.4, kappa=-0.1, tau0=log0.8, k=10"}
    for T in (12, 26):
        for adopt in ("feedback", "random"):
            g, chk = calibrate_gbar(T, adopt, geo)
            res["gbar"]["T%d_%s" % (T, adopt)] = g
            res["gbar_check"]["T%d_%s" % (T, adopt)] = chk
            print("gbar", T, adopt, g, chk, flush=True)
    c, a = calibrate_depletion(res["gbar"]["T26_feedback"], dgp.make_geography(2000))
    res["depletion_pop"] = c
    res["depletion_check_median_max_attack"] = a
    res["target_max_attack"] = TARGET_ATTACK
    res["runtime_sec"] = round(time.time() - t0, 1)
    res["created_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(out_path, "w") as fh:
        json.dump(res, fh, indent=1)
    print("wrote", out_path, res["runtime_sec"], "s")


if __name__ == "__main__":
    main()
