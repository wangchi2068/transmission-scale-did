"""Pilot evidence for the ALPHA_SCALE deviation (see scenarios.py / experiment_design.md).
For alpha-component multipliers c in {1 (spec), 0.5, 0.25}, calibrate gbar crudely (median Y/scale = 20, 4 reps
per bisection step) and report quantiles of county-week counts per unit population scale (8 reps, N = 500).
Output: logs/pilot_alpha_tails.txt.  Not a paper result; design diagnostics only."""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dgp

geo = dgp.make_geography(500)

def stats_for(T, adopt, gbar, c, R=8):
    p = dgp.DGPParams(N=500, T=T, adopt=adopt, gbar=gbar, alpha_state_sd=0.2 * c, alpha_unit_sd=0.2 * c, alpha_pi=0.3 * c)
    Ys, sh = [], []
    for r in range(R):
        s = dgp.simulate(p, geo, [2, r, 0], [2, r, 1])
        Ys.append((s["Yall"][2:] / s["scale"]).ravel()); sh.append(s["Dall"][2:].mean())
    Y = np.concatenate(Ys)
    return np.median(Y), np.quantile(Y, [0.9, 0.99, 0.999]), Y.max(), np.mean(sh)

out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "logs", "pilot_alpha_tails.txt")
with open(out, "w") as fh:
    fh.write("c\tT\tadopt\tgbar\tmedian\tq90\tq99\tq999\tmax\tshare_treated\n")
    for c in (1.0, 0.5, 0.25):
        for T in (12, 26):
            for adopt in ("feedback", "random"):
                lo, hi = -0.6, 0.6
                for _ in range(12):
                    mid = (lo + hi) / 2
                    if stats_for(T, adopt, mid, c, R=4)[0] > 20: hi = mid
                    else: lo = mid
                m, q, mx, sh = stats_for(T, adopt, mid, c)
                fh.write("%.2f\t%d\t%s\t%+.3f\t%.1f\t%.0f\t%.0f\t%.0f\t%.3g\t%.2f\n" % (c, T, adopt, mid, m, *q, mx, sh))
                fh.flush()
print("wrote", out)
