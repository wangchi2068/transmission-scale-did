"""Analysis f (review Comment 7, request 2): tipping point for treatment-linked reporting change.

Reporting model.  Reported cases Y_it = c_i c_t (1 - rho)^{D_it} X_it: a mask order changes the ascertainment
fraction of the adopting county by the factor (1 - rho) for as long as it is in force (rho > 0: reporting falls,
e.g. less testing; rho < 0: reporting rises).  c_i c_t is absorbed by alpha_i, gamma_t (Proposition 1(a));
the treatment-linked part violates Proposition 1's national homogeneity.  Under this model the log reproduction
number of reported cases is shifted by log(1 - rho) * (D_it - weighted D over the generation window), i.e. the
distortion is concentrated in switch weeks and is bounded by |log(1 - rho)| (Proposition 1(b)); through the
neighbour pressure it also enters kappa.
Tipping-point analysis.  For each rho on a grid, undo the reporting change, Y*_it = Y_it / (1 - rho)^{D_it}
(all weeks including the two pre-sample weeks, so L and the neighbour pressure N are rebuilt from Y*), refit M0
(and M0-G / M0-OI J at selected points), and record tau, kappa, tau + kappa with the t_48 sandwich interval.
rho*_0  : rho at which the corrected tau + kappa estimate is 0;
rho*_hi : rho < 0 at which the upper 95% limit reaches 0 (the interval would exclude 0);
rho*_lo : rho > 0 at which the lower 95% limit reaches 0.
Deterministic (no randomness).
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import env  # noqa
import numpy as np
from scipy import optimize
import application as A
import qdgmm
from qdx import tsum

HERE = Path(__file__).resolve().parent


def fit_rho(d, w, rho, tc, extra=False):
    dd = dict(d)
    dd["Y"] = d["Y"] / (1.0 - rho) ** d["D"]
    panel, _ = A.make_panel(dd, w, d["W_adj"])
    m0 = qdgmm.fit(panel, qdgmm.spec_m0(panel))
    if not m0["ok"]:
        return dict(rho=rho, ok=False)
    tot, se = m0["est"]["tau"] + m0["est"]["kappa"], tsum(m0["Vb"], 0, 1)
    r = dict(rho=rho, ok=True, est=m0["est"], se=m0["se"], tpk=tot, tpk_se=se, lo=tot - tc * se, hi=tot + tc * se,
             bound_log_shift=abs(np.log(1 - rho)))
    if extra:
        oi = qdgmm.fit_two_step(panel, qdgmm.spec_m0_oi(panel), m0)
        g = qdgmm.fit_global(panel, selection="oi")
        r.update(J=oi["J"], J_p=oi["J_p"], M0_G=dict(chosen=g.get("chosen"), phi=g["est"]["phi"],
                                                       tpk=g["est"]["tau"] + g["est"]["kappa"]))
    return r


def main():
    t0 = time.time()
    d = A.load_data()
    w = A.gi_weekly(A.GI_SOURCES["ferretti"])["w"]
    G = len(d["states"])
    tc = qdgmm.tcrit(G)
    grid = [round(x, 3) for x in np.arange(-0.6, 0.6001, 0.05)]
    rows = []
    for rho in grid:
        r = fit_rho(d, w, rho, tc, extra=abs(rho * 10 - round(rho * 10)) < 1e-9)
        rows.append(r)
        print("rho %+.2f tau %.3f kappa %.3f phi %.3f tpk %.3f [%.3f, %.3f]%s" % (
            rho, r["est"]["tau"], r["est"]["kappa"], r["est"]["phi"], r["tpk"], r["lo"], r["hi"],
            (" J %.2f G %s %.3f" % (r["J"], r["M0_G"]["chosen"], r["M0_G"]["tpk"])) if "J" in r else ""), flush=True)
    out = dict(grid=rows, tcrit=tc, G=G)
    f = {k: (lambda rho, k=k: fit_rho(d, w, rho, tc)[k]) for k in ("tpk", "hi", "lo")}
    xs = np.array(grid)

    def root(key, sign_region):
        vals = np.array([r[key] for r in rows])
        for a, b, va, vb in zip(xs[:-1], xs[1:], vals[:-1], vals[1:]):
            if sign_region(a, b) and np.sign(va) != np.sign(vb):
                return float(optimize.brentq(f[key], a, b, xtol=1e-5))
        return None
    out["rho_star_tpk_zero"] = root("tpk", lambda a, b: True)
    out["rho_star_upper_limit_zero"] = root("hi", lambda a, b: b <= 0.0 + 1e-12)
    out["rho_star_lower_limit_zero"] = root("lo", lambda a, b: a >= 0.0 - 1e-12)
    for k in ("rho_star_tpk_zero", "rho_star_upper_limit_zero", "rho_star_lower_limit_zero"):
        if out[k] is not None:
            out[k + "_fit"] = fit_rho(d, w, out[k], tc, extra=True)
            out[k + "_log_shift"] = float(np.log(1 - out[k]))
        print(k, out[k])
    out["seconds"] = time.time() - t0
    out["git"] = env.git_sha()
    (HERE / "results.json").write_text(json.dumps(env.clean(out), indent=1))


if __name__ == "__main__":
    main()
