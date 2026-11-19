"""Identification-robust (S / Anderson-Rubin type) confidence sets for phi and tau+kappa (review Comment 3 / W2).

For each phi on a grid, (tau, kappa) and the week effects are solved from their own moments with phi held fixed;
the remaining phi-moment (and optionally the two M0-OI extra moments) is tested with the score-form statistic
    psi_g = g_{t,g} - Gamma_{t,n} Gamma_{n,n}^{-1} g_{n,g}   (n = week, tau, kappa moments / parameters)
    S(phi) = G psibar' (sum_g psi_g psi_g' / G)^{-1} psibar          [psibar = gbar_t since gbar_n = 0]
S(phi) does not use the slope of the phi-moment, so its level does not depend on phi being strongly identified
(Stock & Wright 2000 S-test, profiled over the strongly identified tau, kappa). Critical values: chi2(df) and the
small-sample version S (G - df) / (df (G - 1)) ~ F(df, G - df).
Sets:
  phi:       {phi : S1(phi) <= crit}    (phi-moment only, df = 1)
  phi (AR3): {phi : S3(phi) <= crit}    (phi-moment + 2 M0-OI extra moments, df = 3; also tests the over-identifying restrictions)
  tau+kappa: Bonferroni projection: union over phi in the 97.5% phi-set of the 97.5% t_{G-1} interval for tau+kappa
             with phi fixed (cluster sandwich of the phi-fixed model); overall level >= 95%.
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import env  # noqa
import numpy as np
from scipy import stats
import qdgmm
from qdx import QDX, oi_extras, tsum

HERE = Path(__file__).resolve().parent
GRID = np.r_[0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04, np.arange(0.05, 3.0001, 0.01)]


def s_stats(panel, grid=GRID, phi_bar=0.5):
    exo, _ = oi_extras(panel, phi_bar)
    Xf = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, phi_bar=phi_bar, extras=exo)
    G = Xf.G
    rows = []
    x = np.zeros(2)
    for ph in grid:
        b, ok, Fp = Xf.solve_given_phi(np.array([ph]), x)
        if not ok:
            rows.append(dict(phi=float(ph), ok=False))
            continue
        x = b
        beta = np.r_[b, ph]
        m = Xf.moments(beta, extra=True)
        nW = Xf.nW
        # rows: [week (nW), tau, kappa, phi, extra1, extra2]; cols: [week (nW), tau, kappa, phi]
        n_idx = np.r_[np.arange(nW), nW, nW + 1]
        J = m["J"] / G
        gg = m["gg"]
        Gnn = J[np.ix_(n_idx, n_idx)]
        out = dict(phi=float(ph), ok=True, tau=float(b[0]), kappa=float(b[1]), tpk=float(b[0] + b[1]))
        for lab, t_idx in (("S1", [nW + 2]), ("S3", [nW + 2, nW + 3, nW + 4])):
            Gtn = J[np.ix_(t_idx, n_idx)]
            psi = gg[t_idx] - Gtn @ np.linalg.solve(Gnn, gg[n_idx])
            pbar = psi.sum(1) / G
            S = psi @ psi.T / G
            stat = float(G * pbar @ np.linalg.solve(S, pbar))
            df = len(t_idx)
            out[lab] = stat
            out[lab + "_p_chi2"] = float(stats.chi2.sf(stat, df))
            out[lab + "_p_F"] = float(stats.f.sf(stat * (G - df) / (df * (G - 1)), df, G - df))
        # phi-fixed inference for tau+kappa (week effects, tau, kappa estimated; phi known)
        Xc = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, phi_bar=phi_bar, phi_free=False, phi_fixed=ph)
        inf = Xc.inference(b)
        if inf is not None:
            out["tpk_se_phi_fixed"] = tsum(inf["Vb"], 0, 1)
        rows.append(out)
    return rows


def sets(rows, G, level=0.95):
    def intervals(mask, phis):
        ivs, cur = [], None
        for m, p in zip(mask, phis):
            if m and cur is None:
                cur = [p, p]
            elif m:
                cur[1] = p
            elif cur is not None:
                ivs.append(cur); cur = None
        if cur is not None:
            ivs.append(cur)
        return ivs
    ok = [r for r in rows if r["ok"]]
    phis = [r["phi"] for r in ok]
    res = {}
    for lab in ("S1", "S3"):
        for kind in ("chi2", "F"):
            mask = [r[lab + "_p_" + kind] > 1 - level for r in ok]
            res["%s_%s_%d" % (lab, kind, round(level * 100))] = intervals(mask, phis)
    # Bonferroni projection for tau+kappa: 97.5% phi-set (S1, F) x 97.5% conditional t interval
    tc = float(stats.t.ppf(1 - 0.025 / 2, G - 1))
    inset = [r for r in ok if r["S1_p_F"] > 0.025 and "tpk_se_phi_fixed" in r]
    if inset:
        lo = min(r["tpk"] - tc * r["tpk_se_phi_fixed"] for r in inset)
        hi = max(r["tpk"] + tc * r["tpk_se_phi_fixed"] for r in inset)
        res["tpk_bonferroni_95"] = [lo, hi]
        res["tpk_range_over_phi_set_975"] = [min(r["tpk"] for r in inset), max(r["tpk"] for r in inset)]
        res["phi_set_975_S1_F"] = intervals([r["S1_p_F"] > 0.025 for r in ok], phis)
        comps = []
        for a, b in res["phi_set_975_S1_F"]:
            sub = [r for r in inset if a <= r["phi"] <= b]
            comps.append(dict(phi=[a, b], tpk_bonferroni=[min(r["tpk"] - tc * r["tpk_se_phi_fixed"] for r in sub),
                                                         max(r["tpk"] + tc * r["tpk_se_phi_fixed"] for r in sub)],
                              tpk_range=[min(r["tpk"] for r in sub), max(r["tpk"] for r in sub)]))
        res["tpk_by_component"] = comps
    return res


def main():
    import application as A
    t0 = time.time()
    d = A.load_data()
    w = A.gi_weekly(A.GI_SOURCES["ferretti"])["w"]
    out = dict(grid=dict(n=len(GRID), min=float(GRID.min()), max=float(GRID.max()), step=0.01))
    for phi_bar in (0.5, 0.25, 1.0):
        panel, _ = A.make_panel(d, w, d["W_adj"])
        rows = s_stats(panel, phi_bar=phi_bar)
        G = int(panel["G"])
        out["phi_bar_%g" % phi_bar] = dict(sets=sets(rows, G), rows=rows,
                                            S1_at_0=[r for r in rows if r["phi"] == 0.0][0])
        print(phi_bar, json.dumps(env.clean(out["phi_bar_%g" % phi_bar]["sets"])), flush=True)
    out["seconds"] = time.time() - t0
    out["git"] = env.git_sha()
    (HERE / "results.json").write_text(json.dumps(env.clean(out), indent=1))


if __name__ == "__main__":
    main()
