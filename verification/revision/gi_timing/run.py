"""Analysis d (review Comment 2, request 2): application estimates under generation-interval discretisations that
move the same-week mass into the lagged terms, linked to the paper's gi_w1_* sensitivity runs.

Discretisations of the published daily distribution (p_s = P(lag = s weeks), infector case uniform within its week):
  base      paper: drop p_0 and renormalise (w = p_1:2 / sum)                      [= M0 main fit]
  fold      add the same-week mass to lag 1:  w = (p_0 + p_1, p_2) / sum
  shift_h   add an offset of h days to the generation time before discretising (end-of-week dating of the
            infector's cases / a reporting-delay difference between infector and infectee); h = 1.75, 3.5, 7.
            Same-week mass falls to p_0(h); h = 7 moves all of p_0 into lag 1 and p_1 into lag 2.
for both published sources (Ferretti generation time; Nishiura serial interval).
Each fit reports M0 (tau, kappa, phi, tau+kappa with t_48 interval), Hansen J (M0-OI), M0-G and M0-GI choices
and the profiled S-set for phi (../robust_sets) with the Bonferroni tau+kappa set for its interior component.
The paper's gi_w1_* runs (w = (w1, 1 - w1), w1 = 0.6..1.0) are read from experiment/analysis/results.json.
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "robust_sets"))
import env  # noqa
import numpy as np
import application as A
import qdgmm
import robust
from qdx import tsum

HERE = Path(__file__).resolve().parent


def p_lags(src, h=0.0, S=2, n_u=7000):
    F = A._cdf(src)
    u = (np.arange(n_u) + 0.5) / n_u * 7.0
    return np.array([np.mean(F(np.clip(7 * (s + 1) - u - h, 0, None)) - F(np.clip(7 * s - u - h, 0, None)))
                     for s in range(S + 1)])


def fit_all(d, w, tc):
    panel, _ = A.make_panel(d, w, d["W_adj"])
    sp = qdgmm.spec_m0(panel)
    m0 = qdgmm.fit(panel, sp)
    out = dict(w=list(map(float, w)), M0_ok=bool(m0["ok"]))
    if not m0["ok"]:
        return out
    tot, se = m0["est"]["tau"] + m0["est"]["kappa"], tsum(m0["Vb"], 0, 1)
    out.update(est=m0["est"], se=m0["se"], tau_plus_kappa=dict(est=tot, se=se, ci=[tot - tc * se, tot + tc * se]))
    oi = qdgmm.fit_two_step(panel, qdgmm.spec_m0_oi(panel), m0)
    out.update(J=oi["J"], J_p=oi["J_p"])
    for lab, sel in (("M0_G", "oi"), ("M0_GI", "inv")):
        g = qdgmm.fit_global(panel, selection=sel)
        out[lab] = dict(chosen=g.get("chosen"), Q=g.get("Q"), roots=g.get("roots"),
                        phi=g["est"]["phi"], tau_plus_kappa=g["est"]["tau"] + g["est"]["kappa"])
    rows = robust.s_stats(panel)
    st = robust.sets(rows, int(panel["G"]))
    out["robust"] = dict(phi_set_95_S1_F=st["S1_F_95"], phi_set_975_S1_F=st.get("phi_set_975_S1_F"),
                         S3_F_95=st["S3_F_95"], tpk_by_component=st.get("tpk_by_component"))
    return out


def main():
    t0 = time.time()
    d = A.load_data()
    G = len(np.unique(d["state"]))
    tc = qdgmm.tcrit(G)
    res = dict(variants={}, paper_gi_w1={})
    for sname, src in A.GI_SOURCES.items():
        for lab, h, mode in (("base", 0.0, "drop"), ("fold", 0.0, "fold"), ("shift_1.75", 1.75, "drop"),
                             ("shift_3.5", 3.5, "drop"), ("shift_7", 7.0, "drop")):
            p = p_lags(src, h)
            if mode == "drop":
                w = p[1:] / p[1:].sum()
            else:
                w = np.array([p[0] + p[1], p[2]]) / p.sum()
            key = "%s_%s" % (sname, lab)
            r = fit_all(d, w, tc)
            r.update(p_by_week_lag=p.tolist(), same_week_mass=float(p[0]), offset_days=h, mode=mode)
            res["variants"][key] = r
            print(key, "w", np.round(w, 3), "p0", round(p[0], 3), "tpk", round(r["tau_plus_kappa"]["est"], 3),
                  "ci", np.round(r["tau_plus_kappa"]["ci"], 3), "phi", round(r["est"]["phi"], 3), "J", round(r["J"], 2),
                  "G:", r["M0_G"]["chosen"], round(r["M0_G"]["tau_plus_kappa"], 3), "GI:", r["M0_GI"]["chosen"],
                  round(r["M0_GI"]["tau_plus_kappa"], 3), flush=True)
            (HERE / "results.json").write_text(json.dumps(env.clean(res), indent=1))
    app = json.loads((env.RUN / "experiment" / "analysis" / "results.json").read_text())["application"]["analyses"]
    for k, v in app.items():
        if k.startswith("gi_") or k == "main":
            m = v["M0"]
            res["paper_gi_w1"][k] = dict(tau=m["est"]["tau"], kappa=m["est"]["kappa"], phi=m["est"]["phi"],
                                         tau_plus_kappa=m["est"]["tau"] + m["est"]["kappa"],
                                         M0_G=dict((kk, v["M0_G"].get(kk)) for kk in ("chosen", "est", "Q")) if "M0_G" in v else None)
    res["seconds"] = time.time() - t0
    res["git"] = env.git_sha()
    (HERE / "results.json").write_text(json.dumps(env.clean(res), indent=1))


if __name__ == "__main__":
    main()
