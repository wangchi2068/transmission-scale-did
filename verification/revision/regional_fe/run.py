"""Analysis a (review Comment 1, requests 1-2): refit the application with Census-division x week effects and with
Census-region-specific phi; report tau, kappa, phi, tau+kappa, the over-identification statistic and the root-
selection statistics. Reference distributions under correct specification come from ../calibrated_sim (same
statistics computed on every simulated panel); see ../calibrated_sim/compare.py.

Specifications (all with the paper's instruments, phi_bar = 0.5, Ferretti weekly generation interval):
  N0  national week effects, common phi (= M0; checks the score-form statistics against the paper's J, Q_G, Q_GI)
  R4  Census-region (4) x week effects, common phi
  A1  Census-division (9) x week effects, common phi
  A2  national week effects, region-specific phi (4)
  A3  division x week effects + region-specific phi
J_LM: score-form test of the two M0-OI extra instruments (lagged import share, lagged growth), chi2(2); Hansen's
two-step J is not computable once the week moments outnumber the 49 clusters (Omega singular).
J_LM_inv: same with the M0-GI phi-instruments at reference values 0.05 and 2.0.
Q*_G / Q*_GI: J_LM_oi / J_LM_inv evaluated at each candidate root of a phi-grid search (common-phi specs only);
the chosen root minimises the statistic, as in M0-G / M0-GI.
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import env  # noqa
import numpy as np
from scipy import stats
import application as A
import qdgmm
from qdx import QDX, oi_extras, inv_extras, groups_for_units, tsum, DIV_NAMES, REGION_NAMES

HERE = Path(__file__).resolve().parent


def main():
    t0 = time.time()
    d = A.load_data()
    w = A.gi_weekly(A.GI_SOURCES["ferretti"])["w"]
    panel, _ = A.make_panel(d, w, d["W_adj"])
    div, reg = groups_for_units(d["state"])
    G = int(panel["G"])
    tc = qdgmm.tcrit(G)
    exo, _ = oi_extras(panel)
    exi, _ = inv_extras(panel)
    m0 = qdgmm.fit(panel, qdgmm.spec_m0(panel))
    b0 = np.array([m0["est"][n] for n in ("tau", "kappa", "phi")])
    specs = dict(N0=(None, None), R4=(reg, None), A1=(div, None), A2=(None, reg), A3=(div, reg))
    out = dict(G=G, t_crit=tc, n_units=int(panel["Y"].shape[1]), weeks=int(panel["Y"].shape[0]),
               division_states={DIV_NAMES[k + 1]: int(len(np.unique(d["state"][div == k]))) for k in range(9)},
               region_states={REGION_NAMES[k + 1]: int(len(np.unique(d["state"][reg == k]))) for k in range(4)},
               specs={})
    for lab, (wg, pg) in specs.items():
        X = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, wgroup=wg, phigroup=pg, extras=exo + exi,
                extra_names=["lag_import_share", "lag_growth", "phiinst_0.05", "phiinst_2"])
        P = X.P
        starts = [np.r_[b0[:2], [b0[2]] * P], np.r_[0.0, 0.0, [0.5] * P], np.r_[-0.2, -0.2, [0.2] * P], np.r_[0.0, 0.0, [1.0] * P]]
        x, ok, res = X.fit(starts=starts)
        rec = dict(week_groups=X.R, phi_groups=P, n_moments_week=X.nW, n_params=X.np_, ok=ok, resid=res)
        if ok:
            inf = X.inference(x)
            se = np.sqrt(np.clip(np.diag(inf["Vb"]), 0, None))
            m = X.moments(x, extra=True)
            jo, ji = X.j_lm(x, sel=[0, 1], m=m), X.j_lm(x, sel=[2, 3], m=m)
            tot, stot = float(x[0] + x[1]), tsum(inf["Vb"], 0, 1)
            rec.update(est=dict(zip(X.names, x)), se=dict(zip(X.names, se)),
                       ci={n: [float(v - tc * s), float(v + tc * s)] for n, v, s in zip(X.names, x, se)},
                       tau_plus_kappa=dict(est=tot, se=stot, ci=[tot - tc * stot, tot + tc * stot]),
                       J_LM_oi=jo["J_LM"], J_LM_oi_p=jo["p"], J_LM_inv=ji["J_LM"], J_LM_inv_p=ji["p"],
                       cond_Gamma=float(np.linalg.cond(inf["Gam"])))
            if P > 1:
                # Wald test of equal phi across regions (df P-1), cluster sandwich, F(P-1, G-1) version too
                Vp = inf["Vb"][2:, 2:]
                Rm = np.zeros((P - 1, P)); Rm[:, 0] = 1; Rm[np.arange(P - 1), np.arange(1, P)] = -1
                dv = Rm @ x[2:]
                Wd = float(dv @ np.linalg.solve(Rm @ Vp @ Rm.T, dv))
                rec["equal_phi_wald"] = dict(chi2=Wd, df=P - 1, p_chi2=float(stats.chi2.sf(Wd, P - 1)),
                                             F=Wd / (P - 1), p_F=float(stats.f.sf(Wd / (P - 1), P - 1, G - 1)))
            if X.nW + X.np_ + 2 < G:   # Hansen-type full-Omega Q only where Omega can be non-singular
                Xo = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, wgroup=wg, phigroup=pg, extras=exo)
                Xi = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, wgroup=wg, phigroup=pg, extras=exi)
                rec["Q_full_oi"], rec["Q_full_inv"] = Xo.q_full(x), Xi.q_full(x)
        if P == 1:
            cands, prof = X.global_roots(qdgmm.PHI_GRID)
            cl_ = []
            for kind, beta in cands:
                mm = X.moments(beta, extra=True)
                cl_.append(dict(kind=kind, phi=float(beta[-1]), tau=float(beta[0]), kappa=float(beta[1]),
                                tau_plus_kappa=float(beta[0] + beta[1]),
                                Qstar_G=X.j_lm(beta, sel=[0, 1], m=mm)["J_LM"], Qstar_GI=X.j_lm(beta, sel=[2, 3], m=mm)["J_LM"]))
            rec["candidates"] = cl_
            rec["phi_profile"] = [(float(p), float(F)) for p, _, okp, F in prof if okp]
            if cl_:
                rec["chosen_G"] = min(cl_, key=lambda c: c["Qstar_G"])
                rec["chosen_GI"] = min(cl_, key=lambda c: c["Qstar_GI"])
        out["specs"][lab] = rec
        print(lab, "ok", ok, {k: round(v, 3) for k, v in rec.get("est", {}).items()}, "tpk", rec.get("tau_plus_kappa"),
              "J_LM", rec.get("J_LM_oi"), rec.get("J_LM_inv"), flush=True)
    # paper-code statistics for N0 (validated equality in ../validate_qdx)
    oi = qdgmm.fit_two_step(panel, qdgmm.spec_m0_oi(panel), m0)
    gG = qdgmm.fit_global(panel); gGI = qdgmm.fit_global(panel, selection="inv")
    out["N0_paper_code"] = dict(J=oi["J"], J_p=oi["J_p"], Q_G=gG["Q"], Q_GI=gGI["Q"])
    out["seconds"] = time.time() - t0
    out["git"] = env.git_sha()
    (HERE / "results.json").write_text(json.dumps(env.clean(out), indent=1))


if __name__ == "__main__":
    main()
