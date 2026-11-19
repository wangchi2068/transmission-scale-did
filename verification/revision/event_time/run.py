"""Analysis e (review Comment 4): what does the pooled M0 tau estimate on the M0-ET2 sample average over?

Sample: application panel without the 435 pre-window exposed units (make_panel(..., drop_presample=True)):
2656 counties, 37 state clusters.  Fits: pooled M0 (tau, kappa, phi) and M0-ET2 (tau_e0..tau_e6, kappa, phi).

1. Implicit weights.  Write the M0 estimating equations h(b, kappa, phi) with the linear index sum_e b_e B_e + kappa Dbar
   (M0 instruments: s D, s Dbar, phi instrument; week effects profiled).  M0 imposes b_e = tau for all e.
   Linearising h around a point theta*, the pooled solution satisfies
       tau_pooled ~ sum_e omega_e b_e  + r,   omega = [A^{-1} H_b]_tau-row,  A = [H_b 1, H_kappa, H_phi],
   where H = profiled Jacobian of the M0 moments w.r.t. (b, kappa, phi); sum_e omega_e = 1 by construction.
   r = -[A^{-1} h(theta_ET2)]_tau collects the part of the pooled moments not explained by the bins (it has mean
   zero when the event-time model is correct, because the M0 instruments are predetermined).
   omega is evaluated at the M0-ET2 estimate and at the pooled estimate (b_e = tau_pooled).
2. Contrast pooled tau - sum omega_e b_e (and pooled minus the paper's equal-weight e0..e5 mean, and minus the
   cell-share weighted mean) with a joint cluster sandwich: the two estimators' cluster influence terms
   (-Gamma^{-1} g_c, delta profiled into the full Jacobian) are stacked on the same 37 clusters; t_36 interval.
3. Over-identification on the ET2 sample: two-step Hansen J with the M0-OI extras for M0 and for M0-ET2.
4. phi-profile of the M0-ET2 phi equation over a grid and all roots (QDX with the bins as linear index).
No randomness (deterministic).
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
from qdx import QDX, oi_extras, tsum

HERE = Path(__file__).resolve().parent
GRID = [0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05] + [round(x, 2) for x in np.arange(0.1, 3.0001, 0.05)]


def influence(mod, beta):
    inf = mod.inference(beta)
    Gi = np.linalg.inv(inf["Gam"])
    return (-Gi @ inf["gg"])[mod.T - 1:], inf      # (p, G) parameter rows


def m0_moments_general(panel, bins, b, kappa, phi):
    """M0 instruments, linear index sum b_e B_e + kappa Dbar: profiled moments h (3,) and Jacobian (3, nb+2)."""
    X = QDX(panel, dict(list(bins.items()) + [("kappa", panel["Db"])]))
    # replace the parameter instruments by the M0 ones (s D, s Dbar, phi-instrument)
    s = X.s
    lam_bar = X.L + 0.5 * X.Nb + X.om
    X.Zpar = [s * panel["D"][:-1], s * panel["Db"][:-1], s * 0.5 * X.Nb[:-1] / lam_bar[:-1]]
    beta = np.r_[b, kappa, phi]
    c = X._core(beta)
    F = np.array([(z * c["v"]).sum() for z in X.Zpar]) / X.N
    a, bb, cm, H = X._blocks(c, X.Zpar)
    Jp = (H - np.einsum("mtr,trj->mj", cm, bb / a[:, :, None])) / X.N
    return F, Jp


def weights_at(panel, bins, b, kappa, phi):
    h, Jp = m0_moments_general(panel, bins, b, kappa, phi)
    nb = len(bins)
    Hb, Hk, Hp = Jp[:, :nb], Jp[:, nb], Jp[:, nb + 1]
    Amat = np.column_stack([Hb.sum(1), Hk, Hp])
    sol = np.linalg.solve(Amat, Hb)
    r = -np.linalg.solve(Amat, h)
    return dict(omega=sol[0].tolist(), omega_sum=float(sol[0].sum()), kappa_loading=sol[1].tolist(),
                phi_loading=sol[2].tolist(), remainder_tau=float(r[0]), h=h.tolist())


def main():
    t0 = time.time()
    d = A.load_data()
    w = A.gi_weekly(A.GI_SOURCES["ferretti"])["w"]
    panel, _ = A.make_panel(d, w, d["W_adj"], drop_presample=True)
    G = int(panel["G"])
    tc = qdgmm.tcrit(G)
    out = dict(sample=dict(n_units=int(panel["Y"].shape[1]), G=G, tcrit=tc))
    sp0 = qdgmm.spec_m0(panel)
    m0 = qdgmm.fit(panel, sp0)
    spe = qdgmm.spec_m0_et2(panel)
    me = qdgmm.fit(panel, spe)
    bins = qdgmm.event_bins(panel)
    bn = list(bins)
    out["pooled_M0"] = dict(est=m0["est"], se=m0["se"],
                            tau_plus_kappa=m0["est"]["tau"] + m0["est"]["kappa"], tpk_se=tsum(m0["Vb"], 0, 1))
    out["M0_ET2"] = dict(est=me["est"], se=me["se"], n_cells={k: float(bins[k][:-1].sum()) for k in bn})
    print("pooled", {k: round(v, 3) for k, v in m0["est"].items()}, "se", {k: round(v, 3) for k, v in m0["se"].items()})
    print("ET2", {k: round(v, 3) for k, v in me["est"].items()})
    b_e = np.array([me["est"][k] for k in bn])
    # 1. implicit weights
    W_et = weights_at(panel, bins, b_e, me["est"]["kappa"], me["est"]["phi"])
    W_po = weights_at(panel, bins, np.full(len(bn), m0["est"]["tau"]), m0["est"]["kappa"], m0["est"]["phi"])
    out["implicit_weights"] = dict(bins=bn, at_ET2_estimate=W_et, at_pooled_estimate=W_po)
    for lab, W in (("at ET2", W_et), ("at pooled", W_po)):
        print(lab, "omega", np.round(W["omega"], 3), "sum", round(W["omega_sum"], 4), "remainder", round(W["remainder_tau"], 3))
    # 2. contrasts with joint sandwich
    mod0, mode = qdgmm.QDModel(panel, sp0), qdgmm.QDModel(panel, spe)
    I0, _ = influence(mod0, np.array([m0["est"][k] for k in sp0.names]))
    Ie, _ = influence(mode, np.array([me["est"][k] for k in spe.names]))
    ie_bins = Ie[:len(bn)]
    cells = np.array([bins[k][:-1].sum() for k in bn])
    wsets = {"implicit_omega_at_ET2": np.array(W_et["omega"]), "implicit_omega_at_pooled": np.array(W_po["omega"]),
             "equal_e0_e5": np.r_[np.full(len(bn) - 1, 1.0 / (len(bn) - 1)), 0.0], "cell_share": cells / cells.sum()}
    out["contrasts"] = {}
    for lab, wv in wsets.items():
        avg = float(wv @ b_e)
        diff = m0["est"]["tau"] - avg
        infl = I0[0] - wv @ ie_bins
        se = float(np.sqrt(np.sum(infl ** 2) / G / G * G / (G - 1)))
        se_avg = float(np.sqrt(np.sum((wv @ ie_bins) ** 2) / G / G * G / (G - 1)))
        out["contrasts"][lab] = dict(weights=wv.tolist(), weighted_bins=avg, weighted_bins_se=se_avg,
                                     pooled_minus_weighted=diff, se=se, ci=[diff - tc * se, diff + tc * se],
                                     t=diff / se, p_t=float(2 * stats.t.sf(abs(diff / se), G - 1)))
        print(lab, "avg %.3f (%.3f) diff %.3f se %.3f p %.3f" % (avg, se_avg, diff, se, out["contrasts"][lab]["p_t"]))
    # 3. over-identification on the ET2 sample
    oi0 = qdgmm.fit_two_step(panel, qdgmm.spec_m0_oi(panel), m0)
    speo = qdgmm.spec_m0_et2(panel)
    ex, nm = oi_extras(panel)
    speo.extra_inst, speo.extra_names = ex, nm
    oie = qdgmm.fit_two_step(panel, speo, me)
    out["J"] = dict(M0=dict(J=oi0["J"], df=oi0.get("J_df"), p=oi0["J_p"]),
                    M0_ET2=dict(J=oie["J"], df=oie.get("J_df"), p=oie["J_p"]))
    print("J M0 %.2f p %.2g | J ET2 %.2f p %.2g" % (oi0["J"], oi0["J_p"], oie["J"], oie["J_p"]))
    # 4. phi profile of M0-ET2 and of pooled M0 on this sample
    out["phi_profile"] = {}
    for lab, lin in (("M0_ET2", dict(list(bins.items()) + [("kappa", panel["Db"])])),
                     ("pooled_M0", {"tau": panel["D"], "kappa": panel["Db"]})):
        X = QDX(panel, lin)
        cands, prof = X.global_roots(GRID)
        out["phi_profile"][lab] = dict(grid=[p[0] for p in prof], F_phi=[p[3] for p in prof], ok=[p[2] for p in prof],
                                       roots=[dict(kind=k, est=dict(zip(X.names, x.tolist()))) for k, x in cands])
        print(lab, "roots", [(k, round(float(x[-1]), 4)) for k, x in cands])
    out["seconds"] = time.time() - t0
    out["git"] = env.git_sha()
    (HERE / "results.json").write_text(json.dumps(env.clean(out), indent=1))


if __name__ == "__main__":
    main()
