"""Formal empirical application: German county-level weekly Varicella (2015-2019)
and staggered 16-state summer school holidays.

Reads:  experiment/data/processed/panel.npz (built by experiment/code/data_prep.py)
Writes: experiment/analysis/germany_application_results.json
"""
from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import scipy.sparse as sps
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiment" / "code"))
sys.path.insert(0, str(ROOT / "verification" / "revision" / "common"))
sys.path.insert(0, str(ROOT / "verification" / "application"))

import comparators as cmp
import dgp
import psi as psi_mod
import qdgmm
import wild_bootstrap as WB
from qdx import QDX, inv_extras, oi_extras, tsum

PANEL_PATH = ROOT / "experiment" / "data" / "processed" / "panel.npz"
SUMMARY_PATH = ROOT / "experiment" / "data" / "processed" / "panel_summary.json"
OUT_PATH = ROOT / "experiment" / "analysis" / "germany_application_results.json"

PHI_GRID_FINE = np.r_[0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04, np.arange(0.05, 3.0001, 0.01)]
PHI_GRID_ROOTS = np.r_[0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]


def gi_weekly(mean: float = 14.0, sd: float = 2.5, S: int = 3) -> dict:
    """Discretise a continuous lognormal generation interval (days) onto weekly lags 1..S."""
    sigma2 = np.log(1.0 + (sd / mean) ** 2)
    sigma = float(np.sqrt(sigma2))
    mu = float(np.log(mean) - 0.5 * sigma2)
    from scipy import integrate
    pdf = lambda u: stats.lognorm.pdf(u, s=sigma, scale=np.exp(mu))
    w0, _ = integrate.quad(lambda u: (1.0 - u / 7.0) * pdf(u), 0.0, 7.0)
    raw = []
    for s in range(1, S + 1):
        val, _ = integrate.quad(lambda u: (1.0 - abs(u - 7.0 * s) / 7.0) * pdf(u), 7.0 * (s - 1), 7.0 * (s + 1))
        raw.append(float(val))
    raw_arr = np.array(raw, float)
    w = raw_arr / raw_arr.sum()
    tail = max(0.0, 1.0 - w0 - float(raw_arr.sum()))
    return dict(dist="lognormal", mean=mean, sd=sd, S=S, w0=float(w0), tail_beyond_S=float(tail),
                raw=raw_arr.tolist(), w=w.tolist())


def load_data() -> dict:
    z = np.load(PANEL_PATH, allow_pickle=False)
    W = sps.csr_matrix((z["W_data"], z["W_indices"], z["W_indptr"]), shape=tuple(z["W_shape"]))
    return dict(
        Y=z["Y"],
        D_summer=z["D_summer"],
        D_summer_3d=z["D_summer_3d"],
        D_summer_5d=z["D_summer_5d"],
        Dshare_summer=z["Dshare_summer"],
        D_all=z["D_all"],
        Dshare_all=z["Dshare_all"],
        D_easter=z["D_easter"],
        D_autumn=z["D_autumn"],
        D_winter=z["D_winter"],
        D_pentecost=z["D_pentecost"],
        pop=z["pop"],
        state=z["state"],
        cl=z["cl"],
        states=z["states"],
        region=z["region"],
        nuts2=z["nuts2"],
        nuts2_codes=z["nuts2_codes"],
        ags=z["ags"],
        weeks=z["weeks"],
        edges=z["edges"],
        W=W,
    )


def lag_matrix(M: np.ndarray, lag: int) -> np.ndarray:
    if lag <= 0:
        return M.copy()
    return np.vstack([np.zeros((lag, M.shape[1]), dtype=float), M[:-lag]])


def make_panel(d: dict, w: np.ndarray, D_mat: np.ndarray, lag: int = 1, omega_mult: float = 1.0,
               S: int = 3, t_slice: slice | None = None, Y_override: np.ndarray | None = None) -> dict:
    Y = d["Y"] if Y_override is None else Y_override
    Dl = lag_matrix(D_mat, lag)
    if t_slice is not None:
        Y = Y[t_slice]
        Dl = Dl[t_slice]
    omega = omega_mult * d["pop"] / 1e5
    panel = dgp.build_panel(Y, Dl, d["W"], omega, d["cl"], w=w, S=S)
    panel = qdgmm.attach_presample(panel, Dl[:S])
    panel["Yall"] = Y
    return panel


def score_s3_at_beta(X: QDX, beta: np.ndarray, extra_sel: list[int]) -> dict:
    """Rank-safe 3-moment score statistic S3 over [phi-moment, extra_sel[0], extra_sel[1]]
    orthogonalised against (week, tau, kappa). Valid for any T regardless of G > 3."""
    m = X.moments(beta, extra=True)
    if m is None:
        return dict(S3=float("nan"), p_chi2=float("nan"), p_F=float("nan"))
    G = X.G
    nW = X.nW
    n_idx = np.r_[np.arange(nW), nW, nW + 1]
    t_idx = [nW + 2] + [nW + 3 + k for k in extra_sel]
    J = m["J"] / G
    gg = m["gg"]
    Gnn = J[np.ix_(n_idx, n_idx)]
    Gtn = J[np.ix_(t_idx, n_idx)]
    psi = gg[t_idx] - Gtn @ np.linalg.solve(Gnn, gg[n_idx])
    pbar = psi.sum(1) / G
    Sm = psi @ psi.T / G
    stat = float(G * pbar @ np.linalg.solve(Sm, pbar))
    df = len(t_idx)
    return dict(
        S3=stat,
        p_chi2=float(stats.chi2.sf(stat, df)),
        p_F=float(stats.f.sf(stat * (G - df) / (df * (G - 1)), df, G - df)),
    )


def fit_global_score(panel: dict, phi_bar: float = 0.5, grid=PHI_GRID_ROOTS) -> dict:
    """Global root search over phi in [0, max(grid)] using rank-safe score-form S3 criterion
    for both OI extras (M0-G) and invariance extras (M0-GI)."""
    exo, _ = oi_extras(panel, phi_bar=phi_bar)
    exi, _ = inv_extras(panel)
    X = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, phi_bar=phi_bar, extras=exo + exi)
    cands, prof = X.global_roots(grid)
    interior = [c for c in cands if c[0] == "interior"]
    cand_records = []
    for kind, beta in cands:
        s3_oi = score_s3_at_beta(X, beta, [0, 1])
        s3_inv = score_s3_at_beta(X, beta, [2, 3])
        if kind == "interior":
            inf = X.inference(beta)
            Vb = inf["Vb"]
            se = np.sqrt(np.diag(Vb))
            tpk_se = tsum(Vb, 0, 1)
        else:
            Xc = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, phi_bar=phi_bar, phi_free=False, phi_fixed=0.0)
            inf = Xc.inference(beta[:2])
            Vs = inf["Vb"]
            se = np.r_[np.sqrt(np.diag(Vs)), np.nan]
            tpk_se = tsum(Vs, 0, 1)
        cand_records.append(dict(
            kind=kind,
            est=dict(tau=float(beta[0]), kappa=float(beta[1]), phi=float(beta[2]), tpk=float(beta[0] + beta[1])),
            se=dict(tau=float(se[0]), kappa=float(se[1]), phi=float(se[2]), tpk=float(tpk_se)),
            S3_oi=s3_oi,
            S3_inv=s3_inv,
        ))
    idx_g = int(np.argmin([r["S3_oi"]["S3"] for r in cand_records]))
    idx_gi = int(np.argmin([r["S3_inv"]["S3"] for r in cand_records]))
    return dict(
        n_interior_roots=len(interior),
        interior_phis=[float(c[1][2]) for c in interior],
        profile=[dict(phi=float(p), ok=bool(ok), F_phi=float(fp)) for p, _, ok, fp in prof],
        candidates=cand_records,
        M0_G=cand_records[idx_g],
        M0_GI=cand_records[idx_gi],
    )


def fit_m0_full(panel: dict, phi_bar: float = 0.5, wgroup=None) -> dict:
    """Fit M0 via QDX and compute sandwich inference, J_LM (oi, inv, and 1-df breakdown)."""
    exo, _ = oi_extras(panel, phi_bar=phi_bar)
    exi, _ = inv_extras(panel)
    X = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, wgroup=wgroup, phi_bar=phi_bar, extras=exo + exi)
    x, ok, res = X.fit(starts=[[-0.3, 0.0, 0.3], [-0.2, -0.1, 0.5], [0.0, 0.0, 0.5]])
    if not ok:
        return dict(ok=False)
    inf = X.inference(x)
    Vb = inf["Vb"]
    se = np.sqrt(np.diag(Vb))
    G = X.G
    tc = float(stats.t.ppf(0.975, G - 1))
    tpk = float(x[0] + x[1])
    tpk_se = tsum(Vb, 0, 1)

    mX = X.moments(x, extra=True)
    j_oi = X.j_lm(x, sel=[0, 1], m=mX)
    j_inv = X.j_lm(x, sel=[2, 3], m=mX)
    j_imp = X.j_lm(x, sel=[0], m=mX)
    j_grw = X.j_lm(x, sel=[1], m=mX)

    def f_p(J, df):
        return float(stats.f.sf(J * (G - df) / (df * (G - 1)), df, G - df))

    return dict(
        ok=True,
        G=G,
        N=int(panel["Y"].shape[1]),
        T=int(panel["Y"].shape[0]),
        tcrit=tc,
        est=dict(tau=float(x[0]), kappa=float(x[1]), phi=float(x[2]), tpk=tpk),
        se=dict(tau=float(se[0]), kappa=float(se[1]), phi=float(se[2]), tpk=tpk_se),
        ci_t=dict(
            tau=[float(x[0] - tc * se[0]), float(x[0] + tc * se[0])],
            kappa=[float(x[1] - tc * se[1]), float(x[1] + tc * se[1])],
            phi=[float(x[2] - tc * se[2]), float(x[2] + tc * se[2])],
            tpk=[float(tpk - tc * tpk_se), float(tpk + tc * tpk_se)],
        ),
        J_LM_oi=dict(stat=j_oi["J_LM"], df=2, p_chi2=j_oi["p"], p_F=f_p(j_oi["J_LM"], 2)),
        J_LM_inv=dict(stat=j_inv["J_LM"], df=2, p_chi2=j_inv["p"], p_F=f_p(j_inv["J_LM"], 2)),
        J_LM_oi_components=dict(
            lag_import_share=dict(stat=j_imp["J_LM"], p_chi2=j_imp["p"], p_F=f_p(j_imp["J_LM"], 1)),
            lag_growth=dict(stat=j_grw["J_LM"], p_chi2=j_grw["p"], p_F=f_p(j_grw["J_LM"], 1)),
        ),
        _beta=x,
        _Vb=Vb,
        _inf=inf,
    )


def strip_private(rec: dict) -> dict:
    return {k: v for k, v in rec.items() if not k.startswith("_")}


def jackknife_m0(d: dict, w: np.ndarray, lag: int = 1, phi_bar: float = 0.5) -> dict:
    """Delete-one-state jackknife (re-restricting W to the 15 kept states and row-normalising)."""
    G = len(d["states"])
    ests = []
    oks = []
    by_state = {}
    A_full = (d["W"] > 0).astype(float).tocsr()
    for g in range(G):
        keep = d["cl"] != g
        Ak = A_full[keep][:, keep]
        deg = np.asarray(Ak.sum(1)).ravel()
        # Drop any county left isolated when its enclosing state is removed (e.g. Berlin without Brandenburg)
        sub_keep = deg > 0
        if not sub_keep.all():
            idx_k = np.where(keep)[0][sub_keep]
            keep_mask = np.zeros(len(d["cl"]), dtype=bool)
            keep_mask[idx_k] = True
            Ak = A_full[keep_mask][:, keep_mask]
            deg = np.asarray(Ak.sum(1)).ravel()
        else:
            keep_mask = keep
        Wk = sps.diags(1.0 / deg) @ Ak
        Yk = d["Y"][:, keep_mask]
        Dk = lag_matrix(d["D_summer"][:, keep_mask], lag)
        omk = d["pop"][keep_mask] / 1e5
        clk_raw = d["cl"][keep_mask]
        _, clk = np.unique(clk_raw, return_inverse=True)
        pk = dgp.build_panel(Yk, Dk, Wk.tocsr(), omk, clk, w=w, S=3)
        sp = qdgmm.spec_m0(pk, phi_bar=phi_bar)
        m0 = qdgmm.fit(pk, sp, starts=[[-0.3, 0.0, 0.3], [0.0, 0.0, 0.5]])
        if not m0.get("ok"):
            raise RuntimeError(f"Jackknife fold state={int(d['states'][g])} failed to converge")
        oks.append(bool(m0["ok"]))
        b = np.array([m0["est"]["tau"], m0["est"]["kappa"], m0["est"]["phi"], m0["est"]["tau"] + m0["est"]["kappa"]])
        ests.append(b)
        by_state[int(d["states"][g])] = dict(tau=float(b[0]), kappa=float(b[1]), phi=float(b[2]), tpk=float(b[3]))
    E = np.array(ests)
    ebar = E.mean(0)
    se_jk = np.sqrt((G - 1) / G * ((E - ebar) ** 2).sum(0))
    return dict(
        ok_folds=int(sum(oks)),
        se=dict(tau=float(se_jk[0]), kappa=float(se_jk[1]), phi=float(se_jk[2]), tpk=float(se_jk[3])),
        tpk_range=[float(E[:, 3].min()), float(E[:, 3].max())],
        phi_range=[float(E[:, 2].min()), float(E[:, 2].max())],
        by_dropped_state=by_state,
    )


def stock_wright_sets(panel: dict, phi_bar: float = 0.5, grid=PHI_GRID_FINE) -> dict:
    """Profiled Stock-Wright S1, S3(oi), and S3(inv) confidence sets + Bonferroni 95% projection for tau+kappa."""
    exo, _ = oi_extras(panel, phi_bar=phi_bar)
    exi, _ = inv_extras(panel)
    Xf = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, phi_bar=phi_bar, extras=exo + exi)
    G = Xf.G
    nW = Xf.nW
    n_idx = np.r_[np.arange(nW), nW, nW + 1]
    rows = []
    x = np.array([-0.3, 0.0])
    for ph in grid:
        b, ok, Fp = Xf.solve_given_phi(np.array([ph]), x)
        if not ok:
            rows.append(dict(phi=float(ph), ok=False))
            continue
        x = b
        beta = np.r_[b, ph]
        m = Xf.moments(beta, extra=True)
        J = m["J"] / G
        gg = m["gg"]
        Gnn = J[np.ix_(n_idx, n_idx)]
        rec = dict(phi=float(ph), ok=True, tau=float(b[0]), kappa=float(b[1]), tpk=float(b[0] + b[1]), F_phi=float(Fp[0]))
        for lab, t_idx in (("S1", [nW + 2]), ("S3_oi", [nW + 2, nW + 3, nW + 4]), ("S3_inv", [nW + 2, nW + 5, nW + 6])):
            Gtn = J[np.ix_(t_idx, n_idx)]
            psi = gg[t_idx] - Gtn @ np.linalg.solve(Gnn, gg[n_idx])
            pbar = psi.sum(1) / G
            Sm = psi @ psi.T / G
            stat = float(G * pbar @ np.linalg.solve(Sm, pbar))
            df = len(t_idx)
            rec[lab] = stat
            rec[lab + "_p_chi2"] = float(stats.chi2.sf(stat, df))
            rec[lab + "_p_F"] = float(stats.f.sf(stat * (G - df) / (df * (G - 1)), df, G - df))
        Xc = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, phi_bar=phi_bar, phi_free=False, phi_fixed=ph)
        inf = Xc.inference(b)
        if inf is not None:
            rec["tpk_se_phi_fixed"] = tsum(inf["Vb"], 0, 1)
        rows.append(rec)

    def intervals(mask, phis):
        ivs, cur = [], None
        for m_, p_ in zip(mask, phis):
            if m_ and cur is None:
                cur = [p_, p_]
            elif m_:
                cur[1] = p_
            elif cur is not None:
                ivs.append(cur)
                cur = None
        if cur is not None:
            ivs.append(cur)
        return ivs

    ok_rows = [r for r in rows if r["ok"]]
    phis = [r["phi"] for r in ok_rows]
    sets_out = {}
    for lab in ("S1", "S3_oi", "S3_inv"):
        for kind in ("chi2", "F"):
            sets_out[f"{lab}_{kind}_95"] = intervals([r[f"{lab}_p_{kind}"] > 0.05 for r in ok_rows], phis)
    # Bonferroni 95% projection for tau+kappa using 97.5% S1(F) x 97.5% conditional t_{G-1}
    tc975 = float(stats.t.ppf(1.0 - 0.025 / 2.0, G - 1))
    inset = [r for r in ok_rows if r["S1_p_F"] > 0.025 and "tpk_se_phi_fixed" in r]
    if inset:
        sets_out["phi_set_975_S1_F"] = intervals([r["S1_p_F"] > 0.025 for r in ok_rows], phis)
        sets_out["tpk_bonferroni_95"] = [
            float(min(r["tpk"] - tc975 * r["tpk_se_phi_fixed"] for r in inset)),
            float(max(r["tpk"] + tc975 * r["tpk_se_phi_fixed"] for r in inset)),
        ]
        sets_out["tpk_range_over_975_phi_set"] = [
            float(min(r["tpk"] for r in inset)),
            float(max(r["tpk"] for r in inset)),
        ]
    return dict(sets=sets_out, S1_at_0=ok_rows[0], curve=[
        dict(phi=r["phi"], F_phi=r["F_phi"], S1=r["S1"], S1_p_F=r["S1_p_F"],
             S3_oi=r["S3_oi"], S3_oi_p_F=r["S3_oi_p_F"], S3_inv=r["S3_inv"], S3_inv_p_F=r["S3_inv_p_F"],
             tpk=r["tpk"]) for r in ok_rows
    ])


def event_time_analysis(d: dict, w: np.ndarray, include_leads: bool = True, Y_override: np.ndarray | None = None) -> dict:
    """Event-time bins around Summer holiday onset (e = -2, -1, 0, 1, ..., 7), controlling for Dbar,
    reporting both direct coefficients e_k and total-policy-scale effects tau_k + c_k * kappa."""
    Y = d["Y"] if Y_override is None else Y_override
    D0 = d["D_summer"]
    T_all, N = D0.shape
    ks = list(range(-2, 8)) if include_leads else list(range(0, 8))
    bins = {f"e{k}": np.zeros((T_all, N), dtype=float) for k in ks}
    for i in range(N):
        onsets = [t for t in range(1, T_all) if D0[t, i] == 1.0 and D0[t - 1, i] == 0.0]
        for t0 in onsets:
            for k in ks:
                if 0 <= t0 + k < T_all:
                    bins[f"e{k}"][t0 + k, i] = 1.0
    panel = dgp.build_panel(Y, D0, d["W"], d["pop"] / 1e5, d["cl"], w=w, S=3)
    lin = {k: v[3:] for k, v in bins.items()}
    lin["kappa"] = panel["Db"]
    X = QDX(panel, lin, phi_bar=0.5)
    x0 = np.r_[np.zeros(len(lin)), 0.32]
    x, ok, _ = X.fit(starts=[x0, np.r_[np.zeros(len(lin)), 0.5]])
    inf = X.inference(x)
    Vb = inf["Vb"]
    se = np.sqrt(np.diag(Vb))
    G = X.G
    tc = float(stats.t.ppf(0.975, G - 1))
    k_idx = X.names.index("kappa")
    coefs = {}
    for idx, nm in enumerate(X.names):
        coefs[nm] = dict(est=float(x[idx]), se=float(se[idx]),
                         ci=[float(x[idx] - tc * se[idx]), float(x[idx] + tc * se[idx])])
    total_scale = {}
    for k in ks:
        nm = f"e{k}"
        idx = X.names.index(nm)
        bm = lin[nm] > 0.5
        ck = float(panel["Db"][bm].mean()) if bm.any() else 0.0
        tot = float(x[idx] + ck * x[k_idx])
        var_tot = float(Vb[idx, idx] + (ck ** 2) * Vb[k_idx, k_idx] + 2.0 * ck * Vb[idx, k_idx])
        se_tot = float(np.sqrt(max(var_tot, 1e-30)))
        tstat = tot / se_tot
        total_scale[nm] = dict(
            k=k, ck=ck, est=tot, se=se_tot,
            ci=[float(tot - tc * se_tot), float(tot + tc * se_tot)],
            p_t=float(2.0 * stats.t.sf(abs(tstat), G - 1)),
        )
    out = dict(ok=ok, include_leads=include_leads, coefs=coefs, total_scale=total_scale)
    if include_leads:
        b_pre = x[:2]
        V_pre = Vb[:2, :2]
        W_pre = float(b_pre @ np.linalg.solve(V_pre, b_pre))
        out["pretrend_wald"] = dict(
            stat=W_pre, df=2,
            p_chi2=float(stats.chi2.sf(W_pre, 2)),
            p_F=float(stats.f.sf(W_pre * (G - 2) / (2 * (G - 1)), 2, G - 2)),
        )
    post_idx = [X.names.index(f"e{k}") for k in (1, 2, 3, 4, 5)]
    b_post = x[post_idx]
    V_post = Vb[np.ix_(post_idx, post_idx)]
    W_post = float(b_post @ np.linalg.solve(V_post, b_post))
    out["post_wald_e1_e5"] = dict(
        stat=W_post, df=5,
        p_chi2=float(stats.chi2.sf(W_post, 5)),
        p_F=float(stats.f.sf(W_post * (G - 5) / (5 * (G - 1)), 5, G - 5)),
    )
    # Joint Wald test on total policy scale for post-incubation weeks 1..5 and reopening weeks 6..7
    def wald_total(k_list):
        C = np.zeros((len(k_list), len(X.names)))
        for row_i, kk in enumerate(k_list):
            nm = f"e{kk}"
            C[row_i, X.names.index(nm)] = 1.0
            C[row_i, k_idx] = total_scale[nm]["ck"]
        theta_sub = C @ x[:-1] if C.shape[1] > len(x) else C @ x
        V_sub = C @ Vb @ C.T
        Wt = float(theta_sub @ np.linalg.solve(V_sub, theta_sub))
        df_w = len(k_list)
        return dict(stat=Wt, df=df_w, p_chi2=float(stats.chi2.sf(Wt, df_w)),
                    p_F=float(stats.f.sf(Wt * (G - df_w) / (df_w * (G - 1)), df_w, G - df_w)))

    out["post_total_wald_e1_e5"] = wald_total([1, 2, 3, 4, 5])
    out["reopen_total_wald_e6_e7"] = wald_total([6, 7])
    return out


def remark4_adjacent_dependence(d: dict, panel: dict, beta: np.ndarray) -> dict:
    sp = qdgmm.spec_m0(panel)
    mod = qdgmm.QDModel(panel, sp)
    inf = mod.inference(beta)
    G = mod.G
    gg = inf["gg"]
    edges = d["edges"]
    cl = d["cl"]
    Adj = np.zeros((G, G), dtype=bool)
    a, b = cl[edges[:, 0]], cl[edges[:, 1]]
    Adj[a, b] = True
    Adj[b, a] = True
    np.fill_diagonal(Adj, False)
    Om = gg @ gg.T / G
    cross = sum(np.outer(gg[:, i], gg[:, j]) for i, j in zip(*np.nonzero(Adj))) / G
    Tm1 = mod.T - 1
    names = sp.names
    Gi = np.linalg.inv(inf["Gam"])

    def se_from(Omega):
        Vb = (Gi @ Omega @ Gi.T / G * G / (G - 1))[Tm1:, Tm1:]
        res = {n: float(np.sqrt(Vb[k, k])) if Vb[k, k] > 0 else float("nan") for k, n in enumerate(names)}
        v = Vb[0, 0] + Vb[1, 1] + 2.0 * Vb[0, 1]
        res["tpk"] = float(np.sqrt(v)) if v > 0 else float("nan")
        eigvals = np.linalg.eigvalsh(Vb)
        res["min_eig_Vb"] = float(eigvals.min())
        res["non_psd"] = bool(eigvals.min() <= 0 or any(Vb[k, k] <= 0 for k in range(len(names))))
        return res

    corr = {}
    for k, n in enumerate(names):
        r = gg[Tm1 + k]
        num = sum(r[i] * r[j] for i, j in zip(*np.nonzero(Adj))) / Adj.sum()
        corr[n] = float(num / (r @ r / G))
    se_adj = se_from(Om + cross)
    return dict(
        G=G,
        adjacent_state_pairs=int(Adj.sum() // 2),
        mean_state_neighbours=float(Adj.sum(1).mean()),
        se_cluster=se_from(Om),
        se_with_adjacent=se_adj,
        non_psd=bool(se_adj["non_psd"]),
        adj_pair_score_corr=corr,
    )


def wcr_score_p_tpk0(X: QDX, beta_hat: np.ndarray, tpk0: float, B: int = 999, seed: int = 20261119,
                     weight: str = "webb6") -> float:
    """Kline & Santos (2012) Restricted Score Wild Cluster Bootstrap p-value for H0: tau + kappa = tpk0."""
    G = X.G
    nW = X.nW

    def eq_r(u):
        tau, phi = u[0], u[1]
        F, J = X.profiled_eq(np.array([tau, tpk0 - tau, phi]))
        F_r = np.array([F[0] - F[1], F[2]])
        J_r = np.array([
            [J[0, 0] - J[0, 1] - (J[1, 0] - J[1, 1]), J[0, 2] - J[1, 2]],
            [J[2, 0] - J[2, 1], J[2, 2]],
        ])
        return F_r, J_r

    sol = optimize.root(eq_r, np.array([beta_hat[0], beta_hat[2]]), jac=True, method="hybr")
    if not sol.success:
        return float("nan")
    b_r = np.array([sol.x[0], tpk0 - sol.x[0], sol.x[1]])
    mm = X.moments(b_r, extra=False)
    if mm is None:
        return float("nan")
    J = mm["J"] / G
    gg = mm["gg"]
    Q_m = np.eye(nW + 3)
    Q_m[nW] = 0.0; Q_m[nW, nW] = 1.0; Q_m[nW, nW + 1] = -1.0
    Q_m[nW + 1] = 0.0; Q_m[nW + 1, nW + 2] = 1.0
    Q_m[nW + 2] = 0.0; Q_m[nW + 2, nW] = 1.0; Q_m[nW + 2, nW + 1] = 1.0
    M_c = np.eye(nW + 3)
    M_c[nW:nW + 3, nW:nW + 3] = np.array([[0.5, 0.0, 0.5], [-0.5, 0.0, 0.5], [0.0, 1.0, 0.0]])
    Jt = Q_m @ J @ M_c
    gt = Q_m @ gg
    n_idx = np.arange(nW + 2)
    t_idx = nW + 2
    Gnn = Jt[np.ix_(n_idx, n_idx)]
    Gtn = Jt[t_idx, n_idx]
    psi = gt[t_idx] - Gtn @ np.linalg.solve(Gnn, gt[n_idx])
    se_r = float(np.sqrt(np.sum((psi - psi.mean()) ** 2) / (G * (G - 1))))
    t_obs = abs(float(psi.mean()) / se_r)
    rng_b = np.random.default_rng(seed)
    p1 = (np.sqrt(5.0) + 1.0) / (2.0 * np.sqrt(5.0))
    v1, v2 = (1.0 - np.sqrt(5.0)) / 2.0, (1.0 + np.sqrt(5.0)) / 2.0
    webb = np.array([-np.sqrt(1.5), -1.0, -np.sqrt(0.5), np.sqrt(0.5), 1.0, np.sqrt(1.5)])
    cnt = 0
    for _ in range(B):
        if weight == "webb6":
            wt = rng_b.choice(webb, size=G)
        elif weight == "mammen":
            wt = np.where(rng_b.random(G) < p1, v1, v2)
        else:
            wt = rng_b.choice([-1.0, 1.0], size=G)
        wb = psi * wt
        mb = float(wb.mean())
        seb = float(np.sqrt(np.sum((wb - mb) ** 2) / (G * (G - 1))))
        if abs(mb / seb) >= t_obs:
            cnt += 1
    return (cnt + 1) / (B + 1)


def wcr_score_ci(panel: dict, beta_hat: np.ndarray, B: int = 999, seed: int = 20261119, weight: str = "webb6") -> dict:
    X = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, phi_bar=0.5)
    tpk_hat = float(beta_hat[0] + beta_hat[1])
    grid = np.linspace(tpk_hat - 0.35, tpk_hat + 0.35, 71)
    pvals = [wcr_score_p_tpk0(X, beta_hat, float(g), B=B, seed=seed, weight=weight) for g in grid]
    acc = [float(g) for g, pv in zip(grid, pvals) if np.isfinite(pv) and pv > 0.05]
    p_zero = wcr_score_p_tpk0(X, beta_hat, 0.0, B=B, seed=seed, weight=weight)
    return dict(
        weight=weight, B=B, p_at_zero=float(p_zero),
        ci=[float(min(acc)), float(max(acc))] if acc else [float("nan"), float("nan")],
    )


def ascertainment_tipping_point(d: dict, w: np.ndarray, lag: int = 1) -> dict:
    """Holiday-linked ascertainment tipping point: Y_it = c_i c_t (1 - rho)^{D_it} X_it.
    Reports both lagged-D correction (legacy) and contemporaneous-D0 correction with onset-moment matching."""
    D_lag = lag_matrix(d["D_summer"], lag)
    D0 = d["D_summer"]
    onset_next = ((D0[4:] == 1.0) & (D0[3:-1] == 0.0)).astype(float)

    def fit_rho(rho: float, D_exp: np.ndarray, lg: int, extra: bool = False):
        Y_corr = d["Y"] / ((1.0 - rho) ** D_exp)
        p = make_panel(d, w, D0, lag=lg, Y_override=Y_corr)
        r = fit_m0_full(p)
        if not r["ok"]:
            return dict(rho=float(rho), ok=False)
        X = QDX(p, {"tau": p["D"], "kappa": p["Db"]}, phi_bar=0.5)
        c = X._core(r["_beta"])
        m_on = float((X.s * onset_next * c["v"]).sum() / X.N)
        out = dict(
            rho=float(rho), ok=True,
            est=r["est"], se=r["se"],
            lo=r["ci_t"]["tpk"][0], hi=r["ci_t"]["tpk"][1],
            bound_log_shift=float(abs(np.log(1.0 - rho))),
            m_onset=m_on,
        )
        if extra:
            out["J_LM_oi"] = r["J_LM_oi"]
            out["J_LM_inv"] = r["J_LM_inv"]
        return out

    grid = [round(x, 3) for x in np.arange(-0.40, 0.5001, 0.05)]
    rows = [fit_rho(rho, D_lag, lag, extra=(abs(rho * 10 - round(rho * 10)) < 1e-9)) for rho in grid]

    def find_root(rows_in, D_exp, lg, key: str):
        xs = [r["rho"] for r in rows_in if r["ok"]]
        vs = [r["est"]["tpk"] if key == "tpk" else r[key] for r in rows_in if r["ok"]]
        for a, b, va, vb in zip(xs[:-1], xs[1:], vs[:-1], vs[1:]):
            if np.sign(va) != np.sign(vb):
                fn = lambda rr: (fit_rho(rr, D_exp, lg)["est"]["tpk"] if key == "tpk" else fit_rho(rr, D_exp, lg)[key])
                return float(optimize.brentq(fn, a, b, xtol=1e-4))
        return None

    rho_zero = find_root(rows, D_lag, lag, "tpk")
    rho_hi_zero = find_root(rows, D_lag, lag, "hi")

    # Contemporaneous-D0 correction profiling + onset-moment root m_onset = 0
    c_grid = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    c_rows_l1 = [fit_rho(rho, D0, 1, extra=True) for rho in c_grid]
    c_rows_l2 = [fit_rho(rho, D0, 2, extra=True) for rho in c_grid]
    rho_onset_l1 = find_root(c_rows_l1, D0, 1, "m_onset")
    rho_onset_l2 = find_root(c_rows_l2, D0, 2, "m_onset")
    return dict(
        grid=rows,
        rho_star_tpk_zero=rho_zero,
        rho_star_upper_ci_zero=rho_hi_zero,
        rho_star_tpk_zero_fit=fit_rho(rho_zero, D_lag, lag, extra=True) if rho_zero is not None else None,
        rho_star_upper_ci_zero_fit=fit_rho(rho_hi_zero, D_lag, lag, extra=True) if rho_hi_zero is not None else None,
        contemp_grid_lag1=c_rows_l1,
        contemp_grid_lag2=c_rows_l2,
        rho_star_onset_zero_lag1=rho_onset_l1,
        rho_star_onset_zero_lag2=rho_onset_l2,
        fit_at_onset_zero_lag1=fit_rho(rho_onset_l1, D0, 1, extra=True) if rho_onset_l1 is not None else None,
        fit_at_onset_zero_lag2=fit_rho(rho_onset_l2, D0, 2, extra=True) if rho_onset_l2 is not None else None,
    )


def fit_structured_spec(d: dict, w: np.ndarray, D_mat: np.ndarray | None = None, lag: int = 1,
                        keep_states_exclude: tuple[int, ...] = (), wgroup_mode: str | None = None,
                        trend: bool = False, harmonics: int = 0, year_slice: slice | None = None) -> dict:
    """Flexible QDX specification runner for state exclusions, trends, harmonics, placebos, and time shifts."""
    if D_mat is None:
        D_mat = d["D_summer"]
    keep_mask = ~np.isin(d["state"], list(keep_states_exclude))
    A_full = (d["W"] > 0).astype(float).tocsr()
    Ak = A_full[keep_mask][:, keep_mask]
    deg = np.asarray(Ak.sum(1)).ravel()
    assert (deg > 0).all()
    Wk = sps.diags(1.0 / deg) @ Ak
    Yk = d["Y"][:, keep_mask]
    Dk = lag_matrix(D_mat[:, keep_mask], lag)
    omk = d["pop"][keep_mask] / 1e5
    clk_raw = d["cl"][keep_mask]
    _, clk = np.unique(clk_raw, return_inverse=True)
    regk = d["region"][keep_mask]
    wks = d["weeks"]
    if year_slice is not None:
        Yk = Yk[year_slice]
        Dk = Dk[year_slice]
        wks = wks[year_slice]
    panel = dgp.build_panel(Yk, Dk, Wk.tocsr(), omk, clk, w=w, S=3)
    lin = {"tau": panel["D"], "kappa": panel["Db"]}
    T_p = panel["Y"].shape[0]
    est_wks = wks[3:, 1]
    t_lin = np.linspace(-1.0, 1.0, T_p)[:, None]
    G_sub = int(panel["G"])
    if trend:
        for g in range(1, G_sub):
            mask_g = (clk == g).astype(float)[None, :]
            lin[f"tr_{g}"] = t_lin * mask_g
    if harmonics >= 1:
        ang = 2.0 * np.pi * (est_wks / 52.1775)[:, None]
        for h in range(1, harmonics + 1):
            s_v = np.sin(h * ang)
            c_v = np.cos(h * ang)
            for g in range(1, G_sub):
                mask_g = (clk == g).astype(float)[None, :]
                lin[f"sin{h}_{g}"] = s_v * mask_g
                lin[f"cos{h}_{g}"] = c_v * mask_g
    wg = None
    if wgroup_mode == "ew":
        wg = (regk == 3).astype(int)
    elif wgroup_mode == "reg4":
        _, wg = np.unique(regk, return_inverse=True)
    exo, _ = oi_extras(panel, phi_bar=0.5)
    exi, _ = inv_extras(panel)
    X = QDX(panel, lin, wgroup=wg, phi_bar=0.5, extras=exo + exi)
    x0 = np.r_[[-0.2, 0.0], np.zeros(X.K - 2), [0.32]]
    x, ok, res = X.fit(starts=[x0, np.r_[np.zeros(X.K), 0.5]])
    if not ok:
        return dict(ok=False, res=float(res))
    inf = X.inference(x)
    Vb = inf["Vb"]
    se = np.sqrt(np.diag(Vb))
    G = X.G
    tc = float(stats.t.ppf(0.975, G - 1))
    tpk = float(x[0] + x[1])
    tpk_se = tsum(Vb, 0, 1)
    p_idx = X.np_ - 1
    j_oi = X.j_lm(x, sel=[0, 1])
    j_inv = X.j_lm(x, sel=[2, 3])
    return dict(
        ok=True, N=int(panel["Y"].shape[1]), G=G, T=int(panel["Y"].shape[0]),
        tcrit=tc,
        est=dict(tau=float(x[0]), kappa=float(x[1]), phi=float(x[p_idx]), tpk=tpk),
        se=dict(tau=float(se[0]), kappa=float(se[1]), phi=float(se[p_idx]), tpk=tpk_se),
        ci_t=dict(
            tau=[float(x[0] - tc * se[0]), float(x[0] + tc * se[0])],
            kappa=[float(x[1] - tc * se[1]), float(x[1] + tc * se[1])],
            phi=[float(x[p_idx] - tc * se[p_idx]), float(x[p_idx] + tc * se[p_idx])],
            tpk=[float(tpk - tc * tpk_se), float(tpk + tc * tpk_se)],
        ),
        p_t_tpk=float(2.0 * stats.t.sf(abs(tpk / tpk_se), G - 1)),
        J_LM_oi=dict(stat=j_oi["J_LM"], df=2, p_chi2=j_oi["p"], p_F=float(stats.f.sf(j_oi["J_LM"] * (G - 2) / (2 * (G - 1)), 2, G - 2))),
        J_LM_inv=dict(stat=j_inv["J_LM"], df=2, p_chi2=j_inv["p"], p_F=float(stats.f.sf(j_inv["J_LM"] * (G - 2) / (2 * (G - 1)), 2, G - 2))),
    )


def state_permutation_placebo(d: dict, w: np.ndarray, lag: int = 1, B: int = 200, seed: int = 20261119) -> dict:
    """Permute the 16 federal states' summer holiday schedules across states (randomization inference)."""
    rng = np.random.default_rng(seed + lag)
    st_cols = [np.where(d["cl"] == g)[0][0] for g in range(len(d["states"]))]
    D_st = d["D_summer"][:, st_cols]
    p0 = make_panel(d, w, d["D_summer"], lag=lag)
    r0 = qdgmm.fit(p0, qdgmm.spec_m0(p0), starts=[[-0.3, 0.0, 0.3], [0.0, 0.0, 0.5]])
    obs = float(r0["est"]["tau"] + r0["est"]["kappa"])
    tpks = []
    for _ in range(B):
        perm = rng.permutation(len(st_cols))
        Dp = D_st[:, perm][:, d["cl"]]
        pp = make_panel(d, w, Dp, lag=lag)
        rp = qdgmm.fit(pp, qdgmm.spec_m0(pp), starts=[[-0.2, 0.0, 0.3], [0.0, 0.0, 0.5]])
        if rp["ok"]:
            tpks.append(float(rp["est"]["tau"] + rp["est"]["kappa"]))
    arr = np.asarray(tpks, float)
    return dict(
        lag=lag, B=B, n_ok=int(len(arr)), obs_tpk=obs,
        mean=float(arr.mean()), sd=float(arr.std()),
        q05=float(np.quantile(arr, 0.05)), q50=float(np.median(arr)), q95=float(np.quantile(arr, 0.95)),
        p_one_sided=float((arr <= obs).mean()), p_two_sided=float((np.abs(arr) >= abs(obs)).mean()),
    )


def nuts2_aggregation(d: dict, w: np.ndarray, lag: int = 1) -> dict:
    """Aggregate 400 counties to 38 NUTS-2 regions (G=16 states, low-zero regime)."""
    n2 = d["nuts2"]
    K = int(n2.max()) + 1
    Y2 = np.stack([d["Y"][:, n2 == k].sum(1) for k in range(K)], 1)
    D2 = np.stack([d["D_summer"][:, n2 == k].mean(1) for k in range(K)], 1)
    D2 = (D2 >= 0.5).astype(float)
    pop2 = np.array([d["pop"][n2 == k].sum() for k in range(K)])
    cl2 = np.array([int(d["cl"][n2 == k][0]) for k in range(K)], dtype=int)
    A2 = np.zeros((K, K), float)
    for i, j in d["edges"]:
        if n2[i] != n2[j]:
            A2[n2[i], n2[j]] = 1.0
            A2[n2[j], n2[i]] = 1.0
    deg2 = A2.sum(1)
    W2 = sps.csr_matrix(A2 / np.maximum(deg2[:, None], 1.0))

    D2l = lag_matrix(D2, lag)
    p = dgp.build_panel(Y2, D2l, W2, pop2 / 1e5, cl2, w=w, S=3)
    m0 = strip_private(fit_m0_full(p))
    m0p = qdgmm.fit(p, qdgmm.spec_m0prime(p))
    comps = {}
    for nm, fn in [("M1", cmp.poisson_fe), ("M2", cmp.m2_logr_did), ("M3", cmp.m3_logy_twfe), ("M4", cmp.poisson_qmle)]:
        r = fn(p)
        comps[nm] = dict(ok=r["ok"], est=r.get("est"), se=r.get("se"))
    return dict(
        N=K, G=int(p["G"]), T=int(p["Y"].shape[0]),
        median_Y=float(np.median(p["Y"])),
        mean_Y=float(p["Y"].mean()),
        zero_share=float((p["Y"] == 0).mean()),
        M0=m0,
        M0_prime=dict(est=m0p["est"], se=m0p["se"]),
        comparators=comps,
    )


def main():
    t0 = time.time()
    d = load_data()
    summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    gi_main = gi_weekly(14.0, 2.5, S=3)
    w = np.array(gi_main["w"])
    print(f"Varicella GI (14.0d, SD 2.5d, S=3): w={np.round(w, 4)}, w0={gi_main['w0']:.2e}, tail={gi_main['tail_beyond_S']:.4f}", flush=True)

    # 1. Primary panel & M0
    panel = make_panel(d, w, d["D_summer"], lag=1, omega_mult=1.0, S=3)
    m0_raw = fit_m0_full(panel, phi_bar=0.5)
    beta0 = m0_raw["_beta"]
    Vb0 = m0_raw["_Vb"]
    inf0 = m0_raw["_inf"]
    m0 = strip_private(m0_raw)

    # Wild bootstrap-t (B=999, Mammen & Webb6) + WCR-score + Nonlinear pairs cluster moment bootstrap-t (B=999)
    sp0 = qdgmm.spec_m0(panel, phi_bar=0.5)
    mod0 = qdgmm.QDModel(panel, sp0)
    phi_g = inf0["infl_b"]
    wild_mam = WB.score_wild_bootstrap(phi_g, 999, 20261119, "mammen")
    wild_web = WB.score_wild_bootstrap(phi_g, 999, 20261120, "webb6")
    wcr_web = wcr_score_ci(panel, beta0, 999, 20261119, "webb6")
    wcr_rad = wcr_score_ci(panel, beta0, 999, 20261120, "rademacher")
    print("Running nonlinear pairs cluster bootstrap (B=999)...", flush=True)
    pairs = WB.pairs_moment_bootstrap(mod0, beta0, 999, 20261119)
    tpk = m0["est"]["tpk"]
    tpk_se = m0["se"]["tpk"]
    ci_wild_mam = [float(tpk - wild_mam["q95"][3] * tpk_se), float(tpk + wild_mam["q95"][3] * tpk_se)]
    ci_wild_web = [float(tpk - wild_web["q95"][3] * tpk_se), float(tpk + wild_web["q95"][3] * tpk_se)]
    ci_pairs = [float(tpk - pairs["q95"][3] * tpk_se), float(tpk + pairs["q95"][3] * tpk_se)]
    print(f"M0 primary: tpk={tpk:+.3f} (SE {tpk_se:.3f}) | t15={np.round(m0['ci_t']['tpk'], 3)} | WCR(Webb6)={np.round(wcr_web['ci'], 3)} (p0={wcr_web['p_at_zero']:.4f}) | pairs({pairs['ok']}/999)={np.round(ci_pairs, 3)}", flush=True)

    # Delete-one-state jackknife
    jk = jackknife_m0(d, w, lag=1, phi_bar=0.5)
    tc = m0["tcrit"]
    ci_jk = [float(tpk - tc * jk["se"]["tpk"]), float(tpk + tc * jk["se"]["tpk"])]
    print(f"Jackknife: SE(tpk)={jk['se']['tpk']:.3f}, CI={np.round(ci_jk, 3)}, LOO range={np.round(jk['tpk_range'], 3)}", flush=True)

    # Cases averted psi^dagger (Theorem 4) & omega-scaling sensitivity of psi^dagger
    def idx_fn(b):
        return b[0] * panel["D"] + b[1] * panel["Db"], float(b[2])
    psi_res = psi_mod.psi_with_se(panel["Yall"], panel, beta0, Vb0, idx_fn)
    Y_treated_sum = float(panel["Y"][panel["D"] == 1].sum())
    Y_total_sum = float(panel["Y"].sum())
    psi_res["treated"]["obs_cases"] = Y_treated_sum
    psi_res["treated"]["pct_of_counterfactual"] = float(100.0 * psi_res["treated"]["est"] / (Y_treated_sum + psi_res["treated"]["est"]))
    psi_res["network"]["obs_cases"] = Y_total_sum
    psi_res["network"]["pct_of_counterfactual"] = float(100.0 * psi_res["network"]["est"] / (Y_total_sum + psi_res["network"]["est"]))
    print(f"Cases averted: treated={psi_res['treated']['est']:.0f} (SE {psi_res['treated']['se']:.0f}), network={psi_res['network']['est']:.0f} (SE {psi_res['network']['se']:.0f})", flush=True)

    # 2. Global root search (M0-G, M0-GI)
    glob = fit_global_score(panel, phi_bar=0.5)
    print(f"Global roots: n_interior={glob['n_interior_roots']}, phis={glob['interior_phis']}, M0-G={glob['M0_G']['kind']} ({glob['M0_G']['est']['phi']:.3f}), M0-GI={glob['M0_GI']['kind']} ({glob['M0_GI']['est']['phi']:.3f})", flush=True)

    # 3. Comparators (M0', M1, M2, M3, M4)
    m0p = qdgmm.fit(panel, qdgmm.spec_m0prime(panel))
    comparators = dict(M0_prime=dict(ok=m0p["ok"], est=m0p["est"], se=m0p["se"]))
    for nm, fn in [("M1", cmp.poisson_fe), ("M2", cmp.m2_logr_did), ("M3", cmp.m3_logy_twfe), ("M4", cmp.poisson_qmle)]:
        r = fn(panel)
        rec = dict(ok=r["ok"], est=r.get("est"), se=r.get("se"))
        if r["ok"] and "kappa" in r["est"] and r.get("Vb") is not None:
            rec["est"]["tpk"] = float(r["est"]["tau"] + r["est"]["kappa"])
            rec["se"]["tpk"] = tsum(r["Vb"], 0, 1)
            rec["ci_t_tpk"] = [float(rec["est"]["tpk"] - tc * rec["se"]["tpk"]),
                               float(rec["est"]["tpk"] + tc * rec["se"]["tpk"])]
        elif r["ok"]:
            rec["ci_t_tau"] = [float(r["est"]["tau"] - tc * r["se"]["tau"]),
                               float(r["est"]["tau"] + tc * r["se"]["tau"])]
        comparators[nm] = rec
        print(f"  {nm}: est={r.get('est')}, se={r.get('se')}", flush=True)
    m2_priors = {}
    for lbl, a_pr, b_pr in [("default_a1_b5", 1.0, 5.0), ("moderate_a0.25_b5", 0.25, 5.0),
                            ("weak_a0.1_b5", 0.1, 5.0), ("ultralow_a0.01_b100", 0.01, 100.0)]:
        rp = cmp.twfe_ols(cmp.cori_log_r(panel, a=a_pr, b=b_pr), panel["D"], panel["cl"], int(panel["G"]))
        m2_priors[lbl] = dict(a=a_pr, b=b_pr, tau=rp["est"]["tau"], se=rp["se"]["tau"])
    comparators["M2_prior_sensitivity"] = m2_priors

    # 4. Remark 4 adjacent-state score dependence
    rem4 = remark4_adjacent_dependence(d, panel, beta0)
    print(f"Remark 4 adj score corr: {rem4['adj_pair_score_corr']} | SE(tpk) {rem4['se_cluster']['tpk']:.3f} -> {rem4['se_with_adjacent']['tpk']:.3f} (non_psd={rem4['non_psd']})", flush=True)

    # 5. Stock-Wright robust sets (phi_bar = 0.5, 0.25, 1.0)
    print("Computing Stock-Wright robust sets (304-pt grid x 3 phi_bar)...", flush=True)
    sw = {}
    for pb in (0.5, 0.25, 1.0):
        sw[f"phi_bar_{pb:g}"] = stock_wright_sets(panel, phi_bar=pb)
        print(f"  phi_bar={pb}: {sw[f'phi_bar_{pb:g}']['sets']}", flush=True)

    # 6. Event-time profiles (with and without pre-trend leads, direct e_k + policy-total tau_k+c_k*kappa)
    ev_leads = event_time_analysis(d, w, include_leads=True)
    ev_noleads = event_time_analysis(d, w, include_leads=False)
    print(f"Event-time (with leads): tot_e0={ev_leads['total_scale']['e0']['est']:+.3f}(p={ev_leads['total_scale']['e0']['p_t']:.3f}), tot_e1={ev_leads['total_scale']['e1']['est']:+.3f}, tot_e3={ev_leads['total_scale']['e3']['est']:+.3f}, reopen_pF={ev_leads['reopen_total_wald_e6_e7']['p_F']:.3f}", flush=True)

    # 7. Sensitivity battery & omega psi sensitivity
    print("Running sensitivity battery...", flush=True)
    sens = {}
    omega_psi = {}
    for pb in (0.25, 1.0):
        sens[f"phi_bar_{pb:g}"] = strip_private(fit_m0_full(panel, phi_bar=pb))
    for om in (0.25, 4.0):
        p_om = make_panel(d, w, d["D_summer"], lag=1, omega_mult=om)
        r_om = fit_m0_full(p_om)
        sens[f"omega_x_{om:g}"] = strip_private(r_om)
        def _idx_om(b, p_o=p_om):
            return b[0] * p_o["D"] + b[1] * p_o["Db"], float(b[2])
        omega_psi[f"omega_x_{om:g}"] = psi_mod.psi_with_se(p_om["Yall"], p_om, r_om["_beta"], r_om["_Vb"], _idx_om)
    for lg in (0, 2, 3, 4):
        p_lg = make_panel(d, w, d["D_summer"], lag=lg)
        sens[f"summer_lag_{lg}"] = strip_private(fit_m0_full(p_lg))
    sens["summer_Dshare_lag1"] = strip_private(fit_m0_full(make_panel(d, w, d["Dshare_summer"], lag=1)))
    sens["summer_Dshare_lag2"] = strip_private(fit_m0_full(make_panel(d, w, d["Dshare_summer"], lag=2)))
    sens["summer_min3d_lag1"] = strip_private(fit_m0_full(make_panel(d, w, d["D_summer_3d"], lag=1)))
    sens["summer_min3d_lag2"] = strip_private(fit_m0_full(make_panel(d, w, d["D_summer_3d"], lag=2)))
    sens["summer_min5d_lag1"] = strip_private(fit_m0_full(make_panel(d, w, d["D_summer_5d"], lag=1)))
    sens["summer_min5d_lag2"] = strip_private(fit_m0_full(make_panel(d, w, d["D_summer_5d"], lag=2)))

    # Generation interval perturbations
    gi_sens = {}
    for nm, mn, sd_ in [
        ("Vink_main_14d_sd2.5", 14.0, 2.5),
        ("mean_12d_sd2.5", 12.0, 2.5),
        ("mean_16d_sd2.5", 16.0, 2.5),
        ("mean_14d_sd4.0", 14.0, 4.0),
        ("mean_14d_sd6.0_report", 14.0, 6.0),
    ]:
        g_i = gi_weekly(mn, sd_, S=3)
        p_gi = make_panel(d, np.array(g_i["w"]), d["D_summer"], lag=1)
        gi_sens[nm] = dict(gi=g_i, fit=strip_private(fit_m0_full(p_gi)))

    # Regional x week FE, non-rotating large-state exclusions (BW & BY), state trends & harmonics, 2016-2019
    wks = d["weeks"]
    sl_16_19 = slice(int(np.where(wks[:, 0] == 2016)[0][0]) - 3, None)
    regional = dict(
        four_region_x_week_lag1=strip_private(fit_m0_full(panel, wgroup=d["region"])),
        four_region_x_week_lag2=strip_private(fit_m0_full(make_panel(d, w, d["D_summer"], lag=2), wgroup=d["region"])),
        east_west_x_week_lag1=strip_private(fit_m0_full(panel, wgroup=(d["region"] == 3).astype(int))),
        drop_bw_by_lag1=fit_structured_spec(d, w, keep_states_exclude=(8, 9), lag=1),
        drop_bw_by_lag2=fit_structured_spec(d, w, keep_states_exclude=(8, 9), lag=2),
        reg4_drop_bw_by_lag1=fit_structured_spec(d, w, keep_states_exclude=(8, 9), wgroup_mode="reg4", lag=1),
        state_trends_lag1=fit_structured_spec(d, w, trend=True, lag=1),
        state_trends_lag2=fit_structured_spec(d, w, trend=True, lag=2),
        state_harmonics1_lag1=fit_structured_spec(d, w, harmonics=1, lag=1),
        state_harmonics1_lag2=fit_structured_spec(d, w, harmonics=1, lag=2),
        state_harm1_trends_drop_bw_by_lag1=fit_structured_spec(d, w, keep_states_exclude=(8, 9), trend=True, harmonics=1, lag=1),
        y2016_2019_lag1=fit_structured_spec(d, w, year_slice=sl_16_19, lag=1),
        y2016_2019_drop_bw_by_lag1=fit_structured_spec(d, w, keep_states_exclude=(8, 9), year_slice=sl_16_19, lag=1),
    )

    # Individual calendar years (2015..2019)
    by_year = {}
    for yr in range(2015, 2020):
        idx_yr = np.where(wks[:, 0] == yr)[0]
        sl = slice(int(idx_yr[0]) - 3, int(idx_yr[-1]) + 1)
        p_yr = make_panel(d, w, d["D_summer"], lag=1, t_slice=sl)
        by_year[str(yr)] = strip_private(fit_m0_full(p_yr))

    # All school holidays (2015-2019 lag 0, lag 1, lag 2; 2017-2018 lag 0)
    all_hol = dict(
        all_2015_2019_lag0=strip_private(fit_m0_full(make_panel(d, w, d["D_all"], lag=0))),
        all_2015_2019_lag1=strip_private(fit_m0_full(make_panel(d, w, d["D_all"], lag=1))),
        all_2015_2019_lag2=strip_private(fit_m0_full(make_panel(d, w, d["D_all"], lag=2))),
    )
    idx_17 = np.where(wks[:, 0] == 2017)[0][0]
    idx_18_end = np.where(wks[:, 0] == 2018)[0][-1]
    all_hol["all_2017_2018_lag0"] = strip_private(fit_m0_full(make_panel(d, w, d["D_all"], lag=0, t_slice=slice(int(idx_17) - 3, int(idx_18_end) + 1))))

    # Placebos: other holidays, time-shifted summer, cross-state schedule permutations
    print("Running placebo and permutation suite...", flush=True)
    other_hol_placebos = {}
    for nm, Dmat in [("easter", d["D_easter"]), ("autumn", d["D_autumn"]), ("winter", d["D_winter"]), ("pentecost", d["D_pentecost"])]:
        other_hol_placebos[nm] = {f"lag_{lg}": strip_private(fit_m0_full(make_panel(d, w, Dmat, lag=lg))) for lg in (0, 1, 2)}
    time_shifts = {}
    for sh in (-4, -3, -2, -1, 0, 1, 2, 3, 4):
        if sh < 0:
            k = -sh
            Ds = np.vstack([d["D_summer"][k:], np.zeros((k, d["D_summer"].shape[1]))])
            lg_eff = 1
        else:
            Ds = d["D_summer"]
            lg_eff = 1 + sh
        time_shifts[f"shift_{sh:+d}"] = strip_private(fit_m0_full(make_panel(d, w, Ds, lag=lg_eff)))
    perm_placebo = dict(
        lag1=state_permutation_placebo(d, w, lag=1, B=200, seed=20261119),
        lag2=state_permutation_placebo(d, w, lag=2, B=200, seed=20261119),
    )

    # NUTS-2 spatial aggregation
    nuts2_res = nuts2_aggregation(d, w, lag=1)
    print(f"NUTS-2 (N=38, zero%={100*nuts2_res['zero_share']:.1f}%): M0 tpk={nuts2_res['M0']['est']['tpk']:+.3f} {np.round(nuts2_res['M0']['ci_t']['tpk'], 3)}, M2 tau={nuts2_res['comparators']['M2']['est']['tau']:+.3f}", flush=True)

    # Holiday ascertainment tipping point & contemporaneous reporting profile
    print("Running holiday ascertainment tipping point...", flush=True)
    asc = ascertainment_tipping_point(d, w, lag=1)
    print(f"Ascertainment tipping: rho*(tpk=0)={asc['rho_star_tpk_zero']}, rho*(onset=0, lag1)={asc['rho_star_onset_zero_lag1']}", flush=True)

    out = dict(
        panel_summary=summary,
        gi_main=gi_main,
        primary=dict(
            M0=m0,
            wild_mammen=dict(q95=wild_mam["q95"], ci_tpk=ci_wild_mam, B=999),
            wild_webb6=dict(q95=wild_web["q95"], ci_tpk=ci_wild_web, B=999),
            wcr_webb6=wcr_web,
            wcr_rademacher=wcr_rad,
            pairs_bootstrap=dict(q95=pairs["q95"], ci_tpk=ci_pairs, ok=pairs["ok"], B=999),
            jackknife=dict(**jk, ci_tpk=ci_jk),
            psi=psi_res,
            omega_psi_sensitivity=omega_psi,
        ),
        global_roots=glob,
        comparators=comparators,
        remark4=rem4,
        stock_wright=sw,
        event_time_with_leads=ev_leads,
        event_time_no_leads=ev_noleads,
        sensitivities=sens,
        gi_sensitivity=gi_sens,
        regional_fe=regional,
        by_year=by_year,
        all_holidays=all_hol,
        placebos=dict(other_holidays=other_hol_placebos, time_shifts=time_shifts, state_permutations=perm_placebo),
        nuts2=nuts2_res,
        ascertainment=asc,
        seconds=round(time.time() - t0, 2),
    )
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nSaved {OUT_PATH} in {out['seconds']}s")


if __name__ == "__main__":
    main()
