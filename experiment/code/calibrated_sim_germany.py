"""Calibrated simulation on the exact German Varicella 400-county x 261-week design (R = 500).

Evaluates finite-sample inference under G = 16 state clusters (sandwich t_15, score wild
bootstrap-t, WCR-score Webb6/Rademacher, and delete-one-state jackknife), comparator bias
under known truth (M2 EpiEstim log-R at county vs. NUTS-2 aggregation, and M3 TWFE), root
multiplicity, and calibrated null reference distributions for J_LM and S3.

Reads:  experiment/data/processed/panel.npz, experiment/analysis/germany_application_results.json
Writes: experiment/analysis/germany_calibrated_sim.json
"""
from __future__ import annotations

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import scipy.sparse as sps
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiment" / "code"))
sys.path.insert(0, str(ROOT / "verification" / "revision" / "common"))
sys.path.insert(0, str(ROOT / "verification" / "application"))

import application as app
import comparators as cmp
import dgp
import wild_bootstrap as WB
from qdx import QDX, inv_extras, oi_extras, tsum

OUT_PATH = ROOT / "experiment" / "analysis" / "germany_calibrated_sim.json"

_G_STATE = {}


def _init_worker(state_dict):
    _G_STATE.clear()
    _G_STATE.update(state_dict)


def _run_one(args):
    r, seed_base, k_nb, sig_sw = args
    st = _G_STATE
    Y_pre = st["Y_pre"]
    D_full = st["D_full"]
    W = st["W"]
    omega = st["omega"]
    cl = st["cl"]
    w = st["w"]
    ea = st["ea"]
    eg = st["eg"]
    beta_true = st["beta_true"]
    G = st["G"]
    T_sim = st["T_sim"]
    N = st["N"]
    tc = st["tc"]
    jk_pre = st["jk_pre"]
    nuts_agg = st["nuts_agg"]

    rng = np.random.default_rng(seed_base + r)
    Ys = np.zeros((3 + T_sim, N), dtype=float)
    Ys[:3] = Y_pre
    for j in range(T_sim):
        t = 3 + j
        L = w[0] * Ys[t - 1] + w[1] * Ys[t - 2] + w[2] * Ys[t - 3]
        lam = L + beta_true[2] * (W @ L) + omega
        mu = ea * eg[j] * np.exp(beta_true[0] * D_full[t] + beta_true[1] * (W @ D_full[t])) * lam
        if sig_sw > 0:
            eps = sig_sw * rng.standard_normal(G)
            mu = mu * np.exp(eps - 0.5 * sig_sw ** 2)[cl]
        Ys[t] = dgp._draw_nb(rng, mu, k_nb)

    ps = dgp.build_panel(Ys, D_full, W, omega, cl, w=w, S=3)
    exo, _ = oi_extras(ps, phi_bar=0.5)
    exi, _ = inv_extras(ps)
    X = QDX(ps, {"tau": ps["D"], "kappa": ps["Db"]}, phi_bar=0.5, extras=exo + exi)
    eh, ok, _ = X.fit(starts=[beta_true, [-0.2, 0.0, 0.3], [0.0, 0.0, 0.5]])
    if not ok:
        return dict(ok=False)
    inf = X.inference(eh)
    if inf is None:
        return dict(ok=False)
    Vb = inf["Vb"]
    se = np.sqrt(np.diag(Vb))
    tpk = float(eh[0] + eh[1])
    tpk_se = tsum(Vb, 0, 1)
    tpk_true = float(beta_true[0] + beta_true[1])

    # Score wild bootstrap-t (B=399) + WCR-score (Webb6 & Rademacher)
    wm = WB.score_wild_bootstrap(inf["infl_b"], 399, seed_base + 100000 + r, "mammen")
    ww = WB.score_wild_bootstrap(inf["infl_b"], 399, seed_base + 200000 + r, "webb6")
    pw_w6 = app.wcr_score_p_tpk0(X, eh, tpk_true, 399, seed_base + 300000 + r, "webb6")
    pw_rad = app.wcr_score_p_tpk0(X, eh, tpk_true, 399, seed_base + 400000 + r, "rademacher")

    # Delete-one-state jackknife (G=16) for tpk
    jk_tpks = []
    for keep, Dg, Wg, om_g, cl_inv in jk_pre:
        ps_g = dgp.build_panel(Ys[:, keep], Dg, Wg, om_g, cl_inv, w=w, S=3)
        Xg = QDX(ps_g, {"tau": ps_g["D"], "kappa": ps_g["Db"]}, phi_bar=0.5)
        eh_g, ok_g, _ = Xg.fit(starts=[eh, beta_true])
        if ok_g:
            jk_tpks.append(float(eh_g[0] + eh_g[1]))
    if len(jk_tpks) == G:
        jk_arr = np.array(jk_tpks)
        jk_se = float(np.sqrt(((G - 1.0) / G) * np.sum((jk_arr - jk_arr.mean()) ** 2)))
        cov_jk = bool(abs(tpk - tpk_true) <= tc * jk_se)
    else:
        jk_se = np.nan
        cov_jk = False

    # Comparators under known truth: county M2, county M3, and NUTS-2 M0 & M2
    m2_c = cmp.m2_logr_did(ps)
    m3_c = cmp.m3_logy_twfe(ps)
    m2_c_tau = float(m2_c["est"]["tau"]) if m2_c.get("ok") else np.nan
    m3_c_tau = float(m3_c["est"]["tau"]) if m3_c.get("ok") else np.nan
    c_zero_share = float((ps["Y"] == 0).mean())

    Y2s = Ys @ nuts_agg["A_sum"]
    p2 = dgp.build_panel(Y2s, nuts_agg["D2l"], nuts_agg["W2"], nuts_agg["om2"], nuts_agg["cl2"], w=w, S=3)
    X2 = QDX(p2, {"tau": p2["D"], "kappa": p2["Db"]}, phi_bar=0.5)
    eh2, ok2, _ = X2.fit(starts=[eh, beta_true])
    m0_n2_tpk = float(eh2[0] + eh2[1]) if ok2 else np.nan
    m0_n2_phi = float(eh2[2]) if ok2 else np.nan
    m2_n2 = cmp.m2_logr_did(p2)
    m2_n2_tau = float(m2_n2["est"]["tau"]) if m2_n2.get("ok") else np.nan
    n2_zero_share = float((p2["Y"] == 0).mean())

    # Over-identification statistics
    mX = X.moments(eh, extra=True)
    j_oi = X.j_lm(eh, sel=[0, 1], m=mX)["J_LM"]
    j_inv = X.j_lm(eh, sel=[2, 3], m=mX)["J_LM"]
    s3_oi = app.score_s3_at_beta(X, eh, [0, 1])["S3"]
    s3_inv = app.score_s3_at_beta(X, eh, [2, 3])["S3"]

    # Check root multiplicity on compact grid
    grid_chk = np.array([0.0, 0.05, 0.15, 0.3, 0.5, 0.8, 1.2, 2.0, 3.0])
    cands, _ = X.global_roots(grid_chk)
    n_int = sum(1 for c in cands if c[0] == "interior")

    vals = np.array([eh[0], eh[1], eh[2], tpk])
    trues = np.array([beta_true[0], beta_true[1], beta_true[2], tpk_true])
    ses = np.array([se[0], se[1], se[2], tpk_se])
    err = np.abs(vals - trues)
    return dict(
        ok=True,
        vals=vals.tolist(),
        ses=ses.tolist(),
        cov_t=(err <= tc * ses).tolist(),
        cov_wild_mam=(err <= np.asarray(wm["q95"]) * ses).tolist(),
        cov_wild_web=(err <= np.asarray(ww["q95"]) * ses).tolist(),
        cov_wcr_web=bool(pw_w6 > 0.05) if np.isfinite(pw_w6) else False,
        cov_wcr_rad=bool(pw_rad > 0.05) if np.isfinite(pw_rad) else False,
        jk_se=float(jk_se),
        cov_jackknife=cov_jk,
        jk_ok=(len(jk_tpks) == G),
        c_zero_share=c_zero_share,
        n2_zero_share=n2_zero_share,
        m2_county=m2_c_tau,
        m3_county=m3_c_tau,
        m0_nuts2_tpk=m0_n2_tpk,
        m0_nuts2_phi=m0_n2_phi,
        m2_nuts2=m2_n2_tau,
        j_oi=float(j_oi),
        j_inv=float(j_inv),
        s3_oi=float(s3_oi),
        s3_inv=float(s3_inv),
        n_int=int(n_int),
    )


def calibrate_dgp(panel: dict, beta: np.ndarray) -> dict:
    Y = panel["Y"]
    D = panel["D"]
    Db = panel["Db"]
    cl = panel["cl"]
    G = int(panel["G"])
    lam = panel["L"] + beta[2] * panel["Nb"] + panel["omega"]
    B = np.exp(beta[0] * D + beta[1] * Db) * lam
    ea = np.ones(Y.shape[1])
    eg = np.ones(Y.shape[0])
    for _ in range(300):
        eg = Y.sum(1) / np.maximum((ea[None, :] * B).sum(1), 1e-12)
        ea = Y.sum(0) / np.maximum((eg[:, None] * B).sum(0), 1e-12)
    mu = ea[None, :] * eg[:, None] * B
    inv_k = float(((Y - mu) ** 2 - mu).sum() / (mu ** 2).sum())
    k_nb = 1.0 / inv_k if inv_k > 0 else 1e6
    Yg = np.stack([np.bincount(cl, Y[t], G) for t in range(Y.shape[0])])
    mg = np.stack([np.bincount(cl, mu[t], G) for t in range(Y.shape[0])])
    m2 = np.stack([np.bincount(cl, mu[t] ** 2, G) for t in range(Y.shape[0])])
    ok = Yg > 0
    lr = np.log(Yg[ok] / mg[ok])
    nb_var = ((mg + m2 / k_nb) / mg ** 2)[ok]
    sig_sw = float(np.sqrt(max(lr.var() - nb_var.mean(), 0.0)))
    return dict(ea=ea, eg=eg, k_nb=k_nb, sig_sw=sig_sw)


def _summ_est(arr: np.ndarray, truth: float) -> dict:
    v = arr[np.isfinite(arr)]
    b = v - truth
    return dict(
        n=int(len(v)),
        mean=float(v.mean()),
        bias=float(b.mean()),
        bias_mcse=float(b.std(ddof=1) / np.sqrt(len(v))),
        rmse=float(np.sqrt((b ** 2).mean())),
    )


def summarise_runs(res_list: list[dict], beta_true: np.ndarray, obs_stats: dict) -> dict:
    ok_res = [r for r in res_list if r["ok"]]
    R_ok = len(ok_res)
    vals = np.array([r["vals"] for r in ok_res])
    ses = np.array([r["ses"] for r in ok_res])
    cov_t = np.array([r["cov_t"] for r in ok_res])
    cov_wm = np.array([r["cov_wild_mam"] for r in ok_res])
    cov_ww = np.array([r["cov_wild_web"] for r in ok_res])
    cov_wcr_w = float(np.mean([r["cov_wcr_web"] for r in ok_res]))
    cov_wcr_r = float(np.mean([r["cov_wcr_rad"] for r in ok_res]))
    jk_ok = [r for r in ok_res if r.get("jk_ok")]
    cov_jk = float(np.mean([r["cov_jackknife"] for r in jk_ok]))
    mean_jk_se = float(np.mean([r["jk_se"] for r in jk_ok]))
    tpk_true = float(beta_true[0] + beta_true[1])
    trues = np.array([beta_true[0], beta_true[1], beta_true[2], tpk_true])
    names = ["tau", "kappa", "phi", "tpk"]
    params = {}
    for k, nm in enumerate(names):
        b = vals[:, k] - trues[k]
        rec_p = dict(
            true=float(trues[k]),
            bias=float(b.mean()),
            rmse=float(np.sqrt((b ** 2).mean())),
            sd=float(vals[:, k].std(ddof=1)),
            mean_se=float(ses[:, k].mean()),
            se_ratio=float(ses[:, k].mean() / vals[:, k].std(ddof=1)),
            cov_t15=float(cov_t[:, k].mean()),
            cov_wild_mam=float(cov_wm[:, k].mean()),
            cov_wild_web=float(cov_ww[:, k].mean()),
        )
        if nm == "tpk":
            rec_p["cov_wcr_webb6"] = cov_wcr_w
            rec_p["cov_wcr_rademacher"] = cov_wcr_r
            rec_p["cov_jackknife"] = cov_jk
            rec_p["mean_jk_se"] = mean_jk_se
            rec_p["jk_se_ratio"] = float(mean_jk_se / vals[:, k].std(ddof=1))
            rec_p["jk_ok"] = int(len(jk_ok))
        params[nm] = rec_p

    comparators = dict(
        truth_tpk=tpk_true,
        truth_phi=float(beta_true[2]),
        county_zero_share=float(np.mean([r["c_zero_share"] for r in ok_res])),
        nuts2_zero_share=float(np.mean([r["n2_zero_share"] for r in ok_res])),
        M0_county=_summ_est(vals[:, 3], tpk_true),
        M2_county=_summ_est(np.array([r["m2_county"] for r in ok_res]), tpk_true),
        M3_county=_summ_est(np.array([r["m3_county"] for r in ok_res]), tpk_true),
        M0_nuts2_tpk=_summ_est(np.array([r["m0_nuts2_tpk"] for r in ok_res]), tpk_true),
        M0_nuts2_phi=_summ_est(np.array([r["m0_nuts2_phi"] for r in ok_res]), float(beta_true[2])),
        M2_nuts2=_summ_est(np.array([r["m2_nuts2"] for r in ok_res]), tpk_true),
        M2_county_pos_share=float(np.mean(np.array([r["m2_county"] for r in ok_res]) > 0)),
        M2_nuts2_neg_share=float(np.mean(np.array([r["m2_nuts2"] for r in ok_res]) < 0)),
    )

    diag = {}
    for key in ("j_oi", "j_inv", "s3_oi", "s3_inv"):
        arr = np.array([r[key] for r in ok_res if np.isfinite(r[key])])
        obs_val = obs_stats[key]
        diag[key] = dict(
            n=len(arr),
            obs=float(obs_val),
            mean=float(arr.mean()),
            median=float(np.median(arr)),
            q90=float(np.quantile(arr, 0.90)),
            q95=float(np.quantile(arr, 0.95)),
            q99=float(np.quantile(arr, 0.99)),
            max=float(arr.max()),
            exceed_count=int((arr >= obs_val).sum()),
            p_calibrated=float(((arr >= obs_val).sum() + 1) / (len(arr) + 1)),
        )
    n_ints = np.array([r["n_int"] for r in ok_res])
    return dict(
        R_total=len(res_list),
        R_ok=R_ok,
        params=params,
        comparators=comparators,
        diagnostics=diag,
        root_counts=dict(
            zero=int((n_ints == 0).sum()),
            one=int((n_ints == 1).sum()),
            two_plus=int((n_ints >= 2).sum()),
            share_unique=float((n_ints == 1).mean()),
        ),
    )


def main():
    t0 = time.time()
    d = app.load_data()
    app_res = json.loads((ROOT / "experiment" / "analysis" / "germany_application_results.json").read_text(encoding="utf-8"))
    w = np.array(app_res["gi_main"]["w"])
    panel = app.make_panel(d, w, d["D_summer"], lag=1)
    m0 = app_res["primary"]["M0"]
    beta_true = np.array([m0["est"]["tau"], m0["est"]["kappa"], m0["est"]["phi"]])

    cal = calibrate_dgp(panel, beta_true)
    print(f"Calibrated DGP: k_NB={cal['k_nb']:.3f}, sig_sw={cal['sig_sw']:.4f}", flush=True)

    obs_stats = dict(
        j_oi=m0["J_LM_oi"]["stat"],
        j_inv=m0["J_LM_inv"]["stat"],
        s3_oi=app_res["global_roots"]["M0_G"]["S3_oi"]["S3"],
        s3_inv=app_res["global_roots"]["M0_GI"]["S3_inv"]["S3"],
    )

    D_full = app.lag_matrix(d["D_summer"], 1)
    G = int(panel["G"])
    cl = d["cl"]
    W = d["W"]
    omega = d["pop"] / 1e5
    jk_pre = []
    for g in range(G):
        keep = cl != g
        Wg = W[keep][:, keep]
        rs = np.asarray(Wg.sum(1)).ravel()
        Wg = sps.diags(np.where(rs > 0, 1.0 / np.maximum(rs, 1.0), 0.0)) @ Wg
        _, cl_inv = np.unique(cl[keep], return_inverse=True)
        jk_pre.append((keep, D_full[:, keep], Wg, omega[keep], cl_inv))

    n2 = d["nuts2"]
    K = int(n2.max()) + 1
    A_sum = np.zeros((d["Y"].shape[1], K), float)
    for k in range(K):
        A_sum[n2 == k, k] = 1.0
    D2 = np.stack([d["D_summer"][:, n2 == k].mean(1) for k in range(K)], 1)
    D2 = (D2 >= 0.5).astype(float)
    D2l = app.lag_matrix(D2, 1)
    pop2 = np.array([d["pop"][n2 == k].sum() for k in range(K)])
    cl2 = np.array([int(cl[n2 == k][0]) for k in range(K)], dtype=int)
    A2 = np.zeros((K, K), float)
    for i, j in d["edges"]:
        if n2[i] != n2[j]:
            A2[n2[i], n2[j]] = 1.0
            A2[n2[j], n2[i]] = 1.0
    W2 = sps.csr_matrix(A2 / np.maximum(A2.sum(1)[:, None], 1.0))
    nuts_agg = dict(A_sum=A_sum, D2l=D2l, W2=W2, om2=pop2 / 1e5, cl2=cl2)

    state_dict = dict(
        Y_pre=d["Y"][:3],
        D_full=D_full,
        W=W,
        omega=omega,
        cl=cl,
        w=w,
        ea=cal["ea"],
        eg=cal["eg"],
        beta_true=beta_true,
        G=G,
        T_sim=int(panel["Y"].shape[0]),
        N=int(panel["Y"].shape[1]),
        tc=float(stats.t.ppf(0.975, G - 1)),
        jk_pre=jk_pre,
        nuts_agg=nuts_agg,
    )

    R = 500
    n_workers = max(1, min(os.cpu_count() or 4, 4))
    print(f"Running R={R} base calibrated + R={R} high-noise stress draws on {n_workers} workers...", flush=True)

    out = dict(
        R=R,
        calibration=dict(k_nb=cal["k_nb"], sig_sw=cal["sig_sw"], beta_true=beta_true.tolist()),
        obs_stats=obs_stats,
    )
    for tag, seed_base, k_mult, s_mult in [
        ("calibrated", 2026111900, 1.0, 1.0),
        ("high_noise", 2026121900, 0.5, 2.0),
    ]:
        t1 = time.time()
        tasks = [(r, seed_base, cal["k_nb"] * k_mult, cal["sig_sw"] * s_mult) for r in range(R)]
        with ProcessPoolExecutor(max_workers=n_workers, initializer=_init_worker, initargs=(state_dict,)) as ex:
            results = list(ex.map(_run_one, tasks, chunksize=8))
        summ = summarise_runs(results, beta_true, obs_stats)
        summ["k_nb"] = float(cal["k_nb"] * k_mult)
        summ["sig_sw"] = float(cal["sig_sw"] * s_mult)
        summ["seconds"] = round(time.time() - t1, 2)
        out[tag] = summ
        p_tpk = summ["params"]["tpk"]
        cc = summ["comparators"]
        print(
            f"  [{tag}] {summ['R_ok']}/{R} ok in {summ['seconds']}s | tpk: bias={p_tpk['bias']:+.4f}, "
            f"t15={p_tpk['cov_t15']:.3f}, JK={p_tpk['cov_jackknife']:.3f}, WCR_W6={p_tpk['cov_wcr_webb6']:.3f}, "
            f"WCR_Rad={p_tpk['cov_wcr_rademacher']:.3f} | zero={cc['county_zero_share']:.3f}->{cc['nuts2_zero_share']:.3f} | "
            f"M2_c bias={cc['M2_county']['bias']:+.3f} (pos={cc['M2_county_pos_share']:.2f}), "
            f"M2_n2 bias={cc['M2_nuts2']['bias']:+.3f}, M0_n2 bias={cc['M0_nuts2_tpk']['bias']:+.3f}",
            flush=True,
        )

    out["total_seconds"] = round(time.time() - t0, 2)
    OUT_PATH.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Saved {OUT_PATH} in {out['total_seconds']}s")


if __name__ == "__main__":
    main()
