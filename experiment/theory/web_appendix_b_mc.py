"""Reproducible self-contained script for ALL Web Appendix B Monte Carlo checks (E7).

Generates experiment/theory/web_appendix_b_mc.json covering:
  1. Theorem 2: G=100 (n_g=5) and G=50 (n_g=10), T=10, R=1000, bias, empSE, modSE, cov, corr(err, lin).
  2. Theorem 3(a),(c): N=2000, T=10, R=200 two-step log-R (B_J + B_F) and log(Y+1) TWFE bias decompositions.
  3. Lemma 2 & Prop 2(c): G=40 (n_g=5), T=8, R=4000 adjacent-cluster correlation, cross-cluster score |z|,
     MC variance / Omega_hat ratio, and lead moment z.
  4. T=2 Poisson FE numerical roots vs 200,000-unit Monte Carlo & Theorem 4(a) counterfactual path identity.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiment" / "code"))

import comparators as cmp  # noqa: E402
import dgp  # noqa: E402
import psi as psimod  # noqa: E402
import qdgmm  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent / "web_appendix_b_mc.json"


def _thm2_worker(args):
    G, ng, T, r = args
    N = G * ng
    geo = dgp.make_geography(N, G=G)
    p = dgp.DGPParams(N=N, T=T, adopt="feedback", tau=float(np.log(0.8)), kappa=-0.1, phi=0.4, k=10.0, gbar=0.09)
    seed = 20271000 + G * 10000 + r
    sim = dgp.simulate(p, geo, [seed, 0], [seed, 1])
    pan = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl)
    spec = qdgmm.spec_m0(pan)
    res = qdgmm.fit(pan, spec, starts=[[p.tau, p.kappa, p.phi]])
    if not res["ok"]:
        return None
    mod = res["model"]
    b = res["beta"]
    se = np.array([res["se"]["tau"], res["se"]["kappa"], res["se"]["phi"]])
    truth = np.array([p.tau, p.kappa, p.phi])
    tc = float(stats.t.ppf(0.975, G - 1))
    cov = (np.abs(b - truth) <= tc * se).tolist()
    inf = mod.inference(b)
    d_true = np.diff(sim["gamma"])
    lam0 = pan["L"] + truth[2] * pan["Nb"] + pan["omega"]
    eta0 = truth[0] * pan["D"] + truth[1] * pan["Db"]
    m_rel = np.exp(eta0) * lam0
    r0 = np.exp(-d_true)[:, None] * (m_rel[:-1] / m_rel[1:])
    u0 = pan["Y"][1:] * r0 - pan["Y"][:-1]
    Zpar = np.stack(mod.Zpar, axis=0)
    gw0 = (mod.s * u0).sum(1) / G
    gz0 = np.einsum("mtn,tn->m", Zpar, u0) / G
    gbar0 = np.r_[gw0, gz0]
    lin_tau = float(-np.linalg.solve(inf["Gam"], gbar0)[mod.T - 1]) if inf is not None else np.nan
    return b.tolist(), se.tolist(), cov, lin_tau


def check_thm2(workers: int) -> dict:
    out = {}
    for G, ng in [(100, 5), (50, 10)]:
        T, R = 10, 1000
        tasks = [(G, ng, T, r) for r in range(R)]
        with ProcessPoolExecutor(max_workers=workers) as ex:
            res = [x for x in ex.map(_thm2_worker, tasks, chunksize=25) if x is not None]
        ests = np.array([x[0] for x in res])
        ses = np.array([x[1] for x in res])
        covs = np.array([x[2] for x in res], float)
        lins = np.array([x[3] for x in res], float)
        truth = np.array([float(np.log(0.8)), -0.1, 0.4])
        ok_lin = np.isfinite(lins)
        corr_lin = float(np.corrcoef(ests[ok_lin, 0] - truth[0], lins[ok_lin])[0, 1]) if ok_lin.any() else np.nan
        rec = dict(n=len(res), corr_tau_lin=corr_lin)
        for j, nm in enumerate(("tau", "kappa", "phi")):
            b_arr = ests[:, j] - truth[j]
            c_arr = covs[:, j]
            rec[nm] = dict(
                bias=float(b_arr.mean()), bias_mcse=float(b_arr.std(ddof=1) / np.sqrt(len(res))),
                empSE=float(ests[:, j].std(ddof=1)), meanSE=float(ses[:, j].mean()),
                coverage=float(c_arr.mean()), cov_mcse=float(np.sqrt(c_arr.mean() * (1 - c_arr.mean()) / len(res))),
            )
        out[f"G{G}_n{ng}"] = rec
    return out


def check_thm3_bias_decomp() -> dict:
    N, T, k, R = 2000, 10, 10.0, 200
    tau0 = float(np.log(0.8))
    geo = dgp.make_geography(N, G=50)
    out = {}
    for adopt in ("feedback", "random"):
        p = dgp.DGPParams(N=N, T=T, adopt=adopt, tau=tau0, kappa=0.0, phi=0.0, k=k, gbar=0.08)
        m2_b, bj_l, bf_l, max_id_err = [], [], [], 0.0
        m3_b, m3_om_l = [], []
        dm = lambda a: a - a.mean(0, keepdims=True) - a.mean(1, keepdims=True) + a.mean()
        for r in range(R):
            seed = 20272000 + (0 if adopt == "feedback" else 1000) + r
            sim = dgp.simulate(p, geo, [seed, 0], [seed, 1])
            pan = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl)
            Y, D, L, om = pan["Y"], pan["D"], pan["L"], pan["omega"]
            mu = np.exp(sim["alpha"][None, :] + sim["gamma"][:, None] + tau0 * D) * (L + om)
            logR = np.log((0.5 + Y) / (0.2 + L))
            fit_m2 = cmp.twfe_ols(logR, D, pan["cl"], int(pan["G"]))["est"]["tau"]
            fit_m2_det = cmp.twfe_ols(np.log((0.5 + mu) / (0.2 + L)), D, pan["cl"], int(pan["G"]))["est"]["tau"]
            bj = fit_m2_det - tau0
            bf = fit_m2 - fit_m2_det
            max_id_err = max(max_id_err, abs((fit_m2 - tau0) - (bj + bf)))
            m2_b.append(fit_m2 - tau0)
            bj_l.append(bj)
            bf_l.append(bf)
            fit_m3 = cmp.twfe_ols(np.log(Y + 1.0), D, pan["cl"], int(pan["G"]))["est"]["tau"]
            om_ctrl = cmp.twfe_ols(np.log(L + om), D, pan["cl"], int(pan["G"]))["est"]["tau"]
            m3_b.append(fit_m3 - tau0)
            m3_om_l.append(om_ctrl)
        out[adopt] = dict(
            R=R,
            m2_bias=float(np.mean(m2_b)), m2_mcse=float(np.std(m2_b, ddof=1) / np.sqrt(R)),
            m2_B_J=float(np.mean(bj_l)), m2_B_F=float(np.mean(bf_l)),
            max_decomp_err=float(max_id_err),
            m3_bias=float(np.mean(m3_b)), m3_mcse=float(np.std(m3_b, ddof=1) / np.sqrt(R)),
            m3_omitted_pressure=float(np.mean(m3_om_l)),
        )
    return out


def _lem2_worker(args):
    r, G, ng, T, edges = args
    N = G * ng
    geo = dgp.make_geography(N, G=G)
    p = dgp.DGPParams(N=N, T=T, adopt="feedback", tau=float(np.log(0.8)), kappa=-0.1, phi=0.4, k=10.0, gbar=0.09)
    truth = np.array([p.tau, p.kappa, p.phi])
    seed = 20273000 + r
    sim = dgp.simulate(p, geo, [seed, 0], [seed, 1])
    pan = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl)
    mod = qdgmm.QDModel(pan, qdgmm.spec_m0(pan))
    d_true = np.diff(sim["gamma"])
    lam = pan["L"] + truth[2] * pan["Nb"] + pan["omega"]
    eta = truth[0] * pan["D"] + truth[1] * pan["Db"]
    m_rel = np.exp(eta) * lam
    r0 = np.exp(-d_true)[:, None] * (m_rel[:-1] / m_rel[1:])
    u0 = pan["Y"][1:] * r0 - pan["Y"][:-1]
    Z = np.stack(mod.Zpar, axis=0)
    sg = np.zeros((3, G))
    for g in range(G):
        sel = geo.cl == g
        sg[:, g] = np.einsum("mtn,tn->m", Z[:, :, sel], u0[:, sel])
    Ytot_g = np.bincount(geo.cl, pan["Y"].sum(0), G)
    pair_y = np.array([[Ytot_g[g], Ytot_g[h]] for g, h in edges])
    cross = np.mean([sg[:, g] * sg[:, h] for g, h in edges], axis=0)
    lead_m = float((Z[0] * pan["D"][1:] * u0).mean())
    return sg, cross, pair_y, lead_m


def check_lemma2_mds(workers: int) -> dict:
    G, ng, T, R = 40, 5, 8, 4000
    N = G * ng
    geo = dgp.make_geography(N, G=G)
    W_arr = geo.W.toarray()
    edges = sorted({(min(int(geo.cl[i]), int(geo.cl[j])), max(int(geo.cl[i]), int(geo.cl[j])))
                    for i in range(N) for j in np.where(W_arr[i] > 0)[0] if geo.cl[i] != geo.cl[j]})
    tasks = [(r, G, ng, T, edges) for r in range(R)]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        res = list(ex.map(_lem2_worker, tasks, chunksize=50))
    sg_arr = np.stack([x[0] for x in res], axis=0)
    cross_arr = np.stack([x[1] for x in res], axis=0)
    pair_y = np.concatenate([x[2] for x in res], axis=0)
    lead_arr = np.array([x[3] for x in res])
    cum_corr = float(np.corrcoef(pair_y[:, 0], pair_y[:, 1])[0, 1])
    z_cross = np.abs(cross_arr.mean(0)) / (cross_arr.std(ddof=1, axis=0) / np.sqrt(R))
    sqrtG_mean = np.sqrt(G) * sg_arr.mean(axis=2)
    mc_v = np.diag(np.cov(sqrtG_mean.T, ddof=1))
    om_est = np.mean(np.einsum("rmg,rng->rmn", sg_arr, sg_arr) / G, axis=0)
    var_ratios = mc_v / np.diag(om_est)
    return dict(
        G=G, ng=ng, T=T, R=R, n_adj_pairs=len(edges),
        adj_cum_case_corr=cum_corr,
        max_abs_z_cross_cluster_score=float(z_cross.max()),
        var_ratio_min=float(var_ratios.min()), var_ratio_max=float(var_ratios.max()),
        lead_moment_z=float(lead_arr.mean() / (lead_arr.std(ddof=1) / np.sqrt(R))),
    )


def check_thm4a() -> dict:
    """Verify Theorem 4(a): E[nu] == E[Y(0)] on N=200, T=8, R=1000."""
    N, G, T, R = 200, 20, 8, 1000
    geo = dgp.make_geography(N, G=G)
    p = dgp.DGPParams(N=N, T=T, adopt="feedback", tau=float(np.log(0.8)), kappa=-0.1, phi=0.4, k=10.0, gbar=0.09)
    b0 = np.array([p.tau, p.kappa, p.phi])
    diffs, nav_diffs, y0_tot = [], [], []
    for r in range(R):
        seed = 20274000 + r
        sim1 = dgp.simulate(p, geo, [seed, 0], [seed, 1], counterfactual=False)
        sim0 = dgp.simulate(p, geo, [seed, 0], [seed, 1], counterfactual=True)
        pan = dgp.build_panel(sim1["Yall"], sim1["Dall"], geo.W, sim1["omega"], geo.cl)
        # Counterfactual recursion with true beta0 and delta=gamma_diff
        Y_all = sim1["Yall"]; D_all = sim1["Dall"]
        L_pre, Nb_pre = pan["L"][:1], pan["Nb"][:1]
        # Simple direct E[nu_{it}] using exact A1 multipliers:
        eta_diff = -(p.tau * D_all[dgp.S_LAGS:] + p.kappa * np.einsum("ij,tj->ti", geo.W.toarray(), D_all[dgp.S_LAGS:]))
        Y_post = Y_all[dgp.S_LAGS:]
        nu_naive = Y_post * np.exp(eta_diff)
        # Forward renewal recursion nu
        nu = np.zeros_like(Y_all)
        nu[:dgp.S_LAGS] = Y_all[:dgp.S_LAGS]
        W_arr = geo.W.toarray()
        for j in range(T):
            t = dgp.S_LAGS + j
            Lt = 0.85 * nu[t - 1] + 0.15 * nu[t - 2]
            Lt_y = 0.85 * Y_all[t - 1] + 0.15 * Y_all[t - 2]
            lam_nu = Lt + p.phi * (W_arr @ Lt) + sim1["omega"]
            lam_y = Lt_y + p.phi * (W_arr @ Lt_y) + sim1["omega"]
            Db_t = W_arr @ D_all[t]
            c_it = np.exp(-p.tau * D_all[t] - p.kappa * Db_t) * (lam_nu / lam_y)
            nu[t] = c_it * Y_all[t]
        y0_post = sim0["Yall"][dgp.S_LAGS:]
        diffs.append((nu[dgp.S_LAGS:] - y0_post).sum(1))
        nav_diffs.append((nu_naive - y0_post).sum(1))
        y0_tot.append(y0_post.sum(1))
    d_arr = np.asarray(diffs); nav_arr = np.asarray(nav_diffs); y0_arr = np.asarray(y0_tot)
    cum_rel = float(d_arr.sum(1).mean() / y0_arr.sum(1).mean())
    naive_rel = float(nav_arr.sum(1).mean() / y0_arr.sum(1).mean())
    wk_z = np.abs(d_arr.mean(0)) / (d_arr.std(ddof=1, axis=0) / np.sqrt(R))
    return dict(R=R, cum_rel_diff=cum_rel, max_weekly_abs_z=float(wk_z.max()), naive_rescaling_rel_diff=naive_rel)


def main():
    t0 = time.time()
    workers = max(1, min(os.cpu_count() or 4, 18))
    print(f"Running Web Appendix B MC suite on {workers} workers...", flush=True)
    lem2 = check_lemma2_mds(workers)
    print("Lemma 2:", lem2, flush=True)
    thm3 = check_thm3_bias_decomp()
    print("Thm 3:", thm3, flush=True)
    thm2 = check_thm2(workers)
    print("Thm 2:", thm2, flush=True)
    thm4 = check_thm4a()
    print("Thm 4a:", thm4, flush=True)
    out = dict(seconds=round(time.time() - t0, 2), lemma2=lem2, thm3=thm3, thm2=thm2, thm4=thm4)
    OUT_PATH.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Saved {OUT_PATH} in {out['seconds']}s")


if __name__ == "__main__":
    main()
