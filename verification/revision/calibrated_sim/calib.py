"""Application-calibrated DGP for the few-cluster simulation (review Comment 5).

Truth = the application's M0 estimates (tau, kappa, phi) with the actual geography (3,091 counties, 49 state
clusters, row-normalised adjacency W), the actual mask-order path D (exogenous in the simulation), the Ferretti
weekly generation-interval weights, omega = pop/1e5 and the observed two pre-sample weeks.
alpha_i and gamma_t are the two-way Poisson fixed effects of mu_it = exp(alpha_i + gamma_t) B_it with
B_it = exp(tau D + kappa Dbar) Lambda_it(phi) evaluated on the observed lags (iterated proportional fitting);
the NB2 shape k is the moment estimator 1/k = sum[(Y-mu)^2 - mu] / sum mu^2 on those fitted means.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import env  # noqa
import numpy as np
import application as A, qdgmm, dgp


def calibrate():
    d = A.load_data()
    w = A.gi_weekly(A.GI_SOURCES["ferretti"])["w"]
    panel, Yall = A.make_panel(d, w, d["W_adj"])
    m0 = qdgmm.fit(panel, qdgmm.spec_m0(panel))
    tau, kappa, phi = (m0["est"][n] for n in ("tau", "kappa", "phi"))
    lam = panel["L"] + phi * panel["Nb"] + panel["omega"]
    B = np.exp(tau * panel["D"] + kappa * panel["Db"]) * lam
    Y = panel["Y"]
    ea = np.ones(Y.shape[1]); eg = np.ones(Y.shape[0])
    for _ in range(500):
        eg_new = Y.sum(1) / (ea[None, :] * B).sum(1)
        ea_new = Y.sum(0) / (eg_new[:, None] * B).sum(0)
        done = np.max(np.abs(eg_new - eg) / eg) < 1e-12
        eg, ea = eg_new, ea_new
        if done:
            break
    mu = ea[None, :] * eg[:, None] * B
    inv_k = float(((Y - mu) ** 2 - mu).sum() / (mu ** 2).sum())
    k = 1.0 / inv_k if inv_k > 0 else np.inf
    return dict(d=d, w=np.asarray(w, float), panel=panel, Yall=Yall, tau=tau, kappa=kappa, phi=phi,
                ea=ea, eg=eg, k=k, D_all=d["D"].astype(float), omega=panel["omega"], W=d["W_adj"], cl=d["cl"])


def recalibrate_gamma(cal, n=20, iters=5, seed0=900000):
    """The PML gamma_t reproduce the one-step-ahead mean given observed lags; under simulated feedback the
    national path drifts (Jensen / compounding). Shift gamma_t so that the mean simulated national weekly
    total matches the observed one (geometric mean over n draws) (deterministic seeds; recorded in calibration.json)."""
    S = A.S
    obs = cal["Yall"][S:].sum(1)
    hist = []
    for it in range(iters):
        sims = np.array([simulate(cal, seed0 + r)[1][S:].sum(1) for r in range(n)])
        lr = np.log(sims).mean(0) - np.log(obs)      # geometric mean: stable under heavy-tailed feedback
        hist.append(float(np.max(np.abs(lr))))
        inc = lr - np.r_[0.0, lr[:-1]]
        cal["eg"] = cal["eg"] * np.exp(-inc)
    sims = np.array([simulate(cal, seed0 + r)[1][S:].sum(1) for r in range(n)])
    hist.append(float(np.max(np.abs(np.log(sims).mean(0) - np.log(obs)))))
    cal["recal_hist"] = hist
    return cal


SIGMA_SW = 0.23   # SD of log state x week shocks: var of log(observed/fitted state-week totals) net of NB noise


def state_week_sigma(cal):
    """sqrt(var log(Y_gt / mu_gt) - mean NB noise), over state-weeks with Y_gt > 0, at the PML fit."""
    p = cal["panel"]; Y = p["Y"]; G = int(p["G"]); cl = p["cl"]
    lam = p["L"] + cal["phi"] * p["Nb"] + p["omega"]
    mu = cal["ea"][None, :] * cal["eg"][:, None] * np.exp(cal["tau"] * p["D"] + cal["kappa"] * p["Db"]) * lam
    Yg = np.stack([np.bincount(cl, Y[t], G) for t in range(Y.shape[0])])
    mg = np.stack([np.bincount(cl, mu[t], G) for t in range(Y.shape[0])])
    m2 = np.stack([np.bincount(cl, mu[t] ** 2, G) for t in range(Y.shape[0])])
    ok = Yg > 0
    lr = np.log(Yg[ok] / mg[ok])
    noise = ((mg + m2 / cal["k"]) / mg ** 2)[ok]
    return float(np.sqrt(max(lr.var() - noise.mean(), 0.0))), int((~ok).sum())


RHO_TARGET = 0.21   # mean correlation over weeks of log(Y_gt/mu_gt) between adjacent states (non-adjacent: 0.03)


def state_graph(cal):
    cl = cal["cl"]; G = int(cl.max()) + 1
    import scipy.sparse as sps
    coo = sps.coo_matrix(cal["W"])
    S = np.zeros((G, G))
    m = cl[coo.row] != cl[coo.col]
    S[cl[coo.row[m]], cl[coo.col[m]]] = 1.0
    return np.maximum(S, S.T)


def adjacent_corr(cal):
    p = cal["panel"]; Y = p["Y"]; G = int(p["G"]); cl = p["cl"]
    lam = p["L"] + cal["phi"] * p["Nb"] + p["omega"]
    mu = cal["ea"][None, :] * cal["eg"][:, None] * np.exp(cal["tau"] * p["D"] + cal["kappa"] * p["Db"]) * lam
    Yg = np.stack([np.bincount(cl, Y[t], G) for t in range(Y.shape[0])])
    mg = np.stack([np.bincount(cl, mu[t], G) for t in range(Y.shape[0])])
    lr = np.log(np.maximum(Yg, 0.5) / mg)
    S = state_graph(cal)
    adj = np.argwhere(np.triu(S) > 0)
    non = np.argwhere(np.triu(1 - S - np.eye(G)) > 0)
    c = lambda prs: float(np.mean([np.corrcoef(lr[:, a], lr[:, b])[0, 1] for a, b in prs]))
    return c(adj), c(non)


def shock_matrix(S, rho):
    Ws = S / np.maximum(S.sum(1, keepdims=True), 1)
    M = np.eye(S.shape[0]) + rho * Ws
    M = M / np.sqrt((M ** 2).sum(1, keepdims=True))
    C = M @ M.T
    iu = np.argwhere(np.triu(S) > 0)
    return M, float(np.mean([C[a, b] for a, b in iu]))


def calibrated(sigma=None, rho_target=None):
    cal = calibrate()
    est, nzero = state_week_sigma(cal)
    cal["sigma_sw_est"], cal["zero_state_weeks"] = est, nzero
    cal["sigma_sw"] = SIGMA_SW if sigma is None else sigma
    cal["adj_corr_data"], cal["nonadj_corr_data"] = adjacent_corr(cal)
    tgt = RHO_TARGET if rho_target is None else rho_target
    S = state_graph(cal)
    from scipy.optimize import brentq
    rho = brentq(lambda r: shock_matrix(S, r)[1] - tgt, 0.0, 50.0) if tgt > 0 else 0.0
    cal["shock_M"], cal["shock_adj_corr"] = shock_matrix(S, rho)
    cal["shock_rho"] = rho
    return recalibrate_gamma(cal, n=40)


def simulate(cal, seed):
    """One panel from the calibrated DGP; returns the estimation panel (same builder as the application)."""
    rng = np.random.default_rng(seed)
    S = A.S
    Yobs = cal["Yall"]
    T = Yobs.shape[0] - S
    Y = np.zeros_like(Yobs, dtype=float)
    Y[:S] = Yobs[:S]
    D = cal["D_all"]; W = cal["W"]; w = cal["w"]; om = cal["omega"]
    sig = cal.get("sigma_sw", 0.0)
    G = int(cal["cl"].max()) + 1
    for j in range(T):
        t = S + j
        L = w[0] * Y[t - 1] + w[1] * Y[t - 2]
        lam = L + cal["phi"] * (W @ L) + om
        mu = cal["ea"] * cal["eg"][j] * np.exp(cal["tau"] * D[t] + cal["kappa"] * (W @ D[t])) * lam
        if sig > 0:      # mean-one state x week shocks (common within cluster, spatially correlated across
            # adjacent states, independent of the past): eps = sig * M e, M = rows of (I + rho Ws) scaled to unit variance
            eps = sig * (cal["shock_M"] @ rng.standard_normal(G))
            mu = mu * np.exp(eps - 0.5 * sig ** 2)[cal["cl"]]
        Y[t] = dgp._draw_nb(rng, mu, cal["k"])
    panel = dgp.build_panel(Y, D, W, om, cal["cl"], w=w, S=S)
    return qdgmm.attach_presample(panel, D[:S]), Y


if __name__ == "__main__":
    import json, time
    t0 = time.time()
    cal = calibrate()
    gamma_pml = np.log(cal["eg"]).copy()
    cal = recalibrate_gamma(cal)
    obs = cal["Yall"][A.S:].sum(1)
    sims = np.array([simulate(cal, 1000 + r)[1][A.S:].sum(1) for r in range(20)])
    out = dict(tau=cal["tau"], kappa=cal["kappa"], phi=cal["phi"], k_nb2=cal["k"], w=cal["w"],
               gamma=np.log(cal["eg"]), gamma_pml=gamma_pml, recal_max_abs_log_ratio=cal["recal_hist"], n_units=int(cal["Yall"].shape[1]), G=int(cal["panel"]["G"]),
               n_zero_units=int((cal["Yall"][A.S:].sum(0) == 0).sum()),
               switching_states=int(len(np.unique(cal["cl"][cal["D_all"].any(0)]))),
               observed_weekly_total=obs, simulated_weekly_total_median=np.median(sims, 0),
               simulated_weekly_total_p05=np.quantile(sims, 0.05, 0), simulated_weekly_total_p95=np.quantile(sims, 0.95, 0),
               seconds=time.time() - t0)
    Path(__file__).with_name("calibration.json").write_text(json.dumps(env.clean(out), indent=1))
    print(json.dumps(env.clean(out), indent=1))
