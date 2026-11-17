"""Delete-one-state jackknife, score wild bootstrap-t, and pairs moment cluster bootstrap-t.

Why three few-cluster diagnostics are reported
----------------------------------------------
1. Delete-one-state jackknife: removes each state in turn, rebuilds row-normalised W and refits
   M0. Because Texas contributes 254 of 3091 counties and removing a state alters its neighbours'
   pressure, the jackknife measures sensitivity to dropping whole clusters rather than resampling
   variability.
2. Score wild cluster bootstrap-t (Kline & Santos 2012): perturbs the cluster influence terms
   phi_g = (-Gamma^{-1} g_g)_beta with iid Mammen (1993) two-point weights (and Webb 6-point as a
   check), holding Gamma fixed and studentising each draw by its recentered bootstrap variance so the
   denominator fluctuates across draws (unlike uncentered Rademacher weights, where w_g^2 == 1).
3. Nonlinear pairs cluster bootstrap-t on cluster moments: draws multinomial cluster multiplicities
   m ~ Multinomial(G; 1/G, ..., 1/G), solves sum_g m_g g_g(theta^*) = 0 keeping the spatial graph W
   and lagged pressures intact, and recomputes both the estimate and the cluster sandwich V^*_b in
   every bootstrap draw, capturing both score and Jacobian variability across unequal states.

Run: cd experiment/code && python -I ../../verification/application/wild_bootstrap.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import optimize

HERE = Path(__file__).resolve().parent
CODE = HERE.parents[1] / "experiment" / "code"
sys.path.insert(0, str(CODE))

import application as A  # noqa: E402
import qdgmm  # noqa: E402


def delete_one(d, w, names, est):
    """Refit with each state removed; returns (deviations, sizes) both indexed by state."""
    theta = np.array([est[n] for n in names])
    states = np.unique(d["state"])
    dev = np.full((len(names), states.size), np.nan)
    sizes = np.array([(d["state"] == s).sum() for s in states])
    for k, s in enumerate(states):
        keep = np.where(d["state"] != s)[0]
        pnl, _ = A.make_panel(d, w, d["W_adj"], units=keep)
        f = qdgmm.fit(pnl, qdgmm.spec_m0(pnl), starts=[theta])
        if f.get("ok"):
            dev[:, k] = np.array([f["est"][n] for n in names]) - theta
    return dev, sizes


def score_wild_bootstrap(phi_g: np.ndarray, B: int, seed: int, weight: str = "mammen") -> dict:
    """Kline & Santos (2012) score wild cluster bootstrap-t with recentered bootstrap sandwich SE."""
    rng = np.random.default_rng(seed)
    G = phi_g.shape[1]
    p1 = (np.sqrt(5.0) + 1.0) / (2.0 * np.sqrt(5.0))
    v1, v2 = (1.0 - np.sqrt(5.0)) / 2.0, (1.0 + np.sqrt(5.0)) / 2.0
    webb = np.array([-np.sqrt(1.5), -1.0, -np.sqrt(0.5), np.sqrt(0.5), 1.0, np.sqrt(1.5)])

    t_mat = np.zeros((B, phi_g.shape[0] + 1))
    for b in range(B):
        if weight == "mammen":
            wt = np.where(rng.random(G) < p1, v1, v2)
        elif weight == "webb6":
            wt = rng.choice(webb, size=G)
        else:
            raise ValueError(weight)
        ph_b = phi_g * wt[None, :]
        d_b = ph_b.mean(1)
        res_b = ph_b - d_b[:, None]
        Vb_b = (res_b @ res_b.T) / (G * (G - 1))
        se_b = np.sqrt(np.maximum(np.diag(Vb_b), 1e-300))
        t_mat[b, :-1] = d_b / se_b
        se_tot_b = np.sqrt(max(Vb_b[0, 0] + Vb_b[1, 1] + 2.0 * Vb_b[0, 1], 1e-300))
        t_mat[b, -1] = (d_b[0] + d_b[1]) / se_tot_b
    q = np.quantile(np.abs(t_mat), [0.9, 0.95, 0.975], axis=0)
    return dict(q90=q[0].tolist(), q95=q[1].tolist(), q975=q[2].tolist(), B=B, weight=weight)


def pairs_moment_bootstrap(mod: qdgmm.QDModel, beta_hat: np.ndarray, B: int, seed: int) -> dict:
    """Pairs cluster bootstrap-t solving sum_g m_g g_g(theta^*) = 0 with m ~ Multinomial(G; 1/G)."""
    rng = np.random.default_rng(seed)
    G = mod.G
    Tm1 = mod.T - 1
    ts = []
    fails = 0
    for _ in range(B):
        m_g = rng.multinomial(G, np.full(G, 1.0 / G)).astype(float)
        w_unit = m_g[mod.cl]

        def prof_w(beta):
            b_, ph_ = mod._split(beta)
            lam = mod.L + ph_ * mod.Nb + mod.om
            if np.any(lam <= 0) or not np.all(np.isfinite(beta)):
                return np.full(mod.np_, 1e6 * mod.fscale), np.eye(mod.np_)
            eta = np.zeros_like(mod.Y)
            for bk, x in zip(b_, mod.X):
                if bk != 0.0:
                    eta = eta + bk * x
            q_ = np.exp(np.clip(eta[:-1] - eta[1:], -50, 50)) * lam[:-1] / lam[1:]
            Yq = mod.Y[1:] * q_
            den = (w_unit[None, :] * mod.s * mod.Y[:-1]).sum(1)
            num = (w_unit[None, :] * mod.s * Yq).sum(1)
            if np.any(den <= 0) or np.any(num <= 0):
                return np.full(mod.np_, 1e6 * mod.fscale), np.eye(mod.np_)
            Ed = num / den
            Yr = Yq / Ed[:, None]
            v = Yr - mod.Y[:-1]
            dx = list(mod.dX)
            if mod.spec.phi_free:
                r = mod.Nb / lam
                dx.append(r[:-1] - r[1:])
            F = np.array([(w_unit[None, :] * z * v).sum() for z in mod.Zpar]) / mod.N
            sYr = w_unit[None, :] * mod.s * Yr
            a = sYr.sum(1)
            bb = np.stack([(sYr * d_).sum(1) for d_ in dx], 1)
            cm = np.stack([(w_unit[None, :] * z * Yr).sum(1) for z in mod.Zpar], 0)
            H = np.array([[(w_unit[None, :] * z * Yr * d_).sum() for d_ in dx] for z in mod.Zpar])
            Jp = (H - cm @ (bb / a[:, None])) / mod.N
            return F, Jp

        sol = optimize.root(prof_w, beta_hat, jac=True, method="hybr", options={"xtol": 1e-10, "maxfev": 400})
        F_end, _ = prof_w(sol.x)
        if not sol.success or np.max(np.abs(F_end)) > 1e-7 * mod.fscale or (mod.spec.phi_free and sol.x[mod.K] < 0):
            fails += 1
            continue
        bx = sol.x
        b_, ph_ = mod._split(bx)
        lam = mod.L + ph_ * mod.Nb + mod.om
        eta = np.zeros_like(mod.Y)
        for bk, x in zip(b_, mod.X):
            if bk != 0.0:
                eta = eta + bk * x
        q_ = np.exp(np.clip(eta[:-1] - eta[1:], -50, 50)) * lam[:-1] / lam[1:]
        Yq = mod.Y[1:] * q_
        Ed = (w_unit[None, :] * mod.s * Yq).sum(1) / (w_unit[None, :] * mod.s * mod.Y[:-1]).sum(1)
        Yr = Yq / Ed[:, None]
        v = Yr - mod.Y[:-1]
        dx = list(mod.dX) + [(mod.Nb / lam)[:-1] - (mod.Nb / lam)[1:]]
        gi = np.vstack([mod.s * v, np.stack([(z * v).sum(0) for z in mod.Zpar], 0)])
        gg = mod.cluster_sum(gi)
        sYr = w_unit[None, :] * mod.s * Yr
        a = sYr.sum(1)
        bb = np.stack([(sYr * d_).sum(1) for d_ in dx], 1)
        cm = np.stack([(w_unit[None, :] * z * Yr).sum(1) for z in mod.Zpar], 0)
        H = np.array([[(w_unit[None, :] * z * Yr * d_).sum() for d_ in dx] for z in mod.Zpar])
        J = np.zeros((Tm1 + mod.np_, Tm1 + mod.np_))
        J[np.arange(Tm1), np.arange(Tm1)] = -a
        J[:Tm1, Tm1:] = bb
        J[Tm1:, :Tm1] = -cm
        J[Tm1:, Tm1:] = H
        Gam_b = J / G
        Om_b = (gg * m_g[None, :]) @ gg.T / G
        try:
            Gi_b = np.linalg.inv(Gam_b)
        except np.linalg.LinAlgError:
            fails += 1
            continue
        Vb_b = (Gi_b @ Om_b @ Gi_b.T / G * G / (G - 1))[Tm1:, Tm1:]
        if not np.all(np.isfinite(Vb_b)) or np.any(np.diag(Vb_b) <= 0):
            fails += 1
            continue
        se_b = np.sqrt(np.diag(Vb_b))
        se_tot_b = np.sqrt(max(Vb_b[0, 0] + Vb_b[1, 1] + 2.0 * Vb_b[0, 1], 1e-300))
        db = bx - beta_hat
        ts.append([db[0] / se_b[0], db[1] / se_b[1], db[2] / se_b[2], (db[0] + db[1]) / se_tot_b])

    ts = np.asarray(ts, float)
    q = np.quantile(np.abs(ts), [0.9, 0.95, 0.975], axis=0)
    return dict(q90=q[0].tolist(), q95=q[1].tolist(), q975=q[2].tolist(), ok=int(len(ts)), fails=int(fails), B=B)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=999)
    a = ap.parse_args()
    d = A.load_data()
    A.verify(d)
    gi = A.gi_weekly(A.GI_SOURCES["ferretti"])
    w = gi["w"] if isinstance(gi, dict) else gi
    panel, _ = A.make_panel(d, w, d["W_adj"])
    spec = qdgmm.spec_m0(panel)
    names = list(spec.names)
    mod = qdgmm.QDModel(panel, spec)
    m0 = qdgmm.fit(panel, spec)
    est = np.array([m0["est"][n] for n in names])
    inf = mod.inference(est)
    Vb = inf["Vb"]
    se_s = np.sqrt(np.diag(Vb))
    tot = float(est[0] + est[1])
    se_t = float(np.sqrt(Vb[0, 0] + Vb[1, 1] + 2.0 * Vb[0, 1]))

    # 1. Delete-one-state jackknife
    dev, sizes = delete_one(d, w, names, m0["est"])
    ok = np.all(np.isfinite(dev), axis=0)
    X, sz = dev[:, ok], sizes[ok]
    G = X.shape[1]
    se_j = np.sqrt((G - 1) / G * ((X - X.mean(1, keepdims=True)) ** 2).sum(1))

    # 2. Score wild cluster bootstrap-t (Kline & Santos 2012) on cluster influence vectors
    Tm1 = mod.T - 1
    phi_g = (-np.linalg.inv(inf["Gam"]) @ inf["gg"])[Tm1:, :]
    wild_mammen = score_wild_bootstrap(phi_g, a.B, 20261124, "mammen")
    wild_webb6 = score_wild_bootstrap(phi_g, a.B, 20261125, "webb6")

    # 3. Nonlinear pairs cluster bootstrap-t on cluster moments
    pairs = pairs_moment_bootstrap(mod, est, a.B, 20261126)

    q_w = np.array(wild_mammen["q95"])
    q_p = np.array(pairs["q95"])
    res = dict(
        names=names,
        est=est.tolist(),
        se_sandwich=se_s.tolist(),
        se_jack=se_j.tolist(),
        largest_state=int(sz.max()),
        n_states=int(G),
        B=a.B,
        score_wild_mammen=dict(
            wild_mammen,
            lo=(est - q_w[:3] * se_s).tolist(),
            hi=(est + q_w[:3] * se_s).tolist(),
            tau_plus_kappa=dict(est=tot, se=se_t, q95=float(q_w[3]),
                                lo=float(tot - q_w[3] * se_t), hi=float(tot + q_w[3] * se_t)),
        ),
        score_wild_webb6=wild_webb6,
        pairs_moment_bootstrap=dict(
            pairs,
            lo=(est - q_p[:3] * se_s).tolist(),
            hi=(est + q_p[:3] * se_s).tolist(),
            tau_plus_kappa=dict(est=tot, se=se_t, q95=float(q_p[3]),
                                lo=float(tot - q_p[3] * se_t), hi=float(tot + q_p[3] * se_t)),
        ),
        # keep top-level tau_plus_kappa summary using the nonlinear pairs bootstrap-t
        tau_plus_kappa=dict(
            est=tot,
            se=se_t,
            q95_wild_mammen=float(q_w[3]),
            ci_wild_mammen=[float(tot - q_w[3] * se_t), float(tot + q_w[3] * se_t)],
            q95_pairs=float(q_p[3]),
            ci_pairs=[float(tot - q_p[3] * se_t), float(tot + q_p[3] * se_t)],
            ci_t48=[float(tot - 2.010634757624232 * se_t), float(tot + 2.010634757624232 * se_t)],
        ),
    )

    print("states", G, "largest", int(sz.max()))
    print("sandwich se", se_s.round(4), "| jackknife se", se_j.round(4))
    for i, n in enumerate(names):
        print("  %-6s est %8.4f | wild(Mammen) q95 %5.3f [%7.3f,%7.3f] | pairs q95 %5.3f [%7.3f,%7.3f]"
              % (n, est[i], q_w[i], est[i] - q_w[i] * se_s[i], est[i] + q_w[i] * se_s[i],
                 q_p[i], est[i] - q_p[i] * se_s[i], est[i] + q_p[i] * se_s[i]))
    print("  tau+k  est %8.4f | wild(Mammen) q95 %5.3f [%7.3f,%7.3f] | pairs q95 %5.3f [%7.3f,%7.3f] | t48 [%7.3f,%7.3f]"
          % (tot, q_w[3], tot - q_w[3] * se_t, tot + q_w[3] * se_t,
             q_p[3], tot - q_p[3] * se_t, tot + q_p[3] * se_t,
             tot - 2.010634757624232 * se_t, tot + 2.010634757624232 * se_t))
    (HERE / "wild_bootstrap.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print("written", HERE / "wild_bootstrap.json")


if __name__ == "__main__":
    main()
