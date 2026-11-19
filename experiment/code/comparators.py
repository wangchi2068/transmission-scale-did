"""Comparator estimators (theory/simulation_spec.md, section M).

M1     Poisson fixed effects = conditional multinomial likelihood with offset log Lambda(phi), unit and week
       effects; regressors D, Dbar and phi (through the offset). Exact Newton; cluster-robust sandwich.
M1-BC  Split-panel jackknife (Dhaene & Jochmans 2015): 2*theta_hat - mean(half-panel estimates); SE of M1.
M2     Two-step log R-hat DiD: Cori-type weekly posterior mean R_hat = (a + Y)/(1/b + L), gamma prior
       (shape a = 1, scale b = 5), local pressure only; TWFE OLS of log R_hat on D; CR1 cluster SE.
M3     TWFE OLS of log(Y + 1) on D; CR1 cluster SE.
M4     Pooled renewal Poisson QMLE without unit effects (week effects, tau, kappa, phi); cluster sandwich.
Oracle Poisson QMLE with the true alpha_i as offset (infeasible benchmark); cluster sandwich.
All confidence intervals use t_{G-1}; cluster sandwiches carry the G/(G-1) factor.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

NAMES3 = ("tau", "kappa", "phi")


def _nan_result(names, reason):
    return dict(ok=False, est={n: float("nan") for n in names}, se={n: float("nan") for n in names}, reason=reason)


def _cluster_sandwich(H, scores_unit, cl, G):
    """V = G/(G-1) H^{-1} B H^{-1}, B = sum_g S_g S_g'. H is the (negative-definite) Hessian sum."""
    Sg = np.vstack([np.bincount(cl, s, G) for s in scores_unit])        # (p, G)
    B = Sg @ Sg.T
    Hi = np.linalg.inv(H)
    return Hi @ B @ Hi.T * G / (G - 1)


def _newton_max(fgh, th0, max_iter=200, tol=1e-12):
    """Damped modified-Newton maximiser. fgh(th) -> (ll, grad, Hessian) or None (infeasible).
    The direction solves (-H + mu*diag) d = g with mu raised until the matrix is positive definite;
    Armijo backtracking keeps the iterate feasible and ascending. Convergence: Newton decrement
    g'd < tol*|ll| + 1e-10, or the step falls below 1e-10. Returns (th, ok, iters, reason)."""
    th = np.asarray(th0, float)
    cur = fgh(th)
    if cur is None:
        return th, False, 0, "infeasible start"
    for it in range(1, max_iter + 1):
        ll, g, H = cur
        A = -H
        dA = np.abs(np.diag(A)) + 1e-12
        mu = 0.0
        for _ in range(60):
            try:
                Lc = np.linalg.cholesky(A + mu * np.diag(dA))
                break
            except np.linalg.LinAlgError:
                mu = 1e-8 if mu == 0.0 else mu * 10.0
        else:
            return th, False, it, "no positive-definite modification"
        d = np.linalg.solve(Lc.T, np.linalg.solve(Lc, g))
        dec = float(g @ d)
        if not np.isfinite(dec):
            return th, False, it, "non-finite step"
        if dec < tol * abs(ll) + 1e-10 and mu == 0.0:
            return th, True, it, ""
        t = 1.0
        acc = False
        for _ in range(40):
            cand = th + t * d
            new = fgh(cand)
            if new is not None and new[0] >= ll + 1e-4 * t * dec - 1e-13 * abs(ll):
                acc = True
                break
            t *= 0.5
        if not acc:
            # cannot ascend further at floating-point precision: converged if the decrement is tiny
            ok = dec < 1e-8 * abs(ll) + 1e-6
            return th, ok, it, "" if ok else "line search failed"
        th, cur = cand, new
        if np.max(np.abs(t * d)) < 1e-10:
            return th, True, it, ""
    return th, False, max_iter, "max_iter"


# ------------------------------------------------------------------------------------------------
# M1: Poisson FE (conditional multinomial likelihood)
# ------------------------------------------------------------------------------------------------
def _fe_design(panel, t_idx, with_kappa=True):
    """Return regressor stack X (T', N, K) for weeks t_idx: week dummies 2..T' (first dropped), D[, Dbar]."""
    D, Db = panel["D"][t_idx], panel["Db"][t_idx]
    Tn, N = D.shape
    K = (Tn - 1) + 1 + int(with_kappa)
    X = np.zeros((Tn, N, K))
    for t in range(1, Tn):
        X[t, :, t - 1] = 1.0
    X[:, :, Tn - 1] = D
    if with_kappa:
        X[:, :, Tn] = Db
    return X


def poisson_fe(panel: dict, t_idx=None, phi_start: float = 0.3, max_iter: int = 200,
               with_se: bool = True, with_kappa: bool = True, phi_fixed: Optional[float] = None) -> dict:
    """Conditional (multinomial) Poisson FE with offset log Lambda(phi); returns tau, kappa, phi.
    with_kappa=False drops Dbar (kappa reported as 0); phi_fixed fixes phi (reported, SE NaN)."""
    T = panel["Y"].shape[0]
    t_idx = np.arange(T) if t_idx is None else np.asarray(t_idx)
    Y = panel["Y"][t_idx]
    L, Nb, om = panel["L"][t_idx], panel["Nb"][t_idx], panel["omega"]
    if np.ndim(om) == 2:
        om = om[t_idx]
    X = _fe_design(panel, t_idx, with_kappa)
    Tn, N, K = X.shape
    n = Y.sum(0)
    keep = n > 0
    Y, L, Nb, X, n = Y[:, keep], L[:, keep], Nb[:, keep], X[:, keep], n[keep]
    om = om[..., keep] if np.ndim(om) else om                     # omega is per unit (N,) or (T, N)
    cl = panel["cl"][keep]
    G = int(panel["G"])
    # identification: D must vary within some unit
    Dk = X[:, :, Tn - 1]
    if not np.any(Dk.max(0) != Dk.min(0)):
        return _nan_result(NAMES3, "no within-unit variation in D")
    phi_free = phi_fixed is None
    P = K + int(phi_free)
    th = np.zeros(P)
    if phi_free:
        th[-1] = phi_start
    ok = False

    def parts(th):
        b = th[:K]
        phi = th[K] if phi_free else phi_fixed
        lam = L + phi * Nb + om
        if np.any(lam <= 0):
            return None
        eta = X @ b + np.log(lam)
        eta = eta - eta.max(0, keepdims=True)
        e = np.exp(eta)
        p = e / e.sum(0, keepdims=True)
        r = Nb / lam
        g = np.concatenate([X, r[:, :, None]], 2) if phi_free else X   # d eta / d theta  (T,N,P)
        resid = Y - n[None, :] * p
        ll = float((Y * np.log(np.clip(p, 1e-300, None))).sum())
        score_u = np.einsum("tn,tnp->pn", resid, g)                    # per unit
        pg = np.einsum("tn,tnp->np", p, g)                             # (N,P)
        g2 = g.reshape(-1, g.shape[2])                                  # (T*N, P); avoids a T*N*P*P temporary
        wts = (n[None, :] * p).reshape(-1)
        Hm = -(g2 * wts[:, None]).T @ g2 + (pg * n[:, None]).T @ pg
        if phi_free:
            Hm[K, K] += float((resid * (-(r ** 2))).sum())
        return ll, score_u, Hm

    def fgh(t_):
        pr_ = parts(t_)
        return None if pr_ is None else (pr_[0], pr_[1].sum(1), pr_[2])
    th, ok, it, why = _newton_max(fgh, th, max_iter=max_iter)
    pr = parts(th)
    if pr is None:
        return _nan_result(NAMES3, "Lambda<=0")
    ll, su, Hm = pr
    idx = {"tau": Tn - 1, "kappa": Tn if with_kappa else None, "phi": K if phi_free else None}
    est = {"tau": float(th[Tn - 1]), "kappa": float(th[Tn]) if with_kappa else 0.0,
           "phi": float(th[K]) if phi_free else float(phi_fixed)}
    nan_se = {n_: float("nan") for n_ in NAMES3}
    if not with_se:
        return dict(ok=ok, est=est, se=nan_se, iters=it, reason=why)
    try:
        V = _cluster_sandwich(Hm, su, cl, G)
    except np.linalg.LinAlgError:
        return dict(ok=False, est=est, se=nan_se, reason="singular Hessian at end")
    se = {nm: (float(np.sqrt(max(V[i, i], 0.0))) if i is not None else float("nan")) for nm, i in idx.items()}
    Vb = V[np.ix_([Tn - 1, Tn], [Tn - 1, Tn])] if with_kappa else None
    return dict(ok=ok, est=est, se=se, Vb=Vb, iters=it, reason=why)


def poisson_fe_bc(panel: dict, full: dict) -> dict:
    """Split-panel jackknife: halves are the first and last floor(T/2)/ceil(T/2) weeks."""
    T = panel["Y"].shape[0]
    h = T // 2
    if not full.get("ok"):
        return _nan_result(NAMES3, "M1 failed")
    halves = [poisson_fe(panel, np.arange(0, h), with_se=False), poisson_fe(panel, np.arange(h, T), with_se=False)]
    if not all(x.get("ok") for x in halves):
        return _nan_result(NAMES3, "half-panel fit failed: " + "; ".join(x.get("reason", "nonconv") for x in halves
                                                                        if not x.get("ok")))
    est = {nm: 2 * full["est"][nm] - 0.5 * (halves[0]["est"][nm] + halves[1]["est"][nm]) for nm in NAMES3}
    return dict(ok=True, est=est, se=dict(full["se"]))


# ------------------------------------------------------------------------------------------------
# M2 / M3: TWFE OLS with cluster-robust SE
# ------------------------------------------------------------------------------------------------
def twfe_ols(y: np.ndarray, D: np.ndarray, cl: np.ndarray, G: int) -> dict:
    """Balanced-panel two-way FE slope on D (T x N arrays); CR1 cluster SE with factor G/(G-1)."""
    def dm(a):
        return a - a.mean(0, keepdims=True) - a.mean(1, keepdims=True) + a.mean()
    Dd = dm(D)
    den = float((Dd ** 2).sum())
    if den <= 1e-12:
        return _nan_result(("tau",), "no variation in D")
    yd = dm(y)
    tau = float((Dd * yd).sum() / den)
    e = yd - tau * Dd
    su = (Dd * e).sum(0)
    Sg = np.bincount(cl, su, G)
    var = float((Sg ** 2).sum() / den ** 2 * G / (G - 1))
    return dict(ok=True, est={"tau": tau}, se={"tau": float(np.sqrt(var))})


def cori_log_r(panel: dict, a: float = 1.0, b: float = 5.0) -> np.ndarray:
    """Weekly Cori-type posterior mean of R with Gamma(shape a, scale b) prior and local pressure L."""
    return np.log((a + panel["Y"]) / (1.0 / b + panel["L"]))


def m2_logr_did(panel: dict) -> dict:
    return twfe_ols(cori_log_r(panel), panel["D"], panel["cl"], int(panel["G"]))


def m3_logy_twfe(panel: dict) -> dict:
    return twfe_ols(np.log(panel["Y"] + 1.0), panel["D"], panel["cl"], int(panel["G"]))


# ------------------------------------------------------------------------------------------------
# M4 / Oracle: Poisson QMLE with explicit week effects
# ------------------------------------------------------------------------------------------------
def poisson_qmle(panel: dict, offset: Optional[np.ndarray] = None, max_iter: int = 200,
                 phi_start: float = 0.3) -> dict:
    """mu = exp(offset_i + gamma_t + tau D + kappa Dbar) * Lambda(phi); week effects gamma_1..gamma_T
    (no common intercept). offset = None gives the pooled model (M4); offset = alpha_i gives the Oracle."""
    Y, D, Db = panel["Y"], panel["D"], panel["Db"]
    L, Nb, om = panel["L"], panel["Nb"], panel["omega"]
    T, N = Y.shape
    off = np.zeros(N) if offset is None else np.asarray(offset, float)
    cl, G = panel["cl"], int(panel["G"])
    P = T + 3
    # starting values: week log-ratios of totals
    lam0 = L + phi_start * Nb + om
    g0 = np.log(np.clip(Y.sum(1), 0.5, None) / (np.exp(off)[None, :] * lam0).sum(1))
    th = np.concatenate([g0, [0.0, 0.0, phi_start]])
    ok = False

    def parts(th):
        gam, tau, kap, phi = th[:T], th[T], th[T + 1], th[T + 2]
        lam = L + phi * Nb + om
        if np.any(lam <= 0):
            return None
        eta = off[None, :] + gam[:, None] + tau * D + kap * Db
        mu = np.exp(np.clip(eta, -50, 50)) * lam
        r = Nb / lam
        resid = Y - mu
        ll = float((Y * np.log(np.clip(mu, 1e-300, None)) - mu).sum())
        return mu, r, resid, ll

    def fgh(t_):
        pr_ = parts(t_)
        if pr_ is None:
            return None
        mu_, r_, res_, ll_ = pr_
        g_ = np.concatenate([res_.sum(1), [(res_ * D).sum(), (res_ * Db).sum(), (res_ * r_).sum()]])
        return ll_, g_, _qmle_hessian(mu_, r_, res_, D, Db)
    th, ok, it, why = _newton_max(fgh, th, max_iter=max_iter)
    pr = parts(th)
    if pr is None:
        return _nan_result(NAMES3, "Lambda<=0")
    mu, r, resid, ll = pr
    Hm = _qmle_hessian(mu, r, resid, D, Db)
    su = np.concatenate([resid, (resid * D).sum(0)[None], (resid * Db).sum(0)[None], (resid * r).sum(0)[None]], 0)
    # note: week-effect scores per unit are resid[t, i] (row t), parameter scores are sums over t
    est = {"tau": float(th[T]), "kappa": float(th[T + 1]), "phi": float(th[T + 2])}
    try:
        V = _cluster_sandwich(Hm, su, cl, G)
    except np.linalg.LinAlgError:
        return dict(ok=False, est=est, se={n_: float("nan") for n_ in NAMES3}, reason="singular Hessian at end")
    se = {nm: float(np.sqrt(max(V[T + j, T + j], 0.0))) for j, nm in enumerate(NAMES3)}
    return dict(ok=ok, est=est, se=se, Vb=V[T:T + 3, T:T + 3], iters=it, reason=why)


def _qmle_hessian(mu, r, resid, D, Db):
    T = mu.shape[0]
    P = T + 3
    H = np.zeros((P, P))
    H[np.arange(T), np.arange(T)] = -mu.sum(1)
    cols = [D, Db, r]
    for j, x in enumerate(cols):
        H[:T, T + j] = H[T + j, :T] = -(mu * x).sum(1)
        for k, z in enumerate(cols):
            H[T + j, T + k] = -(mu * x * z).sum()
    H[T + 2, T + 2] += float((resid * (-(r ** 2))).sum())
    return H
