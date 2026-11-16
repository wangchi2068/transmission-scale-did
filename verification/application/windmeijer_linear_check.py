"""Self-contained check of the Windmeijer (2005) correction on a small linear GMM.

The application panel cannot be bootstrapped cheaply (the network and the cluster index must be
rebuilt on every resample), so the correction is validated where it can be: on a linear
over-identified GMM in which the first-step weight matrix is estimated, the second-step estimate is
computed from the same sample, and the finite-sample bias of the two-step standard error is known to
be negative. The script

1. simulates a linear homoskedastic model with q > p and iid clusters;
2. computes the two-step GMM estimate and three variance estimators: the plain sandwich, the
   sandwich times G/(G-1), and the Windmeijer correction;
3. repeats over many replications and compares each estimator's mean with the Monte Carlo standard
   deviation of the estimate.

Run: python -I verification/application/windmeijer_linear_check.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R = 600
G = 50
N_PER_G = 8
P = 2
Q = 4
SEED = 20261116


def one_rep(rng):
    n = G * N_PER_G
    cl = np.repeat(np.arange(G), N_PER_G)
    # cluster-correlated regressors and errors
    zg = rng.normal(size=(G, Q))
    Z = np.repeat(zg, N_PER_G, axis=0) + 0.3 * rng.normal(size=(n, Q))
    X = Z[:, :P].copy()
    u_g = rng.normal(size=G)
    u = np.repeat(u_g, N_PER_G) + rng.normal(size=n)
    beta0 = np.array([1.0, -0.5])
    y = X @ beta0 + u
    Zc = Z - Z.reshape(G, N_PER_G, Q).mean(1).repeat(N_PER_G, axis=0)

    def gsum(b):
        r = y - X @ b
        gi = np.einsum("ij,i->j", Zc, r)
        return np.vstack([np.bincount(cl, gi * (cl == gg), minlength=G) for gg in range(1)])[0]

    def gmat(b):
        r = y - X @ b
        return np.einsum("ij,i->j", Zc, r)

    # first step: W = I
    g1 = gmat(beta0)
    # one-step estimate with W = I
    A1 = Zc.T @ X
    b1 = np.linalg.solve(A1.T @ A1, A1.T @ (Zc.T @ y))
    # first-step weight matrix from the one-step residual
    r1 = y - X @ b1
    gi1 = np.einsum("ij,i->j", Zc, r1)
    Om1 = np.zeros((Q, Q))
    for gg in range(G):
        sel = cl == gg
        v = np.einsum("ij,i->j", Zc[sel], r1[sel])
        Om1 += np.outer(v, v)
    Om1 /= G
    W = np.linalg.inv(Om1)
    A = Zc.T @ X / G
    b2 = np.linalg.solve(A.T @ W @ A, A.T @ W @ (Zc.T @ y / G))

    r2 = y - X @ b2
    V = np.zeros((Q, Q))
    for gg in range(G):
        sel = cl == gg
        v = np.einsum("ij,i->j", Zc[sel], r2[sel])
        V += np.outer(v, v)
    Om2 = V / G
    Ainv = np.linalg.inv(A.T @ W @ A)
    base = Ainv @ (A.T @ W @ Om2 @ W @ A) @ Ainv / G

    # Windmeijer correction: G_n depends on W_hat through the first-step residuals
    # dG_j = d/dbeta_j [ sum_i z_i (y_i - x_i'beta) ] with W_hat held fixed is just -A per column,
    # the correction term is -Ainv * sum_j tr(...) computed numerically on the first-step weight.
    # We approximate d(W_hat)/d beta at beta2 with central differences on the one-step residuals.
    h = 1e-5
    dOm = np.zeros((Q, Q, P))
    for j in range(P):
        bp, bm = b2.copy(), b2.copy()
        bp[j] += h
        bm[j] -= h
        Op = _om(Zc, X, y, bp, cl, G)
        Om_ = _om(Zc, X, y, bm, cl, G)
        dOm[:, :, j] = (Op - Om_) / (2 * h)
    dW = np.zeros((Q, Q, P))
    for j in range(P):
        dW[:, :, j] = -W @ dOm[:, :, j] @ W
    # the moment Jacobian does not depend on the weight matrix
    corr = np.zeros((P, P))
    for j in range(P):
        corr += Ainv @ (A.T @ dW[:, :, j] @ A) @ Ainv
    Vw = base + corr / G
    return b2, base, base * G / (G - 1), Vw


def _om(Zc, X, y, b, cl, G):
    r = y - X @ b
    Om = np.zeros((Zc.shape[1], Zc.shape[1]))
    for gg in range(G):
        sel = cl == gg
        v = np.einsum("ij,i->j", Zc[sel], r[sel])
        Om += np.outer(v, v)
    return Om / G


def main():
    rng = np.random.default_rng(SEED)
    est, vb, vss, vw = [], [], [], []
    for _ in range(R):
        b, a, s, wv = one_rep(rng)
        est.append(b)
        vb.append(np.diag(a))
        vss.append(np.diag(s))
        vw.append(np.diag(wv))
    est = np.asarray(est)
    mc = est.std(0, ddof=1)
    out = dict(
        R=R,
        mc_sd=mc.tolist(),
        mean_plain=np.sqrt(np.asarray(vb)).mean(0).tolist(),
        mean_ss=np.sqrt(np.asarray(vss)).mean(0).tolist(),
        mean_wind=np.sqrt(np.mean(vw, 0)).tolist(),
    )
    print(json.dumps(out, indent=1))
    (HERE / "windmeijer_linear_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
