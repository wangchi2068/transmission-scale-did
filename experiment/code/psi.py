"""Counterfactual recursion nu(theta) and cases averted psi (methods.tex eq:nu, Theorem 4).

nu_it = Y_it for pre-sample weeks; for t >= 1, jointly over units,
    nu_it = c_it Y_it,  c_it = exp(-index_it) * Lambda_it(nu; phi) / Lambda_it(Y; phi),
where index_it is the policy part of the linear predictor (tau D + kappa Dbar, or event-time effects).
psi_P = sum_{(i,t) in P} (nu_it - Y_it). Delta-method SE uses a numerical gradient with respect to beta.
"""
from __future__ import annotations

from typing import Callable, Dict, Optional, Sequence

import numpy as np


def nu_recursion(Yall: np.ndarray, index: np.ndarray, W, omega: np.ndarray, w: Sequence[float], phi: float,
                 S: int = 2) -> np.ndarray:
    """Yall: (S+T) x N observed path including S pre-sample weeks; index: T x N policy index.
    Returns nu for the estimation window (T x N)."""
    w = np.asarray(w, float)
    T = Yall.shape[0] - S
    nu = np.array(Yall, float, copy=True)
    for j in range(T):
        t = S + j
        ly = sum(w[s - 1] * Yall[t - s] for s in range(1, len(w) + 1))
        ln = sum(w[s - 1] * nu[t - s] for s in range(1, len(w) + 1))
        lamY = ly + phi * (W @ ly) + omega
        lamN = ln + phi * (W @ ln) + omega
        nu[t] = np.exp(-index[j]) * lamN / lamY * Yall[t]
    return nu[S:]


def psi_sets(D: np.ndarray) -> Dict[str, np.ndarray]:
    """Target sets: treated = post-adoption county-weeks (D = 1); network = all county-weeks
    (nu - Y = 0 before any policy, so this equals all post-first-adoption weeks)."""
    return {"treated": D == 1, "network": np.ones_like(D, dtype=bool)}


def psi_values(Yall, index, W, omega, w, phi, D, S=2) -> Dict[str, float]:
    nu = nu_recursion(Yall, index, W, omega, w, phi, S)
    diff = nu - Yall[S:]
    return {k: float(diff[m].sum()) for k, m in psi_sets(D).items()}


def psi_with_se(Yall, panel: dict, beta: np.ndarray, Vb: Optional[np.ndarray],
                index_fn: Callable[[np.ndarray], tuple], h: float = 1e-5) -> dict:
    """index_fn(beta) -> (index T x N, phi). Delta method: se = sqrt(Psi Vb Psi')."""
    W, om, w, S, D = panel["W"], panel["omega"], panel["w"], panel["S"], panel["D"]

    def f(b):
        idx, phi = index_fn(b)
        return psi_values(Yall, idx, W, om, w, phi, D, S)

    base = f(beta)
    keys = list(base)
    out = {k: {"est": base[k], "se": float("nan")} for k in keys}
    if Vb is None or not np.all(np.isfinite(Vb)):
        return out
    grad = np.zeros((len(keys), len(beta)))
    for j in range(len(beta)):
        e = np.zeros(len(beta)); e[j] = h
        fp, fm = f(beta + e), f(beta - e)
        for a, k in enumerate(keys):
            grad[a, j] = (fp[k] - fm[k]) / (2 * h)
    for a, k in enumerate(keys):
        v = float(grad[a] @ Vb @ grad[a])
        out[k]["se"] = float(np.sqrt(v)) if v >= 0 else float("nan")
    return out
