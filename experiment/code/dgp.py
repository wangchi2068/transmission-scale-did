"""Data-generating mechanisms of the Monte Carlo study (theory/simulation_spec.md, section D).

All randomness is drawn from two independent streams per replication:
  * the *structural* stream (unit effects, populations, pre-sample counts, adoption propensities and
    random adoption weeks, heterogeneous effects, reporting draws are NOT here), and
  * the *innovation* stream (weekly outcome draws and binomial reporting draws).
Using the same two seeds with ``counterfactual=True`` reproduces the same units and the same innovation
stream under the policy array d = 0 (common random numbers), which is how the population cases-averted
target psi is computed.

Geometry is fixed per N (drawn once from a fixed seed): N points uniform in the unit square, symmetrised
k-nearest-neighbour graph (k = 6), G = 50 contiguous equal-size "states" from a 5 x 10 strip partition.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional, Sequence

import numpy as np
import scipy.sparse as sps
from scipy.spatial import cKDTree

S_LAGS = 2                      # pre-sample weeks / generation-interval length
WGI_DEFAULT = (0.85, 0.15)      # weekly generation-interval weights used for estimation


# ------------------------------------------------------------------------------------------------
# Geography
# ------------------------------------------------------------------------------------------------
@dataclass
class Geography:
    W: sps.csr_matrix           # row-normalised adjacency, zero diagonal
    cl: np.ndarray              # cluster (state) index of each unit, 0..G-1
    coords: np.ndarray          # (N, 2)
    G: int

    @property
    def N(self) -> int:
        return self.W.shape[0]


def make_geography(N: int, G: int = 50, k: int = 6, seed: int = 7, cols: int = 5) -> Geography:
    """Random planar k-NN graph with G contiguous equal-size states (cols x G/cols strips)."""
    if G % cols or N % G:
        raise ValueError("N must be divisible by G and G by cols")
    rng = np.random.default_rng([seed, N])
    xy = rng.random((N, 2))
    tree = cKDTree(xy)
    _, nb = tree.query(xy, k + 1)
    rows = np.repeat(np.arange(N), k)
    colsi = nb[:, 1:].ravel()
    A = sps.csr_matrix((np.ones(rows.size), (rows, colsi)), shape=(N, N))
    A = ((A + A.T) > 0).astype(float).tocsr()
    A.setdiag(0.0)
    A.eliminate_zeros()
    deg = np.asarray(A.sum(1)).ravel()
    W = (sps.diags(1.0 / deg) @ A).tocsr()
    # contiguous equal-size states: split by x into `cols` strips, then each strip by y
    rows_per_col = G // cols
    per_state = N // G
    cl = np.empty(N, dtype=np.int64)
    order_x = np.argsort(xy[:, 0], kind="stable")
    strip = N // cols
    g = 0
    for c in range(cols):
        idx = order_x[c * strip:(c + 1) * strip]
        idx = idx[np.argsort(xy[idx, 1], kind="stable")]
        for r in range(rows_per_col):
            cl[idx[r * per_state:(r + 1) * per_state]] = g
            g += 1
    return Geography(W=W, cl=cl, coords=xy, G=G)


# ------------------------------------------------------------------------------------------------
# Scenario parameters
# ------------------------------------------------------------------------------------------------
@dataclass
class DGPParams:
    N: int = 500
    T: int = 12
    adopt: str = "feedback"                 # "feedback" | "random"
    tau: float = float(np.log(0.8))
    kappa: float = -0.1
    phi: float = 0.4
    k: float = 10.0                         # NB shape; np.inf = Poisson
    gbar: float = 0.0                       # calibrated constant in gamma_t
    w_true: Sequence[float] = WGI_DEFAULT   # generation interval of the DGP
    thr: float = 30.0                       # state threshold (lagged per-capita incidence)
    thr_pi: float = 0.25                    # threshold_g = thr * exp(-thr_pi * pi_g)
    ex_thr: float = 12.0                    # county exemption cut-off (lagged per-capita cases)
    alpha_state_sd: float = 0.2
    alpha_unit_sd: float = 0.2
    alpha_pi: float = 0.3                   # confounding by level: alpha_i += alpha_pi * pi_g
    omega_scale: float = 0.5                # omega_i = omega_scale * scale_i
    pre_mean: float = 20.0                  # pre-sample Y ~ Poisson(pre_mean * scale_i)
    # misspecification switches
    report_p: Optional[float] = None        # X2: Y ~ Binomial(I, p_it)
    report_drop: float = 0.85               # X2: p multiplied by this after state adoption
    tau_sd: float = 0.0                     # X3: tau_i ~ N(tau, tau_sd^2)
    depletion_pop: Optional[float] = None   # X4: P_i = depletion_pop * scale_i
    dyn: Optional[Sequence[float]] = None   # X5: event-time effects (last value carried forward)

    def to_dict(self):
        d = asdict(self)
        d["k"] = "inf" if np.isinf(self.k) else self.k
        d["w_true"] = list(self.w_true)
        d["dyn"] = None if self.dyn is None else list(self.dyn)
        return d


def gamma_path(p: DGPParams) -> np.ndarray:
    t = np.arange(p.T)
    return -np.log1p(p.phi) + p.gbar + 0.1 * np.sin(t / 3.0)


def _draw_nb(rng, mu, k):
    mu = np.maximum(mu, 0.0)
    if np.isinf(k):
        return rng.poisson(mu).astype(float)
    return rng.negative_binomial(k, k / (k + mu)).astype(float)


# ------------------------------------------------------------------------------------------------
# Simulation
# ------------------------------------------------------------------------------------------------
def simulate(p: DGPParams, geo: Geography, seed_struct, seed_innov, counterfactual: bool = False) -> dict:
    """Simulate one weekly panel (S = 2 pre-sample weeks plus T estimation weeks).

    Returns a dict with full-length arrays (length S + T along axis 0):
      Yall (reported cases), Iall (infections; = Yall unless report_p is set), Dall (policy),
    and unit-level quantities alpha, scale, omega (infection scale), tau_i, ev (adoption week
    index within the estimation window, -1 if never), state_adopt_week, eff (true effect array, T x N).
    """
    rs = np.random.default_rng(seed_struct)
    ri = np.random.default_rng(seed_innov)
    N, T, S = geo.N, p.T, S_LAGS
    G, cl, W = geo.G, geo.cl, geo.W
    w = np.asarray(p.w_true, float)
    if len(w) != S:
        raise ValueError("generation interval must have length S = 2")
    # structural draws (always in the same order so that counterfactual runs share them)
    scale = np.exp(rs.normal(0.0, 0.5, N))
    pi_g = rs.normal(0.0, 1.0, G)
    pi_g = (pi_g - pi_g.mean()) / pi_g.std()
    alpha = rs.normal(0.0, p.alpha_state_sd, G)[cl] + rs.normal(0.0, p.alpha_unit_sd, N) + p.alpha_pi * pi_g[cl]
    adopt_week_random = rs.integers(2, T + 3, G)            # 1-based week; >= T+1 means never
    tau_i = p.tau + (rs.normal(0.0, p.tau_sd, N) if p.tau_sd > 0 else np.zeros(N))
    Ipre = rs.poisson(p.pre_mean * scale[None, :] * np.ones((S, 1))).astype(float)
    omega = p.omega_scale * scale
    gam = gamma_path(p)
    thr_g = p.thr * np.exp(-p.thr_pi * pi_g)
    dyn = None if p.dyn is None else np.asarray(p.dyn, float)
    pop = None if p.depletion_pop is None else p.depletion_pop * scale

    I = np.zeros((S + T, N)); Y = np.zeros((S + T, N)); D = np.zeros((S + T, N))
    I[:S] = Ipre
    if p.report_p is not None:
        Y[:S] = ri.binomial(Ipre.astype(np.int64), p.report_p).astype(float)
    else:
        Y[:S] = Ipre
    cum = I[:S].sum(0)
    state_adopted = np.zeros(G, bool)
    state_week = np.full(G, -1)
    ev = np.full(N, -1)
    eff = np.zeros((T, N))
    for j in range(T):                       # j = 0..T-1 is week t = j+1
        t = S + j
        if not counterfactual:
            if p.adopt == "feedback":
                inc = np.bincount(cl, Y[t - 1] / scale, G) / np.bincount(cl, None, G)
                new_state = (~state_adopted) & (inc > thr_g)
                state_week[new_state] = j
                state_adopted |= new_state
                D[t] = np.maximum(D[t - 1], state_adopted[cl] & (Y[t - 1] / scale >= p.ex_thr))
            elif p.adopt == "random":
                on = (j + 1) >= adopt_week_random
                state_week[(state_week < 0) & on] = j
                state_adopted = on
                D[t] = on[cl].astype(float)
            else:
                raise ValueError(p.adopt)
        newly = (D[t] == 1) & (ev < 0)
        ev[newly] = j
        if dyn is not None:
            e = np.clip(j - ev, 0, len(dyn) - 1)
            eff[j] = np.where(D[t] == 1, dyn[e], 0.0)
        else:
            eff[j] = tau_i * D[t]
        Lloc = w[0] * I[t - 1] + w[1] * I[t - 2]
        lam = Lloc + p.phi * (W @ Lloc) + omega
        mu = np.exp(alpha + gam[j] + eff[j] + p.kappa * (W @ D[t])) * lam
        if pop is not None:
            mu = mu * np.clip(1.0 - cum / pop, 0.0, 1.0)
        I[t] = _draw_nb(ri, mu, p.k)
        cum += I[t]
        if p.report_p is not None:
            pr = np.full(N, p.report_p)
            if not counterfactual:
                pr = pr * np.where(state_adopted[cl], p.report_drop, 1.0)
            Y[t] = ri.binomial(I[t].astype(np.int64), pr).astype(float)
        else:
            Y[t] = I[t]
    return dict(Yall=Y, Iall=I, Dall=D, alpha=alpha, scale=scale, omega=omega, tau_i=tau_i, ev=ev,
                state_week=state_week, eff=eff, gamma=gam, pi_g=pi_g, pop=pop)


# ------------------------------------------------------------------------------------------------
# Estimation panel
# ------------------------------------------------------------------------------------------------
def build_panel(Yall: np.ndarray, Dall: np.ndarray, W: sps.csr_matrix, omega: np.ndarray, cl: np.ndarray,
                w: Sequence[float] = WGI_DEFAULT, S: int = S_LAGS) -> dict:
    """Estimation-window arrays (T x N) used by every estimator.

    L = sum_s w_s Y_{t-s} (local pressure), Nb = W L (imported pressure), Db = W D (neighbour policy share),
    ev = adoption week index within the window (-1 if never; absorbing policies only)."""
    w = np.asarray(w, float)
    Tt = Yall.shape[0] - S
    L = np.zeros((Tt, Yall.shape[1]))
    for s in range(1, len(w) + 1):
        L += w[s - 1] * Yall[S - s:S - s + Tt]
    Nb = np.asarray((W @ L.T).T)
    D = Dall[S:].astype(float)
    Db = np.asarray((W @ D.T).T)
    first = np.where(D.any(0), D.argmax(0), -1)
    return dict(Y=Yall[S:].astype(float), D=D, Db=Db, L=L, Nb=Nb, omega=np.asarray(omega, float),
                cl=np.asarray(cl), G=int(np.max(cl)) + 1, ev=first, Ypre=Yall[:S].astype(float),
                w=w, W=W, S=S)


def panel_descriptives(panel: dict) -> dict:
    Y, D = panel["Y"], panel["D"]
    sw = ((D[1:] == 1) & (D[:-1] == 0)).sum()
    return dict(median_Y=float(np.median(Y)), mean_Y=float(Y.mean()), frac_zero=float((Y == 0).mean()),
                share_treated=float(D.mean()), n_switch=int(sw),
                n_ever=int((panel["ev"] >= 0).sum()), n_switch_at_1=int(((panel["ev"] == 0)).sum()))
