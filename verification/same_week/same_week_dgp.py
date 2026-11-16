"""Within-week ("same-week") propagation DGP for the M2 sensitivity (verification/same_week).

Weekly skeleton identical to ``experiment/code/dgp.simulate``; the weekly generation interval keeps the
lag-1/lag-2 ratio of the estimation weights for every w0:

    w1 = 0.85 (1 - w0),   w2 = 0.15 (1 - w0),   same-week share w0.

Weeks with within-week generation (engine weeks, weeks 2-12 of the estimation window; the first two weeks
of the window, whose lags 1-2 reach into the pre-sample weeks, use the weekly draw):

  * day d of week t:  Y_d ~ NB(k, k / (k + mu_d)),
        mu_d = exp(alpha + gamma_t + eff + kappa Dbar) * [(Lp + omega + r_inf) / 7 + c S_<d],
    Lp = (1 - w0) [0.85 (Y_{t-1} + phi W Y_{t-1}) + 0.15 (Y_{t-2} + phi W Y_{t-2})]   (lagged, own + neighbours),
    S_<d = Y_0 + ... + Y_{d-1} of week t,   r_inf = w0 (Lp + omega)                  (same-week mass),
    c = (w0 / 7) g                                                        (same-week re-entry rate),
  with g a calibration constant (pretest.py): the re-entry c = (w0 / 7) g makes the expected week-d
    total increase geometrically inside the week, c (1 + c + c^2 + ...) = w0 g, so week t's expected total
    is exp(...) (Lp + omega) (1 + w0 g), and the lag-1 W content (0.85) is carried at (1 - w0) 0.85 as in the
    weekly model.  g = 1 is the plain within-week re-entry (A = w0); g is fitted on the panel check.

The re-entry of a same-week infection therefore enters its own unit's later days of the same week with rate c
(no within-week cross-unit transfer); at the aggregate this is a within-week renewal 1 / (1 - c) that the
weekly model cannot express, which is the misspecification of interest.  g is calibrated in pretest.py by the
panel check (engine / weekly means of the main core scenario).

w0 = 0 uses exactly the weekly draw of dgp.simulate (bit-identical, no extra RNG draws), so the benchmark arm
is the main simulation.
"""
from __future__ import annotations

import json
import os

import numpy as np

from dgp import DGPParams, Geography, S_LAGS, gamma_path

W_DAYS = 7
ENGINE_FROM = 2                  # window index j of the first engine week (j = 2 -> t = 4)
CAL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "calib.json")
_G = 0.5                         # fallback; overwritten by calib.json when present
if os.path.exists(CAL_PATH):
    with open(CAL_PATH) as _fh:
        _G = float(json.load(_fh).get("g", _G))
CAL = {"g": _G}


def weekly_weights(w0: float) -> tuple:
    return 0.85 * (1.0 - w0), 0.15 * (1.0 - w0)


def reentry_rate(w0: float, g: float = None) -> float:
    g = CAL["g"] if g is None else g
    return w0 / W_DAYS * g


def _draw_nb(rng, mu, k):
    mu = np.maximum(mu, 0.0)
    if np.isinf(k):
        return rng.poisson(mu).astype(float)
    return rng.negative_binomial(k, k / (k + mu)).astype(float)


def simulate_same_week(p: DGPParams, geo: Geography, seed_struct, seed_innov, w0: float,
                       counterfactual: bool = False) -> dict:
    """One weekly panel; same dict as dgp.simulate plus ``engine`` (weeks generated within-week)."""
    rs = np.random.default_rng(seed_struct)
    ri = np.random.default_rng(seed_innov)
    N, T, S = geo.N, p.T, S_LAGS
    G, cl, W = geo.G, geo.cl, geo.W
    if S != 2:
        raise ValueError("S = 2 only")
    w1, w2 = weekly_weights(w0)
    c = reentry_rate(w0)
    # ---- structural draws: same order as dgp.simulate
    scale = np.exp(rs.normal(0.0, 0.5, N))
    pi_g = rs.normal(0.0, 1.0, G)
    pi_g = (pi_g - pi_g.mean()) / pi_g.std()
    alpha = rs.normal(0.0, p.alpha_state_sd, G)[cl] + rs.normal(0.0, p.alpha_unit_sd, N) + p.alpha_pi * pi_g[cl]
    adopt_week_random = rs.integers(2, T + 3, G)
    tau_i = p.tau + (rs.normal(0.0, p.tau_sd, N) if p.tau_sd > 0 else np.zeros(N))
    Ipre = rs.poisson(p.pre_mean * scale[None, :] * np.ones((S, 1))).astype(float)
    omega = p.omega_scale * scale
    gam = gamma_path(p)
    thr_g = p.thr * np.exp(-p.thr_pi * pi_g)
    dyn = None if p.dyn is None else np.asarray(p.dyn, float)
    pop = None if p.depletion_pop is None else p.depletion_pop * scale
    I = np.zeros((S + T, N))
    Y = np.zeros((S + T, N))
    D = np.zeros((S + T, N))
    I[:S] = Ipre
    Y[:S] = Ipre
    cum = I[:S].sum(0)
    state_adopted = np.zeros(G, bool)
    state_week = np.full(G, -1)
    ev = np.full(N, -1)
    eff = np.zeros((T, N))
    engine = np.zeros(T, bool)
    for j in range(T):
        t = S + j
        if not counterfactual:
            if p.adopt == "feedback":
                inc = np.bincount(cl, np.clip(Y[t - 1] / scale, 0.0, None), G) / np.bincount(cl, None, G)
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
        base = np.exp(alpha + gam[j] + eff[j] + p.kappa * np.asarray(W @ D[t]))
        Lp = w1 * (Y[t - 1] + p.phi * np.asarray(W @ Y[t - 1])) + w2 * (Y[t - 2] + p.phi * np.asarray(W @ Y[t - 2]))
        if not ((w0 > 0) and (j >= ENGINE_FROM)):
            mu = base * (Lp + omega)
            if pop is not None:
                mu = mu * np.clip(1.0 - cum / pop, 0.0, 1.0)
            draw = _draw_nb(ri, mu, p.k)
            Y[t] = draw
            I[t] = draw
        else:
            engine[j] = True
            bg = (Lp + omega + w0 * (Lp + omega)) / W_DAYS     # per-day outside + same-week-seed mass
            cur = np.zeros((W_DAYS, N))
            esc = np.zeros(N)
            for d in range(W_DAYS):
                mu = base * (bg + c * esc)
                if pop is not None:
                    mu = mu * np.clip(1.0 - cum / pop, 0.0, 1.0)
                cur[d] = _draw_nb(ri, mu, p.k)
                esc = esc + cur[d]
            I[t] = cur.sum(0)
            Y[t] = I[t]
        cum += I[t]
    return dict(Yall=Y, Iall=I, Dall=D, alpha=alpha, scale=scale, omega=omega, tau_i=tau_i, ev=ev,
                state_week=state_week, eff=eff, gamma=gam, pi_g=pi_g, pop=pop, engine=engine)
