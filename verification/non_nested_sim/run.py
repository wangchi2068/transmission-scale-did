"""Non-nested daily stochastic epidemic DGP aggregated to weeks (Major Comment 1(b)-(c)).

Data are NOT generated from the working QD-GMM weekly equation:
  - Daily (sub-weekly) stochastic branching renewal across d = 1..7*(T+S) days using a continuous
    Gamma(mean=12d, sd=3.5d) generation interval discretised on days 1..28.
  - True background importation fluctuates weakly over time & deviates from unit population.
  - Reporting delay: each day-d infection is reported at day d + Delta_d with Delta_d drawn from a
    discrete multinomial Gamma(mean=delay_mean, sd=1.5d) delay distribution across days 0..14.
  - Analyst does NOT know true generation weights w or true omega: the analyst uses noisy empirical
    serial-interval weights w_hat (perturbed mean 13.5d, sd 4.0d -> lags 1..3) and pop-scaled omega_hat.
  - Scenarios (R = 500 each):
      NN1_daily_gamma_gi     : daily renewal, 0d delay, mis-measured w_hat & omega_hat
      NN2_delay_2d           : + stochastic reporting delay (mean 2 days)
      NN3_delay_4d           : + stronger reporting delay (mean 4 days)
      NN4_seir_depletion     : + finite susceptible depletion (5%) & post-adoption day-10% reporting jitter
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiment" / "code"))

import comparators as cmp  # noqa: E402
import dgp  # noqa: E402
import qdgmm  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent / "results.json"
_G = {}


def _daily_gamma_kernel(mean: float, sd: float, max_d: int) -> np.ndarray:
    shape = (mean / sd) ** 2
    scale = (sd ** 2) / mean
    edges = np.arange(0.5, max_d + 1.5, 1.0)
    cdf = stats.gamma.cdf(edges, a=shape, scale=scale)
    k = np.diff(cdf)
    # Zero out days 0..3 so multi-week structure has negligible <4d generation
    return k / k.sum()


def _init_worker():
    N, G, T, S = 500, 50, 16, 3
    geo = dgp.make_geography(N, G=G)
    # True daily generation kernel (mean 12d, sd 3.5d, days 1..28)
    k_gi = _daily_gamma_kernel(12.0, 3.5, 28)
    # Working weekly weights estimated imperfectly from noisy literature (mean 13.5d, sd 4.0d)
    days = np.arange(1, 22)
    shape_e = (13.5 / 4.0) ** 2
    scale_e = (4.0 ** 2) / 13.5
    p_d = stats.gamma.cdf(days + 0.5, a=shape_e, scale=scale_e) - stats.gamma.cdf(days - 0.5, a=shape_e, scale=scale_e)
    w_raw = np.array([p_d[:7].sum(), p_d[7:14].sum(), p_d[14:21].sum()])
    w_hat = w_raw / w_raw.sum()
    _G.update(N=N, G=G, T=T, S=S, geo=geo, k_gi=k_gi, w_hat=tuple(w_hat.tolist()))


def _sim_daily(seed: int, delay_mean: float, use_depl: bool, rep_jitter: float):
    N, G, T, S = _G["N"], _G["G"], _G["T"], _G["S"]
    geo, k_gi = _G["geo"], _G["k_gi"]
    W_arr = geo.W.toarray()
    cl = geo.cl
    rs = np.random.default_rng(seed)
    ri = np.random.default_rng(seed + 50_000_000)

    tau0 = float(np.log(0.8))
    kappa0 = -0.10
    phi0 = 0.35
    k_nb = 12.0
    D_gi = len(k_gi)
    tot_weeks = S + T + 4
    n_days = tot_weeks * 7
    pre_days = 4 * 7

    # Structural draws
    pop_scale = np.exp(rs.normal(0.0, 0.45, N))
    pi_g = rs.normal(0.0, 1.0, G); pi_g = (pi_g - pi_g.mean()) / pi_g.std()
    alpha = -0.15 + 0.25 * (rs.normal(0.0, 0.2, G)[cl] + rs.normal(0.0, 0.2, N) + 0.3 * pi_g[cl])
    omega_true = 1.8 * pop_scale * np.exp(rs.normal(0.0, 0.15, N))  # daily omega*7; analyst sees pop_scale only
    omega_hat = 1.6 * pop_scale  # misspecified background scaling used by analyst
    gam_day = np.repeat(0.08 * np.sin(np.linspace(0, 1.6 * np.pi, S + T)), 7)
    thr_g = 38.0 * np.exp(-0.22 * pi_g)

    # Delay kernel (0..14 days)
    if delay_mean <= 0:
        d_probs = np.array([1.0])
    else:
        max_dl = 14
        sh = (delay_mean / max(0.5, delay_mean * 0.5)) ** 2
        sc = (max(0.5, delay_mean * 0.5) ** 2) / delay_mean
        ed = np.arange(-0.5, max_dl + 1.5, 1.0)
        d_probs = np.diff(stats.gamma.cdf(ed, a=sh, scale=sc))
        d_probs[0] += max(0.0, 1.0 - d_probs.sum())
        d_probs /= d_probs.sum()

    I_day = np.zeros((n_days, N))
    Y_rep_day = np.zeros((n_days + 15, N))
    I_day[:pre_days] = rs.poisson(2.2 * pop_scale[None, :], size=(pre_days, N)).astype(float)
    Y_rep_day[:pre_days] = I_day[:pre_days]
    D_week = np.zeros((S + T, N))
    state_on = np.zeros(G, bool)
    cum0 = I_day[:pre_days].sum(0)
    cum = cum0.copy()
    pop_cap = (1500.0 * pop_scale) if use_depl else None

    for w_idx in range(S + T):
        # Weekly feedback decision at start of window week w_idx >= S
        if w_idx >= S:
            p0 = (w_idx + 4 - 1) * 7
            p1 = (w_idx + 4) * 7
            prev_inc = np.bincount(cl, Y_rep_day[p0:p1].sum(0) / np.maximum(pop_scale, 0.1), G) / np.bincount(cl, None, G)
            state_on |= (prev_inc > thr_g)
            D_week[w_idx] = state_on[cl].astype(float)
        Db_w = W_arr @ D_week[w_idx]
        for dd in range(7):
            d = (w_idx + 4) * 7 + dd
            if d >= n_days:
                continue
            hist = I_day[d - D_gi : d][::-1]
            L_d = np.einsum("k,kn->n", k_gi, hist)
            Nb_d = L_d @ W_arr.T
            depl = np.clip(1.0 - (cum - cum0) / pop_cap, 0.2, 1.0) if pop_cap is not None else 1.0
            mu_d = np.exp(alpha + gam_day[d - pre_days] + tau0 * D_week[w_idx] + kappa0 * Db_w) * (L_d + phi0 * Nb_d + omega_true / 7.0) * depl
            I_new = ri.negative_binomial(k_nb, k_nb / (k_nb + np.maximum(mu_d, 1e-12))).astype(float)
            I_day[d] = I_new
            cum += I_new
            tot_c = I_new.astype(int)
            if rep_jitter > 0:
                tot_c = ri.binomial(tot_c, 1.0 - rep_jitter * D_week[w_idx])
            if len(d_probs) == 1:
                Y_rep_day[d] += tot_c
            else:
                nz = np.where(tot_c > 0)[0]
                for i in nz:
                    alloc = ri.multinomial(int(tot_c[i]), d_probs)
                    Y_rep_day[d : d + len(d_probs), i] += alloc

    # Aggregate into weekly counts over the S+T weeks
    Y_week = np.zeros((S + T, N))
    for w_idx in range(S + T):
        d0 = (w_idx + 4) * 7
        Y_week[w_idx] = Y_rep_day[d0 : d0 + 7].sum(0)

    pan = dgp.build_panel(Y_week, D_week, geo.W, omega_hat, geo.cl, w=_G["w_hat"], S=S)
    pan["D_mean"] = float(pan["D"].mean())
    return pan, Y_week, D_week, tau0, float(delay_mean)


def _run_one(args):
    r, delay_mean, use_depl, rep_jitter, seed_off = args
    pan, Y_week, D_week, tau0, dl = _sim_daily(20278000 + seed_off + r, delay_mean, use_depl, rep_jitter)
    G = int(pan["G"])
    tc = float(stats.t.ppf(0.975, G - 1))
    out = {"D_mean": pan["D_mean"]}

    try:
        m0 = qdgmm.fit(pan, qdgmm.spec_m0(pan, phi_bar=0.5))
    except np.linalg.LinAlgError:
        m0 = {"ok": False}
    out["M0"] = (float(m0["est"]["tau"]), float(m0["se"]["tau"]), abs(m0["est"]["tau"] - tau0) <= tc * m0["se"]["tau"]) if m0["ok"] else None

    if dl > 0:
        frac = min(1.0, dl / 7.0)
        D_adj = np.zeros_like(D_week)
        D_adj[1:] = (1.0 - frac) * D_week[1:] + frac * D_week[:-1]
        pan_adj = dgp.build_panel(Y_week, D_adj, pan["W"], pan["omega"], pan["cl"], w=pan["w"], S=3)
        try:
            m0a = qdgmm.fit(pan_adj, qdgmm.spec_m0(pan_adj, phi_bar=0.5))
        except np.linalg.LinAlgError:
            m0a = {"ok": False}
        out["M0_delay_adj"] = (float(m0a["est"]["tau"]), float(m0a["se"]["tau"]), abs(m0a["est"]["tau"] - tau0) <= tc * m0a["se"]["tau"]) if m0a["ok"] else None
    else:
        out["M0_delay_adj"] = out["M0"]

    m1 = cmp.poisson_fe(pan)
    out["M1"] = (float(m1["est"]["tau"]), float(m1["se"]["tau"]), abs(m1["est"]["tau"] - tau0) <= tc * m1["se"]["tau"]) if m1["ok"] else None

    m2_w = cmp.twfe_ols(cmp.cori_log_r(pan, a=0.1, b=5.0), pan["D"], pan["cl"], G)
    out["M2_weak"] = (float(m2_w["est"]["tau"]), float(m2_w["se"]["tau"]), abs(m2_w["est"]["tau"] - tau0) <= tc * m2_w["se"]["tau"]) if m2_w["ok"] else None

    m2_d = cmp.m2_logr_did(pan)
    out["M2_def"] = (float(m2_d["est"]["tau"]), float(m2_d["se"]["tau"]), abs(m2_d["est"]["tau"] - tau0) <= tc * m2_d["se"]["tau"]) if m2_d["ok"] else None

    m3 = cmp.m3_logy_twfe(pan)
    out["M3"] = (float(m3["est"]["tau"]), float(m3["se"]["tau"]), False) if m3["ok"] else None

    m4 = cmp.poisson_qmle(pan)
    out["M4"] = (float(m4["est"]["tau"]), float(m4["se"]["tau"]), abs(m4["est"]["tau"] - tau0) <= tc * m4["se"]["tau"]) if m4["ok"] else None
    return out


def _summ(rows, tau0):
    ok = [x for x in rows if x is not None]
    n = len(ok)
    if n == 0:
        return None
    arr = np.array(ok, float)
    b = arr[:, 0] - tau0
    return dict(
        n=n, bias=float(b.mean()), bias_mcse=float(b.std(ddof=1) / np.sqrt(n)),
        empSE=float(arr[:, 0].std(ddof=1)), meanSE=float(arr[:, 1].mean()),
        rmse=float(np.sqrt((b ** 2).mean())), coverage=float(arr[:, 2].mean()),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=500)
    ap.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 4, 18)))
    args = ap.parse_args()

    t0 = time.time()
    tau0 = float(np.log(0.8))
    specs = [
        ("NN1_daily_gamma_gi", 0.0, False, 0.0, 0),
        ("NN2_delay_2d", 2.0, False, 0.0, 100000),
        ("NN3_delay_4d", 4.0, False, 0.0, 200000),
        ("NN4_seir_depletion", 2.0, True, 0.05, 300000),
    ]
    out_sc = {}
    for name, dl, dp, rj, off in specs:
        tasks = [(r, dl, dp, rj, off) for r in range(args.R)]
        with ProcessPoolExecutor(max_workers=args.workers, initializer=_init_worker) as ex:
            res = list(ex.map(_run_one, tasks, chunksize=5))
        block = {"mean_D": float(np.mean([x["D_mean"] for x in res]))}
        for m in ("M0", "M0_delay_adj", "M1", "M2_weak", "M2_def", "M3", "M4"):
            block[m] = _summ([x[m] for x in res], tau0)
        out_sc[name] = block
        def _fmt(k):
            b = block[k]
            return f"{b['bias']:+.3f}({100*b['coverage']:4.1f}%)" if b else "   NA      "
        print(
            f"{name:<20} (D={block['mean_D']:.2f}) | M0={_fmt('M0')} | M0adj={_fmt('M0_delay_adj')} "
            f"| M1={_fmt('M1')} | M2w={_fmt('M2_weak')} | M2d={_fmt('M2_def')} | M4={_fmt('M4')}",
            flush=True,
        )

    out = dict(R=args.R, tau_true=tau0, seconds=round(time.time() - t0, 2), scenarios=out_sc)
    OUT_PATH.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Saved {OUT_PATH} in {out['seconds']}s")


if __name__ == "__main__":
    main()
