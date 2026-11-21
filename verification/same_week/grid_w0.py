"""Fine grid of within-week transmission w0 in {0.0, 0.05, 0.10, 0.20, 0.30} plus exact share
decomposition and treated-share-matched control (Major Comment 2 & E6).
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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiment" / "code"))

import dgp  # noqa: E402
import qdgmm  # noqa: E402
import same_week_dgp as sw  # noqa: E402
import scenarios as sc  # noqa: E402

OUT_PATH = HERE / "grid_results.json"
_G = {}


def _init_worker():
    p = sc.build_scenarios()[1].params
    _G["p"] = p
    _G["geo"] = dgp.make_geography(p.N)
    _G["w_est"] = sc.build_scenarios()[1].w_est
    _G["cal_g"] = float(json.loads((HERE / "calib.json").read_text(encoding="utf-8"))["g"])


def _run_fixed_D_same_week(p, geo, seed_struct, seed_innov, w0: float, D_override: np.ndarray, g_val: float):
    """Exact copy of simulate_same_week forcing Dall = D_override to equalise treated share."""
    rs = np.random.default_rng(seed_struct)
    ri = np.random.default_rng(seed_innov)
    N, T, S = geo.N, p.T, dgp.S_LAGS
    G, cl, W = geo.G, geo.cl, geo.W
    w1, w2 = sw.weekly_weights(w0)
    c = w0 / sw.W_DAYS * g_val
    scale = np.exp(rs.normal(0.0, 0.5, N))
    pi_g = rs.normal(0.0, 1.0, G)
    pi_g = (pi_g - pi_g.mean()) / pi_g.std()
    alpha = rs.normal(0.0, p.alpha_state_sd, G)[cl] + rs.normal(0.0, p.alpha_unit_sd, N) + p.alpha_pi * pi_g[cl]
    _ = rs.integers(2, T + 3, G)
    Ipre = rs.poisson(p.pre_mean * scale[None, :] * np.ones((S, 1))).astype(float)
    omega = p.omega_scale * scale
    gam = dgp.gamma_path(p)
    Y = np.zeros((S + T, N))
    Y[:S] = Ipre
    D = D_override
    reinf_tot = 0.0
    seed_tot = 0.0
    all_tot = 0.0
    for j in range(T):
        t = S + j
        base = np.exp(alpha + gam[j] + p.tau * D[t] + p.kappa * np.asarray(W @ D[t]))
        Lp = w1 * (Y[t - 1] + p.phi * np.asarray(W @ Y[t - 1])) + w2 * (Y[t - 2] + p.phi * np.asarray(W @ Y[t - 2]))
        if not ((w0 > 0) and (j >= sw.ENGINE_FROM)):
            Y[t] = sw._draw_nb(ri, base * (Lp + omega), p.k)
        else:
            lag_only = base * (Lp + omega)
            seed_only = base * w0 * (Lp + omega)
            bg = (Lp + omega + w0 * (Lp + omega)) / sw.W_DAYS
            cur = np.zeros((sw.W_DAYS, N))
            esc = np.zeros(N)
            reinf_week = np.zeros(N)
            for d_day in range(sw.W_DAYS):
                mu_day = base * (bg + c * esc)
                reinf_week += base * c * esc
                cur[d_day] = sw._draw_nb(ri, mu_day, p.k)
                esc = esc + cur[d_day]
            Y[t] = cur.sum(0)
            reinf_tot += float(reinf_week.sum())
            seed_tot += float(seed_only.sum())
            all_tot += float((lag_only + seed_only + reinf_week).sum())
    return dict(Yall=Y, Dall=D, omega=omega, reinf_share=(reinf_tot / all_tot if all_tot > 0 else 0.0),
                seeded_inc_share=((reinf_tot + seed_tot) / all_tot if all_tot > 0 else 0.0))


def _fit_panel(Yall, Dall, geo, omega, w_est, params):
    pan = dgp.build_panel(Yall, Dall, geo.W, omega, geo.cl, w=w_est)
    spec = qdgmm.spec_m0(pan)
    res = qdgmm.fit(pan, spec, starts=[[params.tau, params.kappa, params.phi], [0.0, 0.0, 0.4]])
    if not res["ok"]:
        return None
    b = res["beta"]
    se = np.array([res["se"]["tau"], res["se"]["kappa"], res["se"]["phi"]])
    V = res["Vb"]
    tpk = float(b[0] + b[1])
    tpk_se = float(np.sqrt(max(V[0, 0] + V[1, 1] + 2 * V[0, 1], 0.0)))
    tc = float(stats.t.ppf(0.975, int(pan["G"]) - 1))
    jres = qdgmm.fit_two_step(pan, qdgmm.spec_m0_oi(pan), res)
    return dict(
        tau=float(b[0]), kappa=float(b[1]), phi=float(b[2]), tpk=tpk,
        tau_se=float(se[0]), tpk_se=tpk_se,
        cov_tau=(abs(b[0] - params.tau) <= tc * se[0]),
        cov_tpk=(abs(tpk - (params.tau + params.kappa)) <= tc * tpk_se),
        j_rej=bool(jres.get("J_p", np.nan) < 0.05) if np.isfinite(jres.get("J_p", np.nan)) else False,
        treated_share=float(pan["D"].mean()),
    )


def _run_one(args):
    r, w0 = args
    p, geo, w_est, g_val = _G["p"], _G["geo"], _G["w_est"], _G["cal_g"]
    sw.CAL["g"] = g_val
    seed = sc.seed_for(1, r) + 100_000
    sim_bench = dgp.simulate(p, geo, [seed, 0], [seed, 1])
    if w0 == 0.0:
        out_free = _fit_panel(sim_bench["Yall"], sim_bench["Dall"], geo, sim_bench["omega"], w_est, p)
        return dict(w0=w0, free=out_free, matched=out_free, reinf_share=0.0, seeded_inc_share=0.0)

    # Free-feedback arm (matching aggregate.py / same_week_dgp.py)
    sim_sw = sw.simulate_same_week(p, geo, [seed, 0], [seed, 1], w0)
    out_free = _fit_panel(sim_sw["Yall"], sim_sw["Dall"], geo, sim_sw["omega"], w_est, p)

    # Treated-share matched arm using exact D from w0=0
    sim_m = _run_fixed_D_same_week(p, geo, [seed, 0], [seed, 1], w0, sim_bench["Dall"], g_val)
    out_matched = _fit_panel(sim_m["Yall"], sim_m["Dall"], geo, sim_m["omega"], w_est, p)
    return dict(
        w0=w0, free=out_free, matched=out_matched,
        reinf_share=sim_m["reinf_share"], seeded_inc_share=sim_m["seeded_inc_share"],
    )


def _summ(rows, truth_tau, truth_tpk):
    ok = [x for x in rows if x is not None]
    n = len(ok)
    taus = np.array([x["tau"] for x in ok])
    tpks = np.array([x["tpk"] for x in ok])
    t_b = taus - truth_tau
    k_b = tpks - truth_tpk
    return dict(
        n=n,
        treated_share=float(np.mean([x["treated_share"] for x in ok])),
        tau_bias=float(t_b.mean()), tau_bias_mcse=float(t_b.std(ddof=1) / np.sqrt(n)),
        tau_cov=float(np.mean([x["cov_tau"] for x in ok])),
        tpk_bias=float(k_b.mean()), tpk_bias_mcse=float(k_b.std(ddof=1) / np.sqrt(n)),
        tpk_cov=float(np.mean([x["cov_tpk"] for x in ok])),
        j_rej=float(np.mean([x["j_rej"] for x in ok])),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 4, 18)))
    args = ap.parse_args()

    t0 = time.time()
    p = sc.build_scenarios()[1].params
    grid_w0 = [0.0, 0.05, 0.10, 0.20, 0.30]
    summary = {}
    for w0 in grid_w0:
        tasks = [(r, w0) for r in range(args.R)]
        with ProcessPoolExecutor(max_workers=args.workers, initializer=_init_worker) as ex:
            res = list(ex.map(_run_one, tasks, chunksize=5))
        summary[f"w0_{w0:.2f}"] = dict(
            w0=w0,
            free=_summ([x["free"] for x in res], p.tau, p.tau + p.kappa),
            matched_D=_summ([x["matched"] for x in res], p.tau, p.tau + p.kappa),
            mean_reinf_share=float(np.mean([x["reinf_share"] for x in res])),
            mean_seeded_inc_share=float(np.mean([x["seeded_inc_share"] for x in res])),
        )
        s_f = summary[f"w0_{w0:.2f}"]["free"]
        s_m = summary[f"w0_{w0:.2f}"]["matched_D"]
        print(
            f"w0={w0:.2f} | reinf={summary[f'w0_{w0:.2f}']['mean_reinf_share']:.3f} "
            f"(+seed={summary[f'w0_{w0:.2f}']['mean_seeded_inc_share']:.3f}) | "
            f"free(D={s_f['treated_share']:.3f}): tpk_bias={s_f['tpk_bias']:+.4f}, cov={s_f['tpk_cov']:.3f}, J={s_f['j_rej']:.3f} | "
            f"matchedD(D={s_m['treated_share']:.3f}): tpk_bias={s_m['tpk_bias']:+.4f}, cov={s_m['tpk_cov']:.3f}",
            flush=True,
        )

    out = dict(R=args.R, seconds=round(time.time() - t0, 2), grid=summary)
    OUT_PATH.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Saved {OUT_PATH} in {out['seconds']}s")


if __name__ == "__main__":
    main()
