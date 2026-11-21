"""Controlled few-cluster & unbalanced-cluster size grid (Major Comment 3 & Minor Comment 4).

Grid A (pure G effect at fixed N = 2000, T = 12, feedback, spillover):
    G in {10, 16, 20, 30, 50, 100} balanced clusters -> isolates G cleanly from N!
Grid B (cluster skewness x T decomposition at N = 400):
    (balanced G=16 vs German-skewed G=16 [1..96/400]) x (T=12 vs T=60) -> quantifies exact source
    of the ~15% sandwich SE understatement (leverage skewness vs G vs T week profiling).
Compares: t_{G-1}, Hwang-style t_{G-2} FE-df adj, CV3/delete-one-cluster jackknife, WCR-Webb6, WCR-Rademacher.
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
import scipy.sparse as sps
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiment" / "code"))
sys.path.insert(0, str(ROOT / "verification" / "revision" / "common"))
sys.path.insert(0, str(ROOT / "verification" / "application"))

import dgp  # noqa: E402
import qdgmm  # noqa: E402
from qdx import QDX, tsum  # noqa: E402
import application as app  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent / "cluster_grid_results.json"
DE_SIZES = [1, 1, 3, 5, 6, 8, 11, 14, 15, 18, 26, 36, 44, 45, 71, 96]  # sum=400 German federal state county counts


def _make_skewed_geo(sizes: list[int]) -> dgp.Geography:
    N = sum(sizes)
    G = len(sizes)
    # Assign skewed cluster sizes contiguously along sorted x-coordinate strips so clusters are spatially connected
    geo0 = dgp.make_geography(N, G=G, cols=4)
    order = np.argsort(geo0.coords[:, 0], kind="stable")
    cl = np.empty(N, dtype=int)
    pos = 0
    for g, sz in enumerate(sizes):
        cl[order[pos : pos + sz]] = g
        pos += sz
    return dgp.Geography(W=geo0.W, cl=cl, coords=geo0.coords, G=G)


def _run_cell_worker(args):
    cell_name, N, G, T, skewed, adopt_mode, seed_base, r = args
    cols = 4 if G == 16 else (10 if G == 100 else 5)
    geo = _make_skewed_geo(DE_SIZES) if skewed else dgp.make_geography(N, G=G, cols=cols)
    gbar = 0.0 if adopt_mode == "random" else -0.02
    p = dgp.DGPParams(N=N, T=T, adopt=adopt_mode, tau=float(np.log(0.8)), kappa=-0.1, phi=0.4, k=10.0, gbar=gbar, alpha_state_sd=0.05, alpha_unit_sd=0.05, alpha_pi=0.075, thr=22.0)
    seed = seed_base + r
    sim = dgp.simulate(p, geo, [seed, 0], [seed, 1])
    pan = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl)
    tpk_true = float(p.tau + p.kappa)

    X = QDX(pan, {"tau": pan["D"], "kappa": pan["Db"]}, phi_bar=0.5)
    eh, ok, _ = X.fit(starts=[[-0.2, -0.1, 0.4], [0.0, 0.0, 0.4]])
    if not ok:
        return None
    inf = X.inference(eh)
    if inf is None:
        return None
    Vb = inf["Vb"]
    tpk = float(eh[0] + eh[1])
    se_tpk = float(tsum(Vb, 0, 1))
    tc_g1 = float(stats.t.ppf(0.975, G - 1))
    tc_hw = float(stats.t.ppf(0.975, max(2, G - 2)))  # fixed-G t(G - r) style adjustment

    # WCR-score Webb6 & Rademacher (B=399)
    pw_w6 = app.wcr_score_p_tpk0(X, eh, tpk_true, 399, seed + 100_000, "webb6")
    pw_rad = app.wcr_score_p_tpk0(X, eh, tpk_true, 399, seed + 200_000, "rademacher")

    # Delete-one-cluster jackknife (unified dropping any isolated unit)
    A_full = (geo.W > 0).astype(float).tocsr()
    jk_vals = []
    for g in range(G):
        keep = geo.cl != g
        Ak = A_full[keep][:, keep]
        deg = np.asarray(Ak.sum(1)).ravel()
        sub = deg > 0
        if not sub.all():
            idx = np.where(keep)[0][sub]
            km = np.zeros(N, dtype=bool); km[idx] = True
            Ak = A_full[km][:, km]; deg = np.asarray(Ak.sum(1)).ravel()
        else:
            km = keep
        Wk = sps.diags(1.0 / deg) @ Ak
        _, cl_inv = np.unique(geo.cl[km], return_inverse=True)
        ps_g = dgp.build_panel(sim["Yall"][:, km], sim["Dall"][:, km], Wk.tocsr(), sim["omega"][km], cl_inv)
        Xg = QDX(ps_g, {"tau": ps_g["D"], "kappa": ps_g["Db"]}, phi_bar=0.5)
        eh_g, ok_g, _ = Xg.fit(starts=[[-0.2, -0.1, 0.4], [0.0, 0.0, 0.4]])
        if ok_g:
            jk_vals.append(float(eh_g[0] + eh_g[1]))
    if len(jk_vals) == G:
        jk_arr = np.array(jk_vals)
        jk_se = float(np.sqrt(((G - 1.0) / G) * np.sum((jk_arr - jk_arr.mean()) ** 2)))
        cov_jk = bool(abs(tpk - tpk_true) <= tc_g1 * jk_se)
    else:
        jk_se = np.nan
        cov_jk = False

    return dict(
        tpk=tpk, se=se_tpk, jk_se=jk_se,
        cov_tG1=bool(abs(tpk - tpk_true) <= tc_g1 * se_tpk),
        cov_hwG2=bool(abs(tpk - tpk_true) <= tc_hw * se_tpk),
        cov_jk=cov_jk, jk_ok=(len(jk_vals) == G),
        cov_wcr_w6=bool(pw_w6 > 0.05) if np.isfinite(pw_w6) else False,
        cov_wcr_rad=bool(pw_rad > 0.05) if np.isfinite(pw_rad) else False,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=400)
    ap.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 4, 18)))
    args = ap.parse_args()

    t0 = time.time()
    cells = [
        # Grid A: fixed N=2400, T=12, vary G purely across all divisors (feedback adoption)
        ("G10_N2400_T12_bal", 2400, 10, 12, False, "feedback", 2028010000),
        ("G16_N2400_T12_bal", 2400, 16, 12, False, "feedback", 2028020000),
        ("G20_N2400_T12_bal", 2400, 20, 12, False, "feedback", 2028030000),
        ("G30_N2400_T12_bal", 2400, 30, 12, False, "feedback", 2028040000),
        ("G50_N2400_T12_bal", 2400, 50, 12, False, "feedback", 2028050000),
        ("G100_N2400_T12_bal", 2400, 100, 12, False, "feedback", 2028060000),
        # Grid B: fixed N=400, G=16 -> balanced vs German-skewed x (T=12 vs T=60) under staggered random adoption
        ("G16_N400_T12_bal", 400, 16, 12, False, "random", 2028110000),
        ("G16_N400_T60_bal", 400, 16, 60, False, "random", 2028120000),
        ("G16_N400_T12_skew", 400, 16, 12, True, "random", 2028130000),
        ("G16_N400_T60_skew", 400, 16, 60, True, "random", 2028140000),
    ]
    out_cells = {}
    tpk_true = float(np.log(0.8)) - 0.10
    for nm, N, G, T, sk, adm, sb in cells:
        tasks = [(nm, N, G, T, sk, adm, sb, r) for r in range(args.R)]
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            res = [x for x in ex.map(_run_cell_worker, tasks, chunksize=5) if x is not None]
        tpks = np.array([x["tpk"] for x in res])
        ses = np.array([x["se"] for x in res])
        jk_ses = np.array([x["jk_se"] for x in res if np.isfinite(x["jk_se"])])
        sd = float(tpks.std(ddof=1))
        out_cells[nm] = dict(
            N=N, G=G, T=T, skewed=sk, adopt=adm, n_ok=len(res),
            bias=float((tpks - tpk_true).mean()),
            empSD=sd, meanSE=float(ses.mean()), medSE=float(np.median(ses)),
            iqrSD=float((np.quantile(tpks, 0.75) - np.quantile(tpks, 0.25)) / 1.349),
            se_ratio=float(ses.mean() / sd),
            jk_se_ratio=float(jk_ses.mean() / sd) if len(jk_ses) else np.nan,
            cov_tG1=float(np.mean([x["cov_tG1"] for x in res])),
            cov_hwG2=float(np.mean([x["cov_hwG2"] for x in res])),
            cov_jk=float(np.mean([x["cov_jk"] for x in res if x["jk_ok"]])),
            cov_wcr_w6=float(np.mean([x["cov_wcr_w6"] for x in res])),
            cov_wcr_rad=float(np.mean([x["cov_wcr_rad"] for x in res])),
        )
        c = out_cells[nm]
        print(
            f"{nm:<20} | n={c['n_ok']} | se_ratio={c['se_ratio']:.3f} | t(G-1)={c['cov_tG1']:.3f} "
            f"| t(G-2)={c['cov_hwG2']:.3f} | JK={c['cov_jk']:.3f}({c['jk_se_ratio']:.3f}) "
            f"| WCR_W6={c['cov_wcr_w6']:.3f} | WCR_Rad={c['cov_wcr_rad']:.3f}",
            flush=True,
        )

    out = dict(R=args.R, seconds=round(time.time() - t0, 2), cells=out_cells)
    OUT_PATH.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Saved {OUT_PATH} in {out['seconds']}s")


if __name__ == "__main__":
    main()
