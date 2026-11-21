"""Weak-prior and ultra-low-prior EpiEstim log-R comparator (M2) across all 25 scenarios (R = 8000).

Evaluates M2 with Gamma(a=0.1, b=5) and Gamma(a=0.01, b=100) priors on the identical
DGP seeds used in experiment/code/simulate.py (Major Comment 1(d)).
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
import scenarios as sc  # noqa: E402

OUT_PATH = Path(__file__).resolve().parent / "results.json"
_GEO = {}
_SCEN = None


def _init_worker():
    global _SCEN
    _SCEN = sc.build_scenarios()
    for n in (500, 2000):
        _GEO[n] = dgp.make_geography(n)


def _run_block(args):
    sid, reps = args
    scen = _SCEN[sid]
    geo = _GEO[scen.params.N]
    t0 = scen.params.tau
    rows = []
    for r in reps:
        seed = sc.seed_for(sid, r)
        sim = dgp.simulate(scen.params, geo, [seed, 0], [seed, 1])
        pan = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl, w=scen.w_est)
        tc = float(stats.t.ppf(0.975, int(pan["G"]) - 1))
        out_r = {}
        for lbl, a_val, b_val in [("M2_def", 1.0, 5.0), ("M2_weak", 0.1, 5.0), ("M2_ultra", 0.01, 100.0)]:
            fit = cmp.twfe_ols(cmp.cori_log_r(pan, a=a_val, b=b_val), pan["D"], pan["cl"], int(pan["G"]))
            if fit.get("ok"):
                est = float(fit["est"]["tau"])
                se = float(fit["se"]["tau"])
                out_r[lbl] = (est, se, abs(est - t0) <= tc * se)
            else:
                out_r[lbl] = None
        rows.append(out_r)
    return sid, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=8000)
    ap.add_argument("--workers", type=int, default=max(1, min(os.cpu_count() or 4, 18)))
    args = ap.parse_args()

    t0 = time.time()
    scen = sc.build_scenarios()
    _init_worker()
    blocks = []
    for sid in sorted(scen):
        for start in range(0, args.R, 100):
            blocks.append((sid, list(range(start, min(start + 100, args.R)))))

    agg = {sid: {m: [] for m in ("M2_def", "M2_weak", "M2_ultra")} for sid in scen}
    with ProcessPoolExecutor(max_workers=args.workers, initializer=_init_worker) as ex:
        for sid, rows in ex.map(_run_block, blocks, chunksize=1):
            for r in rows:
                for m, val in r.items():
                    if val is not None:
                        agg[sid][m].append(val)

    summary = {}
    for sid in sorted(scen):
        s = scen[sid]
        t_true = float(np.mean(s.params.dyn)) if s.params.dyn is not None else float(s.params.tau)
        rec = dict(sid=sid, scenario=s.name, block=s.block, N=s.params.N, T=s.params.T,
                   adopt=s.params.adopt, spill=(s.params.phi > 0), truth=t_true)
        for m in ("M2_def", "M2_weak", "M2_ultra"):
            arr = np.array(agg[sid][m], float) if agg[sid][m] else np.zeros((0, 3))
            if len(arr):
                ests, ses, covs = arr[:, 0], arr[:, 1], arr[:, 2]
                b = ests - t_true
                rec[m] = dict(
                    n=int(len(arr)),
                    bias=float(b.mean()),
                    bias_mcse=float(b.std(ddof=1) / np.sqrt(len(arr))),
                    empSE=float(ests.std(ddof=1)),
                    modelSE=float(ses.mean()),
                    rmse=float(np.sqrt((b ** 2).mean())),
                    coverage=float(covs.mean()),
                    cov_mcse=float(np.sqrt(covs.mean() * (1 - covs.mean()) / len(arr))),
                )
            else:
                rec[m] = None
        summary[str(sid)] = rec

    out = dict(R=args.R, seconds=round(time.time() - t0, 2), summary=summary)
    OUT_PATH.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Saved {OUT_PATH} in {out['seconds']}s")


if __name__ == "__main__":
    main()
