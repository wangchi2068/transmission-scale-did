"""M2 simulation: same-week propagation (w0 = 0.3) against its matched benchmark (w0 = 0).

Both arms use the same seeds, the same DGP parameters (core scenario N=500, T=12, G=50, phi=0.4, kappa=-0.1,
tau=log 0.8, NB k=10, feedback adoption, alpha scaling as in experiment/code/scenarios.py), the same
estimator (M0, phi free, phi_bar = 0.5) and the same J test (M0-OI two-step, Hansen J).  At w0 = 0 the DGP is
bit-identical to experiment/code/dgp.simulate, so the benchmark arm is the main simulation.

Rows: one JSON line per replication with the estimates, standard errors, convergence flags and J for both
arms (arms = dict).  Resumable: existing (rep) rows are skipped.

Usage (from this directory):
  python -I sim.py --R 500 --workers 16            # full run  -> rows.jsonl
  python -I sim.py --R 3 --workers 2 --out smoke.jsonl   # smoke test
"""
from __future__ import annotations

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    __import__("os").environ.setdefault(_v, "1")

import argparse
import json
import math
import os
import sys
import time
import traceback
from multiprocessing import get_context

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, "..", "..", "experiment", "code"))
sys.path.insert(0, HERE)
sys.path.insert(0, CODE)

import dgp  # noqa: E402
import qdgmm  # noqa: E402
import same_week_dgp as sw  # noqa: E402
import scenarios as sc  # noqa: E402

SID = 1                         # core_N500_T12_fb_spill: seeds 20261112 + 1e5*1 + rep as in the main study
ARMS = (("benchmark", 0.0), ("same_week", 0.3))
_GEO = {}


def geography():
    if "geo" not in _GEO:
        p = sc.build_scenarios()[SID].params
        _GEO["geo"] = dgp.make_geography(p.N)
        _GEO["p"] = p
    return _GEO["geo"], _GEO["p"]


def _num(x):
    x = float(x)
    return x if math.isfinite(x) else None


def run_rep(rep: int) -> dict:
    geo, p = geography()
    seed = sc.seed_for(SID, rep)
    row = dict(rep=rep, seed=seed, timestamp=time.time())
    for arm, w0 in ARMS:
        a = dict(w0=w0)
        t0 = time.perf_counter()
        try:
            sim = sw.simulate_same_week(p, geo, [seed, 0], [seed, 1], w0)
            panel = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl, w=(0.85, 0.15))
            a["desc"] = dgp.panel_descriptives(panel)
            a["desc"]["mean_Y_last"] = float(panel["Y"][-1].mean())
            m0 = qdgmm.fit(panel, qdgmm.spec_m0(panel, phi_bar=0.5))
            a["ok"] = bool(m0.get("ok"))
            a["beta"] = {k: _num(v) for k, v in m0.get("est", {}).items()}
            a["se"] = {k: _num(v) for k, v in m0.get("se", {}).items()}
            if m0.get("Vb") is not None:
                Vb = np.asarray(m0["Vb"], float)
                a["Vb"] = [[_num(v) for v in r_] for r_ in Vb]
            try:
                moi = qdgmm.fit_two_step(panel, qdgmm.spec_m0_oi(panel, 0.5), m0)
                a["J"] = _num(moi.get("J"))
                a["J_df"] = int(moi.get("J_df", 0))
                a["J_p"] = _num(moi.get("J_p"))
                a["J_ok"] = bool(moi.get("ok"))
            except Exception:  # noqa: BLE001 - J failure is recorded per replication
                a["J"] = a["J_p"] = None
                a["J_ok"] = False
        except Exception as e:  # noqa: BLE001 - a failed replication is recorded, never fatal
            a.update(ok=False, beta={}, se={}, J=None, J_p=None, J_ok=False,
                     error="%s: %s" % (type(e).__name__, str(e)[:200]))
        a["sec"] = time.perf_counter() - t0
        row[arm] = a
    return row


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        o = float(o)
        return o if math.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=500)
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--out", default=os.path.join(HERE, "rows.jsonl"))
    ap.add_argument("--block", type=int, default=25)
    a = ap.parse_args()
    done = set()
    if os.path.exists(a.out):
        with open(a.out, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(int(json.loads(line)["rep"]))
                except (ValueError, KeyError):
                    pass
    todo = [r for r in range(a.R) if r not in done]
    print("R=%d done=%d todo=%d workers=%d g=%.3f" % (a.R, len(done), len(todo), a.workers, sw.CAL["g"]),
          flush=True)
    t0 = time.time()
    with get_context("spawn").Pool(a.workers) as pool, open(a.out, "a", encoding="utf-8") as fh:
        for i, row in enumerate(pool.imap_unordered(run_rep, todo, chunksize=max(1, a.block // a.workers)), 1):
            fh.write(json.dumps(_clean(row)) + "\n")
            if i % a.block == 0 or i == len(todo):
                fh.flush()
                print("%d/%d  %.0fs" % (i, len(todo), time.time() - t0), flush=True)
    print("done %.0fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
