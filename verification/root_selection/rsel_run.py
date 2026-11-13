"""根选择规则验证的最小试运行。

对 rsel_config.CELLS 中的每个格子，生成数据、拟合 M0（局部解）与 M0-GI（全局 + inv 选择），
记录 τ 的估计与区间、φ 的估计、是否出现多根、M0-GI 是否选到边界。

用法：
  python rsel_run.py --R 200 --workers 8 --code <repo>/experiment/code [--cells 0,1,2]
"""
from __future__ import annotations

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    __import__("os").environ.setdefault(_v, "1")

import argparse
import json
import os
import sys
import time
from multiprocessing import get_context

import numpy as np
import rsel_config as RC


def worker(args):
    cell_i, rep, code_dir = args
    sys.path.insert(0, code_dir)
    import dgp
    import qdgmm
    import simulate as S

    cell = RC.build()[101 + cell_i]
    p = cell["params"]
    geo = dgp.make_geography(p.N, G=int(50 * cell["g_scale"]))
    seed = RC.seed_for(cell["sid"], rep)
    sim = dgp.simulate(p, geo, [seed, 0], [seed, 1])
    panel = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl, w=(0.85, 0.15))
    out = dict(cell=cell_i, sid=cell["sid"], rep=rep)
    try:
        m0 = qdgmm.fit(panel, qdgmm.spec_m0(panel, cell["phi_bar"]))
        out.update(ok0=True, tau0=float(m0["est"]["tau"]), se0=float(m0["se"]["tau"]), phi0=float(m0["est"]["phi"]))
    except Exception as e:
        out.update(ok0=False, err0=type(e).__name__)
    try:
        mg = qdgmm.fit_global(panel, phi_bar=cell["phi_bar"], selection="inv")
        out.update(okgi=True, tau_gi=float(mg["est"]["tau"]), se_gi=float(mg["se"]["tau"]),
                   phi_gi=float(mg["est"]["phi"]), at_boundary=bool(mg.get("at_boundary")),
                   n_roots=mg.get("n_roots"))
    except Exception as e:
        out.update(okgi=False, errgi=type(e).__name__)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=200)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--code", required=True)
    ap.add_argument("--cells", default="")
    ap.add_argument("--out", default="rsel_rows.jsonl")
    a = ap.parse_args()
    cells = [int(x) for x in a.cells.split(",")] if a.cells else list(range(len(RC.CELLS)))
    tasks = [(c, r, os.path.abspath(a.code)) for c in cells for r in range(a.R)]
    t0 = time.time()
    with get_context("spawn").Pool(a.workers) as pool, open(a.out, "w", encoding="utf-8") as fh:
        for i, row in enumerate(pool.imap_unordered(worker, tasks, chunksize=8), 1):
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")
            if i % 200 == 0:
                print(f"{i}/{len(tasks)}  {time.time()-t0:.0f}s", flush=True)
    print("done", len(tasks), f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
