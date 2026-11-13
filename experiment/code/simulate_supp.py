"""Supplementary Monte Carlo runner for M0-GI (post-hoc addition in run-sim-application; see
experiment_design.md, deviation 7). M0-GI = M0-G with the root-selection criterion built from the
phi-instruments at other reference values (qdgmm.fit_global(selection="inv")).

The data sets are *identical* to the main run (src/simulate.py): replication r of scenario sid is regenerated
from the same seeds (scenarios.seed_for, streams [seed, 0] / [seed, 1]) and the same panel construction, and
only M0-GI (+ cases averted psi) is fitted. Rows (kind = "supp") are written to sim/rows_m0gi.jsonl, resumable,
and src/aggregate.py merges them into the main replication rows by (sid, rep); only replications present in
both files enter the M0-GI performance measures (n_ok is reported).

Usage:  python src/simulate_supp.py --R 1000 --workers 8
"""
from __future__ import annotations

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    __import__("os").environ.setdefault(_v, "1")

import argparse
import json
import os
import sys
import time
import traceback
from multiprocessing import get_context

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dgp  # noqa: E402
import psi as psimod  # noqa: E402
import qdgmm  # noqa: E402
import scenarios as sc  # noqa: E402
import simulate as S  # noqa: E402

RUN = S.RUN
VERSION = "1.0.0-supp-m0gi"


def run_rep_supp(sid: int, rep: int) -> dict:
    scen = S.scenario_table()[sid]
    p = scen.params
    geo = S.geography(p.N)
    seed = sc.seed_for(sid, rep)
    sim = dgp.simulate(p, geo, [seed, 0], [seed, 1])
    panel = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl, w=scen.w_est)
    row = dict(kind="supp", sid=sid, rep=rep, seed=seed, version=VERSION,
               max_Y=float(panel["Y"].max()))   # checked against the main row by aggregate.py (same data)

    def idx0(b):
        return b[0] * panel["D"] + b[1] * panel["Db"], b[2]
    mg, t = S._timed(qdgmm.fit_global, panel, selection="inv")
    res = S._slim(mg)
    if mg.get("ok"):
        ps, _ = S._timed(psimod.psi_with_se, sim["Yall"], panel, mg["beta"], mg["Vb"], idx0)
        res["psi"] = ps
    row["res"] = {"M0_GI": res}
    row["time"] = {"M0_GI": t}
    return S._clean(row)


def _work(task):
    t0 = time.perf_counter()
    try:
        _, sid, reps = task
        return task, [run_rep_supp(sid, r) for r in reps], time.perf_counter() - t0, None
    except Exception:  # noqa: BLE001
        return task, [], time.perf_counter() - t0, traceback.format_exc()[-2000:]


def read_done(path):
    done = set()
    if not os.path.exists(path):
        return done
    good, bad = [], 0
    with open(path) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                bad += 1
                continue
            good.append(line if line.endswith("\n") else line + "\n")
            done.add((r["sid"], r["rep"]))
    if bad:
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            fh.writelines(good)
            fh.flush(); os.fsync(fh.fileno())
        os.replace(tmp, path)
    return done


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--R", type=int, default=1000)
    ap.add_argument("--first-wave", type=int, default=200)
    ap.add_argument("--block", type=int, default=10)
    ap.add_argument("--scenarios", default="all")
    ap.add_argument("--out", default=os.path.join(RUN, "sim", "rows_m0gi.jsonl"))
    ap.add_argument("--workers", type=int, default=sc.default_workers())
    ap.add_argument("--max-minutes", type=float, default=float("inf"))
    a = ap.parse_args(argv)
    scen = S.scenario_table()
    order = sc.run_order(scen)
    sids = order if a.scenarios == "all" else [s for s in order if s in {int(x) for x in a.scenarios.split(",")}]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    done = read_done(a.out)
    tasks = [t for t in S.build_tasks(sids, a.R, a.first_wave, a.block, 0, 1, set(), set())]
    tasks = [("rep", sid, [r for r in reps if (sid, r) not in done]) for _, sid, reps in tasks]
    tasks = [t for t in tasks if t[2]]
    with open(os.path.splitext(a.out)[0] + "_meta.json", "w") as fh:
        json.dump(dict(version=VERSION, R=a.R, first_wave=a.first_wave, block=a.block, scenarios=sids,
                       selection="inv", refs=list(qdgmm.SEL_PHI_REFS), base_seed=sc.BASE_SEED), fh, indent=1)
    print("[supp] %d reps done; %d tasks to run on %d workers" % (len(done), len(tasks), a.workers), flush=True)
    if not tasks:
        return 0
    t_start, n_done, n_err = time.time(), 0, 0
    ctx = get_context("fork" if sys.platform.startswith("linux") else "spawn")
    with open(a.out, "a") as fh, ctx.Pool(a.workers, maxtasksperchild=200) as pool:
        for task, rows, secs, err in pool.imap_unordered(_work, tasks, chunksize=1):
            if err:
                n_err += 1
                sys.stderr.write("[supp] task %s failed:\n%s\n" % (str(task[:2]), err))
            for r in rows:
                fh.write(json.dumps(r) + "\n")
            fh.flush(); os.fsync(fh.fileno())
            n_done += 1
            if n_done % 20 == 0 or n_done == len(tasks):
                el = time.time() - t_start
                print("[supp] %d/%d tasks, %.1f min elapsed, eta %.1f min" %
                      (n_done, len(tasks), el / 60, el / n_done * (len(tasks) - n_done) / 60), flush=True)
            if (time.time() - t_start) / 60 > a.max_minutes:
                print("[supp] --max-minutes reached; stopping", flush=True)
                pool.terminate()
                break
    print("[supp] finished: %d tasks, %d failed, %.1f min" % (n_done, n_err, (time.time() - t_start) / 60), flush=True)
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main())
