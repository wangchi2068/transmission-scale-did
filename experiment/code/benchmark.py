"""Timing / memory benchmark for the Monte Carlo study and the choice of R (writes run_plan.json).

For every scenario a fresh worker process runs --reps replications (rep indices 99_000+, disjoint from the
study's reps; output discarded - these are NOT results) and one psi-truth block of --psi-size pairs, and
reports seconds per replication, seconds per psi pair and the worker's peak RSS. Workers run --workers at a
time so that timings include the memory-bandwidth contention of a loaded machine.

R is then chosen for the target machine (default 16 cores x 4 h) as the largest value on the grid
1000, 1500, ..., 8000 whose projected wall time (core-seconds / cores / efficiency / speed factor) is at most
--headroom (default 0.75) of the budget; the remaining quarter absorbs slower target cores and rep-to-rep
variation in cost (the benchmark times only a few reps). The grid stops at 8000 because aggregate.py holds
all rows in memory (~18 kB per replication row, measured: 12,000 rows -> 313 MB peak), i.e. ~3.6 GB at
R = 8000 x 25 scenarios. If even R = 1000 does not fit, the largest R in steps of 100 that fits is
reported together with the coverage MCSE it implies, so the shortfall is explicit rather than hidden.

Usage: python src/benchmark.py [--reps 3] [--psi-size 10] [--workers 4] [--out run_plan.json]
"""
from __future__ import annotations

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    __import__("os").environ.setdefault(_v, "1")

import argparse
import json
import math
import os
import platform
import resource
import subprocess
import sys
import time
from multiprocessing import get_context

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simulate  # noqa: E402
import scenarios as sc  # noqa: E402

RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REP0 = 99_000


def _bench(args):
    sid, reps, psi_size = args
    t0 = time.perf_counter()
    simulate.run_rep(sid, REP0 - 1)                     # warm-up: geography / imports (not timed)
    t_warm = time.perf_counter() - t0
    times, fails = [], 0
    for r in range(reps):
        t = time.perf_counter()
        row = simulate.run_rep(sid, REP0 + r)
        times.append(time.perf_counter() - t)
        fails += sum(1 for v in row["res"].values() if not v.get("ok"))
    t = time.perf_counter()
    simulate.run_psi_block(sid, 9_999, psi_size)
    t_psi = (time.perf_counter() - t) / psi_size
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0     # Linux: kB -> MB
    return dict(sid=sid, sec_per_rep=sum(times) / len(times), sec_rep_max=max(times), sec_warmup=t_warm,
                sec_per_psi_pair=t_psi, peak_rss_mb=rss, method_failures_in_bench=fails, reps_timed=reps)


AGG_KB_PER_ROW = 18.0


def choose_R(per_scen, psi_pairs, cores, hours, eff, speed, headroom=0.75, grid=tuple(range(1000, 8001, 500))):
    budget = cores * hours * 3600.0 * eff * headroom
    rep_cost = sum(s["sec_per_rep"] for s in per_scen) / speed          # core-s per replication of all scen
    psi_cost = sum(s["sec_per_psi_pair"] for s in per_scen) * psi_pairs / speed

    def wall(R):
        return (R * rep_cost + psi_cost) / (cores * eff)
    fits = [R for R in grid if wall(R) <= headroom * hours * 3600]
    if fits:
        R = max(fits)
        why = ("largest R on the grid 1000..8000 (step 500) whose projected wall time is <= %.0f%% of %.1f h on "
               "%d cores" % (100 * headroom, hours, cores))
    else:
        R = max(100, int((budget - psi_cost) / rep_cost // 100 * 100))
        why = ("R = 1000 does not fit; largest multiple of 100 that fits. Coverage MCSE at 95%% nominal = %.4f "
               "(vs 0.0069 at R = 1000)" % math.sqrt(0.95 * 0.05 / R))
    return R, why, dict(core_sec_per_rep_all_scenarios=rep_cost, core_sec_psi_truth=psi_cost,
                        projected_wall_hours=wall(R) / 3600, usable_core_hours=budget / 3600,
                        projected_aggregate_mem_gib=R * len(per_scen) * AGG_KB_PER_ROW / 2 ** 20 + 0.2)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--psi-size", type=int, default=10)
    ap.add_argument("--workers", type=int, default=sc.default_workers())
    ap.add_argument("--scenarios", default="all")
    ap.add_argument("--target-cores", type=int, default=16)
    ap.add_argument("--target-hours", type=float, default=4.0)
    ap.add_argument("--target-mem-gib", type=float, default=32.0)
    ap.add_argument("--efficiency", type=float, default=0.85, help="parallel efficiency assumed on target")
    ap.add_argument("--speed", type=float, default=1.0, help="target core speed relative to this machine")
    ap.add_argument("--headroom", type=float, default=0.75, help="fraction of the time budget to plan for")
    ap.add_argument("--psi-pairs", type=int, default=5000)
    ap.add_argument("--out", default=os.path.join(RUN, "run_plan.json"))
    a = ap.parse_args(argv)
    scen = simulate.scenario_table()
    order = sc.run_order(scen)
    sids = order if a.scenarios == "all" else [s for s in order if s in {int(x) for x in a.scenarios.split(",")}]
    t0 = time.time()
    ctx = get_context("fork" if sys.platform.startswith("linux") else "spawn")
    with ctx.Pool(a.workers, maxtasksperchild=1) as pool:               # fresh process -> per-scenario RSS
        res = []
        for r in pool.imap_unordered(_bench, [(s, a.reps, a.psi_size) for s in sids]):
            res.append(r)
            print("[bench] sid %2d  %.2f s/rep (max %.2f)  %.3f s/psi pair  %.0f MB" %
                  (r["sid"], r["sec_per_rep"], r["sec_rep_max"], r["sec_per_psi_pair"], r["peak_rss_mb"]),
                  flush=True)
    res.sort(key=lambda r: r["sid"])
    for r in res:
        r["name"] = scen[r["sid"]].name
    out = dict(task="epi-novel-method", created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               note="benchmark timings only - not results", bench_wall_sec=round(time.time() - t0, 1),
               machine=dict(cpus=os.cpu_count(), workers=a.workers, python=platform.python_version(),
                            processor=platform.processor() or platform.machine()),
               per_scenario=res)
    if a.scenarios == "all":
        R, why, cost = choose_R(res, a.psi_pairs, a.target_cores, a.target_hours, a.efficiency, a.speed,
                                a.headroom)
        peak = max(r["peak_rss_mb"] for r in res)
        mem = peak * a.target_cores / 1024 + 0.5                         # + parent / OS headroom
        out["plan"] = dict(
            target=dict(cores=a.target_cores, hours=a.target_hours, mem_gib=a.target_mem_gib,
                        efficiency=a.efficiency, speed_factor=a.speed, headroom=a.headroom),
            R=R, reason=why, first_wave=min(200, R), block=10, psi_pairs=a.psi_pairs, psi_block=250,
            mcse_coverage_at_R=math.sqrt(0.95 * 0.05 / R), projected_peak_mem_gib=mem,
            memory_ok=bool(max(mem, cost["projected_aggregate_mem_gib"]) < 0.8 * a.target_mem_gib), workers=a.target_cores, **cost,
            command=("QDG_WORKERS=%d bash run_all.sh   # == python src/simulate.py --R %d --first-wave %d "
                     "--block 10 --psi-pairs %d --psi-block 250" % (a.target_cores, R, min(200, R), a.psi_pairs)),
            priority="waves: every scenario to %d reps + psi truth first (run_order: core feedback/spill -> "
                     "core no-spill -> core random -> null -> misspec), then all to 1000, then all to R" % min(200, R),
            resume="rerun the same command; finished (scenario, rep) and psi blocks in sim/rows.jsonl are skipped")
    try:
        out["git_sha"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=RUN, capture_output=True,
                                        text=True).stdout.strip() or None
    except OSError:
        out["git_sha"] = None
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(out, fh, indent=1)
    os.replace(tmp, a.out)
    if "plan" in out:
        p = out["plan"]
        print("[bench] plan: R = %d (%s); projected %.2f h, peak memory %.1f GiB" %
              (p["R"], p["reason"], p["projected_wall_hours"], p["projected_peak_mem_gib"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
