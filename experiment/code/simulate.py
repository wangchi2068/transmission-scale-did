"""Monte Carlo runner (theory/simulation_spec.md).

Usage (from the experiment directory):
    python src/simulate.py --R 1000                      # all 25 scenarios, resumable
    python src/simulate.py --R 4 --scenarios 1,17 --out smoke/rows.jsonl --psi-pairs 20   # smoke test

Design
  * Work unit = (scenario, rep-block) of --block replications (default 10), or a psi-truth block of
    --psi-block counterfactual pairs. Workers compute; the parent process is the single writer and appends
    one JSON line per replication (kind = "rep") / per psi block (kind = "psi_truth") immediately after
    each work unit finishes, followed by flush + fsync.
  * Resume: on start the existing rows file is read (a truncated last line from a kill is ignored and
    rewritten), finished (sid, rep) and (sid, psi block) pairs are skipped.
  * Priority: waves. Wave 1 brings every scenario (in scenarios.run_order priority) to --first-wave reps
    (default 200) plus its psi truth; wave 2 brings all scenarios to 1000 reps (the spec minimum); wave 3
    completes all to --R (run_plan.json). A run cut short therefore has every scenario at >= 200 reps
    (coverage MCSE <= 0.016), then at >= 1000 (<= 0.0069), before any scenario goes beyond.
  * Seeds: rep r of scenario sid uses seed = 20261112 + sid*1e5 + r with two streams [seed, 0]
    (structure) and [seed, 1] (innovations); psi truth pair j uses seed + 9e8 analogously. Results do not
    depend on the number of workers or the order of execution.
  * Parallelism: QDG_WORKERS (default: usable cores, see scenarios.default_workers) processes, OMP/MKL/OpenBLAS threads = 1.
    Memory per worker < 300 MB at N = 2000, T = 26 (measured by benchmark.py), i.e. < 5 GB on 16 workers.
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comparators as cmp  # noqa: E402
import dgp  # noqa: E402
import psi as psimod  # noqa: E402
import qdgmm  # noqa: E402
import scenarios as sc  # noqa: E402

RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
METHODS = ("M0", "M0_G", "M0_ET", "M0_OI", "M0p", "M1", "M1_BC", "M2", "M3", "M4", "Oracle")
VERSION = "1.0.0"
MIN_R = 1000                       # spec minimum; an intermediate wave so every scenario reaches it first

_GEO = {}
_SCEN = None


def geography(N):
    if N not in _GEO:
        _GEO[N] = dgp.make_geography(N)
    return _GEO[N]


def scenario_table():
    global _SCEN
    if _SCEN is None:
        _SCEN = sc.build_scenarios()
    return _SCEN


def _clean(x):
    """JSON-safe floats (NaN/inf -> None)."""
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, (np.floating, float)):
        x = float(x)
        return x if math.isfinite(x) else None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, np.ndarray):
        return _clean(x.tolist())
    return x


def _slim(res: dict, keep=("ok", "est", "se", "reason", "start_used", "resid", "cond", "iters", "J", "J_df", "J_p",
                           "n_roots", "roots", "at_boundary", "chosen", "Q")):
    return {k: res[k] for k in keep if k in res}


def _timed(fn, *a, **k):
    t0 = time.perf_counter()
    try:
        out = fn(*a, **k)
    except Exception as e:  # noqa: BLE001 - a failure of one method is recorded, never fatal
        out = dict(ok=False, est={}, se={}, reason="exception: %s" % (repr(e)[:200]))
    return out, time.perf_counter() - t0


# ------------------------------------------------------------------------------------------------
# One replication
# ------------------------------------------------------------------------------------------------
def true_index(sim, p: dgp.DGPParams, panel_true: dict) -> np.ndarray:
    return sim["eff"] + p.kappa * panel_true["Db"]


def x3_weights(sim, p, panel_est, panel_true, phi_bar=0.5):
    """Per-replication sums entering Proposition 3's first-order pseudo-true value (theory.md section 7)."""
    D = panel_est["D"]
    lam_bar = panel_est["L"] + phi_bar * panel_est["Nb"] + panel_est["omega"]
    s = 1.0 / (lam_bar[:-1] + 1.0)
    lam0 = panel_true["L"] + p.phi * panel_true["Nb"] + panel_true["omega"]
    mu = np.exp(sim["alpha"][None, :] + sim["gamma"][:, None] + sim["eff"] + p.kappa * panel_true["Db"]) * lam0
    sw = D[1:] - D[:-1]
    ti = sim["tau_i"][None, :]
    A = (s * D[:-1] * mu[:-1]).sum(1)
    M = (s * mu[:-1]).sum(1)
    K = (s * mu[:-1] * sw).sum(1)
    H = (s * mu[:-1] * sw * ti).sum(1)
    switchers = panel_est["ev"] >= 1
    return dict(A=A, M=M, K=K, H=H, tau_switch_mean=float(sim["tau_i"][switchers].mean()) if switchers.any() else None,
                n_switchers=int(switchers.sum()))


def run_rep(sid: int, rep: int) -> dict:
    scen = scenario_table()[sid]
    p = scen.params
    geo = geography(p.N)
    seed = sc.seed_for(sid, rep)
    sim = dgp.simulate(p, geo, [seed, 0], [seed, 1])
    panel = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl, w=scen.w_est)
    panel_true = dgp.build_panel(sim["Iall"], sim["Dall"], geo.W, sim["omega"], geo.cl, w=p.w_true)
    T, N = panel["Y"].shape
    row = dict(kind="rep", sid=sid, rep=rep, seed=seed, version=VERSION)
    desc = dgp.panel_descriptives(panel)
    desc["max_Y"] = float(panel["Y"].max())
    row["desc"] = desc
    res, tm = {}, {}

    # M0 and its variants
    sp0 = qdgmm.spec_m0(panel)
    m0, tm["M0"] = _timed(qdgmm.fit, panel, sp0)
    res["M0"] = _slim(m0)
    def idx0(b):
        return b[0] * panel["D"] + b[1] * panel["Db"], b[2]
    if m0.get("ok"):
        t0 = time.perf_counter()
        ps, _ = _timed(psimod.psi_with_se, sim["Yall"], panel, m0["beta"], m0["Vb"], idx0)
        res["M0"]["psi"] = ps
        tm["M0_psi"] = time.perf_counter() - t0
    mg, tm["M0_G"] = _timed(qdgmm.fit_global, panel)
    res["M0_G"] = _slim(mg)
    if mg.get("ok"):
        ps, _ = _timed(psimod.psi_with_se, sim["Yall"], panel, mg["beta"], mg["Vb"], idx0)
        res["M0_G"]["psi"] = ps
    spet = qdgmm.spec_m0_et(panel)
    met, tm["M0_ET"] = _timed(qdgmm.fit, panel, spet)
    res["M0_ET"] = _slim(met)
    if met.get("ok") and p.dyn is not None:
        bins = [spet.lin[nm] for nm in spet.names if nm.startswith("tau_e")]

        def idxet(b):
            ix = sum(bk * B for bk, B in zip(b[:len(bins)], bins)) + b[len(bins)] * panel["Db"]
            return ix, b[len(bins) + 1]
        ps, _ = _timed(psimod.psi_with_se, sim["Yall"], panel, met["beta"], met["Vb"], idxet)
        res["M0_ET"]["psi"] = ps
    moi, tm["M0_OI"] = _timed(qdgmm.fit_two_step, panel, qdgmm.spec_m0_oi(panel), m0)
    res["M0_OI"] = _slim(moi)
    mp, tm["M0p"] = _timed(qdgmm.fit, panel, qdgmm.spec_m0prime(panel))
    res["M0p"] = _slim(mp)
    # comparators
    m1, tm["M1"] = _timed(cmp.poisson_fe, panel)
    res["M1"] = _slim(m1)
    m1b, tm["M1_BC"] = _timed(cmp.poisson_fe_bc, panel, m1)
    res["M1_BC"] = _slim(m1b)
    m2, tm["M2"] = _timed(cmp.m2_logr_did, panel)
    res["M2"] = _slim(m2)
    m3, tm["M3"] = _timed(cmp.m3_logy_twfe, panel)
    res["M3"] = _slim(m3)
    m4, tm["M4"] = _timed(cmp.poisson_qmle, panel)
    res["M4"] = _slim(m4)
    mo, tm["Oracle"] = _timed(cmp.poisson_qmle, panel, offset=sim["alpha"])
    res["Oracle"] = _slim(mo)
    row["res"] = res
    row["time"] = tm

    # replication-level truths
    truth = dict(tau=p.tau, kappa=p.kappa, phi=p.phi)
    ix = true_index(sim, p, panel_true)
    truth["psi_dagger"] = psimod.psi_values(sim["Yall"], ix, geo.W, sim["omega"], p.w_true, p.phi, panel["D"],
                                            panel["S"])
    if p.tau_sd > 0:
        truth["x3"] = x3_weights(sim, p, panel, panel_true)
    row["truth"] = truth
    return _clean(row)


# ------------------------------------------------------------------------------------------------
# Population psi truth (common random numbers, d = 0 counterfactual)
# ------------------------------------------------------------------------------------------------
def run_psi_block(sid: int, block: int, size: int) -> dict:
    scen = scenario_table()[sid]
    p = scen.params
    geo = geography(p.N)
    S = dgp.S_LAGS
    acc = {k: [0.0, 0.0] for k in ("treated", "network")}
    n = 0
    for j in range(block * size, (block + 1) * size):
        seed = sc.seed_for(sid, j) + sc.PSI_SEED_OFFSET
        s1 = dgp.simulate(p, geo, [seed, 0], [seed, 1])
        s0 = dgp.simulate(p, geo, [seed, 0], [seed, 1], counterfactual=True)
        diff = s0["Yall"][S:] - s1["Yall"][S:]
        for k, m in psimod.psi_sets(s1["Dall"][S:]).items():
            v = float(diff[m].sum())
            acc[k][0] += v
            acc[k][1] += v * v
        n += 1
    return _clean(dict(kind="psi_truth", sid=sid, block=block, n=n, seed_offset=sc.PSI_SEED_OFFSET,
                       sum={k: v[0] for k, v in acc.items()}, sumsq={k: v[1] for k, v in acc.items()},
                       version=VERSION))


def _work(task):
    kind = task[0]
    t0 = time.perf_counter()
    try:
        if kind == "rep":
            _, sid, reps = task
            rows = [run_rep(sid, r) for r in reps]
        else:
            _, sid, block, size = task
            rows = [run_psi_block(sid, block, size)]
        return task, rows, time.perf_counter() - t0, None
    except Exception:  # noqa: BLE001 - reported by the parent, the task is retried on resume
        return task, [], time.perf_counter() - t0, traceback.format_exc()[-2000:]


# ------------------------------------------------------------------------------------------------
# Scheduling, resume, writing
# ------------------------------------------------------------------------------------------------
def read_done(path):
    done_rep, done_psi = set(), set()
    if not os.path.exists(path):
        return done_rep, done_psi, 0
    good, bad = [], 0
    with open(path) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                bad += 1
                continue
            good.append(line if line.endswith("\n") else line + "\n")
            if r.get("kind") == "rep":
                done_rep.add((r["sid"], r["rep"]))
            elif r.get("kind") == "psi_truth":
                done_psi.add((r["sid"], r["block"]))
    if bad:   # a kill during a write leaves a truncated last line: rewrite the file without it
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            fh.writelines(good)
            fh.flush(); os.fsync(fh.fileno())
        os.replace(tmp, path)
    return done_rep, done_psi, bad


def build_tasks(sids, R, first_wave, block, psi_pairs, psi_block, done_rep, done_psi):
    tasks = []
    # wave limits: first_wave (all scenarios + psi truth), then the R = 1000 floor of the spec, then R
    waves = sorted({min(x, R) for x in (first_wave or R, MIN_R, R)})
    lo = 0
    for wi, hi in enumerate(waves):
        for sid in sids:
            for b0 in range(lo, hi, block):
                reps = [r for r in range(b0, min(b0 + block, hi)) if (sid, r) not in done_rep]
                if reps:
                    tasks.append(("rep", sid, reps))
            if wi == 0 and psi_pairs > 0:
                for b in range(int(math.ceil(psi_pairs / psi_block))):
                    if (sid, b) not in done_psi:
                        tasks.append(("psi", sid, b, psi_block))
        lo = hi
    return tasks


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--R", type=int, default=1000)
    ap.add_argument("--first-wave", type=int, default=200)
    ap.add_argument("--block", type=int, default=10)
    ap.add_argument("--scenarios", default="all", help="'all' or comma-separated scenario ids")
    ap.add_argument("--out", default=os.path.join(RUN, "sim", "rows.jsonl"))
    ap.add_argument("--psi-pairs", type=int, default=5000)
    ap.add_argument("--psi-block", type=int, default=250)
    ap.add_argument("--workers", type=int, default=sc.default_workers())
    ap.add_argument("--max-minutes", type=float, default=float("inf"),
                    help="stop scheduling new work after this wall time (finished blocks are kept)")
    a = ap.parse_args(argv)
    scen = scenario_table()
    order = sc.run_order(scen)
    sids = order if a.scenarios == "all" else [s for s in order if s in {int(x) for x in a.scenarios.split(",")}]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    done_rep, done_psi, bad = read_done(a.out)
    tasks = build_tasks(sids, a.R, a.first_wave, a.block, a.psi_pairs, a.psi_block, done_rep, done_psi)
    meta_path = os.path.splitext(a.out)[0] + "_meta.json"
    meta = dict(version=VERSION, R=a.R, first_wave=a.first_wave, block=a.block, scenarios=sids,
                psi_pairs=a.psi_pairs, psi_block=a.psi_block, base_seed=sc.BASE_SEED,
                psi_seed_offset=sc.PSI_SEED_OFFSET, alpha_scale=sc.ALPHA_SCALE,
                scenario_table={s: scen[s].to_dict() for s in sids}, calibration=sc.load_calibration())
    with open(meta_path, "w") as fh:
        json.dump(_clean(meta), fh, indent=1)
    print("[sim] %d reps and %d psi blocks already done (%d truncated lines dropped); %d tasks to run on %d workers"
          % (len(done_rep), len(done_psi), bad, len(tasks), a.workers), flush=True)
    if not tasks:
        return 0
    t_start = time.time()
    n_done, n_err = 0, 0
    ctx = get_context("fork" if sys.platform.startswith("linux") else "spawn")
    with open(a.out, "a") as fh, ctx.Pool(a.workers, maxtasksperchild=200) as pool:
        it = pool.imap_unordered(_work, tasks, chunksize=1)
        for task, rows, secs, err in it:
            if err:
                n_err += 1
                sys.stderr.write("[sim] task %s failed:\n%s\n" % (str(task[:2]), err))
            for r in rows:
                fh.write(json.dumps(r) + "\n")
            fh.flush(); os.fsync(fh.fileno())
            n_done += 1
            if n_done % 20 == 0 or n_done == len(tasks):
                el = time.time() - t_start
                print("[sim] %d/%d tasks, %.1f min elapsed, eta %.1f min" %
                      (n_done, len(tasks), el / 60, el / n_done * (len(tasks) - n_done) / 60), flush=True)
            if (time.time() - t_start) / 60 > a.max_minutes:
                print("[sim] --max-minutes reached; stopping (resume by re-running the same command)", flush=True)
                pool.terminate()
                break
    print("[sim] finished: %d tasks, %d failed, %.1f min" % (n_done, n_err, (time.time() - t_start) / 60), flush=True)
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main())
