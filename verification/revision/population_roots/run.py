"""Analysis h (review Comment 3): is the second ("spurious") root of the profiled phi-equation a population
feature in all zero-contagion designs, not only scenario 7?

For each zero-contagion design (phi0 = kappa0 = 0; scenarios 5, 6, 13, 14 = T 12/26 x feedback/random adoption;
N is overridden, so scenarios 7, 8, 15, 16 differ from these only in N) the panel is simulated with the paper's
DGP (experiment/code/dgp.py, scenarios.py, seeds scenarios.seed_for(sid, rep)) at large N (G = 50 clusters as in
the DGP), and the profiled phi-equation F_phi(phi) (kappa, tau solved given phi; experiment/code/qdgmm.py
_solve_given_phi) is evaluated on a fine phi grid. Sign changes away from phi = 0 that persist as N grows are
population roots, i.e. A6 global identification fails. At each sign change the linearly interpolated (tau, kappa)
is recorded (true tau = log 0.8 = -0.223, kappa = 0).

Run (from this folder):  OMP_NUM_THREADS=1 python3 -I run.py --workers 2 > log.txt 2>&1
"""
import argparse
import dataclasses
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import env  # noqa: E402
import numpy as np  # noqa: E402
import dgp  # noqa: E402
import qdgmm  # noqa: E402
import scenarios as sc  # noqa: E402

HERE = Path(__file__).resolve().parent
SCEN = sc.build_scenarios(sc.load_calibration())
GRID = np.r_[np.arange(-0.10, 0.0, 0.02), np.arange(0.0, 0.2, 0.02), np.arange(0.2, 1.5, 0.05), np.arange(1.5, 3.01, 0.25)]
_GEO = {}


def job(args):
    sid, N, rep, phi_bar = args
    t0 = time.time()
    scen = SCEN[sid]
    p = dataclasses.replace(scen.params, N=N)
    if N not in _GEO:
        _GEO[N] = dgp.make_geography(N)
    geo = _GEO[N]
    seed = sc.seed_for(sid, rep)
    sim = dgp.simulate(p, geo, [seed, 0], [seed, 1])
    panel = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl, w=scen.w_est)
    mod = qdgmm.QDModel(panel, qdgmm.spec_m0(panel, phi_bar))
    x = np.zeros(mod.K)
    prof = []
    for ph in GRID:
        b, ok, Fp = qdgmm._solve_given_phi(mod, float(ph), x)
        if ok:
            x = b
        prof.append((float(ph), float(Fp) if ok else None, [float(v) for v in b] if ok else None))
    roots = []
    for (p0, f0, b0), (p1, f1, b1) in zip(prof[:-1], prof[1:]):
        if f0 is not None and f1 is not None and np.sign(f0) != np.sign(f1):
            w = abs(f0) / (abs(f0) + abs(f1))
            roots.append(dict(phi=(1 - w) * p0 + w * p1, bracket=[p0, p1],
                              coef=[(1 - w) * u + w * v for u, v in zip(b0, b1)]))
    f0 = [f for (ph, f, _) in prof if abs(ph) < 1e-12][0]
    return dict(sid=sid, scenario=scen.name, N=N, rep=rep, seed=seed, phi_bar=phi_bar, G=int(panel["G"]),
                truth=dict(tau=float(p.tau), kappa=float(p.kappa), phi=float(p.phi)), F_at_0=f0, roots=roots,
                spurious=[r for r in roots if r["phi"] > 0.05], secs=round(time.time() - t0, 1),
                profile=[[ph, f] for ph, f, _ in prof])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--sids", default="5,6,13,14")
    ap.add_argument("--reps", type=int, default=3)
    a = ap.parse_args()
    sids = [int(s) for s in a.sids.split(",")]
    jobs = [(s, 20000, r, pb) for s in sids for pb in (0.25, 0.5, 1.0) for r in range(a.reps)]
    jobs += [(s, 50000, 0, 0.5) for s in sids]
    print("jobs", len(jobs), "git", env.git_sha(), flush=True)
    if a.workers > 1:
        import multiprocessing as mp
        it = mp.get_context("fork").Pool(a.workers).imap_unordered(job, jobs)
    else:
        it = map(job, jobs)
    out = []
    for rec in it:
        out.append(rec)
        print(time.strftime("%H:%M:%S", time.gmtime()), rec["sid"], rec["N"], rec["phi_bar"], rec["rep"], rec["secs"], "s",
              "roots", [(round(r["phi"], 3), [round(c, 3) for c in r["coef"]]) for r in rec["roots"]], flush=True)
    out.sort(key=lambda r: (r["sid"], r["N"], r["phi_bar"], r["rep"]))
    summ = []
    for s in sids:
        for N in (20000, 50000):
            for pb in (0.25, 0.5, 1.0):
                d = [r for r in out if r["sid"] == s and r["N"] == N and r["phi_bar"] == pb]
                if not d:
                    continue
                sp = [r["spurious"] for r in d]
                summ.append(dict(sid=s, scenario=d[0]["scenario"], N=N, phi_bar=pb, draws=len(d),
                                 draws_with_spurious=sum(1 for x in sp if x),
                                 spurious_phi=[round(x["phi"], 3) for y in sp for x in y],
                                 spurious_tau=[round(x["coef"][0], 3) for y in sp for x in y],
                                 spurious_kappa=[round(x["coef"][1], 3) for y in sp for x in y] if len(d[0]["roots"]) and
                                 len(d[0]["roots"][0]["coef"]) > 1 else None))
    res = dict(git=env.git_sha(), grid=GRID.tolist(), seeds="scenarios.seed_for(sid, rep) = BASE_SEED + sid*100000 + rep",
               summary=summ, draws=out)
    (HERE / "results.json").write_text(json.dumps(env.clean(res), indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
