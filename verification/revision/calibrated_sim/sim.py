"""Application-calibrated few-cluster simulation (review Comment 5; also the reference distributions for
Comment 1 / analysis b: J, J_LM, Q_G, Q_GI and the grouped-specification J_LM under correct specification).

Per replication r (seed SEED_DATA + r), written as one JSON line to reps.jsonl immediately:
  * M0 estimate, cluster-sandwich SEs (t_{G-1} = t_48 intervals), tau+kappa and its SE;
  * score wild cluster bootstrap-t (Kline-Santos, Mammen weights; B_WILD draws, seed SEED_WILD + r);
  * pairs moment cluster bootstrap-t (B_PAIRS draws, seed SEED_PAIRS + r) -- the application's three methods;
  * M0-OI Hansen J (two-step), score-form J_LM, M0-G (Q_G) and M0-GI (Q_GI) chosen root and Q at it;
  * grouped specifications of analysis a: A1 division x week effects, A2 region-specific phi,
    A3 division x week + region-specific phi: estimates, SEs, J_LM with the M0-OI extras and with the M0-GI extras.
Resumable: replications already in reps.jsonl are skipped.

Designs: "cal" -- state x week shock SD estimated from the application (calib.SIGMA_SW = 0.23), spatially correlated
across adjacent states (mean adjacent correlation 0.21 as in the data); "stress" -- shock SD 0.50, chosen so that
the simulated sandwich SE of tau+kappa approaches the application's 0.311 (the "cal" design gives about 0.13).
Run (from this folder):  OMP_NUM_THREADS=1 nohup python3 -I sim.py --design cal --reps 600 --workers 4 > run_cal.log 2>&1 &
"""
import argparse
import json
import os
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import env  # noqa
import numpy as np
import qdgmm
import application as A
import wild_bootstrap as WB
from qdx import QDX, oi_extras, inv_extras, groups_for_units, tsum
import calib

HERE = Path(__file__).resolve().parent
SEED_DATA, SEED_WILD, SEED_PAIRS = 20261118000, 20261118500000, 20261118900000
NAMES = ["tau", "kappa", "phi"]
_CAL = None
SIGMA = None


def _cal():
    global _CAL
    if _CAL is None:
        _CAL = calib.calibrated(SIGMA)
        div, reg = groups_for_units(_CAL["d"]["state"])
        _CAL["div"], _CAL["reg"] = div, reg
    return _CAL


def grouped(panel, truth, wgroup, phigroup, P):
    exo, _ = oi_extras(panel)
    exi, _ = inv_extras(panel)
    X = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, wgroup=wgroup, phigroup=phigroup, extras=exo + exi)
    x0 = np.r_[truth[:2], [truth[2]] * P]
    x, ok, res = X.fit(starts=[x0, np.r_[0.0, 0.0, [0.5] * P]])
    out = dict(ok=ok, est=x.tolist())
    if not ok:
        return out
    inf = X.inference(x)
    if inf is None:
        out["ok"] = False
        return out
    se = np.sqrt(np.clip(np.diag(inf["Vb"]), 0, None))
    m = X.moments(x, extra=True)
    out.update(se=se.tolist(), tpk=float(x[0] + x[1]), tpk_se=tsum(inf["Vb"], 0, 1),
               J_LM_oi=X.j_lm(x, sel=[0, 1], m=m)["J_LM"], J_LM_inv=X.j_lm(x, sel=[2, 3], m=m)["J_LM"])
    return out


def one_rep(r, B_wild, B_pairs):
    cal = _cal()
    t0 = time.time()
    panel, _ = calib.simulate(cal, SEED_DATA + r)
    truth = np.array([cal["tau"], cal["kappa"], cal["phi"]])
    rec = dict(rep=r, seed=SEED_DATA + r)
    sp = qdgmm.spec_m0(panel)
    m0 = qdgmm.fit(panel, sp, starts=[truth, [0.0, 0.0, 0.5]])
    rec["m0_ok"] = bool(m0["ok"])
    if not m0["ok"]:
        rec["seconds"] = time.time() - t0
        return rec
    est = np.array([m0["est"][n] for n in NAMES])
    mod = qdgmm.QDModel(panel, sp)
    inf = mod.inference(est)
    Vb = inf["Vb"]
    se = np.sqrt(np.diag(Vb))
    rec.update(est=est.tolist(), se=se.tolist(), tpk=float(est[0] + est[1]), tpk_se=tsum(Vb, 0, 1))
    Tm1 = mod.T - 1
    phi_g = (-np.linalg.inv(inf["Gam"]) @ inf["gg"])[Tm1:, :]
    wild = WB.score_wild_bootstrap(phi_g, B_wild, SEED_WILD + r, "mammen")
    rec["wild_q95"] = wild["q95"]
    pairs = WB.pairs_moment_bootstrap(mod, est, B_pairs, SEED_PAIRS + r)
    rec["pairs_q95"] = pairs["q95"]
    rec["pairs_ok"] = pairs["ok"]
    # over-identification and root selection
    oi = qdgmm.fit_two_step(panel, qdgmm.spec_m0_oi(panel), m0)
    rec["J"] = oi["J"]
    exo, _ = oi_extras(panel)
    exi, _ = inv_extras(panel)
    X = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, extras=exo + exi)
    mX = X.moments(est, extra=True)
    rec["J_LM_oi"] = X.j_lm(est, sel=[0, 1], m=mX)["J_LM"]
    rec["J_LM_inv"] = X.j_lm(est, sel=[2, 3], m=mX)["J_LM"]
    for lab, sel in (("G", "oi"), ("GI", "inv")):
        g = qdgmm.fit_global(panel, selection=sel)
        if g.get("Q"):
            k = int(np.argmin(g["Q"]))
            rec["Q_" + lab] = float(g["Q"][k])
            rec["Q_%s_all" % lab] = [float(q) for q in g["Q"]]
            rec["%s_chosen" % lab] = g.get("chosen")
            rec["%s_phi" % lab] = float(g["est"]["phi"])
            rec["%s_tpk" % lab] = float(g["est"]["tau"] + g["est"]["kappa"])
            rec["%s_nroots" % lab] = int(g["n_roots"])
    # grouped specifications (analysis a reference distributions)
    rec["A1_divweek"] = grouped(panel, truth, cal["div"], None, 1)
    rec["A2_regphi"] = grouped(panel, truth, None, cal["reg"], 4)
    rec["A3_divweek_regphi"] = grouped(panel, truth, cal["div"], cal["reg"], 4)
    rec["seconds"] = time.time() - t0
    return rec


def _work(args):
    r, B_wild, B_pairs = args
    try:
        return env.clean(one_rep(r, B_wild, B_pairs))
    except Exception:  # noqa: BLE001 - record the failure, keep the run going
        return dict(rep=r, error=traceback.format_exc()[-2000:])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=800)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--B_wild", type=int, default=999)
    ap.add_argument("--B_pairs", type=int, default=199)
    ap.add_argument("--sigma", type=float, default=None, help="SD of log state-week shocks (default calib.SIGMA_SW)")
    ap.add_argument("--design", default="cal")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    global SIGMA
    SIGMA = a.sigma
    cal = _cal()
    if a.out is None:
        a.out = str(HERE / ("reps_%s.jsonl" % a.design))
    out = Path(a.out)
    done = set()
    if out.exists():
        for line in out.read_text().splitlines():
            try:
                done.add(json.loads(line)["rep"])
            except Exception:  # noqa: BLE001
                pass
    todo = [(r, a.B_wild, a.B_pairs) for r in range(a.start, a.start + a.reps) if r not in done]
    meta = dict(started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), reps=a.reps, start=a.start, todo=len(todo),
                workers=a.workers, B_wild=a.B_wild, B_pairs=a.B_pairs, seeds=dict(data=SEED_DATA, wild=SEED_WILD, pairs=SEED_PAIRS),
                git=env.git_sha(), design=a.design, sigma_sw=cal["sigma_sw"], sigma_sw_est=cal["sigma_sw_est"],
                shock_rho=cal["shock_rho"], shock_adj_corr=cal["shock_adj_corr"], adj_corr_data=cal["adj_corr_data"],
                nonadj_corr_data=cal["nonadj_corr_data"], k_nb2=cal["k"], recal_hist=cal["recal_hist"],
                truth=dict(tau=cal["tau"], kappa=cal["kappa"], phi=cal["phi"]), gamma=np.log(cal["eg"]).tolist())
    (HERE / ("run_meta_%s.json" % a.design)).write_text(json.dumps(env.clean(meta), indent=1))
    print(json.dumps(meta), flush=True)
    if a.workers <= 1:
        it = map(_work, todo)
    else:
        import multiprocessing as mp
        pool = mp.get_context("fork").Pool(a.workers, initializer=_cal)
        it = pool.imap_unordered(_work, todo)
    with out.open("a") as fh:
        for k, rec in enumerate(it, 1):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            os.fsync(fh.fileno())
            if k % 10 == 0 or "error" in rec:
                print(time.strftime("%H:%M:%S", time.gmtime()), "done", k, "/", len(todo), "rep", rec["rep"],
                      "err" if "error" in rec else "%.1fs" % rec.get("seconds", -1), flush=True)


if __name__ == "__main__":
    main()
