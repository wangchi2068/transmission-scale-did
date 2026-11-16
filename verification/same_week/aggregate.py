"""Aggregate rows.jsonl into results.json (M2: same-week propagation vs matched benchmark).

Reported per arm (w0 = 0, w0 = 0.3) for the M0 estimates of tau, kappa, phi and tau+kappa:
  bias, empirical SD, mean model SE, 95% coverage (t_{G-1} intervals), J rejection rate at 5%,
  convergence rate, each with its Monte Carlo standard error; J power = rejection(w0=0.3) - rejection(w0=0)
  with its MCSE (independent arms: same seeds, different data.

Usage:  python -I aggregate.py [--rows rows.jsonl] [--out results.json]
"""
from __future__ import annotations

import argparse
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TAU = float(np.log(0.8))
TRUTH = {"tau": TAU, "kappa": -0.1, "phi": 0.4, "tau+kappa": TAU - 0.1, "tau_simple": TAU}
GCL = 50
TCRIT = 2.009575       # t_{G-1, 0.975} with G = 50


def arm_metrics(rows, arm):
    r = [x[arm] for x in rows]
    R = len(r)
    ok = np.array([bool(x.get("ok")) for x in r])
    conv = float(ok.mean())
    out = {"R": R, "conv": conv, "conv_mcse": math.sqrt(conv * (1 - conv) / R)}
    for name in ("tau", "kappa", "phi", "tau+kappa"):
        if name == "tau+kappa":
            est = np.array([(x["beta"]["tau"] + x["beta"]["kappa"]) if x.get("ok") and x.get("beta") else np.nan for x in r])
            se = np.array([math.sqrt(max(x["Vb"][0][0] + x["Vb"][1][1] + 2 * x["Vb"][0][1], 0.0))
                           if x.get("ok") and x.get("Vb") else np.nan for x in r])
        else:
            est = np.array([x["beta"][name] if x.get("ok") and name in x.get("beta", {}) else np.nan for x in r])
            se = np.array([x["se"][name] if x.get("ok") and name in x.get("se", {}) else np.nan for x in r])
        good = np.isfinite(est) & np.isfinite(se) & (se > 0)
        n = int(good.sum())
        rec = {"n": n}
        if n:
            e, s = est[good], se[good]
            err = e - TRUTH[name]
            rec["mean_est"] = float(e.mean())
            rec["bias"] = float(err.mean())
            rec["empSD"] = float(err.std(ddof=1))
            rec["bias_mcse"] = rec["empSD"] / math.sqrt(n)
            rec["modelSE"] = float(s.mean())
            rec["modelSE_mcse"] = float(s.std(ddof=1) / math.sqrt(n))
            cov = float((np.abs(err) <= TCRIT * s).mean())
            rec["coverage"] = cov
            rec["coverage_mcse"] = math.sqrt(cov * (1 - cov) / n)
            rej = float((np.abs(e / s) > TCRIT).mean())
            rec["reject0"] = rej
            rec["reject0_mcse"] = math.sqrt(rej * (1 - rej) / n)
            rec["median_est"] = float(np.median(e))
        out[name] = rec
    jp = np.array([x.get("J_p") if x.get("J_ok") and x.get("J_p") is not None else np.nan for x in r])
    j = np.isfinite(jp)
    if j.any():
        rej = float((jp[j] < 0.05).mean())
        out["J_reject_05"] = {"p": rej, "mcse": math.sqrt(rej * (1 - rej) / int(j.sum())), "n": int(j.sum()),
                              "J_df": int(r[0].get("J_df") or 0)}
    else:
        out["J_reject_05"] = {"p": None, "mcse": None, "n": 0}
    tr = np.array([x["desc"]["share_treated"] for x in r if x.get("desc")])
    ya = np.array([x["desc"]["mean_Y_last"] for x in r if x.get("desc")])
    out["design"] = {"share_treated": float(tr.mean()), "share_treated_mcse": float(tr.std(ddof=1) / math.sqrt(tr.size)),
                     "mean_Y_last": float(ya.mean())}
    return out


def power(a0, a1):
    p0r, p1r = a0["J_reject_05"], a1["J_reject_05"]
    if p0r["p"] is None or p1r["p"] is None:
        return None
    p0, p1 = p0r["p"], p1r["p"]
    d = p1 - p0
    return {"w0_0": p0, "w0_0_mcse": p0r["mcse"], "w0_03": p1, "w0_03_mcse": p1r["mcse"], "diff": d,
            "diff_mcse": math.sqrt(p0r["mcse"] ** 2 + p1r["mcse"] ** 2),
            "z": d / math.sqrt(p0r["mcse"] ** 2 + p1r["mcse"] ** 2)}


def diagnostics(rows):
    """Same seeds in both arms, so the paired differences show what the same-week term does beyond the
    shared innovation noise: means of Delta_tau, Delta_kappa, paired t statistic, treated-share gap, and the
    share of replications in which the two arms select the same parameter signs."""
    b = [x["benchmark"]["beta"] for x in rows if x["benchmark"].get("ok")]
    s = [x["same_week"]["beta"] for x in rows if x["same_week"].get("ok")]
    n = min(len(b), len(s))
    out = {"n_pairs": n}
    for nm in ("tau", "kappa", "phi"):
        d = np.array([s[i][nm] - b[i][nm] for i in range(n)])
        out["d_" + nm] = {"mean": float(d.mean()), "mcse": float(d.std(ddof=1) / math.sqrt(n)),
                          "sd": float(d.std(ddof=1)), "corr": float(np.corrcoef([s[i][nm] for i in range(n)],
                                                                                [b[i][nm] for i in range(n)])[0, 1])}
    tb = np.array([x["benchmark"]["desc"]["share_treated"] for x in rows if x["benchmark"].get("desc")])
    ts = np.array([x["same_week"]["desc"]["share_treated"] for x in rows if x["same_week"].get("desc")])
    out["share_treated"] = {"w0_0": float(tb.mean()), "w0_03": float(ts.mean()),
                            "gap": float(ts.mean() - tb.mean())}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", default=os.path.join(HERE, "rows.jsonl"))
    ap.add_argument("--out", default=os.path.join(HERE, "results.json"))
    a = ap.parse_args()
    with open(a.rows, encoding="utf-8") as fh:
        rows = [json.loads(line) for line in fh if line.strip()]
    rows.sort(key=lambda r: r["rep"])
    res = {"task": "M2 same-week propagation sensitivity", "scenario": "core_N500_T12_fb_spill (sid 1)",
           "R": len(rows), "seeds": rows[0]["seed"] if rows else None, "truth": TRUTH,
           "arms": {"benchmark": {"w0": 0.0}, "same_week": {"w0": 0.3}}, "tcrit": TCRIT, "G": GCL}
    for arm in ("benchmark", "same_week"):
        res["arms"][arm].update(arm_metrics(rows, arm))
    res["J_power"] = power(res["arms"]["benchmark"], res["arms"]["same_week"])
    res["diagnostics"] = diagnostics(rows)
    b, s = res["arms"]["benchmark"], res["arms"]["same_week"]
    res["diff"] = {name: {"bias_same_week_minus_benchmark": s[name]["bias"] - b[name]["bias"],
                          "mcse": math.sqrt(b[name]["bias_mcse"] ** 2 + s[name]["bias_mcse"] ** 2)}
                   for name in ("tau", "kappa", "phi", "tau+kappa")}
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    # console table
    print("%-10s %-8s %8s %8s %8s %8s %8s %8s %8s %8s" % ("arm", "par", "bias", "mcse", "empSD", "modelSE",
                                                          "cover", "mcse", "Jrej", "conv"))
    for arm in ("benchmark", "same_week"):
        A = res["arms"][arm]
        for name in ("tau", "kappa", "phi", "tau+kappa"):
            m = A[name]
            print("%-10s %-8s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f" % (
                arm, name, m["bias"], m["bias_mcse"], m["empSD"], m["modelSE"], m["coverage"],
                m["coverage_mcse"], A["J_reject_05"]["p"], A["conv"]))
    print("J power:", res["J_power"])


if __name__ == "__main__":
    main()
