"""Aggregate Monte Carlo rows into performance measures with Monte Carlo standard errors.

Reads whatever rows exist in sim/rows.jsonl (partial runs, truncated last line and missing scenarios are
fine) and writes
  * results.json                 03-results-protocol schema (task/dataset/metrics/seeds/statistics/provenance/
                                 summary) plus the long-format table `metrics_table` and scenario-level blocks;
  * sim/metrics.csv              one row per scenario x method x parameter x target (simulation_spec.md, Output).
If application/application_results.json exists, it is copied into results.json["application"].

Performance measures (Morris, White & Crowther 2019, Stat Med 38:2074), per scenario x method x parameter,
over the n replications in which the method converged and returned a finite estimate:
  bias       mean(est) - truth                          MCSE empSE / sqrt(n)
  empSE      sd(est, ddof=1)                            MCSE empSE / sqrt(2(n-1))
  modelSE    mean(se)
  relerr_SE  modelSE / empSE - 1                        MCSE (modelSE/empSE) sqrt(var(se^2)/(4 n modelSE^4) + 1/(2(n-1)))
  RMSE       sqrt(mean((est - truth)^2))                MCSE jackknife
  coverage   P(|est - truth| <= t_{G-1} se)             MCSE sqrt(p(1-p)/n)
  ci_width   mean(2 t_{G-1} se)                         MCSE sd / sqrt(n)
  reject0    P(|est / se| > t_{G-1})  (H0: parameter = 0; type-I error in the null block for tau)
  conv       n_ok / R  (denominator R = replications present for the scenario)
Replication-specific truths (psi-dagger, X3 switcher mean) use the per-replication value in est - truth.
Run:  python src/aggregate.py [--rows sim/rows.jsonl] [--out results.json]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time
from collections import defaultdict

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import scenarios as sc  # noqa: E402

METHODS = ("M0", "M0_G", "M0_GI", "M0_ET", "M0_OI", "M0p", "M1", "M1_BC", "M2", "M3", "M4", "Oracle")
# Methods fitted by the supplementary runner (src/simulate_supp.py, post-hoc, experiment_design.md deviation 7):
# merged into the main rows by (sid, rep) and evaluated over the replications on which they were run.
SUPP_METHODS = ("M0_GI",)
G_CLUSTERS = 50                    # dgp.make_geography: 50 "states" for every N
STD_DDOF = 1
METRIC_KEYS = ("bias", "bias_mcse", "empSE", "empSE_mcse", "modelSE", "relerr_SE", "relerr_SE_mcse", "RMSE",
               "RMSE_mcse", "coverage", "coverage_mcse", "ci_width", "ci_width_mcse", "reject0", "reject0_mcse",
               "n_ok", "n_se", "R", "conv")


# ------------------------------------------------------------------------------------------------
# I/O
# ------------------------------------------------------------------------------------------------
def load_rows(path):
    reps, psis, bad = defaultdict(list), defaultdict(list), 0
    if not os.path.exists(path):
        return reps, psis, bad
    seen = set()
    with open(path) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                bad += 1
                continue
            if r.get("kind") == "rep":
                key = (r["sid"], r["rep"])
                if key in seen:          # duplicates cannot arise from the runner, but be safe
                    continue
                seen.add(key)
                reps[r["sid"]].append(r)
            elif r.get("kind") == "psi_truth":
                key = ("psi", r["sid"], r["block"])
                if key in seen:
                    continue
                seen.add(key)
                psis[r["sid"]].append(r)
    for sid in reps:
        reps[sid].sort(key=lambda r: r["rep"])
    return reps, psis, bad


def merge_supp(reps, path):
    """Attach supplementary-method results (kind = "supp") to the main replication rows by (sid, rep). The
    supplementary runner regenerates the same data set from the same seeds; max_Y is checked as a fingerprint
    and a mismatch is counted and not merged."""
    info = dict(path=os.path.relpath(path, RUN) if path else None, n_rows=0, merged=0, mismatch=0, orphan=0,
                bad_lines=0)
    if not path or not os.path.exists(path):
        return info
    idx = {(r["sid"], r["rep"]): r for v in reps.values() for r in v}
    seen = set()
    with open(path) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                info["bad_lines"] += 1
                continue
            key = (r["sid"], r["rep"])
            if r.get("kind") != "supp" or key in seen:
                continue
            seen.add(key)
            info["n_rows"] += 1
            main = idx.get(key)
            if main is None:
                info["orphan"] += 1
                continue
            if main["seed"] != r["seed"] or _f(main["desc"].get("max_Y")) != _f(r.get("max_Y")):
                info["mismatch"] += 1
                continue
            main["res"].update(r["res"])
            main.setdefault("time", {}).update(r.get("time", {}))
            info["merged"] += 1
    return info


def git_sha():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=RUN, capture_output=True, text=True,
                              timeout=10).stdout.strip() or None
    except Exception:  # noqa: BLE001
        return None


def _f(x):
    """float or nan"""
    return float(x) if x is not None else float("nan")


def _num(x, nd=None):
    if x is None:
        return None
    if isinstance(x, (float, np.floating)):
        x = float(x)
        if not math.isfinite(x):
            return None
        return round(x, nd) if nd else x
    if isinstance(x, (np.integer,)):
        return int(x)
    return x


# ------------------------------------------------------------------------------------------------
# Performance measures
# ------------------------------------------------------------------------------------------------
def perf(est, se, truth, R, tcrit):
    """est, se, truth: arrays over replications (nan = failed / missing). Returns dict of measures."""
    est, se, truth = (np.asarray(a, float) for a in (est, se, truth))
    ok = np.isfinite(est) & np.isfinite(truth)
    n = int(ok.sum())
    out = dict(R=int(R), n_ok=n, conv=n / R if R else None)
    if n == 0:
        return out
    e, th = est[ok], truth[ok]
    err = e - th
    out["bias"] = float(err.mean())
    out["mean_est"] = float(e.mean())
    out["truth_mean"] = float(th.mean())
    if n >= 2:
        # empSE is the sd of the estimates; with a per-replication truth it is the sd of est - truth
        const_truth = np.ptp(th) == 0
        es = float(np.std(e if const_truth else err, ddof=1))
        out["empSE"] = es
        out["bias_mcse"] = es / math.sqrt(n)
        out["empSE_mcse"] = es / math.sqrt(2 * (n - 1))
        mse = float((err ** 2).mean())
        out["RMSE"] = math.sqrt(mse)
        loo = np.sqrt(np.maximum((mse * n - err ** 2) / (n - 1), 0.0))
        out["RMSE_mcse"] = float(math.sqrt((n - 1) / n * ((loo - loo.mean()) ** 2).sum()))
    s_ok = ok & np.isfinite(se) & (se > 0)
    ns = int(s_ok.sum())
    out["n_se"] = ns
    if ns >= 1:
        s, ee, tt = se[s_ok], est[s_ok], truth[s_ok]
        out["modelSE"] = float(s.mean())
        cov = (np.abs(ee - tt) <= tcrit * s).astype(float)
        p = float(cov.mean())
        out["coverage"] = p
        out["coverage_mcse"] = math.sqrt(p * (1 - p) / ns)
        w = 2 * tcrit * s
        out["ci_width"] = float(w.mean())
        out["ci_width_mcse"] = float(w.std(ddof=1) / math.sqrt(ns)) if ns >= 2 else None
        rej = (np.abs(ee / s) > tcrit).astype(float)
        q = float(rej.mean())
        out["reject0"] = q
        out["reject0_mcse"] = math.sqrt(q * (1 - q) / ns)
        if ns >= 2 and out.get("empSE"):
            es = out["empSE"]
            m = out["modelSE"]
            out["relerr_SE"] = m / es - 1.0
            v2 = float(np.var(s ** 2, ddof=1))
            out["relerr_SE_mcse"] = (m / es) * math.sqrt(v2 / (4 * ns * m ** 4) + 1.0 / (2 * (n - 1)))
    return out


def prop(x):
    x = np.asarray([v for v in x if v is not None], float)
    n = len(x)
    if n == 0:
        return dict(p=None, mcse=None, n=0)
    p = float(x.mean())
    return dict(p=p, mcse=math.sqrt(p * (1 - p) / n), n=n)


# ------------------------------------------------------------------------------------------------
# Scenario-level truths
# ------------------------------------------------------------------------------------------------
def param_truth(p, name, rows):
    """Truth vector over replications for parameter `name` of the generating model."""
    R = len(rows)
    if name in ("tau", "kappa", "phi"):
        return np.full(R, getattr(p, name))
    if name.startswith("tau_e"):
        k = int(name[5:])
        if p.dyn is not None:
            return np.full(R, p.dyn[min(k, len(p.dyn) - 1)])
        return np.full(R, p.tau)                 # static effect: every event-time bin equals tau
    return np.full(R, np.nan)


def x3_targets(rows):
    """X3: per-replication switcher mean of tau_i and the pooled first-order pseudo-true value of
    Proposition 3, tau_tilde = sum_t (A_t/M_t) H_t / sum_t (A_t/M_t) K_t, with the population limits
    A_t, M_t, H_t, K_t replaced by Monte Carlo means over the available replications."""
    sw = np.array([_f(r["truth"]["x3"].get("tau_switch_mean")) for r in rows])
    A = np.mean([r["truth"]["x3"]["A"] for r in rows], 0)
    M = np.mean([r["truth"]["x3"]["M"] for r in rows], 0)
    K = np.mean([r["truth"]["x3"]["K"] for r in rows], 0)
    H = np.mean([r["truth"]["x3"]["H"] for r in rows], 0)
    wgt = np.where(M > 0, A / np.where(M > 0, M, 1.0), 0.0)
    den = float((wgt * K).sum())
    first = float((wgt * H).sum() / den) if den != 0 else float("nan")
    return sw, first


def psi_population(psis):
    """Population psi from common-random-number pairs: mean and its Monte Carlo SE."""
    out = {}
    if not psis:
        return out
    n = sum(b["n"] for b in psis)
    for k in psis[0]["sum"]:
        s = sum(b["sum"][k] for b in psis)
        ss = sum(b["sumsq"][k] for b in psis)
        m = s / n
        var = max(ss / n - m * m, 0.0) * n / (n - 1) if n > 1 else float("nan")
        out[k] = dict(mean=m, mcse=math.sqrt(var / n) if n > 1 else None, sd=math.sqrt(var) if n > 1 else None,
                      n_pairs=n)
    return out


# ------------------------------------------------------------------------------------------------
# Main aggregation
# ------------------------------------------------------------------------------------------------
def aggregate(reps, psis, scen):
    tcrit = float(stats.t.ppf(0.975, G_CLUSTERS - 1))
    table, scen_info = [], {}
    for sid in sorted(reps):
        rows = reps[sid]
        s = scen[sid]
        p = s.params
        R = len(rows)
        info = dict(sid=sid, name=s.name, block=s.block, R=R, reps=[rows[0]["rep"], rows[-1]["rep"]],
                    seeds=[rows[0]["seed"], rows[-1]["seed"]], params=s.params.to_dict(), w_est=list(s.w_est))
        popsi = psi_population(psis.get(sid, []))
        info["psi_population"] = popsi
        x3sw, x3first = (None, None)
        if p.tau_sd > 0 and all("x3" in r["truth"] for r in rows):
            x3sw, x3first = x3_targets(rows)
            info["x3_first_order_tau"] = x3first
            info["x3_switch_mean_tau_avg"] = float(np.nanmean(x3sw))
        # descriptives
        dk = rows[0]["desc"].keys()
        info["desc_median_over_reps"] = {k: float(np.median([_f(r["desc"].get(k)) for r in rows])) for k in dk}
        info["desc_max_over_reps"] = {k: float(np.max([_f(r["desc"].get(k)) for r in rows])) for k in dk}
        info["time_mean_sec"] = {m: float(np.mean([_f(r["time"].get(m)) for r in rows if m in r["time"]]))
                                 for m in rows[0]["time"]}
        info["time_total_mean_sec"] = float(np.mean([sum(v for v in r["time"].values() if v) for r in rows]))
        rows_all, x3sw_all = rows, x3sw
        for m in METHODS:
            rows, x3sw = rows_all, x3sw_all
            if m in SUPP_METHODS:         # evaluate over the replications the supplementary run reached
                keep = [i for i, r in enumerate(rows_all) if m in r["res"]]
                info.setdefault("supp_R", {})[m] = len(keep)
                if not keep:
                    continue
                rows = [rows_all[i] for i in keep]
                if x3sw_all is not None:
                    x3sw = x3sw_all[keep]
            Rm = len(rows)
            res = [r["res"].get(m, {}) for r in rows]
            names = []
            for x in res:
                for k in x.get("est", {}):
                    if k not in names:
                        names.append(k)
            conv = prop([1.0 if x.get("ok") else 0.0 for x in res])
            info.setdefault("convergence", {})[m] = conv
            for nm in names:
                est = np.array([_f(x.get("est", {}).get(nm)) if x.get("ok") else np.nan for x in res])
                se = np.array([_f(x.get("se", {}).get(nm)) if x.get("ok") else np.nan for x in res])
                targets = [("model", param_truth(p, nm, rows))]
                if nm == "tau" and x3sw is not None:
                    targets += [("switch_mean", x3sw), ("first_order", np.full(Rm, x3first))]
                for tname, tr in targets:
                    rec = dict(sid=sid, scenario=s.name, block=s.block, method=m, param=nm, target=tname)
                    rec.update(perf(est, se, tr, Rm, tcrit))
                    table.append(rec)
            # cases averted
            for k in ("treated", "network"):
                pe = np.array([_f(x.get("psi", {}).get(k, {}).get("est")) if x.get("ok") and x.get("psi") else np.nan
                               for x in res])
                if not np.isfinite(pe).any():
                    continue
                ps = np.array([_f(x.get("psi", {}).get(k, {}).get("se")) if x.get("ok") and x.get("psi") else np.nan
                               for x in res])
                dag = np.array([_f(r["truth"]["psi_dagger"].get(k)) for r in rows])
                rec = dict(sid=sid, scenario=s.name, block=s.block, method=m, param="psi_" + k, target="dagger")
                rec.update(perf(pe, ps, dag, Rm, tcrit))
                ok = np.isfinite(pe) & np.isfinite(dag) & (dag != 0)
                rec["relbias"] = float(np.mean((pe[ok] - dag[ok]) / np.abs(dag[ok]))) if ok.any() else None
                rec["relbias_median"] = float(np.median((pe[ok] - dag[ok]) / np.abs(dag[ok]))) if ok.any() else None
                table.append(rec)
                if k in popsi:
                    rec = dict(sid=sid, scenario=s.name, block=s.block, method=m, param="psi_" + k,
                               target="population")
                    rec.update(perf(pe, ps, np.full(Rm, popsi[k]["mean"]), Rm, tcrit))
                    rec["truth_mcse"] = popsi[k]["mcse"]
                    rec["relbias"] = (rec["bias"] / abs(popsi[k]["mean"])) if rec.get("bias") is not None and \
                        popsi[k]["mean"] else None
                    table.append(rec)
        rows, x3sw = rows_all, x3sw_all
        gi = [r["res"]["M0_GI"] for r in rows if "M0_GI" in r["res"]]
        if gi:
            info["M0GI_at_boundary"] = prop([float(bool(x.get("at_boundary"))) if x.get("ok") else None for x in gi])
            info["M0GI_differs_from_M0G"] = prop([
                float(max(abs(_f(r["res"]["M0_GI"]["est"].get(k)) - _f(r["res"]["M0_G"]["est"].get(k)))
                          for k in ("tau", "kappa", "phi")) > 1e-6)
                if r["res"]["M0_GI"].get("ok") and r["res"].get("M0_G", {}).get("ok") else None
                for r in rows if "M0_GI" in r["res"]])
        # J test (M0-OI), M0-G root diagnostics
        oi = [r["res"].get("M0_OI", {}) for r in rows]
        info["J_reject_05"] = prop([float(x["J_p"] < 0.05) if x.get("ok") and x.get("J_p") is not None else None
                                    for x in oi])
        info["J_df"] = next((x.get("J_df") for x in oi if x.get("J_df") is not None), None)
        g = [r["res"].get("M0_G", {}) for r in rows]
        m0 = [r["res"].get("M0", {}) for r in rows]
        info["M0G_multiple_roots"] = prop([float(x.get("n_roots", 0) > 1) if x.get("ok") else None for x in g])
        info["M0G_no_interior_root"] = prop([float(x.get("n_roots", 0) == 0) if x.get("ok") else None for x in g])
        info["M0G_at_boundary"] = prop([float(bool(x.get("at_boundary"))) if x.get("ok") else None for x in g])
        info["M0G_differs_from_M0"] = prop([
            float(max(abs(_f(a["est"].get(k)) - _f(b["est"].get(k))) for k in ("tau", "kappa", "phi")) > 1e-6)
            if a.get("ok") and b.get("ok") else None for a, b in zip(g, m0)])
        info["M1_BC_failure_reasons"] = _reasons([r["res"].get("M1_BC", {}) for r in rows])
        info["failure_reasons"] = {m: _reasons([r["res"].get(m, {}) for r in rows]) for m in METHODS}
        scen_info[sid] = info
    return table, scen_info, tcrit


def _reasons(res):
    c = defaultdict(int)
    for x in res:
        if not x.get("ok"):
            c[str(x.get("reason", "not ok"))[:80]] += 1
    return dict(c)


def to_summary(table):
    """protocol summary: method -> metric -> cell(scenario[|param|target]) -> {mean, std, n_seeds, mcse}.
    For bias, `mean` is the bias and `std` the sd of est - truth (ddof 1); for coverage / rejection, `mean`
    is the proportion and `std` the sd of the 0/1 indicator; n_seeds = replications entering the measure."""
    out = defaultdict(lambda: defaultdict(dict))
    for r in table:
        cell = "%s|%s" % (r["scenario"], r["param"]) + ("" if r["target"] == "model" else "@" + r["target"])
        n = r.get("n_ok", 0)
        if r.get("bias") is not None:
            sd = r.get("empSE")
            out[r["method"]]["bias"][cell] = dict(mean=r["bias"], std=sd, n_seeds=n, mcse=r.get("bias_mcse"))
        if r.get("RMSE") is not None:
            out[r["method"]]["RMSE"][cell] = dict(mean=r["RMSE"], std=None, n_seeds=n, mcse=r.get("RMSE_mcse"))
        if r.get("empSE") is not None:
            out[r["method"]]["empSE"][cell] = dict(mean=r["empSE"], std=None, n_seeds=n, mcse=r.get("empSE_mcse"))
        for key in ("coverage", "reject0"):
            if r.get(key) is not None:
                p, ns = r[key], r["n_se"]
                out[r["method"]][key][cell] = dict(mean=p, std=math.sqrt(p * (1 - p) * ns / (ns - 1)) if ns > 1 else None,
                                                   n_seeds=ns, mcse=r.get(key + "_mcse"))
        if r.get("relerr_SE") is not None:
            out[r["method"]]["relerr_SE"][cell] = dict(mean=r["relerr_SE"], std=None, n_seeds=r["n_se"],
                                                       mcse=r.get("relerr_SE_mcse"))
        if r.get("ci_width") is not None:
            out[r["method"]]["ci_width"][cell] = dict(mean=r["ci_width"], std=None, n_seeds=r["n_se"],
                                                      mcse=r.get("ci_width_mcse"))
        if r["target"] == "model" and r["param"] in ("tau", "tau_e0"):
            out[r["method"]]["convergence"][r["scenario"]] = dict(mean=r["conv"], std=None, n_seeds=r["R"])
    return {m: {k: dict(v) for k, v in d.items()} for m, d in out.items()}


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    return _num(x)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rows", default=os.path.join(RUN, "sim", "rows.jsonl"))
    ap.add_argument("--out", default=os.path.join(RUN, "results.json"))
    ap.add_argument("--csv", default=None, help="default: <rows dir>/metrics.csv")
    ap.add_argument("--application", default=os.path.join(RUN, "application", "application_results.json"))
    ap.add_argument("--supp", default=os.path.join(RUN, "sim", "rows_m0gi.jsonl"),
                    help="supplementary rows (src/simulate_supp.py, M0-GI); merged by (sid, rep) if present")
    ap.add_argument("--label", default=None, help="free-text label, e.g. 'smoke test - not paper numbers'")
    a = ap.parse_args(argv)
    reps, psis, bad = load_rows(a.rows)
    supp_info = merge_supp(reps, a.supp)
    meta_path = os.path.splitext(a.rows)[0] + "_meta.json"
    meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {}
    scen = sc.build_scenarios()
    table, scen_info, tcrit = aggregate(reps, psis, scen)
    n_rep = sum(len(v) for v in reps.values())
    R_target = meta.get("R")
    complete = bool(R_target) and all(len(reps.get(s, [])) >= R_target for s in meta.get("scenarios", []))
    res = dict(
        task="Monte Carlo evaluation of QD-GMM (quasi-differenced GMM for policy effects on a transmission scale) "
             "and comparators; theory/simulation_spec.md (ADEMP)",
        dataset="simulated county-week epidemic panels (25 pre-specified scenarios; src/dgp.py, src/scenarios.py)"
                + ("; application: US county COVID-19 cases and state mask orders (data_contract.md)"
                   if os.path.exists(a.application) else ""),
        metrics=["bias", "empSE", "relerr_SE", "RMSE", "coverage", "ci_width", "reject0", "convergence",
                 "J_reject_05", "M0G_at_boundary"],
        seeds=dict(rule="seed = %d + scenario_id * 100000 + rep; streams [seed, 0] structure, [seed, 1] innovations; "
                        "population psi pairs use seed + %d" % (sc.BASE_SEED, sc.PSI_SEED_OFFSET),
                   base_seed=sc.BASE_SEED, reps_per_scenario={str(s): len(v) for s, v in reps.items()}),
        statistics=dict(std_ddof=STD_DDOF, ci="95%% t_{G-1} intervals, G = %d clusters, t = %.4f" % (G_CLUSTERS, tcrit),
                        mcse="Morris, White & Crowther (2019); `std` in summary cells is the across-replication sd "
                             "of the error (bias), of the 0/1 indicator (coverage, reject0) or null where the "
                             "measure is itself a dispersion; `mcse` is given in every cell"),
        simulated=False,
        simulated_note="false = every number was produced by actually running this code (Monte Carlo study and "
                       "real-data application); the experiment-suite 'simulated' flag denotes placeholder numbers, "
                       "which are not used anywhere. The Monte Carlo data sets are synthetic by design (src/dgp.py).",
        provenance=dict(mode="measured",
                        source="Monte Carlo run of src/simulate.py, rows %s, aggregated by src/aggregate.py"
                               % os.path.relpath(a.rows, RUN),
                        dataset_contract="./data_contract.md",
                        timestamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), git_sha=git_sha(),
                        simulate_version=meta.get("version"), R_target=R_target, n_rep_rows=n_rep,
                        n_truncated_lines_ignored=bad, run_complete=complete, label=a.label,
                        supplementary=dict(supp_info, methods=list(SUPP_METHODS),
                                           note="M0_GI is a post-hoc addition (run-sim-application, "
                                                "experiment_design.md deviation 7), fitted by src/simulate_supp.py "
                                                "on the same seeds/data sets; its measures use the replications "
                                                "it reached (scenarios[sid].supp_R)")),
        summary=to_summary(table),
        scenarios={str(k): v for k, v in scen_info.items()},
        metrics_table=table,
        notes="Partial runs are aggregated as they stand: each scenario reports its own R and each measure its own "
              "denominator (n_ok, n_se). M0_G inference for phi is undefined when the selected solution is the "
              "boundary phi = 0 (se = null); its phi coverage is computed over interior selections only (n_se). "
              "psi values are sums over county-weeks and can be very large in scenarios with explosive growth.",
    )
    if os.path.exists(a.application):
        with open(a.application) as fh:
            res["application"] = json.load(fh)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(clean(res), fh, indent=1)
    os.replace(tmp, a.out)
    csv_path = a.csv or os.path.join(os.path.dirname(os.path.abspath(a.rows)), "metrics.csv")
    cols = ["sid", "scenario", "block", "method", "param", "target", "truth_mean", "mean_est"] + list(METRIC_KEYS) \
        + ["relbias", "relbias_median", "truth_mcse"]
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in table:
            w.writerow({k: _num(r.get(k)) for k in cols})
    print("[aggregate] %d rep rows over %d scenarios (%d truncated lines ignored); %d metric rows -> %s, %s"
          % (n_rep, len(reps), bad, len(table), os.path.relpath(a.out, RUN), os.path.relpath(csv_path, RUN)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
