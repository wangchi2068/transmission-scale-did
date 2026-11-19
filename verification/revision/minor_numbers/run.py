"""Numbers for the minor corrections m1-m4, m14 and t2 (review section 4.2/4.3), read from the archived results only.

Sources (no re-estimation):
  experiment/analysis/metrics.csv            simulation summaries (R = 8000 per scenario)
  experiment/analysis/results.json           application table (main + sensitivity runs, M0 / M0-G / M0-GI ...)
  verification/application/wild_bootstrap.json  jackknife and sandwich SEs
Output: results.json and log.txt in this folder.
Run (from this folder): python3 -I run.py > log.txt
"""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import env  # noqa: E402

RUN = Path(__file__).resolve().parents[3]
rows = list(csv.DictReader(open(RUN / "experiment/analysis/metrics.csv")))
app = json.load(open(RUN / "experiment/analysis/results.json"))["application"]
wb = json.load(open(RUN / "verification/application/wild_bootstrap.json"))


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def pick(method, param, target=None, sids=None):
    out = []
    for r in rows:
        if r["method"] == method and r["param"] == param and (target is None or r["target"] == target):
            if sids is None or int(r["sid"]) in sids:
                out.append(r)
    return out


res = {"git": env.git_sha()}

# m1: jackknife vs sandwich
ratio = [a / b - 1 for a, b in zip(wb["se_jack"], wb["se_sandwich"])]
res["m1_jackknife_excess_pct"] = dict(values=[round(100 * x, 1) for x in ratio], min=round(100 * min(ratio), 1),
                                      max=round(100 * max(ratio), 1), se_jack=wb["se_jack"], se_sandwich=wb["se_sandwich"])

# m2 / m14: M0 vs M0-G vs M0-GI across application analyses
tab = {}
for r in app["table"]:
    tab.setdefault(r["analysis"], {})[r["method"]] = r
dep = []
allruns = []
for an, d in tab.items():
    m0, g, gi = d.get("M0"), d.get("M0_G"), d.get("M0_GI")
    if not m0 or not m0.get("ok") or "tau_plus_kappa" not in m0 or "phi" not in m0:
        continue   # event-time fits etc. have no pooled tau+kappa
    rec = dict(analysis=an, M0_tau=m0["tau"], M0_tpk=m0["tau_plus_kappa"], M0_tpk_lo=m0.get("tau_plus_kappa_lo"),
               M0_tpk_hi=m0.get("tau_plus_kappa_hi"), M0_phi=m0["phi"])
    for lab, x in (("G", g), ("GI", gi)):
        if x and x.get("ok"):
            rec[lab + "_phi"], rec[lab + "_tpk"], rec[lab + "_tau"] = x["phi"], x["tau_plus_kappa"], x["tau"]
            rec[lab + "_departs"] = bool(abs(x["phi"] - m0["phi"]) > 1e-4)
            rec[lab + "_kind"] = ("boundary" if x["phi"] <= 1e-8 else "interior") if rec[lab + "_departs"] else "same"
            rec[lab + "_sign_reversal"] = bool(rec[lab + "_departs"] and (x["tau_plus_kappa"] > 0) != (m0["tau_plus_kappa"] > 0))
    allruns.append(rec)
    if rec.get("G_departs"):
        dep.append(rec)
res["m2_M0G_departures"] = dep
res["m2_n_departures"] = len(dep)
res["m2_n_sign_reversals"] = sum(r["G_sign_reversal"] for r in dep)
res["application_runs"] = allruns
excl0 = [r["analysis"] for r in allruns if r["M0_tpk_hi"] is not None and r["M0_tpk_hi"] < 0]
res["m2_M0_ci_excludes_zero"] = excl0
tpk = [r["M0_tpk"] for r in allruns]
res["m14_M0_tpk_range"] = [min(tpk), max(tpk)]
res["m14_M0_tau_positive_runs"] = [r["analysis"] for r in allruns if r["M0_tau"] > 0]

# m3: 2-MCSE flags for phi (M0, model target), all scenarios 1-20
flags = []
for r in pick("M0", "phi", "model", set(range(1, 21))):
    b, s = f(r["bias"]), f(r["bias_mcse"])
    if b is None or s is None or s == 0:
        continue
    flags.append(dict(sid=int(r["sid"]), scenario=r["scenario"], bias=b, mcse=s, z=b / s, flag=abs(b) > 2 * s,
                      coverage=f(r["coverage"])))
res["m3_phi_2mcse"] = flags
res["m3_flagged_sids"] = [x["sid"] for x in flags if x["flag"]]
SPILL, NULL, NOSPILL = {1, 2, 3, 4, 9, 10, 11, 12}, {17, 18, 19, 20}, {5, 6, 7, 8, 13, 14, 15, 16}
res["m3_flagged_spill"] = [x["sid"] for x in flags if x["flag"] and x["sid"] in SPILL]
res["m3_flagged_null"] = [x["sid"] for x in flags if x["flag"] and x["sid"] in NULL]
res["m3_null_bias_range"] = [min(x["bias"] for x in flags if x["sid"] in NULL), max(x["bias"] for x in flags if x["sid"] in NULL)]
res["m3_note"] = "phi0 = 0 (no-spillover) scenarios 5-8, 13-16 are reported separately in the paper (spurious-root block)"

# m4: psi coverage, dagger vs population (M0)
psi = []
for r in rows:
    if r["method"] == "M0" and r["param"] in ("psi_treated", "psi_network") and r["target"] in ("dagger", "population"):
        psi.append(dict(sid=int(r["sid"]), scenario=r["scenario"], param=r["param"], target=r["target"],
                        coverage=f(r["coverage"]), relbias_median=f(r["relbias_median"]), n_se=f(r["n_se"])))
res["m4_psi"] = psi
for tg in ("dagger", "population"):
    cv = [x["coverage"] for x in psi if x["target"] == tg and x["coverage"] is not None and x["sid"] <= 20]
    res["m4_%s_coverage_range_core_null" % tg] = [min(cv), max(cv)] if cv else None
for tg in ("dagger", "population"):
    cv = [x["coverage"] for x in psi if x["target"] == tg and x["sid"] in SPILL]
    res["m4_%s_coverage_range_spill" % tg] = [min(cv), max(cv)]
x4 = [x for x in psi if x["sid"] == 24]
res["m4_X4"] = x4
res["m4_application_4week"] = app["psi"]["M0"]["by_horizon"][0]

# t2: M1-BC convergence and coverage per scenario
bc = []
for r in pick("M1_BC", "tau", "model"):
    bc.append(dict(sid=int(r["sid"]), scenario=r["scenario"], conv=f(r["conv"]), n_ok=f(r["n_ok"]), n_se=f(r["n_se"]),
                   R=f(r["R"]), coverage=f(r["coverage"])))
res["t2_M1BC_tau"] = bc
core = [x for x in bc if x["sid"] <= 16]
res["t2_core_failure_pct"] = {x["sid"]: round(100 * (1 - x["conv"]), 1) for x in core if x["conv"] < 1}
res["t2_core_coverage_range"] = [min(x["coverage"] for x in core), max(x["coverage"] for x in core)]
res["t2_note"] = "coverage is computed over converged replications (n_se) of R = 8000"
print(json.dumps(env.clean(res), indent=1))
(Path(__file__).resolve().parent / "results.json").write_text(json.dumps(env.clean(res), indent=1))
