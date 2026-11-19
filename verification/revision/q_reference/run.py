"""Analysis b (review Comment 1, request 2; Comment 3): reference distributions for the application's
over-identification and root-selection statistics.

Sources (no new simulation here; numbers are read from script outputs):
  1. ../calibrated_sim/results.json  -- application-calibrated DGP (49 state clusters, actual policy path, W, GI,
     NB2 k = 15.9, state x week shocks), correct specification of M0; designs "cal" (shock SD 0.23) and
     "stress" (0.50).  Same moment set and cluster count as the application, so this is the like-for-like reference.
  2. ../../../../rerun/out/qcal_s1_0_500.json, qcal_s12_0_300.json -- the paper's simulation scenarios 1 and 12
     (correct specification, T = 12 weeks: 11 week moments instead of 25), quoted for continuity with the paper's
     earlier "outside the simulated range" statement.
Output: for each statistic, application value, chi2 reference (where it applies), simulated q95/q99/max and
the exceedance share with a Clopper-Pearson upper 95% bound.
"""
import json
from pathlib import Path

from scipy import stats

HERE = Path(__file__).resolve().parent
REV = HERE.parent
RERUN = REV.parents[2] / "rerun" / "out"


def main():
    cs = json.loads((REV / "calibrated_sim" / "results.json").read_text())
    app = cs["application_values"]
    out = dict(application_values=app, calibrated={}, paper_scenarios={})
    for design, D in cs["designs"].items():
        rows = dict(D["reference"])
        for k, v in D["grouped"].items():
            rows[k + "_J_LM_oi"] = v["J_LM_oi"]
            rows[k + "_J_LM_inv"] = v["J_LM_inv"]
        out["calibrated"][design] = dict(R_ok=D["R_ok"], rows=rows)
    for fn in ("qcal_s1_0_500.json", "qcal_s12_0_300.json"):
        q = json.loads((RERUN / fn).read_text())
        out["paper_scenarios"]["scenario_%d" % q["sid"]] = dict(n=q["n"], Q_G=q["Q_G"], Q_GI=q["Q_GI"], source="rerun/out/" + fn)
    out["chi2_2"] = dict(q95=float(stats.chi2.ppf(0.95, 2)), q99=float(stats.chi2.ppf(0.99, 2)),
                         p_J=float(stats.chi2.sf(app["J"], 2)))
    lines = []
    for design, C in out["calibrated"].items():
        for k, v in C["rows"].items():
            lines.append("%-7s %-14s app %7.2f | q95 %6.2f q99 %6.2f max %6.2f | exceed %d/%d (upper95 %.3f)" % (
                design, k, v["application_value"], v["q95"], v["q99"], v["max"], v["exceed_count"], v["R"], v["exceed_upper95"]))
    print("\n".join(lines))
    (HERE / "results.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
