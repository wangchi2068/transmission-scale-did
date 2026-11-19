"""Summarise the calibrated simulation (reps_<design>.jsonl) -> results.json.

For each design (cal, stress):
  * M0 tau, kappa, phi, tau+kappa: bias, empirical SD, mean sandwich SE, SE ratio, RMSE;
    coverage of nominal 95% intervals est +/- q * se for q = t_48 (sandwich), the score wild cluster bootstrap-t
    |t| 95% quantile (Mammen) and the pairs moment cluster bootstrap-t |t| 95% quantile; mean interval widths;
    Monte Carlo SE of coverage sqrt(p(1-p)/R).
  * Grouped specifications A1/A2/A3 (division x week effects, region phi): tau+kappa bias and t_48 coverage.
  * M0-G / M0-GI: share of replications where the chosen root differs from the M0 root, bias of tau+kappa at the
    chosen root.
  * Reference distributions under correct specification of J (Hansen two-step, chi2_2), J_LM (M0-OI extras,
    chi2_2), J_LM with the M0-GI extras, Q_G and Q_GI at the chosen root, and the grouped-spec J_LM:
    quantiles, rejection rates at the chi2 5% point, and the exceedance share of the application's value
    (with a Clopper-Pearson upper 95% bound when no replication exceeds it).
Application values are read from ../regional_fe/results.json (computed on the observed panel).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
NAMES = ["tau", "kappa", "phi", "tau_plus_kappa"]


def load(design):
    recs = []
    p = HERE / ("reps_%s.jsonl" % design)
    for line in p.read_text().splitlines():
        try:
            recs.append(json.loads(line))
        except Exception:  # noqa: BLE001
            pass
    meta = json.loads((HERE / ("run_meta_%s.json" % design)).read_text())
    return recs, meta


def cov_block(est, se, truth, qs, G):
    est, se = np.asarray(est), np.asarray(se)
    R = len(est)
    out = dict(R=R, bias=float(np.mean(est - truth)), emp_sd=float(np.std(est, ddof=1)),
               mean_se=float(np.mean(se)), median_se=float(np.median(se)),
               se_ratio=float(np.mean(se) / np.std(est, ddof=1)), rmse=float(np.sqrt(np.mean((est - truth) ** 2))),
               bias_mcse=float(np.std(est, ddof=1) / np.sqrt(R)))
    for lab, q in qs.items():
        q = np.asarray(q)
        ok = np.isfinite(q)
        c = np.abs(est[ok] - truth) <= q[ok] * se[ok]
        p = float(np.mean(c))
        out["cov_" + lab] = p
        out["cov_%s_mcse" % lab] = float(np.sqrt(p * (1 - p) / ok.sum()))
        out["width_" + lab] = float(np.mean(2 * q[ok] * se[ok]))
        out["n_" + lab] = int(ok.sum())
    return out


def ref_block(x, app_value, df=None):
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)], float)
    R = len(x)
    out = dict(R=R, mean=float(np.mean(x)), q50=float(np.quantile(x, 0.5)), q90=float(np.quantile(x, 0.9)),
               q95=float(np.quantile(x, 0.95)), q99=float(np.quantile(x, 0.99)), max=float(np.max(x)))
    if df is not None:
        out["chi2_df"] = df
        out["reject_5pct_chi2"] = float(np.mean(x > stats.chi2.ppf(0.95, df)))
        out["chi2_q95"] = float(stats.chi2.ppf(0.95, df))
    if app_value is not None:
        k = int(np.sum(x >= app_value))
        out["application_value"] = float(app_value)
        out["exceed_share"] = k / R
        out["exceed_count"] = k
        out["exceed_upper95"] = float(stats.beta.ppf(0.95, k + 1, R - k)) if k < R else 1.0
    return out


def main():
    G = 49
    tc = float(stats.t.ppf(0.975, G - 1))
    app = json.loads((HERE.parent / "regional_fe" / "results.json").read_text())
    sp = app["specs"]
    app_vals = dict(J=app["N0_paper_code"]["J"], J_LM_oi=sp["N0"]["J_LM_oi"], J_LM_inv=sp["N0"]["J_LM_inv"],
                    Q_G=min(app["N0_paper_code"]["Q_G"]), Q_GI=min(app["N0_paper_code"]["Q_GI"]))
    for k in ("A1", "A2", "A3"):
        app_vals[k + "_J_LM_oi"] = sp[k]["J_LM_oi"]
        app_vals[k + "_J_LM_inv"] = sp[k]["J_LM_inv"]
    out = dict(application_values=app_vals, tcrit=tc, designs={})
    for design in ("cal", "stress"):
        p = HERE / ("reps_%s.jsonl" % design)
        if not p.exists():
            continue
        recs, meta = load(design)
        ok = [r for r in recs if r.get("m0_ok") and "est" in r]
        truth = meta["truth"]
        tv = [truth["tau"], truth["kappa"], truth["phi"], truth["tau"] + truth["kappa"]]
        D = dict(R_total=len(recs), R_ok=len(ok), errors=sum("error" in r for r in recs),
                 meta={k: meta[k] for k in ("sigma_sw", "shock_adj_corr", "k_nb2", "seeds", "B_wild", "B_pairs", "truth", "git")},
                 M0={})
        est = np.array([r["est"] + [r["tpk"]] for r in ok])
        se = np.array([r["se"] + [r["tpk_se"]] for r in ok])
        wq = np.array([r["wild_q95"] for r in ok], float)
        pq = np.array([r["pairs_q95"] if r.get("pairs_q95") is not None else [np.nan] * 4 for r in ok], float)
        for j, nm in enumerate(NAMES):
            D["M0"][nm] = cov_block(est[:, j], se[:, j], tv[j], dict(t48=np.full(len(ok), tc), wild=wq[:, j], pairs=pq[:, j]), G)
        D["application_se_ratio"] = None
        D["grouped"] = {}
        for k, key in (("A1", "A1_divweek"), ("A2", "A2_regphi"), ("A3", "A3_divweek_regphi")):
            g = [r[key] for r in ok if r.get(key, {}).get("ok")]
            if g:
                D["grouped"][k] = dict(tau_plus_kappa=cov_block([x["tpk"] for x in g], [x["tpk_se"] for x in g], tv[3],
                                                                dict(t48=np.full(len(g), tc)), G),
                                       J_LM_oi=ref_block([x["J_LM_oi"] for x in g], app_vals[k + "_J_LM_oi"], 2),
                                       J_LM_inv=ref_block([x["J_LM_inv"] for x in g], app_vals[k + "_J_LM_inv"], 2),
                                       R_ok=len(g))
        sel = {}
        for lab in ("G", "GI"):
            rr = [r for r in ok if r.get("%s_tpk" % lab) is not None]
            diff = [abs(r["%s_phi" % lab] - r["est"][2]) > 1e-6 for r in rr]
            sel[lab] = dict(R=len(rr), share_root_differs_from_M0=float(np.mean(diff)) if rr else None,
                            share_boundary=float(np.mean([r.get("%s_chosen" % lab) == "boundary" for r in rr])) if rr else None,
                            share_multiple_roots=float(np.mean([r.get("%s_nroots" % lab, 1) > 1 for r in rr])) if rr else None,
                            tpk_bias=float(np.mean([r["%s_tpk" % lab] - tv[3] for r in rr])) if rr else None)
        D["root_selection"] = sel
        D["reference"] = dict(
            J=ref_block([r.get("J") for r in ok], app_vals["J"], 2),
            J_LM_oi=ref_block([r.get("J_LM_oi") for r in ok], app_vals["J_LM_oi"], 2),
            J_LM_inv=ref_block([r.get("J_LM_inv") for r in ok], app_vals["J_LM_inv"], 2),
            Q_G=ref_block([r.get("Q_G") for r in ok], app_vals["Q_G"]),
            Q_GI=ref_block([r.get("Q_GI") for r in ok], app_vals["Q_GI"]))
        out["designs"][design] = D
        m = D["M0"]["tau_plus_kappa"]
        print(design, "R", len(ok), "tpk bias %.3f sd %.3f se %.3f cov t48 %.3f wild %.3f pairs %.3f" % (
            m["bias"], m["emp_sd"], m["mean_se"], m["cov_t48"], m["cov_wild"], m["cov_pairs"]))
        for nm in NAMES[:3]:
            m = D["M0"][nm]
            print("   %-6s bias %.3f sd %.3f se %.3f cov t48 %.3f wild %.3f pairs %.3f" % (
                nm, m["bias"], m["emp_sd"], m["mean_se"], m["cov_t48"], m["cov_wild"], m["cov_pairs"]))
        for k, v in D["reference"].items():
            print("   ref %-8s q95 %.2f q99 %.2f max %.2f app %.2f exceed %d/%d rej5 %s" % (
                k, v["q95"], v["q99"], v["max"], v["application_value"], v["exceed_count"], v["R"], v.get("reject_5pct_chi2")))
        for k, v in D["grouped"].items():
            print("   %s tpk cov %.3f bias %.3f | J_LM_oi q95 %.2f exceed %d | J_LM_inv q95 %.2f rej5 %.3f exceed %d" % (
                k, v["tau_plus_kappa"]["cov_t48"], v["tau_plus_kappa"]["bias"], v["J_LM_oi"]["q95"],
                v["J_LM_oi"]["exceed_count"], v["J_LM_inv"]["q95"], v["J_LM_inv"]["reject_5pct_chi2"], v["J_LM_inv"]["exceed_count"]))
        print("   root sel", sel)
    (HERE / "results.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
