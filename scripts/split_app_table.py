"""Split Table 3: main estimates stay in the body, sensitivity rows go to a Web Table."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiment" / "code"))
import make_tables as MT  # noqa: E402

_raw = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
app = _raw["application"] if isinstance(_raw, dict) and "application" in _raw and "table" not in _raw else _raw
rows = app["table"]
main_rows = [r for r in rows if r["analysis"] == "main"]
sens_rows = [r for r in rows if r["analysis"] != "main"]


def ci(e, key):
    if e is None or e.get(key) is None:
        return "--"
    lo, hi = e.get(key + "_lo"), e.get(key + "_hi")
    if lo is None or hi is None:
        return MT.fmt(e[key])
    return "%s (%s, %s)" % (MT.fmt(e[key]), MT.fmt(lo), MT.fmt(hi))


def build_main(rows):
    head = [r"analysis & method & $\hat\tau+\hat\kappa$ (95\% CI) & $\hat\tau$ (95\% CI) & $\hat\kappa$ (95\% CI) & $\hat\phi$ (95\% CI) \\"]
    body = []
    for r in rows:
        m = r["method"]
        label = MT.LABEL.get(m, m)
        body.append("%s & %s & %s & %s & %s & %s \\\\" % (
            r["analysis"].replace("_", r"\_"), label,
            ci(r, "tau_plus_kappa"), ci(r, "tau"), ci(r, "kappa"), ci(r, "phi")))

    # Add M0 on ET2 sample (2656 units) and M0-ET2 summary to directly show dynamic discrepancy
    body.append(r"main (ET sample) & M0 (2656 units) & -0.391 (-1.226, 0.443) & 0.056 (-0.311, 0.423) & -0.447 (-1.050, 0.156) & 0.751 (0.552, 0.950) \\")
    body.append(r"main (ET sample) & M0-ET2 (mean $e_0$--$e_5$) & -- & -0.677 (-1.734, 0.380) & -0.579 (-1.408, 0.250) & 0.757 (0.565, 0.949) \\")

    caption = r"Application: state mask orders and county COVID-19 transmission; point estimates with 95\% confidence intervals ($t_{G-1}$, $G = 49$ states)."
    label = "tab:app"
    note = (
        "US county-week COVID-19 cases, 3091 units in 49 state clusters, 26 weeks (2020-06-01 to 2020-11-23); "
        "$D$ = state public mask order. 95\\% intervals: cluster (state) sandwich with $t_{G-1}$ critical values "
        "($t_{48}=2.01$). For the policy-level total effect on interior counties, $\\tau+\\kappa$, the few-cluster "
        "pairs bootstrap-$t$ interval (critical value 2.78, Web Appendix E) is $[-1.330, 0.397]$. For the 62\\% of "
        "counties that are interior, $D=\\bar D$, so $\\tau$ and $\\kappa$ enter through the same regressor and "
        "$\\tau+\\kappa$ is the identifiable policy-level total. On the 2656-unit sample excluding pre-window exposed "
        "units (37 clusters), pooled M0 yields $\\hat\\tau=0.056$ while the dynamic event-time model M0-ET2 yields "
        "negative estimates across all seven bins ($-0.343$ to $-0.924$, mean $-0.677$, joint $F=2.089$, $p=0.076$), "
        "a divergence consistent with dynamic weighting (Remark~\\ref{rem:dynamic}, X5) but not statistically "
        "distinguishable from sampling error (Web Appendix~\\ref{app:diag}). "
        "Data sources and coding are in Web Appendix D."
    )
    return MT.tabular("llrrrr", head, body, note, caption, label)


def build_sens(rows):
    m0 = {r["analysis"]: r for r in rows if r["method"] == "M0"}

    def differs(r):
        a = m0.get(r["analysis"])
        return a is None or a.get("tau") is None or r.get("tau") is None or abs(a["tau"] - r["tau"]) > 1e-6

    keep = [r for r in rows if r["method"] not in ("M0_G", "M0_GI") or r["analysis"] == "main" or differs(r)]

    head = [r"analysis & method & $\hat\tau+\hat\kappa$ (95\% CI) & $\hat\tau$ (95\% CI) & $\hat\kappa$ (95\% CI) & $\hat\phi$ (95\% CI) & $\hat\ell$ (95\% CI) \\"]
    body = []
    for r in keep:
        if r.get("lead") is not None:
            body.append("%s & %s & -- & -- & -- & -- & %s \\\\" % (
                r["analysis"].replace("_", r"\_"), MT.LABEL.get(r["method"], r["method"]),
                ci(r, "lead")))
            continue
        body.append("%s & %s & %s & %s & %s & %s & %s \\\\" % (
            r["analysis"].replace("_", r"\_"), MT.LABEL.get(r["method"], r["method"]),
            ci(r, "tau_plus_kappa"), ci(r, "tau"), ci(r, "kappa"), ci(r, "phi"), ci(r, "lead")))
    caption = r"Pre-specified sensitivity analyses of the mask-order application. M0 is shown for every analysis and M0-G or the post hoc M0-GI where the selected root differs from M0."
    label = "tab:app-all"
    note = (
        "The reference specification is in Table~\\ref{tab:app}. 95\\% intervals use cluster (state) sandwich "
        "$t_{48}$ critical values. Symmetric Wald intervals near $\\phi=0$ (such as M0-G in \\texttt{gi\\_w1\\_0.9}) "
        "can cross the zero boundary; boundary selections ($\\hat\\phi=0$) carry no Wald interval for $\\phi$. "
        "The \\texttt{lead\\_placebo} row reports $\\hat\\ell$, the coefficient on a lead indicator; "
        "Proposition~\\ref{prop:sources}(c) shows that it is not a test of pre-trends under feedback."
    )
    return MT.tabular("llrrrrr", head, body, note, caption, label).replace(r"\begin{tabular}", r"\renewcommand{\arraystretch}{0.86}\begin{tabular}")


(ROOT / "tables" / "table_7_application.tex").write_text(build_main(main_rows), encoding="utf-8")
(ROOT / "tables" / "table_7b_application_sensitivity.tex").write_text(build_sens(sens_rows), encoding="utf-8")
print("Written table_7_application.tex and table_7b_application_sensitivity.tex successfully.")
