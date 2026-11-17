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


def build(rows, caption, label, note_extra=""):
    m0 = {r["analysis"]: r for r in rows if r["method"] == "M0"}

    def differs(r):
        a = m0.get(r["analysis"])
        return a is None or a.get("tau") is None or r.get("tau") is None or abs(a["tau"] - r["tau"]) > 1e-6

    keep = [r for r in rows if r["method"] not in ("M0_G", "M0_GI") or r["analysis"] == "main" or differs(r)]

    def ci(e, key):
        if e is None or e.get(key) is None:
            return "--"
        lo, hi = e.get(key + "_lo"), e.get(key + "_hi")
        if lo is None or hi is None:
            return MT.fmt(e[key])
        return "%s (%s, %s)" % (MT.fmt(e[key]), MT.fmt(lo), MT.fmt(hi))

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
    note = ("US county-week COVID-19 cases, %s units in %s state clusters, %s weeks (%s to %s); $D$ = state "
            "public mask order. 95\\%% intervals: cluster (state) sandwich with $t_{G-1}$ critical values. "
            "The first column is the model-conditional total effect of a state order on an interior county, "
            "$\\tau+\\kappa$, taken from the joint covariance; about 62\\%% of counties are interior, so for "
            "them the two coefficients load on the same regressor. The lead\\_placebo row reports $\\hat\\ell$, "
            "the coefficient on a lead indicator from a specification that also contains $D_{it}$ and "
            "$\\bar D_{it}$; Proposition~\\ref{prop:sources}(c) shows that it is not a test of pre-trends under "
            "feedback. Data sources, checksums and the mask-order coding are given in Web Appendix D. %s"
            % (app.get("n_units", "--"), app.get("G", "--"), app.get("T", "--"),
               app.get("first_week", "--"), app.get("last_week", "--"), note_extra))
    return MT.tabular("llrrrrr", head, body, note, caption, label)


(ROOT / "tables" / "table_7_application.tex").write_text(
    build(main_rows, r"Application: state mask orders and county COVID-19 transmission; point estimates with 95\% confidence intervals ($t_{G-1}$, $G = 49$ states).",
          "tab:app"), encoding="utf-8")
(ROOT / "tables" / "table_7b_application_sensitivity.tex").write_text(
    build(sens_rows, r"Pre-specified sensitivity analyses of the mask-order application. M0 is shown for every analysis and M0-G or the post hoc M0-GI where the selected root differs from M0.",
          "tab:app-all", "The reference specification is in Table~\\ref{tab:app}. Symmetric Wald intervals near $\\phi=0$ (such as M0-G in \\texttt{gi\\_w1\\_0.9}) can cross the zero boundary; boundary selections ($\\hat\\phi=0$) carry no Wald interval for $\\phi$.").replace(r"\begin{tabular}", r"\renewcommand{\arraystretch}{0.86}\begin{tabular}"), encoding="utf-8")
print("main rows", len(main_rows), "sensitivity rows", len(sens_rows))
