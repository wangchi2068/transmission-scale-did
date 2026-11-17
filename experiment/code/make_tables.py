"""LaTeX tables generated from results.json (never typed by hand); writes tables/*.tex and tables/manifest.json.

Works with partial runs: a scenario without replications is printed as '--' and every table states the
replications it rests on. Numbers are performance measures (bias, coverage, ...) over Monte Carlo
replications; the largest MCSE of each column block is reported in the table note.
LaTeX requirements: booktabs, adjustbox.
Run:  python src/make_tables.py [--results results.json] [--outdir tables]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import scenarios as sc  # noqa: E402

LABEL = {"M0": "M0", "M0_G": "M0-G", "M0_GI": "M0-GI", "M0_ET": "M0-ET", "M0_OI": "M0-OI", "M0p": "M0$'$", "M1": "M1",
         "M1_BC": "M1-BC", "M2": "M2", "M3": "M3", "M4": "M4", "Oracle": "Oracle"}
PAR = {"tau": r"$\tau$", "kappa": r"$\kappa$", "phi": r"$\phi$"}


class T:
    """Lookup helper over results.json."""

    def __init__(self, res):
        self.res = res
        self.idx = {}
        for r in res.get("metrics_table", []):
            self.idx[(r["sid"], r["method"], r["param"], r["target"])] = r
        self.scen = {int(k): v for k, v in res.get("scenarios", {}).items()}

    def get(self, sid, m, par, key, target="model"):
        r = self.idx.get((sid, m, par, target))
        return None if r is None else r.get(key)

    def R(self, sid):
        return self.scen.get(sid, {}).get("R")


def fmt(x, nd=3, pct=False):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "--"
    if pct:
        return "%.1f" % (100 * x)
    if abs(x) >= 1e4:
        return "%.2e" % x
    return ("%." + str(nd) + "f") % x


def maxof(vals):
    v = [x for x in vals if x is not None]
    return max(v) if v else None


def short_name(s):
    p = s.params
    if s.block == "core":
        return "%d & %d & %s & %s" % (p.N, p.T, "FB" if p.adopt == "feedback" else "RD",
                                       "yes" if p.phi > 0 else "no")
    return s.name


def tabular(colspec, header, body, note, caption, label):
    lines = [r"\begin{table}[!t]", r"\centering", r"\caption{%s}" % caption, r"\label{%s}" % label,
             r"\footnotesize", r"\setlength{\tabcolsep}{3.5pt}",
             r"\begin{adjustbox}{max width=\linewidth}",   # layout-safe: shrink only if wider than the text block
             r"\begin{tabular}{%s}" % colspec, r"\toprule"]
    lines += header + [r"\midrule"] + body + [r"\bottomrule", r"\end{tabular}", r"\end{adjustbox}"]
    lines += [r"\par\smallskip\parbox{\linewidth}{\scriptsize %s}" % note, r"\end{table}", ""]
    return "\n".join(lines)


def r_range(t, sids):
    Rs = [t.R(s) for s in sids if t.R(s)]
    if not Rs:
        return "no replications available", 0
    return ("R = %d" % Rs[0] if min(Rs) == max(Rs) else "R = %d--%d" % (min(Rs), max(Rs))), min(Rs)


# ------------------------------------------------------------------------------------------------
def table_core_tau(t, scen):
    sids = [s for s in sorted(scen) if scen[s].block == "core"]
    meths = ["M0", "M0_G", "M1", "M1_BC", "M2", "M3", "Oracle"]
    head = [r"$N$ & $T$ & adopt & spill & " + " & ".join(r"\multicolumn{2}{c}{%s}" % LABEL[m] for m in meths)
            + r" \\",
            "".join(r"\cmidrule(lr){%d-%d}" % (5 + 2 * i, 6 + 2 * i) for i in range(len(meths))),
            r" & & & & " + " & ".join(["bias & cov"] * len(meths)) + r" \\"]
    body, mb, mc = [], [], []
    for s in sids:
        cells = []
        for m in meths:
            cells += [fmt(t.get(s, m, "tau", "bias")), fmt(t.get(s, m, "tau", "coverage"), pct=True)]
            mb.append(t.get(s, m, "tau", "bias_mcse"))
            mc.append(t.get(s, m, "tau", "coverage_mcse"))
        body.append(short_name(scen[s]) + " & " + " & ".join(cells) + r" \\")
    rr, n = r_range(t, sids)
    note = ("Core scenarios, policy effect $\\tau_0=\\log 0.8$. bias = mean$(\\hat\\tau)-\\tau_0$; cov = coverage (\\%%) "
            "of 95\\%% $t_{G-1}$ intervals. FB = feedback adoption, RD = random adoption; spill = $(\\phi_0,\\kappa_0)"
            "=(0.4,-0.1)$ vs $(0,0)$. %s per scenario; largest MCSE: bias %s, coverage %s (percentage points). "
            "Biases and coverages are computed over converged replications (convergence rates in the supplementary "
            "table). M2 and M3 target $\\tau$ on different scales and are shown against $\\tau_0$ as in the "
            "pre-specified protocol." % (rr, fmt(maxof(mb)), fmt(maxof(mc), pct=True)))
    return tabular("cccc" + "rr" * len(meths), head, body, note,
                   r"Bias (mean of estimate minus true value over replications) and coverage of nominal 95\% confidence intervals for $\tau$ in the 16 core scenarios.", "tab:sim-core-tau"), n


def table_core_spill(t, scen):
    sids = [s for s in sorted(scen) if scen[s].block == "core"]
    meths = ["M0", "M0_G", "M1", "M4", "Oracle"]
    head = [r"$N$ & $T$ & adopt & spill & " + " & ".join(r"\multicolumn{4}{c}{%s}" % LABEL[m] for m in meths)
            + r" \\",
            "".join(r"\cmidrule(lr){%d-%d}" % (5 + 4 * i, 8 + 4 * i) for i in range(len(meths))),
            r" & & & & " + " & ".join([r"\multicolumn{2}{c}{$\kappa$} & \multicolumn{2}{c}{$\phi$}"] * len(meths))
            + r" \\",
            r" & & & & " + " & ".join(["bias & cov & bias & cov"] * len(meths)) + r" \\"]
    body, mb, mc = [], [], []
    for s in sids:
        cells = []
        for m in meths:
            for par in ("kappa", "phi"):
                cells += [fmt(t.get(s, m, par, "bias")), fmt(t.get(s, m, par, "coverage"), pct=True)]
                mb.append(t.get(s, m, par, "bias_mcse"))
                mc.append(t.get(s, m, par, "coverage_mcse"))
        body.append(short_name(scen[s]) + " & " + " & ".join(cells) + r" \\")
    rr, n = r_range(t, sids)
    note = ("Spillover parameters in the core scenarios (truth $(\\kappa_0,\\phi_0)=(-0.1,0.4)$ with spillover, "
            "$(0,0)$ without). For M0-G, $\\phi$ coverage is computed over replications whose selected solution is "
            "interior; boundary selections ($\\hat\\phi=0$) have no Wald interval for $\\phi$ (shares in "
            "Table~\\ref{tab:sim-roots}). %s per scenario; largest MCSE: bias %s, coverage %s (percentage points)."
            % (rr, fmt(maxof(mb)), fmt(maxof(mc), pct=True)))
    return tabular("cccc" + "rrrr" * len(meths), head, body, note,
                   r"Bias (mean of estimate minus true value) and coverage of nominal 95\% confidence intervals for the spillover parameters $\kappa$ and $\phi$ (core scenarios).",
                   "tab:sim-core-spill"), n


def table_null(t, scen):
    sids = [s for s in sorted(scen) if scen[s].block == "null"]
    meths = ["M0", "M0_G", "M0_ET", "M1", "M1_BC", "M2", "M3", "M4", "Oracle"]
    head = [r"scenario & " + " & ".join(LABEL[m] for m in meths) + r" & $J$ (M0-OI) \\"]
    body, mc = [], []
    for s in sids:
        p = scen[s].params
        nm = "%s, $k=%s$" % ("FB" if p.adopt == "feedback" else "RD", "10" if math.isfinite(p.k) else r"\infty")
        cells = []
        for m in meths:
            par = "tau_e0" if m == "M0_ET" else "tau"
            cells.append(fmt(t.get(s, m, par, "reject0"), pct=True))
            mc.append(t.get(s, m, par, "reject0_mcse"))
        J = t.scen.get(s, {}).get("J_reject_05", {}).get("p")
        body.append(nm + " & " + " & ".join(cells) + " & " + fmt(J, pct=True) + r" \\")
    rr, n = r_range(t, sids)
    note = ("Null block ($\\tau_0=0$, $N=500$, $T=12$, spillover present): rejection rate (\\%%) of $H_0{:}\\,\\tau=0$ "
            "at nominal 5\\%% with $t_{G-1}$ critical values (M0-ET: first event-time bin), and rejection rate of the "
            "Hansen $J$ test of M0-OI at 5\\%%. $k$ = negative-binomial dispersion ($\\infty$ = Poisson). %s per "
            "scenario; largest MCSE %s percentage points." % (rr, fmt(maxof(mc), pct=True)))
    return tabular("l" + "r" * (len(meths) + 1), head, body, note,
                   r"Type I error: rejection rate of the 5\%-level $t_{G-1}$ test of $H_0\colon \tau = 0$ when the null holds.", "tab:sim-null"), n


def table_misspec(t, scen):
    sids = [s for s in sorted(scen) if scen[s].block == "misspec"]
    head = [r"scenario & \multicolumn{2}{c}{$\tau$} & \multicolumn{2}{c}{$\kappa$} & \multicolumn{2}{c}{$\phi$} & "
            r"$J$ rej. & M1 $\tau$ bias & M2 $\tau$ bias \\",
            r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
            r" & bias & cov & bias & cov & bias & cov & & & \\"]
    body = []
    for s in sids:
        cells = []
        for par in ("tau", "kappa", "phi"):
            cells += [fmt(t.get(s, "M0", par, "bias")), fmt(t.get(s, "M0", par, "coverage"), pct=True)]
        J = t.scen.get(s, {}).get("J_reject_05", {}).get("p")
        cells += [fmt(J, pct=True), fmt(t.get(s, "M1", "tau", "bias")), fmt(t.get(s, "M2", "tau", "bias"))]
        body.append(scen[s].name.replace("_", r"\_") + " & " + " & ".join(cells) + r" \\")
    # X3 targets and X5 event-time rows
    extra = []
    if 23 in t.scen:
        for tg, lab in (("switch_mean", "switcher mean of $\\tau_i$"), ("first_order", "first-order $\\tilde\\tau$")):
            extra.append(r"\multicolumn{10}{l}{X3, M0 $\tau$ vs %s: bias %s, coverage %s\%%} \\"
                         % (lab, fmt(t.get(23, "M0", "tau", "bias", tg)), fmt(t.get(23, "M0", "tau", "coverage", tg),
                                                                            pct=True)))
    if 25 in t.scen:
        cov = [fmt(t.get(25, "M0_ET", "tau_e%d" % k, "coverage"), pct=True) for k in range(7)]
        bia = [fmt(t.get(25, "M0_ET", "tau_e%d" % k, "bias")) for k in range(7)]
        extra.append(r"\multicolumn{10}{l}{X5, M0-ET bins $e=0,\dots,6^{+}$: bias %s; coverage (\%%) %s} \\"
                     % (", ".join(bia), ", ".join(cov)))
    if extra:
        body += [r"\midrule"] + extra
    rr, n = r_range(t, sids)
    note = ("Misspecification block (base: $N=2000$, $T=26$, feedback adoption, spillover). M0 bias and coverage "
            "(\\%%) against the generating parameters; $J$ rej. = rejection rate of the M0-OI Hansen test at 5\\%% "
            "(power, since every scenario violates the working model). X3 compares $\\hat\\tau$ also with the "
            "replication-specific switcher mean and with the first-order pseudo-true value of Proposition~3; X5 "
            "reports the event-time specification. %s per scenario." % rr)
    return tabular("lrrrrrrrrr", head, body, note, r"Robustness of M0 under misspecification: bias (mean over replications) and coverage of nominal 95\% confidence intervals.", "tab:sim-misspec"), n


def table_psi(t, scen):
    sids = [s for s in sorted(scen) if s in t.scen]
    head = [r"scenario & \multicolumn{3}{c}{M0, treated} & \multicolumn{3}{c}{M0, network} & "
            r"\multicolumn{2}{c}{M0-G, treated} \\",
            r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-9}",
            r" & rel.\ bias & cov$^\dagger$ & cov$^{\rm pop}$ & rel.\ bias & cov$^\dagger$ & cov$^{\rm pop}$ & "
            r"rel.\ bias & cov$^\dagger$ \\"]
    body = []
    for s in sids:
        cells = []
        for m, k in (("M0", "treated"), ("M0", "network")):
            par = "psi_" + k
            cells += [fmt(t.get(s, m, par, "relbias_median", "dagger")),
                      fmt(t.get(s, m, par, "coverage", "dagger"), pct=True),
                      fmt(t.get(s, m, par, "coverage", "population"), pct=True)]
        cells += [fmt(t.get(s, "M0_G", "psi_treated", "relbias_median", "dagger")),
                  fmt(t.get(s, "M0_G", "psi_treated", "coverage", "dagger"), pct=True)]
        body.append(str(s) + " " + scen[s].name.replace("_", r"\_") + " & " + " & ".join(cells) + r" \\")
    rr, n = r_range(t, sids)
    note = ("Cases averted. rel.\\ bias = median over replications of $(\\hat\\psi-\\psi^\\dagger)/|\\psi^\\dagger|$; "
            "cov$^\\dagger$ = coverage (\\%%) of the replication-specific target $\\psi^\\dagger$ (Theorem~4(b)); "
            "cov$^{\\rm pop}$ = coverage of the population $\\psi$ estimated from common-random-number "
            "counterfactual pairs (descriptive; under-coverage is expected with spillover, Remark~3). %s." % rr)
    return tabular("l" + "r" * 8, head, body, note, r"Inference for the number of cases averted $\psi$: relative bias and coverage of nominal 95\% confidence intervals.",
                   "tab:sim-psi"), n


def table_roots(t, scen):
    sids = [s for s in sorted(scen) if s in t.scen]
    head = [r"scenario & $\phi_0$ & conv.\ M0 & conv.\ M0-G & $>1$ root & no interior root & boundary & "
            r"M0-G $\neq$ M0 \\"]
    body = []
    for s in sids:
        d = t.scen[s]
        c = d.get("convergence", {})
        cells = [fmt(scen[s].params.phi, 1), fmt(c.get("M0", {}).get("p"), pct=True),
                 fmt(c.get("M0_G", {}).get("p"), pct=True)]
        cells += [fmt(d.get(k, {}).get("p"), pct=True) for k in ("M0G_multiple_roots", "M0G_no_interior_root",
                                                                  "M0G_at_boundary", "M0G_differs_from_M0")]
        body.append(str(s) + " " + scen[s].name.replace("_", r"\_") + " & " + " & ".join(cells) + r" \\")
    rr, n = r_range(t, sids)
    note = ("Root diagnostics of the profiled QD-GMM equations over the grid $\\phi\\in[0,3]$ (M0-G): share (\\%%) "
            "of replications with more than one interior root, with no interior root, in which the minimum-$Q$ "
            "selection is the boundary $\\phi=0$, and in which M0-G differs from M0 (local root from the "
            "pre-specified start). %s." % rr)
    return tabular("lrrrrrrr", head, body, note, r"Multiple roots and boundary solutions of the QD-GMM equations.",
                   "tab:sim-roots"), n


def table_m0gi(t, scen):
    """Post-hoc M0-GI (run-sim-application, experiment_design.md deviation 7) next to the pre-specified M0 and
    M0-G, all 25 scenarios, tau; M0-GI over the replications the supplementary run reached."""
    sids = [s for s in sorted(scen) if s in t.scen]
    if not any(t.get(s, "M0_GI", "tau", "n_ok") for s in sids):
        return None, 0
    head = [r"scenario & $\phi_0$ & \multicolumn{2}{c}{M0} & \multicolumn{2}{c}{M0-G} & "
            r"\multicolumn{5}{c}{M0-GI (post hoc)} & bound.\ & $R_{\rm GI}$ \\",
            r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}\cmidrule(lr){7-11}",
            r" & & bias & cov & bias & cov & bias & empSE & modSE & RMSE & cov & GI & \\"]
    body, mb, mc, Rs = [], [], [], []
    for s in sids:
        cells = [fmt(scen[s].params.phi, 1)]
        for m in ("M0", "M0_G"):
            cells += [fmt(t.get(s, m, "tau", "bias")), fmt(t.get(s, m, "tau", "coverage"), pct=True)]
        cells += [fmt(t.get(s, "M0_GI", "tau", k)) for k in ("bias", "empSE", "modelSE", "RMSE")]
        cells += [fmt(t.get(s, "M0_GI", "tau", "coverage"), pct=True),
                  fmt((t.scen[s].get("M0GI_at_boundary") or {}).get("p"), pct=True)]
        r = t.get(s, "M0_GI", "tau", "R")
        cells.append(str(r) if r else "--")
        if r:
            Rs.append(r)
        mb.append(t.get(s, "M0_GI", "tau", "bias_mcse"))
        mc.append(t.get(s, "M0_GI", "tau", "coverage_mcse"))
        body.append(str(s) + " " + scen[s].name.replace("_", r"\_") + " & " + " & ".join(cells) + r" \\")
    rr, n = r_range(t, sids)
    note = ("Post-hoc root-selection variant M0-GI (not pre-specified; added after the no-spillover failure of "
            "M0 and M0-G was observed, experiment\\_design.md deviation~7): same candidate roots as M0-G, "
            "selection criterion built from the $\\phi$-instrument at reference values $0.05$ and $2$ instead of "
            "the M0-OI instruments. Same data sets (seeds) as the main run; M0 and M0-G: %s, M0-GI: $R_{\\rm GI}$ "
            "replications (bias, cov over the replications it reached). bound.\\ GI = share (\\%%) of boundary "
            "($\\hat\\phi=0$) selections. Largest MCSE for M0-GI: bias %s, coverage %s (percentage points)."
            % (rr, fmt(maxof(mb)), fmt(maxof(mc), pct=True)))
    return tabular("lr" + "rr" * 2 + "rrrrr" + "rr", head, body, note,
                   r"Post-hoc root selection M0-GI versus the pre-specified M0 and M0-G ($\tau$): bias (mean over replications) and coverage of nominal 95\% confidence intervals.", "tab:sim-m0gi"), \
        (min(Rs) if Rs else 0)


def table_m0_full(t, scen):
    """Every performance measure with its MCSE for the primary estimator M0 (tau), all scenarios."""
    sids = [s for s in sorted(scen) if s in t.scen]
    head = [r"scenario & bias (MCSE) & empSE (MCSE) & modSE & rel.\ err.\ SE (MCSE) & RMSE (MCSE) & "
            r"cov.\ (MCSE) & conv. \\"]
    body = []

    def pm(k, pct=False):
        def f(s):
            v, e = t.get(s, "M0", "tau", k), t.get(s, "M0", "tau", k + "_mcse")
            return "%s (%s)" % (fmt(v, pct=pct), fmt(e, pct=pct))
        return f
    for s in sids:
        cells = [pm("bias")(s), pm("empSE")(s), fmt(t.get(s, "M0", "tau", "modelSE")), pm("relerr_SE")(s),
                 pm("RMSE")(s), pm("coverage", True)(s), fmt(t.get(s, "M0", "tau", "conv"), pct=True)]
        body.append(str(s) + " " + scen[s].name.replace("_", r"\_") + " & " + " & ".join(cells) + r" \\")
    rr, n = r_range(t, sids)
    note = ("Primary estimator M0, $\\tau$: Monte Carlo performance measures with Monte Carlo standard errors "
            "(Morris, White and Crowther, 2019). empSE = sd of $\\hat\\tau$ (ddof 1); modSE = mean cluster-robust "
            "SE; rel.\\ err.\\ SE = modSE/empSE$-1$; cov.\\ (\\%%) of 95\\%% $t_{G-1}$ intervals; conv.\\ = share "
            "(\\%%) of replications in which the solver converged. %s per scenario." % rr)
    return tabular("lrrrrrrr", head, body, note, r"Full performance measures for M0 ($\tau$), all scenarios: bias (mean over replications), empirical and model standard errors, RMSE and coverage of nominal 95\% confidence intervals, with Monte Carlo standard errors.",
                   "tab:sim-m0-full"), n


def table_application(app):
    if not app:
        return None, 0
    rows = app.get("table", [])
    # sensitivity rows: M0-G is listed only where it selects a different solution than M0
    m0 = {r["analysis"]: r for r in rows if r["method"] == "M0"}

    def differs(r):
        a = m0.get(r["analysis"])
        return a is None or a.get("tau") is None or r.get("tau") is None or abs(a["tau"] - r["tau"]) > 1e-6
    rows = [r for r in rows if r["method"] not in ("M0_G", "M0_GI") or r["analysis"] == "main" or differs(r)]
    head = [r"analysis & method & $\hat\tau+\hat\kappa$ (95\% CI) & $\hat\tau$ (95\% CI) & $\hat\kappa$ (95\% CI) & $\hat\phi$ (95\% CI) & $\hat\ell$ (95\% CI) \\"]
    body = []

    def ci(e, key):
        if e is None or e.get(key) is None:
            return "--"
        lo, hi = e.get(key + "_lo"), e.get(key + "_hi")
        if lo is None or hi is None:
            return fmt(e[key])
        return "%s (%s, %s)" % (fmt(e[key]), fmt(lo), fmt(hi))
    for r in rows:
        if r.get("lead") is not None:
            # the lead row reports the lead coefficient only
            body.append("%s & %s & -- & -- & -- & -- & %s \\\\" % (
                r["analysis"].replace("_", r"\_"), LABEL.get(r["method"], r["method"]), ci(r, "lead")))
            continue
        body.append("%s & %s & %s & %s & %s & %s & %s \\\\" % (
            r["analysis"].replace("_", r"\_"), LABEL.get(r["method"], r["method"]),
            ci(r, "tau_plus_kappa"), ci(r, "tau"), ci(r, "kappa"), ci(r, "phi"), ci(r, "lead")))
    note = ("US county-week COVID-19 cases, %s units in %s state clusters, %s weeks (%s to %s); $D$ = state "
            "public mask order. 95\\%% intervals: cluster (state) sandwich with $t_{G-1}$ critical values for "
            "all methods. The first column is the model-conditional total effect of a state order on an "
            "interior county, $\\tau+\\kappa$, taken from the joint covariance of $\\hat\\tau$ and $\\hat\\kappa$; "
            "about 62\\%% of counties are interior, so for them the two coefficients load on the same "
            "regressor. Sensitivity rows show M0, and M0-G / M0-GI (post hoc) only where they select a different "
            "root than M0. The lead\\_placebo row reports $\\hat\\ell$, the coefficient on a lead indicator from a "
            "specification that also contains $D_{it}$ and $\\bar D_{it}$; Proposition~\\ref{prop:sources}(c) shows "
            "that this coefficient is not a test of pre-trends under feedback. Data sources, checksums and the "
            "mask-order coding are given in Web Appendix D."
            % (app.get("n_units", "--"), app.get("G", "--"), app.get("T", "--"), app.get("first_week", "--"),
               app.get("last_week", "--")))
    return tabular("llrrrrr", head, body, note, r"Application: state mask orders and county COVID-19 transmission; point estimates with 95\% confidence intervals ($t_{G-1}$, $G = 49$ states).",
                   "tab:app"), 1


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(RUN, "results.json"))
    ap.add_argument("--outdir", default=os.path.join(RUN, "tables"))
    a = ap.parse_args(argv)
    with open(a.results) as fh:
        res = json.load(fh)
    t = T(res)
    scen = sc.build_scenarios()
    os.makedirs(a.outdir, exist_ok=True)
    specs = [("tab:sim-core-tau", "table_1_core_tau.tex", table_core_tau(t, scen)),
             ("tab:sim-core-spill", "table_2_core_spill.tex", table_core_spill(t, scen)),
             ("tab:sim-null", "table_3_null.tex", table_null(t, scen)),
             ("tab:sim-misspec", "table_4_misspec.tex", table_misspec(t, scen)),
             ("tab:sim-psi", "table_5_psi.tex", table_psi(t, scen)),
             ("tab:sim-roots", "table_6_roots.tex", table_roots(t, scen)),
             ("tab:app", "table_7_application.tex", table_application(res.get("application"))),
             ("tab:sim-m0gi", "table_8_m0gi_posthoc.tex", table_m0gi(t, scen)),
             ("tab:sim-m0-full", "table_9_m0_full.tex", table_m0_full(t, scen))]
    entries = []
    for tid, fn, (tex, n) in specs:
        if tex is None:
            continue
        path = os.path.join(a.outdir, fn)
        with open(path, "w") as fh:
            fh.write(tex)
        with open(path, "rb") as fh:
            dig = hashlib.sha256(fh.read()).hexdigest()
        entries.append(dict(id=tid, path=fn, sha256=dig, generator="src/make_tables.py",
                            data_sources=[os.path.relpath(a.results, RUN)],
                            display="raw", n_seeds=n,
                            std_ddof=res.get("statistics", {}).get("std_ddof", 1),
                            note="cells are Monte Carlo performance measures over n_seeds replications (MCSE in "
                                 "table note)" if tid != "tab:app" else "point estimates with 95% CIs (single "
                                                                        "observed data set)"))
    man = dict(schema_version=1, provenance=res.get("provenance"), tables=entries)
    with open(os.path.join(a.outdir, "manifest.json"), "w") as fh:
        json.dump(man, fh, indent=1)
    print("[tables] wrote %d tables to %s" % (len(entries), os.path.relpath(a.outdir, RUN)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
