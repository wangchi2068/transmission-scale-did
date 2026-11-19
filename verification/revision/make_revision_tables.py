"""Web Tables for the revision analyses (review Comments 1-7), written from the results.json files only.

Writes paper/tables/table_11 ... table_17 and adds/updates their entries in paper/tables/manifest.json
(idempotent: entries are keyed by id; existing entries for other tables are left alone). Also registers
table_7b (scripts/split_app_table.py) and table_10 (verification/root_selection/build_rsel_table.py), which
existed on disk but were missing from the manifest.

NOTE: experiment/code/make_tables.py rewrites manifest.json from scratch; re-run this script after it.
Run (from paper/):  python3 -I verification/revision/make_revision_tables.py
"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN = HERE.parents[1]
TAB = RUN / "tables"
GEN = "verification/revision/make_revision_tables.py"


def load(name):
    return json.load(open(HERE / name / "results.json"))


def f(x, nd=3):
    if x is None:
        return "--"
    return ("%." + str(nd) + "f") % x


def ci(c, nd=3):
    return "(%s, %s)" % (f(c[0], nd), f(c[1], nd))


def pv(p):
    return "$<$0.001" if p < 0.001 else f(p, 3)


def ivs(sets, nd=3):
    if not sets:
        return "empty"
    return " $\\cup$ ".join("[%s, %s]" % (f(a, nd), f(b, nd)) for a, b in sets)


def table(colspec, header, body, note, caption, label):
    if "Source: " in note:
        head, src = note.split("Source: ", 1)
        note = head + "Source: \\texttt{" + src.rstrip(".").replace("_", "\\_") + "}."
    L = ["\\begin{table}[!t]", "\\centering", "\\caption{%s}" % caption, "\\label{%s}" % label, "\\footnotesize",
         "\\setlength{\\tabcolsep}{3.5pt}", "\\begin{adjustbox}{max width=\\linewidth}",
         "\\begin{tabular}{%s}" % colspec, "\\toprule", header + " \\\\", "\\midrule"]
    L += body
    L += ["\\bottomrule", "\\end{tabular}", "\\end{adjustbox}",
          "\\par\\smallskip\\parbox{\\linewidth}{\\scriptsize %s}" % note, "\\end{table}", ""]
    return "\n".join(L)


# ---------------------------------------------------------------- a: regional
def t_regional():
    r = load("regional_fe")
    cs = load("calibrated_sim")["designs"]["cal"]
    names = dict(N0="national week effects (paper M0)", R4="region $\\times$ week", A1="division $\\times$ week",
                 A2="national week, region-specific $\\phi$", A3="division $\\times$ week, region-specific $\\phi$")
    body = []
    for k in ("N0", "R4", "A1", "A2", "A3"):
        s = r["specs"][k]
        phis = [v for kk, v in s["est"].items() if kk.startswith("phi")]
        ph = f(phis[0]) if len(phis) == 1 else "%s--%s" % (f(min(phis), 2), f(max(phis), 2))
        tpk = s["tau_plus_kappa"]
        ref = cs["reference"]["J_LM_oi"] if k == "N0" else (cs["grouped"][k]["J_LM_oi"] if k in cs["grouped"] else None)
        refs = "%s (%d/%d)" % (f(ref["q95"], 2), ref["exceed_count"] if "exceed_count" in ref else 0, ref["R"]) if ref else "--"
        eq = s.get("equal_phi_wald")
        body.append("%s & %s & %s %s & %s (%s) & %s (%s) & %s & %s \\\\" % (
            names[k], ph, f(tpk["est"]), ci(tpk["ci"]), f(s["J_LM_oi"], 2), pv(s["J_LM_oi_p"]),
            f(s["J_LM_inv"], 2), pv(s["J_LM_inv_p"]), refs, f(eq["p_F"], 3) if eq else "--"))
    note = ("Grouped extension of M0 (common/qdx.py). Week effects are profiled within the stated groups; $\\phi$ is common or "
            "region-specific (Northeast, Midwest, South, West). $J_{LM}$: score-form over-identification statistic for the two "
            "imported-pressure moments (oi) or the two invariance moments (inv), $\\chi^2_2$ $p$-values. Reference: 95th percentile "
            "of $J_{LM,oi}$ in the application-calibrated simulation (design \\emph{cal}, Web Table~\\ref{tab:rev-calsim}) and the "
            "number of draws at or above the application value / draws that converged. Equal-$\\phi$: $F$-test $p$-value of equal "
            "regional $\\phi$. 95\\% intervals: cluster (state) sandwich with $t_{48}$. Specifications R4--A3 were added in revision "
            "(post hoc). Source: verification/revision/regional_fe/results.json.")
    return table("lrlllll", "specification & $\\hat\\phi$ & $\\hat\\tau+\\hat\\kappa$ (95\\% CI) & $J_{LM,oi}$ ($p$) & "
                 "$J_{LM,inv}$ ($p$) & ref.\\ q95 (exceed) & equal-$\\phi$ $p$", body, note,
                 "Application under regional week effects and region-specific transmission scale.", "tab:rev-regional"), None, ["verification/revision/regional_fe/results.json", "verification/revision/calibrated_sim/results.json"]


# ---------------------------------------------------------------- g/b: calibrated simulation
def t_calsim():
    d = load("calibrated_sim")
    D = d["designs"]
    body = ["\\multicolumn{9}{l}{\\emph{(a) Interval performance, M0}} \\\\"]
    lab = dict(tau_plus_kappa="$\\tau+\\kappa$", tau="$\\tau$", kappa="$\\kappa$", phi="$\\phi$")
    for dz in ("cal", "stress"):
        for p in ("tau_plus_kappa", "tau", "kappa", "phi"):
            m = D[dz]["M0"][p]
            body.append("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\" % (
                dz if p == "tau_plus_kappa" else "", lab[p], f(m["bias"]), f(m["emp_sd"]), f(m["mean_se"]),
                f(m["cov_t48"]), f(m["cov_wild"]), f(m["cov_pairs"]),
                "%s / %s / %s" % (f(m["width_t48"], 2), f(m["width_wild"], 2), f(m["width_pairs"], 2))))
    body += ["\\midrule", "\\multicolumn{9}{l}{\\emph{(b) Reference distributions of the diagnostics}} \\\\",
             "design & statistic & mean & q95 & q99 & max & application & exceed & 95\\% upper \\\\", "\\midrule"]
    slab = dict(J="$J$", J_LM_inv="$J_{LM,inv}$", Q_G="$Q_G$", Q_GI="$Q_{GI}$")
    for dz in ("cal", "stress"):
        for s in ("J", "J_LM_inv", "Q_G", "Q_GI"):
            m = D[dz]["reference"][s]
            body.append("%s & %s & %s & %s & %s & %s & %s & %d/%d & %s \\\\" % (
                dz if s == "J" else "", slab[s], f(m["mean"], 2), f(m["q95"], 2), f(m["q99"], 2), f(m["max"], 2),
                f(m["application_value"], 2), m["exceed_count"], m["R"], f(m["exceed_upper95"], 3)))
    meta = D["cal"]["meta"]
    note = ("Application-calibrated simulation (added in revision, post hoc): observed 49-state panel, county counts, adoption "
            "dates, network and generation-interval weights; truth = M0 estimate ($\\tau$ %s, $\\kappa$ %s, $\\phi$ %s); NB2 noise "
            "($k$ = %s) and spatially correlated state $\\times$ week shocks, sd %s (\\emph{cal}, R = %d) or %s (\\emph{stress}, "
            "R = %d; chosen so that the standard error of $\\hat\\tau+\\hat\\kappa$ is close to the application's 0.311). Coverage of nominal 95\\%% intervals: cluster "
            "sandwich with $t_{48}$, score wild bootstrap (Mammen, B = %d), pairs cluster bootstrap (B = %d); MCSE about 0.010 "
            "(cal) and 0.013 (stress). Widths for $\\tau+\\kappa$ etc.\\ in the order $t_{48}$ / wild / pairs. Panel (b): "
            "\\emph{exceed} = draws at or above the application value; \\emph{95\\%% upper} = Clopper--Pearson upper limit of "
            "the exceedance probability. Seeds: data 20261118000+r, wild 20261118500000+r, pairs 20261118900000+r. Source: "
            "verification/revision/calibrated_sim/results.json.") % (
        f(meta["truth"]["tau"]), f(meta["truth"]["kappa"]), f(meta["truth"]["phi"]), f(meta["k_nb2"], 1),
        f(D["cal"]["meta"]["sigma_sw"], 2), D["cal"]["R_ok"], f(D["stress"]["meta"]["sigma_sw"], 2), D["stress"]["R_ok"],
        meta["B_wild"], meta["B_pairs"])
    hdr = "design & parameter & bias & emp.\\ SD & mean SE & cov.\\ $t_{48}$ & cov.\\ wild & cov.\\ pairs & width"
    return table("llrrrrrrr", hdr, body, note, "Application-calibrated simulation: interval coverage and reference "
                 "distributions of $J$ and $Q$.", "tab:rev-calsim"), D["cal"]["R_ok"], ["verification/revision/calibrated_sim/results.json"]


# ---------------------------------------------------------------- d: GI timing
def t_gi():
    d = load("gi_timing")["variants"]
    body = []
    for k, v in d.items():
        name = k.replace("_", "\\_")
        g = v["M0_G"]
        gdep = abs(g["phi"] - v["est"]["phi"]) > 1e-4
        body.append("%s & %s & %s & %s %s & %s (%s) & %s \\\\" % (
            name, f(v["same_week_mass"]), f(v["est"]["phi"]), f(v["tau_plus_kappa"]["est"]), ci(v["tau_plus_kappa"]["ci"]),
            f(v["J"], 2), pv(v["J_p"]), ("%s / %s" % (f(g["phi"], 4), f(g["tau_plus_kappa"]))) if gdep else "same as M0"))
    note = ("Generation interval of Ferretti et al.\\ (reference) or Nishiura et al., discretised to weeks. \\emph{base}: paper's "
            "discretisation (same-week mass dropped); \\emph{fold}: same-week mass moved into lag 1; \\emph{shift}: infector "
            "timing shifted by 1.75, 3.5 or 7 days before discretisation, which moves mass out of the same week. Same-week mass is "
            "the share of the generation interval within the same week. M0-G column: selected $\\hat\\phi$ / $\\hat\\tau+\\hat\\kappa$ "
            "where M0-G picks a root different from M0; M0-GI keeps the M0 root in all variants. $J$: Hansen statistic, "
            "$\\chi^2_2$. Variants added in revision (post hoc). Source: verification/revision/gi_timing/results.json.")
    return table("lrrlll", "variant & same-week mass & $\\hat\\phi$ & $\\hat\\tau+\\hat\\kappa$ (95\\% CI) & $J$ ($p$) & "
                 "M0-G $\\hat\\phi$ / $\\hat\\tau+\\hat\\kappa$", body, note,
                 "Application under alternative generation-interval discretisations.", "tab:rev-gi"), None, ["verification/revision/gi_timing/results.json"]


# ---------------------------------------------------------------- e: event-time weights
def t_event():
    d = load("event_time")
    est, se = d["M0_ET2"]["est"], d["M0_ET2"]["se"]
    bins = d["implicit_weights"]["bins"]
    w2 = d["implicit_weights"]["at_ET2_estimate"]["omega"]
    wp = d["implicit_weights"]["at_pooled_estimate"]["omega"]
    body = []
    for b, a, c in zip(bins, w2, wp):
        lab = "$e=%s$" % b[-1] + ("+" if b == bins[-1] else "")
        body.append("%s & %s (%s) & %s & %s \\\\" % (lab, f(est[b]), f(se[b]), f(a), f(c)))
    body += ["\\midrule", "\\multicolumn{4}{l}{\\emph{Pooled minus weighted bins}} \\\\"]
    cn = dict(implicit_omega_at_ET2="implicit weights at M0-ET2", implicit_omega_at_pooled="implicit weights at M0",
              equal_e0_e5="equal weights $e=0,\\dots,5$", cell_share="treated-cell shares")
    for k, lab in cn.items():
        c = d["contrasts"][k]
        body.append("%s & %s (%s) & \\multicolumn{2}{l}{$p$ = %s} \\\\" % (lab, f(c["pooled_minus_weighted"]), f(c["se"]), f(c["p_t"])))
    J = d["J"]
    note = ("Sample: %d counties in the %d states with event-time cells (identical for M0 and M0-ET2). Implicit weights $\\omega$ "
            "of the pooled estimand on the event-time effects, $\\hat\\tau_{pooled}\\approx\\sum_e\\omega_e\\tau_e$ "
            "($\\omega$ = $\\tau$-row of $A^{-1}H_b$, $\\sum_e\\omega_e = 1$; Web Appendix A.10), evaluated at the M0-ET2 and at the "
            "pooled estimates. Contrasts use the joint cluster sandwich of the stacked influence functions with $t_{%d}$. Pooled "
            "$\\hat\\tau$ = %s (SE %s). $J$ ($\\chi^2_2$): M0 %s ($p$ = %s), M0-ET2 %s ($p$ = %s). Added in revision (post hoc). "
            "Source: verification/revision/event_time/results.json.") % (
        d["sample"]["n_units"], d["sample"]["G"], d["sample"]["G"] - 1, f(d["pooled_M0"]["est"]["tau"]),
        f(d["pooled_M0"]["se"]["tau"]), f(J["M0"]["J"], 2), f(J["M0"]["p"], 4), f(J["M0_ET2"]["J"], 2), f(J["M0_ET2"]["p"], 4))
    return table("lrrr", "event time & $\\hat\\tau_e$ (SE) & $\\omega$ at M0-ET2 & $\\omega$ at M0", body, note,
                 "Implicit event-time weights of the pooled estimand and pooled--dynamic contrasts.", "tab:rev-event"), None, ["verification/revision/event_time/results.json"]


# ---------------------------------------------------------------- f: ascertainment
def t_asc():
    d = load("ascertainment")
    body = []
    for g in d["grid"]:
        if g.get("J") is None or g["J"] != g["J"]:
            continue
        body.append("%s & %s %s & %s & %s & %s \\\\" % (
            f(g["rho"], 2), f(g["tpk"]), ci([g["lo"], g["hi"]]), f(g["est"]["phi"]), f(g["J"], 2),
            ("near-zero root" if g["M0_G"]["chosen"] == "interior" and g["M0_G"]["phi"] < 0.01 else g["M0_G"]["chosen"]) + (" (%s)" % f(g["M0_G"]["tpk"]) if abs(g["M0_G"]["phi"] - g["est"]["phi"]) > 1e-4 else "")))
    body.append("\\midrule")
    for k, lab in (("rho_star_upper_limit_zero", "upper limit = 0"), ("rho_star_tpk_zero", "$\\hat\\tau+\\hat\\kappa$ = 0"),
                   ("rho_star_lower_limit_zero", "lower limit = 0")):
        g = d[k + "_fit"]
        body.append("%s (%s) & %s %s & %s & %s & \\\\" % (f(g["rho"], 3), lab, f(g["tpk"]), ci([g["lo"], g["hi"]]),
                                                          f(g["est"]["phi"]), f(g["J"], 2)))
    note = ("Observed counts $Y^{obs}_{it} = c_i c_t (1-\\rho)^{D_{it}} X_{it}$: treated counties report a share $1-\\rho$ of "
            "what they would report untreated ($\\rho>0$: reporting falls after the order). Counts are corrected by "
            "$(1-\\rho)^{-D}$, $L$ and $N$ rebuilt and M0 refitted. M0-G column: selected root type, with "
            "$\\hat\\tau+\\hat\\kappa$ where it differs from M0. Bottom rows: values of $\\rho$ at which the estimate or a 95\\% "
            "limit crosses zero. Added in revision (post hoc). Source: verification/revision/ascertainment/results.json.")
    return table("lllrl", "$\\rho$ & $\\hat\\tau+\\hat\\kappa$ (95\\% CI) & $\\hat\\phi$ & $J$ & M0-G", body, note,
                 "Application: sensitivity to a treatment-linked change in reporting.", "tab:rev-ascert"), None, ["verification/revision/ascertainment/results.json"]


# ---------------------------------------------------------------- c: robust sets
def t_robust():
    d = load("robust_sets")
    body = []
    for pb in ("0.25", "0.5", "1"):
        s = d["phi_bar_" + pb]["sets"]
        comp = s["tpk_by_component"]
        inter = [c for c in comp if c["phi"][0] > 0.05]
        body.append("%s & %s & %s & %s & %s \\\\" % (
            pb, ivs(s["S1_F_95"]), ivs(s["S3_F_95"]), ci(s["tpk_bonferroni_95"], 2),
            ci(inter[0]["tpk_bonferroni"], 2) if inter else "--"))
    s0 = d["phi_bar_0.5"]["S1_at_0"]
    note = ("Profiled Stock--Wright $S$ statistics on a grid of %d values of $\\phi$ in [0, 3] (step 0.01, plus points near 0), "
            "small-sample $F$ reference, $S(G-d)/\\{d(G-1)\\}\\sim F(d, G-d)$ with $d$ moments and $G = 49$. $S_1$: the three just-identifying moments; $S_3$: "
            "adds the two over-identifying (imported-pressure) moments. $\\tau+\\kappa$: Bonferroni projection (97.5\\%% set "
            "for $\\phi$, 97.5\\%% interval for $\\tau+\\kappa$ given $\\phi$), over the whole set and over its interior component. "
            "At $\\phi$ = 0 ($\\bar\\phi$ = 0.5): $\\hat\\tau$ = %s, $\\hat\\kappa$ = %s, $S_1$ $p$ = %s. Added in revision (post "
            "hoc). Source: verification/revision/robust_sets/results.json.") % (
        d["grid"]["n"], f(s0["tau"], 2), f(s0["kappa"], 2), f(s0["S1_p_F"], 2))
    return table("rllll", "$\\bar\\phi$ & $S_1$ 95\\% set for $\\phi$ & $S_3$ 95\\% set & $\\tau+\\kappa$ (all) & "
                 "$\\tau+\\kappa$ (interior)", body, note,
                 "Application: identification-robust confidence sets.", "tab:rev-robust"), None, ["verification/revision/robust_sets/results.json"]


# ---------------------------------------------------------------- h: population roots
def t_roots():
    d = load("population_roots")
    lab = {5: "S5: T 12, feedback", 6: "S6: T 26, feedback", 13: "S13: T 12, random", 14: "S14: T 26, random"}
    body = []
    for s in d["summary"]:
        rng = lambda v: "%s to %s" % (f(min(v)), f(max(v))) if v else "--"
        body.append("%s & %d & %s & %d/%d & %s & %s \\\\" % (
            lab.get(s["sid"], str(s["sid"])), s["N"], f(s["phi_bar"], 2), s["draws_with_spurious"], s["draws"],
            rng(s["spurious_phi"]), rng(s["spurious_tau"])))
    note = ("Zero-contagion designs ($\\phi_0=\\kappa_0=0$, $\\tau_0=\\log 0.8=-0.223$) simulated with the paper's DGP at large "
            "$N$ (seeds of the paper, scenarios.seed\\_for). The profiled equation for $\\phi$ is evaluated on a grid from $-0.1$ "
            "to 3; a sign change at $\\phi>0.05$ is a second root. Its persistence at large $N$ and its movement with "
            "$\\bar\\phi$ indicate a population root, i.e.\\ failure of global identification (A6). $\\hat\\tau$ at the second "
            "root is interpolated. Added in revision (post hoc). Source: verification/revision/population_roots/results.json.")
    return table("lrrrll", "design & $N$ & $\\bar\\phi$ & draws with 2nd root & $\\phi$ at 2nd root & $\\tau$ at 2nd root",
                 body, note, "Second roots of the profiled $\\phi$-equation in zero-contagion designs at large $N$.",
                 "tab:rev-poproots"), sum(s["draws"] for s in d["summary"]), ["verification/revision/population_roots/results.json"]


SPECS = [("tab:rev-regional", "table_11_regional.tex", t_regional, "point estimates, 95% CIs and test statistics (single observed data set); reference column from the calibrated simulation"),
         ("tab:rev-calsim", "table_12_calibrated_sim.tex", t_calsim, "Monte Carlo measures; n_seeds = cal design (stress design R = 300)"),
         ("tab:rev-gi", "table_13_gi_timing.tex", t_gi, "point estimates with 95% CIs (single observed data set)"),
         ("tab:rev-event", "table_14_event_weights.tex", t_event, "point estimates and SEs (single observed data set)"),
         ("tab:rev-ascert", "table_15_ascertainment.tex", t_asc, "point estimates with 95% CIs (single observed data set)"),
         ("tab:rev-robust", "table_16_robust_sets.tex", t_robust, "confidence sets (single observed data set)"),
         ("tab:rev-poproots", "table_17_population_roots.tex", t_roots, "roots per simulated draw; n_seeds = number of draws")]

EXISTING = [("tab:app-all", "table_7b_application_sensitivity.tex", "scripts/split_app_table.py",
             ["experiment/analysis/results.json"], None, "point estimates with 95% CIs (single observed data set); registered in revise-prep (was missing from the manifest)"),
            ("tab:rsel", "table_10_rsel.tex", "verification/root_selection/build_rsel_table.py",
             ["verification/root_selection"], None, "root-selection rerun summary; registered in revise-prep (was missing from the manifest)")]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    man = json.load(open(TAB / "manifest.json"))
    byid = {t["id"]: t for t in man["tables"]}
    for tid, fn, fun, note in SPECS:
        tex, n, src = fun()
        (TAB / fn).write_text(tex)
        byid[tid] = dict(id=tid, path=fn, sha256=sha(TAB / fn), generator=GEN, data_sources=src, display="raw",
                         n_seeds=n, std_ddof=1 if n else None, note=note + " (revision analysis, post hoc)")
        print("wrote", fn)
    for tid, fn, gen, src, n, note in EXISTING:
        if (TAB / fn).exists():
            byid[tid] = dict(id=tid, path=fn, sha256=sha(TAB / fn), generator=gen, data_sources=src, display="raw",
                             n_seeds=n, std_ddof=None, note=note)
    order = [t["id"] for t in man["tables"]]
    order += [k for k in byid if k not in order]
    man["tables"] = [byid[k] for k in order]
    (TAB / "manifest.json").write_text(json.dumps(man, indent=1))
    print("manifest entries", len(man["tables"]))


if __name__ == "__main__":
    main()
