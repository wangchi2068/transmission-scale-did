"""Generate all main-text and Supporting Information LaTeX tables for the
German Varicella + staggered summer school holidays application.

Reads:  experiment/analysis/germany_application_results.json
        experiment/analysis/germany_calibrated_sim.json
Writes: tables/table_7_application.tex
        tables/table_7b_application_sensitivity.tex
        tables/table_11_regional.tex
        tables/table_12_calibrated_sim.tex
        tables/table_13_gi_timing.tex
        tables/table_14_event_weights.tex
        tables/table_15_ascertainment.tex
        tables/table_16_robust_sets.tex
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
APP_PATH = ROOT / "experiment" / "analysis" / "germany_application_results.json"
CAL_PATH = ROOT / "experiment" / "analysis" / "germany_calibrated_sim.json"
TAB_DIR = ROOT / "tables"


def fmt_est_ci(est: float, lo: float, hi: float) -> str:
    return f"{est:.3f} ({lo:.3f}, {hi:.3f})"


def fmt_est_se_ci(est: float, se: float, tc: float) -> str:
    return fmt_est_ci(est, est - tc * se, est + tc * se)


def fmt_ivs(ivs: list) -> str:
    if not ivs:
        return "empty"
    return r" $\cup$ ".join(f"[{a:.2f}, {b:.2f}]" for a, b in ivs)


def write_table_7(app: dict):
    m0 = app["primary"]["M0"]
    tc = m0["tcrit"]
    g_m0g = app["global_roots"]["M0_G"]
    g_m0gi = app["global_roots"]["M0_GI"]
    comps = app["comparators"]
    nuts = app["nuts2"]["M0"]
    wild_m = app["primary"]["wild_mammen"]["ci_tpk"]
    pairs = app["primary"]["pairs_bootstrap"]
    jk = app["primary"]["jackknife"]
    sw05 = app["stock_wright"]["phi_bar_0.5"]["sets"]["tpk_bonferroni_95"]
    psi = app["primary"]["psi"]

    rows = [
        ("Primary (400 counties)", "M0",
         fmt_est_ci(m0["est"]["tpk"], *m0["ci_t"]["tpk"]),
         fmt_est_ci(m0["est"]["tau"], *m0["ci_t"]["tau"]),
         fmt_est_ci(m0["est"]["kappa"], *m0["ci_t"]["kappa"]),
         fmt_est_ci(m0["est"]["phi"], *m0["ci_t"]["phi"])),
        ("Primary (400 counties)", "M0-G",
         fmt_est_se_ci(g_m0g["est"]["tpk"], g_m0g["se"]["tpk"], tc),
         fmt_est_se_ci(g_m0g["est"]["tau"], g_m0g["se"]["tau"], tc),
         fmt_est_se_ci(g_m0g["est"]["kappa"], g_m0g["se"]["kappa"], tc),
         fmt_est_se_ci(g_m0g["est"]["phi"], g_m0g["se"]["phi"], tc)),
        ("Primary (400 counties)", "M0-GI",
         fmt_est_se_ci(g_m0gi["est"]["tpk"], g_m0gi["se"]["tpk"], tc),
         fmt_est_se_ci(g_m0gi["est"]["tau"], g_m0gi["se"]["tau"], tc),
         fmt_est_se_ci(g_m0gi["est"]["kappa"], g_m0gi["se"]["kappa"], tc),
         fmt_est_se_ci(g_m0gi["est"]["phi"], g_m0gi["se"]["phi"], tc)),
        ("Primary (400 counties)", r"M0$'$",
         "--",
         fmt_est_se_ci(comps["M0_prime"]["est"]["tau"], comps["M0_prime"]["se"]["tau"], tc),
         "--", "--"),
        ("Primary (400 counties)", "M1 (Poisson FE)",
         fmt_est_ci(comps["M1"]["est"]["tpk"], *comps["M1"]["ci_t_tpk"]),
         fmt_est_se_ci(comps["M1"]["est"]["tau"], comps["M1"]["se"]["tau"], tc),
         fmt_est_se_ci(comps["M1"]["est"]["kappa"], comps["M1"]["se"]["kappa"], tc),
         fmt_est_se_ci(comps["M1"]["est"]["phi"], comps["M1"]["se"]["phi"], tc)),
        ("Primary (400 counties)", r"M2 (EpiEstim $\log R$)",
         "--",
         fmt_est_ci(comps["M2"]["est"]["tau"], *comps["M2"]["ci_t_tau"]),
         "--", "--"),
        ("Primary (400 counties)", r"M3 ($\log(Y+1)$ TWFE)",
         "--",
         fmt_est_ci(comps["M3"]["est"]["tau"], *comps["M3"]["ci_t_tau"]),
         "--", "--"),
        ("Primary (400 counties)", "M4 (Pooled QMLE)",
         fmt_est_ci(comps["M4"]["est"]["tpk"], *comps["M4"]["ci_t_tpk"]),
         fmt_est_se_ci(comps["M4"]["est"]["tau"], comps["M4"]["se"]["tau"], tc),
         fmt_est_se_ci(comps["M4"]["est"]["kappa"], comps["M4"]["se"]["kappa"], tc),
         fmt_est_se_ci(comps["M4"]["est"]["phi"], comps["M4"]["se"]["phi"], tc)),
        ("NUTS-2 (38 regions)", "M0",
         fmt_est_ci(nuts["est"]["tpk"], *nuts["ci_t"]["tpk"]),
         fmt_est_ci(nuts["est"]["tau"], *nuts["ci_t"]["tau"]),
         fmt_est_ci(nuts["est"]["kappa"], *nuts["ci_t"]["kappa"]),
         fmt_est_ci(nuts["est"]["phi"], *nuts["ci_t"]["phi"])),
    ]
    body = "\n".join(" & ".join(r) + r" \\" for r in rows)
    pct_red = (1.0 - float(__import__("math").exp(m0["est"]["tpk"]))) * 100.0
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Application: staggered German summer school holidays and county-level varicella transmission (2015--2019); point estimates with 95\% confidence intervals ($t_{{15}}$, $G = 16$ federal states).}}
\label{{tab:app}}
\footnotesize
\setlength{{\tabcolsep}}{{3.5pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{llrrrr}}
\toprule
sample & method & $\hat\tau+\hat\kappa$ (95\% CI) & $\hat\tau$ (95\% CI) & $\hat\kappa$ (95\% CI) & $\hat\phi$ (95\% CI) \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Weekly RKI SurvStat varicella cases across 400 German counties in $G = 16$ federal states, 261 estimation weeks (ISO 2015-W01 to 2019-W52, 115{{,}}129 cases, plus three pre-sample lag weeks from late 2014); $D$ = state summer school holiday indicator lagged one week to align with varicella's 14-day generation/incubation interval ($w = (0.141, 0.721, 0.138)$, same-week mass $w_0 = 2.6 \times 10^{{-6}}$). 95\% intervals use the state-clustered sandwich with $t_{{15}} = {tc:.2f}$. Because holidays are set at the state level, $\tau+\kappa$ is the policy-level total for interior counties (primary M0 estimate ${m0['est']['tpk']:.3f}$, a {pct_red:.1f}\% reduction in $R$). Alternative 95\% intervals for M0 $\tau+\kappa$: Mammen score wild bootstrap-$t$ ($B=999$) $[{wild_m[0]:.3f}, {wild_m[1]:.3f}]$; nonlinear pairs cluster bootstrap-$t$ ({pairs['ok']}/{pairs['B']} converged) $[{pairs['ci_tpk'][0]:.3f}, {pairs['ci_tpk'][1]:.3f}]$; delete-one-state jackknife $[{jk['ci_tpk'][0]:.3f}, {jk['ci_tpk'][1]:.3f}]$ (leave-one-state-out range $[{jk['tpk_range'][0]:.3f}, {jk['tpk_range'][1]:.3f}]$); identification-robust Stock--Wright Bonferroni projection $[{sw05[0]:.3f}, {sw05[1]:.3f}]$. Propagated cases averted $\hat\psi^\dagger$: {psi['treated']['est']:,.0f} (SE {psi['treated']['se']:,.0f}) during treated county-weeks and {psi['network']['est']:,.0f} (SE {psi['network']['se']:,.0f}) network-wide. Global search finds a unique interior root at $\hat\phi = {m0['est']['phi']:.3f}$, so M0, M0-G and M0-GI coincide. Data sources and coding are in Web Appendix D.}}
\end{{table}}
"""
    (TAB_DIR / "table_7_application.tex").write_text(tex, encoding="utf-8")


def write_table_7b(app: dict):
    m0 = app["primary"]["M0"]
    s = app["sensitivities"]
    by_yr = app["by_year"]
    ah = app["all_holidays"]
    nuts = app["nuts2"]["M0"]

    specs = [
        (r"Primary (summer, lag 1, $\bar\phi=0.5$)", m0),
        (r"Instrument ref.\ $\bar\phi = 0.25$", s["phi_bar_0.25"]),
        (r"Instrument ref.\ $\bar\phi = 1.00$", s["phi_bar_1"]),
        (r"Background $\omega_i \times 0.25$", s["omega_x_0.25"]),
        (r"Background $\omega_i \times 4.0$", s["omega_x_4"]),
        ("Summer holiday lag 0 (contemp.)", s["summer_lag_0"]),
        ("Summer holiday lag 2", s["summer_lag_2"]),
        ("Summer holiday lag 3", s["summer_lag_3"]),
        ("Summer holiday lag 4", s["summer_lag_4"]),
        ("Continuous holiday share (lag 1)", s["summer_Dshare_lag1"]),
        ("Single year: 2015 ($T=53$)", by_yr["2015"]),
        ("Single year: 2016 ($T=52$)", by_yr["2016"]),
        ("Single year: 2017 ($T=52$)", by_yr["2017"]),
        ("Single year: 2018 ($T=52$)", by_yr["2018"]),
        ("Single year: 2019 ($T=52$)", by_yr["2019"]),
        ("All school holidays 2015--2019 (lag 0)", ah["all_2015_2019_lag0"]),
        ("All school holidays 2015--2019 (lag 1)", ah["all_2015_2019_lag1"]),
        ("All school holidays 2017--2018 (lag 0)", ah["all_2017_2018_lag0"]),
        ("NUTS-2 regions ($N=38$, 2.9\% zeros)", nuts),
    ]
    lines = []
    for lbl, r in specs:
        lines.append(" & ".join([
            lbl,
            fmt_est_ci(r["est"]["tpk"], *r["ci_t"]["tpk"]),
            fmt_est_ci(r["est"]["tau"], *r["ci_t"]["tau"]),
            fmt_est_ci(r["est"]["kappa"], *r["ci_t"]["kappa"]),
            fmt_est_ci(r["est"]["phi"], *r["ci_t"]["phi"]),
            f"{r['J_LM_oi']['stat']:.2f} ({r['J_LM_oi']['p_F']:.3f})",
            f"{r['J_LM_inv']['stat']:.2f} ({r['J_LM_inv']['p_F']:.3f})",
        ]) + r" \\")
    body = "\n".join(lines)
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Sensitivity analyses of the German varicella application (M0, state-clustered $t_{{15}}$ 95\% CIs and score-form over-identification tests $J_{{LM}}$ with $F_{{2,14}}$ $p$-values).}}
\label{{tab:app-all}}
\footnotesize
\setlength{{\tabcolsep}}{{3.0pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\renewcommand{{\arraystretch}}{{0.88}}\begin{{tabular}}{{lrrrrrr}}
\toprule
specification & $\hat\tau+\hat\kappa$ (95\% CI) & $\hat\tau$ (95\% CI) & $\hat\kappa$ (95\% CI) & $\hat\phi$ (95\% CI) & $J_{{LM,\mathrm{{oi}}}}$ ($p_F$) & $J_{{LM,\mathrm{{inv}}}}$ ($p_F$) \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Primary specification is Table~\ref{{tab:app}}. $J_{{LM,\mathrm{{oi}}}}$ tests the two lagged predetermined moments (lagged imported share and lagged one-week log growth); $J_{{LM,\mathrm{{inv}}}}$ tests reference-value invariance ($\phi_r \in \{{0.05, 2.0\}}$); $p_F$ uses the finite-cluster $F_{{2,14}}$ reference $J_{{LM}}(G-2)/\{{\,2(G-1)\,\}}$. Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_7b_application_sensitivity.tex").write_text(tex, encoding="utf-8")


def write_table_11(app: dict):
    m0 = app["primary"]["M0"]
    reg = app["regional_fe"]
    specs = [
        ("National week FE (primary, lag 1)", m0),
        (r"East/West $\times$ week FE (lag 1)", reg["east_west_x_week_lag1"]),
        (r"4-region $\times$ week FE (lag 1)", reg["four_region_x_week_lag1"]),
        (r"4-region $\times$ week FE (lag 2)", reg["four_region_x_week_lag2"]),
    ]
    lines = []
    for lbl, r in specs:
        lines.append(" & ".join([
            lbl,
            fmt_est_ci(r["est"]["tpk"], *r["ci_t"]["tpk"]),
            fmt_est_ci(r["est"]["tau"], *r["ci_t"]["tau"]),
            fmt_est_ci(r["est"]["kappa"], *r["ci_t"]["kappa"]),
            fmt_est_ci(r["est"]["phi"], *r["ci_t"]["phi"]),
            f"{r['J_LM_oi']['stat']:.2f} ({r['J_LM_oi']['p_F']:.3f})",
            f"{r['J_LM_inv']['stat']:.2f} ({r['J_LM_inv']['p_F']:.3f})",
        ]) + r" \\")
    body = "\n".join(lines)
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{German varicella application with regional $\times$ week fixed effects (M0, $G = 16$ federal states, $t_{{15}}$ 95\% CIs).}}
\label{{tab:rev-regional}}
\footnotesize
\setlength{{\tabcolsep}}{{3.2pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{lrrrrrr}}
\toprule
week-effect structure & $\hat\tau+\hat\kappa$ (95\% CI) & $\hat\tau$ (95\% CI) & $\hat\kappa$ (95\% CI) & $\hat\phi$ (95\% CI) & $J_{{LM,\mathrm{{oi}}}}$ ($p_F$) & $J_{{LM,\mathrm{{inv}}}}$ ($p_F$) \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Four German macro-regions: North (SH, HH, NI, HB, MV), West (NW, HE, RP, SL), South (BW, BY) and East (BE, BB, SN, ST, TH); East/West contrasts the five eastern states plus Berlin against the ten western states. Week effects are profiled in closed form within each region. Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_11_regional.tex").write_text(tex, encoding="utf-8")


def write_table_12(cal: dict):
    lines_a = []
    for tag, lbl in [("calibrated", "cal"), ("high_noise", "high-noise")]:
        block = cal[tag]["params"]
        for i, (pk, plbl) in enumerate([("tpk", r"$\tau+\kappa$"), ("tau", r"$\tau$"), ("kappa", r"$\kappa$"), ("phi", r"$\phi$")]):
            p = block[pk]
            lines_a.append(" & ".join([
                lbl if i == 0 else "",
                plbl,
                f"{p['bias']:.3f}",
                f"{p['sd']:.3f}",
                f"{p['mean_se']:.3f}",
                f"{p['cov_t15']:.3f}",
                f"{p['cov_wild_mam']:.3f}",
                f"{p['cov_wild_web']:.3f}",
            ]) + r" \\")
    lines_b = []
    for tag, lbl in [("calibrated", "cal"), ("high_noise", "high-noise")]:
        block = cal[tag]["diagnostics"]
        for i, (dk, dlbl) in enumerate([
            ("j_oi", r"$J_{LM,\mathrm{oi}}$"),
            ("j_inv", r"$J_{LM,\mathrm{inv}}$"),
            ("s3_oi", r"$S_{3,\mathrm{oi}}$"),
            ("s3_inv", r"$S_{3,\mathrm{inv}}$"),
        ]):
            d = block[dk]
            lines_b.append(" & ".join([
                lbl if i == 0 else "",
                dlbl,
                f"{d['mean']:.2f}",
                f"{d['q95']:.2f}",
                f"{d['q99']:.2f}",
                f"{d['max']:.2f}",
                f"{d['obs']:.2f}",
                f"{d['exceed_count']}/{d['n']} ({d['p_calibrated']:.3f})",
            ]) + r" \\")
    c_cal = cal["calibration"]
    r_tot = cal["R"]
    rc_cal = cal["calibrated"]["root_counts"]
    rc_hn = cal["high_noise"]["root_counts"]
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Application-calibrated simulation on the German varicella design ($N=400$, $G=16$ states, $T=261$, $R={r_tot}$): finite-cluster interval coverage and null reference distributions of $J_{{LM}}$ and $S_3$.}}
\label{{tab:rev-calsim}}
\footnotesize
\setlength{{\tabcolsep}}{{3.5pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{llrrrrrr}}
\toprule
design & parameter & bias & emp.\ SD & mean SE & cov.\ $t_{{15}}$ & cov.\ wild (M) & cov.\ wild (W6) \\
\midrule
\multicolumn{{8}}{{l}}{{\emph{{(a) Interval performance, M0}}}} \\
{chr(10).join(lines_a)}
\midrule
\multicolumn{{8}}{{l}}{{\emph{{(b) Finite-sample null distributions of over-identification and root-selection statistics}}}} \\
design & statistic & mean & q95 & q99 & max & application & exceed ($p_{{\mathrm{{cal}}}}$) \\
\midrule
{chr(10).join(lines_b)}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Calibrated DGP uses the observed 400-county German network, 261-week summer holiday schedule, Vink et al.\ (2014) generation weights, truth set to the primary M0 estimate ($\tau={c_cal['beta_true'][0]:.3f}$, $\kappa={c_cal['beta_true'][1]:.3f}$, $\phi={c_cal['beta_true'][2]:.3f}$), NB2 shape $k={c_cal['k_nb']:.2f}$ and state $\times$ week log-normal shock SD $\sigma={c_cal['sig_sw']:.3f}$ (\emph{{cal}}), or $k/2$ and $2\sigma$ (\emph{{high-noise}}). Unique interior root in {100*rc_cal['share_unique']:.1f}\% (\emph{{cal}}) and {100*rc_hn['share_unique']:.1f}\% (\emph{{high-noise}}) of draws. Coverage uses nominal 95\% $t_{{15}}$ sandwich and score wild bootstrap-$t$ with Mammen (M) and Webb 6-point (W6) weights ($B=399$; MCSE $\approx 0.010$). Panel (b): $p_{{\mathrm{{cal}}}} = (\text{{exceed}} + 1)/(R_{{\mathrm{{ok}}}} + 1)$. Source: \texttt{{experiment/analysis/germany\_calibrated\_sim.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_12_calibrated_sim.tex").write_text(tex, encoding="utf-8")


def write_table_13(app: dict):
    gis = app["gi_sensitivity"]
    order = [
        ("Vink_main_14d_sd2.5", "Vink et al.\ (2014) primary (14.0 d, SD 2.5)"),
        ("mean_12d_sd2.5", "Shorter mean (12.0 d, SD 2.5)"),
        ("mean_16d_sd2.5", "Longer mean (16.0 d, SD 2.5)"),
        ("mean_14d_sd4.0", "Wider dispersion (14.0 d, SD 4.0)"),
        ("mean_14d_sd6.0_report", "Reporting-convolved (14.0 d, SD 6.0)"),
    ]
    lines = []
    for k, lbl in order:
        rec = gis[k]
        g = rec["gi"]
        f = rec["fit"]
        w_str = f"({g['w'][0]:.3f}, {g['w'][1]:.3f}, {g['w'][2]:.3f})"
        w0_str = f"{g['w0']:.1e}" if g["w0"] < 1e-3 else f"{g['w0']:.4f}"
        lines.append(" & ".join([
            lbl, w_str, w0_str,
            fmt_est_ci(f["est"]["tpk"], *f["ci_t"]["tpk"]),
            fmt_est_ci(f["est"]["phi"], *f["ci_t"]["phi"]),
            f"{f['J_LM_oi']['stat']:.2f} ({f['J_LM_oi']['p_F']:.3f})",
            f"{f['J_LM_inv']['stat']:.2f} ({f['J_LM_inv']['p_F']:.3f})",
        ]) + r" \\")
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Generation-interval sensitivity for the German varicella application (M0, summer holidays lagged one week).}}
\label{{tab:rev-gi}}
\footnotesize
\setlength{{\tabcolsep}}{{3.0pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{lccrrrr}}
\toprule
generation-interval specification & $(w_1, w_2, w_3)$ & $w_0$ & $\hat\tau+\hat\kappa$ (95\% CI) & $\hat\phi$ (95\% CI) & $J_{{LM,\mathrm{{oi}}}}$ ($p_F$) & $J_{{LM,\mathrm{{inv}}}}$ ($p_F$) \\
\midrule
{chr(10).join(lines)}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Lognormal generation intervals discretised onto weekly lags $s=1,2,3$ under a uniform within-week infection day; $w_0$ is the resulting same-week transmission share. All specifications yield a unique interior root and a significant transmission reduction. Under the broader reporting-convolved dispersion (SD 6.0 d), both $J_{{LM,\mathrm{{oi}}}}$ and $J_{{LM,\mathrm{{inv}}}}$ pass at the 5\% $F_{{2,14}}$ level. Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_13_gi_timing.tex").write_text(tex, encoding="utf-8")


def write_table_14(app: dict):
    evL = app["event_time_with_leads"]
    evN = app["event_time_no_leads"]
    ks = list(range(-2, 8))
    lines = []
    for k in ks:
        nm = f"e{k}"
        cL = evL["coefs"][nm]
        sL = fmt_est_ci(cL["est"], *cL["ci"])
        if nm in evN["coefs"]:
            cN = evN["coefs"][nm]
            sN = fmt_est_ci(cN["est"], *cN["ci"])
        else:
            sN = "--"
        role = "Pre-trend lead" if k < 0 else ("Onset week (incubation)" if k == 0 else ("Holiday weeks (post-incubation)" if k <= 5 else "Reopening / tail"))
        lines.append(f"$e = {k:+d}$ & {role} & {sL} & {sN} \\\\")
    for nm, lbl in [("kappa", r"Neighbour holiday ($\hat\kappa$)"), ("phi", r"Contagion coupling ($\hat\phi$)")]:
        cL = evL["coefs"][nm]
        cN = evN["coefs"][nm]
        lines.append(f"{lbl} & Control / structural & {fmt_est_ci(cL['est'], *cL['ci'])} & {fmt_est_ci(cN['est'], *cN['ci'])} \\\\")
    pw = evL["pretrend_wald"]
    postw = evL["post_wald_e1_e5"]
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Event-time dynamic profile around German summer school holiday onset ($e = -2, \dots, 7$ weeks, controlling for neighbour holiday share $\bar D_{{it}}$, $t_{{15}}$ 95\% CIs).}}
\label{{tab:rev-event}}
\footnotesize
\setlength{{\tabcolsep}}{{4.0pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{llrr}}
\toprule
event week / parameter & biological / timing role & with pre-trend leads & without leads \\
\midrule
{chr(10).join(lines)}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Indicators $\mathbf 1\{{t - t_{{0,i}} = e\}}$ mark weeks relative to each state's summer holiday onset $t_{{0,i}}$ across 2015--2019 (typical holiday duration 6 weeks). Pre-trend joint Wald test ($H_0\colon \tau_{{-2}}=\tau_{{-1}}=0$): $\chi^2_2 = {pw['stat']:.2f}$ ($p = {pw['p_chi2']:.3f}$, $F_{{2,14}}$ $p = {pw['p_F']:.3f}$). Onset week $e=0$ is indistinguishable from zero, ruling out an immediate holiday under-reporting drop, followed by a post-incubation decline across $e=1,\dots,5$ (joint Wald $\chi^2_5 = {postw['stat']:.2f}$, $p = {postw['p_chi2']:.3f}$, $F_{{5,11}}$ $p = {postw['p_F']:.3f}$). Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_14_event_weights.tex").write_text(tex, encoding="utf-8")


def write_table_15(app: dict):
    asc = app["ascertainment"]
    keep_rhos = {-0.20, -0.10, 0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40}
    rows_g = [r for r in asc["grid"] if r["ok"] and round(r["rho"], 2) in keep_rhos]
    lines = []
    for r in rows_g:
        lines.append(" & ".join([
            f"{r['rho']:+.2f}",
            f"{r['bound_log_shift']:.3f}",
            fmt_est_ci(r["est"]["tpk"], r["lo"], r["hi"]),
            f"{r['est']['phi']:.3f} ({r['se']['phi']:.3f})",
        ]) + r" \\")
    r_hi = asc["rho_star_upper_ci_zero"]
    r_0 = asc["rho_star_tpk_zero"]
    if r_hi is not None:
        f_hi = asc["rho_star_upper_ci_zero_fit"]
        lines.append(rf"$\rho^*_{{\mathrm{{CI}}}} = {r_hi:+.3f}$ & {f_hi['bound_log_shift']:.3f} & {fmt_est_ci(f_hi['est']['tpk'], f_hi['lo'], f_hi['hi'])} & {f_hi['est']['phi']:.3f} ({f_hi['se']['phi']:.3f}) \\")
    if r_0 is not None:
        f_0 = asc["rho_star_tpk_zero_fit"]
        lines.append(rf"$\rho^*_0 = {r_0:+.3f}$ & {f_0['bound_log_shift']:.3f} & {fmt_est_ci(f_0['est']['tpk'], f_0['lo'], f_0['hi'])} & {f_0['est']['phi']:.3f} ({f_0['se']['phi']:.3f}) \\")
    note_hi = f"{100*r_hi:.1f}\\%" if r_hi is not None else "beyond grid"
    note_0 = f"{100*r_0:.1f}\\%" if r_0 is not None else "beyond grid"
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Holiday-linked ascertainment tipping-point analysis: M0 refitted after undoing a proportional reporting change $(1-\rho)^{{D_{{it}}}}$ during summer holidays.}}
\label{{tab:rev-ascert}}
\footnotesize
\setlength{{\tabcolsep}}{{4.0pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{rrrr}}
\toprule
$\rho$ & $|\log(1-\rho)|$ & $\hat\tau+\hat\kappa$ (95\% $t_{{15}}$ CI) & $\hat\phi$ (SE) \\
\midrule
{chr(10).join(lines)}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Under $Y_{{it}} = c_i c_t (1-\rho)^{{D_{{it}}}} I_{{it}}$, $\rho > 0$ represents reduced case ascertainment during summer holidays (e.g.\ fewer clinic visits while travelling); counts are corrected via $Y^*_{{it}} = Y_{{it}} / (1-\rho)^{{D_{{it}}}}$ before rebuilding infection pressures and refitting M0. The 95\% upper confidence limit crosses zero only at $\rho^*_{{\mathrm{{CI}}}} =$ {note_hi}, and the point estimate reaches zero only at $\rho^*_0 =$ {note_0}. Combined with the null onset-week estimate ($e=0$, Web Table~\ref{{tab:rev-event}}), plausible holiday reporting fluctuations cannot account for the estimated transmission drop. Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_15_ascertainment.tex").write_text(tex, encoding="utf-8")


def write_table_16(app: dict):
    sw = app["stock_wright"]
    lines = []
    for pb in (0.25, 0.5, 1.0):
        st = sw[f"phi_bar_{pb:g}"]["sets"]
        bp = st["tpk_bonferroni_95"]
        lines.append(" & ".join([
            f"{pb:g}",
            fmt_ivs(st["S1_F_95"]),
            fmt_ivs(st["S1_chi2_95"]),
            fmt_ivs(st["S3_inv_chi2_95"]),
            fmt_ivs(st["S3_oi_F_95"]),
            f"({bp[0]:.3f}, {bp[1]:.3f})",
        ]) + r" \\")
    s0 = sw["phi_bar_0.5"]["S1_at_0"]
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{German varicella application: identification-robust profiled Stock--Wright confidence sets for $\phi$ and Bonferroni 95\% projection for $\tau+\kappa$.}}
\label{{tab:rev-robust}}
\footnotesize
\setlength{{\tabcolsep}}{{3.2pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{rlllll}}
\toprule
$\bar\phi$ & $S_1(F)$ 95\% set & $S_1(\chi^2)$ 95\% set & $S_{{3,\mathrm{{inv}}}}(\chi^2)$ 95\% set & $S_{{3,\mathrm{{oi}}}}(F)$ 95\% set & $\tau+\kappa$ (Bonferroni 95\%) \\
\midrule
{chr(10).join(lines)}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Profiled Stock--Wright $S$ statistics evaluated on 304 values of $\phi \in [0, 3]$ ($G = 16$ states), with week effects and $(\tau, \kappa)$ profiled out via their just-identifying moments (presuming strong identification of $(\delta, \tau, \kappa)$ at the tested $\phi$). Small-sample $F$ reference uses $S(G-d)/\{{\,d(G-1)\,\}} \sim F_{{d, G-d}}$. Across all three reference values $\bar\phi$, $S_1(F)$ and $S_1(\chi^2)$ form a single connected interior interval excluding $\phi = 0$ (at $\phi=0$, $\bar\phi=0.5$: $S_1 = {s0['S1']:.2f}$, $p_F = {s0['S1_p_F']:.3f}$), both $S_{{3,\mathrm{{inv}}}}$ and $S_{{3,\mathrm{{oi}}}}(F)$ are non-empty, and the Bonferroni 95\% projection for $\tau+\kappa$ (union of 97.5\% conditional $t_{{15}}$ intervals over the 97.5\% $S_1(F)$ set for $\phi$) strictly excludes zero. Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_16_robust_sets.tex").write_text(tex, encoding="utf-8")


def main():
    app = json.loads(APP_PATH.read_text(encoding="utf-8"))
    cal = json.loads(CAL_PATH.read_text(encoding="utf-8"))
    write_table_7(app)
    write_table_7b(app)
    write_table_11(app)
    write_table_12(cal)
    write_table_13(app)
    write_table_14(app)
    write_table_15(app)
    write_table_16(app)
    print("Generated all 8 German Varicella application LaTeX tables in tables/.")


if __name__ == "__main__":
    main()
