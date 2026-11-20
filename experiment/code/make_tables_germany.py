"""Generate all main-text and Supporting Information LaTeX tables for the
German Varicella empirical stress test (2015--2019).

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
        tables/table_18_placebos.tex
"""
from __future__ import annotations

import json
import math
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


def _row_m0(lbl: str, r: dict) -> str:
    if not r.get("ok"):
        return f"{lbl} & \\multicolumn{{6}}{{c}}{{did not converge}} \\\\"
    return " & ".join([
        lbl,
        fmt_est_ci(r["est"]["tpk"], *r["ci_t"]["tpk"]),
        fmt_est_ci(r["est"]["tau"], *r["ci_t"]["tau"]),
        fmt_est_ci(r["est"]["kappa"], *r["ci_t"]["kappa"]),
        fmt_est_ci(r["est"]["phi"], *r["ci_t"]["phi"]),
        f"{r['J_LM_oi']['stat']:.2f} ({r['J_LM_oi']['p_F']:.3f})",
        f"{r['J_LM_inv']['stat']:.2f} ({r['J_LM_inv']['p_F']:.3f})",
    ]) + r" \\"


def write_table_7(app: dict):
    m0 = app["primary"]["M0"]
    tc = m0["tcrit"]
    g_m0g = app["global_roots"]["M0_G"]
    g_m0gi = app["global_roots"]["M0_GI"]
    comps = app["comparators"]
    nuts = app["nuts2"]["M0"]
    wcr_w = app["primary"]["wcr_webb6"]
    wcr_r = app["primary"]["wcr_rademacher"]
    wild_m = app["primary"]["wild_mammen"]["ci_tpk"]
    pairs = app["primary"]["pairs_bootstrap"]
    jk = app["primary"]["jackknife"]
    sw05 = app["stock_wright"]["phi_bar_0.5"]["sets"]["tpk_bonferroni_95"]
    psi = app["primary"]["psi"]

    nuts_m2 = app["nuts2"]["comparators"]["M2"]
    nuts_m2_ci = (
        nuts_m2["est"]["tau"] - tc * nuts_m2["se"]["tau"],
        nuts_m2["est"]["tau"] + tc * nuts_m2["se"]["tau"],
    )

    def _phi_cell(rec):
        if rec["kind"] == "boundary":
            return "0.000 (boundary)"
        return fmt_est_se_ci(rec["est"]["phi"], rec["se"]["phi"], tc)

    rows = [
        ("Primary (58.8\\% zeros)", "M0",
         fmt_est_ci(m0["est"]["tpk"], *m0["ci_t"]["tpk"]),
         fmt_est_ci(m0["est"]["tau"], *m0["ci_t"]["tau"]),
         fmt_est_ci(m0["est"]["kappa"], *m0["ci_t"]["kappa"]),
         fmt_est_ci(m0["est"]["phi"], *m0["ci_t"]["phi"])),
        ("Primary (58.8\\% zeros)", "M0-G",
         fmt_est_se_ci(g_m0g["est"]["tpk"], g_m0g["se"]["tpk"], tc),
         fmt_est_se_ci(g_m0g["est"]["tau"], g_m0g["se"]["tau"], tc),
         fmt_est_se_ci(g_m0g["est"]["kappa"], g_m0g["se"]["kappa"], tc),
         _phi_cell(g_m0g)),
        ("Primary (58.8\\% zeros)", "M0-GI",
         fmt_est_se_ci(g_m0gi["est"]["tpk"], g_m0gi["se"]["tpk"], tc),
         fmt_est_se_ci(g_m0gi["est"]["tau"], g_m0gi["se"]["tau"], tc),
         fmt_est_se_ci(g_m0gi["est"]["kappa"], g_m0gi["se"]["kappa"], tc),
         _phi_cell(g_m0gi)),
        ("Primary (58.8\\% zeros)", r"M0$'$",
         "--",
         fmt_est_se_ci(comps["M0_prime"]["est"]["tau"], comps["M0_prime"]["se"]["tau"], tc),
         "--", "--"),
        ("Primary (58.8\\% zeros)", "M1 (Poisson FE)",
         fmt_est_ci(comps["M1"]["est"]["tpk"], *comps["M1"]["ci_t_tpk"]),
         fmt_est_se_ci(comps["M1"]["est"]["tau"], comps["M1"]["se"]["tau"], tc),
         fmt_est_se_ci(comps["M1"]["est"]["kappa"], comps["M1"]["se"]["kappa"], tc),
         fmt_est_se_ci(comps["M1"]["est"]["phi"], comps["M1"]["se"]["phi"], tc)),
        ("Primary (58.8\\% zeros)", r"M2 (EpiEstim $\log R$)",
         "--",
         fmt_est_ci(comps["M2"]["est"]["tau"], *comps["M2"]["ci_t_tau"]),
         "--", "--"),
        ("Primary (58.8\\% zeros)", r"M3 ($\log(Y+1)$ TWFE)",
         "--",
         fmt_est_ci(comps["M3"]["est"]["tau"], *comps["M3"]["ci_t_tau"]),
         "--", "--"),
        ("Primary (58.8\\% zeros)", "M4 (Pooled QMLE)",
         fmt_est_ci(comps["M4"]["est"]["tpk"], *comps["M4"]["ci_t_tpk"]),
         fmt_est_se_ci(comps["M4"]["est"]["tau"], comps["M4"]["se"]["tau"], tc),
         fmt_est_se_ci(comps["M4"]["est"]["kappa"], comps["M4"]["se"]["kappa"], tc),
         fmt_est_se_ci(comps["M4"]["est"]["phi"], comps["M4"]["se"]["phi"], tc)),
        ("NUTS-2 (2.9\\% zeros)", "M0",
         fmt_est_ci(nuts["est"]["tpk"], *nuts["ci_t"]["tpk"]),
         fmt_est_ci(nuts["est"]["tau"], *nuts["ci_t"]["tau"]),
         fmt_est_ci(nuts["est"]["kappa"], *nuts["ci_t"]["kappa"]),
         fmt_est_ci(nuts["est"]["phi"], *nuts["ci_t"]["phi"])),
        ("NUTS-2 (2.9\\% zeros)", r"M2 (EpiEstim $\log R$)",
         "--",
         fmt_est_ci(nuts_m2["est"]["tau"], *nuts_m2_ci),
         "--", "--"),
    ]
    body = "\n".join(" & ".join(r) + r" \\" for r in rows)
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Illustrative estimator comparison in a sparse surveillance panel (German varicella, 2015--2019; $N=400$ counties, $G=16$ states, $T=261$ weeks).}}
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
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Weekly RKI SurvStat varicella notifications across 400 German counties in $G = 16$ federal states, 261 estimation weeks (ISO 2015-W01 to 2019-W52, 115{{,}}129 notifications, plus three pre-sample lag weeks); $D$ = one-week-lagged state summer holiday indicator; $w = (0.141, 0.720, 0.138)$ ($w_0 = 2.6 \times 10^{{-6}}$, log-normal working surrogate fitted to Vink et al.\ 2014 serial interval moments). Working-model finite-cluster 95\% intervals for M0 $\tau+\kappa$: WCR-score Webb6 $[{wcr_w['ci'][0]:.3f}, {wcr_w['ci'][1]:.3f}]$; WCR-score Rademacher $[{wcr_r['ci'][0]:.3f}, {wcr_r['ci'][1]:.3f}]$; unrestricted score wild (Mammen) $[{wild_m[0]:.3f}, {wild_m[1]:.3f}]$; pairs cluster bootstrap-$t$ ({pairs['ok']}/{pairs['B']} converged) $[{pairs['ci_tpk'][0]:.3f}, {pairs['ci_tpk'][1]:.3f}]$; delete-one-state jackknife $[{jk['ci_tpk'][0]:.3f}, {jk['ci_tpk'][1]:.3f}]$; Stock--Wright Bonferroni projection $[{sw05[0]:.3f}, {sw05[1]:.3f}]$. These working-model coefficients illustrate numerical computation, M2 aggregation shifts (${comps['M2']['est']['tau']:+.3f}$ to ${nuts_m2['est']['tau']:+.3f}$, driven in maintained-model simulations by EpiEstim Gamma(1,5) prior pseudo-counts under small expected counts; at county level moderate $a=0.25$ gives ${comps['M2_prior_sensitivity']['moderate_a0.25_b5']['tau']:+.3f}$ and weak $a=0.1$ gives ${comps['M2_prior_sensitivity']['weak_a0.1_b5']['tau']:+.3f}$), and finite-cluster calibration rather than causal school-holiday effects; full calendar, placebo, over-identification, and seasonal diagnostics appear in Web Appendices~D--F.}}
\end{{table}}
"""
    (TAB_DIR / "table_7_application.tex").write_text(tex, encoding="utf-8")


def write_table_7b(app: dict):
    m0 = app["primary"]["M0"]
    s = app["sensitivities"]
    by_yr = app["by_year"]
    reg = app["regional_fe"]
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
        (r"Summer $\ge 3$ days (lag 1)", s["summer_min3d_lag1"]),
        (r"Summer $\ge 3$ days (lag 2)", s["summer_min3d_lag2"]),
        (r"Summer $\ge 5$ days (lag 1)", s["summer_min5d_lag1"]),
        (r"Summer $\ge 5$ days (lag 2)", s["summer_min5d_lag2"]),
        ("Continuous holiday share (lag 1)", s["summer_Dshare_lag1"]),
        ("Continuous holiday share (lag 2)", s["summer_Dshare_lag2"]),
        ("Single year: 2015 ($T=53$)", by_yr["2015"]),
        ("Single year: 2016 ($T=52$)", by_yr["2016"]),
        ("Single year: 2017 ($T=52$)", by_yr["2017"]),
        ("Single year: 2018 ($T=52$)", by_yr["2018"]),
        ("Single year: 2019 ($T=52$)", by_yr["2019"]),
        ("Pooled 2016--2019 (drop 2015, lag 1)", reg["y2016_2019_lag1"]),
        ("All school holidays 2015--2019 (lag 0)", ah["all_2015_2019_lag0"]),
        ("All school holidays 2015--2019 (lag 1)", ah["all_2015_2019_lag1"]),
        ("All school holidays 2015--2019 (lag 2)", ah["all_2015_2019_lag2"]),
        ("All school holidays 2017--2018 (lag 0)", ah["all_2017_2018_lag0"]),
        (r"NUTS-2 regions ($N=38$, 2.9\% zeros)", nuts),
    ]
    body = "\n".join(_row_m0(lbl, r) for lbl, r in specs)
    j_n2 = nuts["J_LM_oi"]
    jg_n2 = nuts["J_LM_oi_components"]["lag_growth"]
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Temporal, threshold, background, and subsample sensitivity specifications for the German varicella panel (M0, state-clustered $t_{{G-1}}$ 95\% CIs and score-form $J_{{LM}}$ tests with $F_{{2,G-2}}$ $p$-values).}}
\label{{tab:app-all}}
\footnotesize
\setlength{{\tabcolsep}}{{3.0pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\renewcommand{{\arraystretch}}{{0.86}}\begin{{tabular}}{{lrrrrrr}}
\toprule
specification & $\hat\tau+\hat\kappa$ (95\% CI) & $\hat\tau$ (95\% CI) & $\hat\kappa$ (95\% CI) & $\hat\phi$ (95\% CI) & $J_{{LM,\mathrm{{oi}}}}$ ($p_F$) & $J_{{LM,\mathrm{{inv}}}}$ ($p_F$) \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Primary specification is Table~\ref{{tab:app}}. $J_{{LM,\mathrm{{oi}}}}$ tests lagged imported share and lagged one-week log growth; $J_{{LM,\mathrm{{inv}}}}$ tests reference-value invariance ($\phi_r \in \{{0.05, 2.0\}}$). Note substantial sensitivity across years (three of five single years cover zero), lag 2 (covers zero), and all school holidays at lag 1 (covers zero). At NUTS-2 aggregation, overall $J_{{LM,\mathrm{{oi}}}}$ is ${j_n2['stat']:.2f}$ ($p_F={j_n2['p_F']:.3f}$) while its 1-df lagged-growth component remains elevated (${jg_n2['stat']:.2f}$, $p_F={jg_n2['p_F']:.3f}$). Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
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
        ("Drop non-rotating BW \\& BY ($N=260, G=14$, lag 1)", reg["drop_bw_by_lag1"]),
        ("Drop non-rotating BW \\& BY ($N=260, G=14$, lag 2)", reg["drop_bw_by_lag2"]),
        (r"4-region $\times$ week + Drop BW \& BY (lag 1)", reg["reg4_drop_bw_by_lag1"]),
        ("Pooled 2016--2019 + Drop BW \\& BY (lag 1)", reg["y2016_2019_drop_bw_by_lag1"]),
        ("State linear trends (16 states, lag 1)", reg["state_trends_lag1"]),
        ("State linear trends (16 states, lag 2)", reg["state_trends_lag2"]),
        ("State 1-harmonic seasonal curves (lag 1)", reg["state_harmonics1_lag1"]),
        ("State 1-harmonic seasonal curves (lag 2)", reg["state_harmonics1_lag2"]),
        ("State 1-harm.\\ + trends + Drop BW \\& BY (lag 1)", reg["state_harm1_trends_drop_bw_by_lag1"]),
    ]
    body = "\n".join(_row_m0(lbl, r) for lbl, r in specs)
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Regional week-effect controls, exclusion of non-rotating southern states (Baden-W\"urttemberg and Bavaria), state trends, and state seasonal harmonics (M0, state-clustered $t_{{G-1}}$ 95\% CIs).}}
\label{{tab:rev-regional}}
\footnotesize
\setlength{{\tabcolsep}}{{3.0pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{lrrrrrr}}
\toprule
specification & $\hat\tau+\hat\kappa$ (95\% CI) & $\hat\tau$ (95\% CI) & $\hat\kappa$ (95\% CI) & $\hat\phi$ (95\% CI) & $J_{{LM,\mathrm{{oi}}}}$ ($p_F$) & $J_{{LM,\mathrm{{inv}}}}$ ($p_F$) \\
\midrule
{body}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Four macro-regions: North (SH, HH, NI, HB, MV), West (NW, HE, RP, SL), South (BW, BY) and East (BE, BB, SN, ST, TH). Baden-W\"urttemberg (08, 44 counties) and Bavaria (09, 96 counties) start summer holidays last (ISO weeks 30--31) in all five years without rotation. Dropping BW and BY leaves 14 rotating states ($N=260$ counties): the lag-1 interval widens to cover zero, and the lag-2 estimate collapses near zero. Four-region $\times$ week fixed effects likewise include zero at both lags. Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_11_regional.tex").write_text(tex, encoding="utf-8")


def write_table_12(cal: dict):
    lines_a = []
    for tag, lbl in [("calibrated", "cal"), ("high_noise", "high-noise")]:
        block = cal[tag]["params"]
        for i, (pk, plbl) in enumerate([("tpk", r"$\tau+\kappa$"), ("tau", r"$\tau$"), ("kappa", r"$\kappa$"), ("phi", r"$\phi$")]):
            p = block[pk]
            w6_c = f"{p['cov_wcr_webb6']:.3f}" if "cov_wcr_webb6" in p else "--"
            wr_c = f"{p['cov_wcr_rademacher']:.3f}" if "cov_wcr_rademacher" in p else "--"
            jk_c = f"{p['cov_jackknife']:.3f}" if "cov_jackknife" in p else "--"
            pr_c = f"{p['cov_pairs']:.3f}" if ("cov_pairs" in p and math.isfinite(p["cov_pairs"])) else "--"
            lines_a.append(" & ".join([
                lbl if i == 0 else "",
                plbl,
                f"{p['bias']:.3f}",
                f"{p['sd']:.3f}",
                f"{p['mean_se']:.3f}",
                f"{p['cov_t15']:.3f}",
                f"{p['cov_wild_mam']:.3f}",
                f"{p['cov_wild_web']:.3f}",
                w6_c,
                wr_c,
                jk_c,
                pr_c,
            ]) + r" \\")
    c_cal = cal["calibration"]
    r_tot = cal["R"]
    r_pairs = cal.get("R_pairs", 200)
    b_pairs = cal.get("B_pairs", 199)
    rc_cal = cal["calibrated"]["root_counts"]
    rc_hn = cal["high_noise"]["root_counts"]
    tpk_cal = cal["calibrated"]["params"]["tpk"]
    tpk_hn = cal["high_noise"]["params"]["tpk"]
    lines_b_fixed = []
    for tag, lbl in [("calibrated", "cal"), ("high_noise", "high-noise")]:
        block = cal[tag]["diagnostics"]
        for i, (dk, dlbl) in enumerate([
            ("j_oi", r"$J_{LM,\mathrm{oi}}$"),
            ("j_inv", r"$J_{LM,\mathrm{inv}}$"),
            ("s3_oi", r"$S_{3,\mathrm{oi}}$"),
            ("s3_inv", r"$S_{3,\mathrm{inv}}$"),
        ]):
            d = block[dk]
            lines_b_fixed.append(" & ".join([
                lbl if i == 0 else "",
                dlbl,
                f"{d['mean']:.2f}",
                f"{d['q95']:.2f}",
                f"{d['q99']:.2f}",
                f"{d['max']:.2f}",
                f"{d['obs']:.2f}",
                rf"\multicolumn{{5}}{{r}}{{{d['exceed_count']}/{d['n']} ({d['p_calibrated']:.3f})}}",
            ]) + r" \\")

    lines_c = []
    for tag, lbl in [("calibrated", "cal"), ("high_noise", "high-noise")]:
        cc = cal[tag]["comparators"]
        c_zs = f"{100 * cc['county_zero_share']:.1f}\\%"
        n2_zs = f"{100 * cc['nuts2_zero_share']:.1f}\\%"
        for i, (m_lbl, lvl, rec, extra) in enumerate([
            ("M0", f"county ({c_zs} zeros)", cc["M0_county"], f"$n={cc['M0_county']['n']}/{r_tot}$"),
            ("M2 (EpiEstim, $a=1,b=5$)", f"county ({c_zs} zeros)", cc["M2_county"], f"pos: {100*cc['M2_county_pos_share']:.1f}\\%"),
            ("M2 (mod.\\ prior $a=0.25$)", f"county ({c_zs} zeros)", cc["M2_county_mod"], "intermediate prior bias"),
            ("M2 (weak prior $a=0.1$)", f"county ({c_zs} zeros)", cc["M2_county_weak"], "halves bias"),
            ("M2 (ultra-low $a=0.01$)", f"county ({c_zs} zeros)", cc["M2_county_ultralow"], "smallest prior bias"),
            ("M2 (noise-free $Y=\\mu$, $a=1$)", f"county (deterministic)", cc["M2_county_noisefree"], f"noise/0 share: {cc['M2_county_noise_zero_component']:+.3f}"),
            (r"M3 ($\log(Y+1)$ TWFE)", f"county ({c_zs} zeros)", cc["M3_county"], "outcome slope"),
            (r"M0 (global $\phi\ge 0$, primary)", f"NUTS-2 ({n2_zs} zeros)", cc["M0_nuts2_global_tpk"], f"$n={cc['M0_nuts2_global_tpk']['n']}$, bnd: {100*cc['M0_nuts2_global_boundary_share']:.1f}\\%, $\\hat\\phi$: {cc['M0_nuts2_global_phi']['mean']:.3f}"),
            (r"M0 (unconstrained, $\hat\phi\ge 0$)", f"NUTS-2 ({n2_zs} zeros)", cc["M0_nuts2_nonneg_tpk"], f"$n={cc['M0_nuts2_nonneg_tpk']['n']}$, $\\hat\\phi$: {cc['M0_nuts2_nonneg_phi']['mean']:.3f}"),
            ("M0 (unconstrained, all conv.)", f"NUTS-2 ({n2_zs} zeros)", cc["M0_nuts2_tpk"], f"$n={cc['M0_nuts2_tpk']['n']}$ ({cc['M0_nuts2_neg_phi_count']} $\\hat\\phi<0$), $\\hat\\phi$: {cc['M0_nuts2_phi']['mean']:.3f}"),
            ("M2 (EpiEstim, $a=1,b=5$)", f"NUTS-2 ({n2_zs} zeros)", cc["M2_nuts2"], f"neg: {100*cc['M2_nuts2_neg_share']:.1f}\\%"),
            ("M2 (weak prior $a=0.1$)", f"NUTS-2 ({n2_zs} zeros)", cc["M2_nuts2_weak"], "within 2 MCSE of 0"),
        ]):
            lines_c.append(
                f"{lbl if i == 0 else ''} & {m_lbl} & \\multicolumn{{2}}{{l}}{{{lvl}}} & "
                f"{rec['mean']:.3f} & {rec['bias']:+.3f} ({rec['bias_mcse']:.3f}) & {rec['rmse']:.3f} & "
                f"\\multicolumn{{5}}{{r}}{{{extra}}} \\\\"
            )
    body_a = "\n".join(lines_a)
    body_b = "\n".join(lines_b_fixed)
    body_c = "\n".join(lines_c)
    cc_cal = cal["calibrated"]["comparators"]
    cc_hn = cal["high_noise"]["comparators"]
    mcse_tot = math.sqrt(0.95 * 0.05 / r_tot)
    mcse_pairs = math.sqrt(0.95 * 0.05 / r_pairs)
    fr_cal = cc_cal["M0_nuts2_fail_reasons"]
    fr_hn = cc_hn["M0_nuts2_fail_reasons"]
    f_tot_cal = fr_cal["lam_nonpos"] + fr_cal["tol_neg_phi"] + fr_cal["tol_other"]
    f_tot_hn = fr_hn["lam_nonpos"] + fr_hn["tol_neg_phi"] + fr_hn["tol_other"]
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Application-calibrated simulation on the German varicella design ($N=400$, $G=16$ states, $T=261$, $R={r_tot}$; nested pairs $R={r_pairs}, B={b_pairs}$): finite-cluster interval coverage, null reference distributions of $J_{{LM}}$ and $S_3$, and comparator distortion under known truth ($\tau+\kappa={tpk_cal['true']:.3f}$).}}
\label{{tab:rev-calsim}}
\footnotesize
\setlength{{\tabcolsep}}{{2.4pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{llrrrrrrrrrr}}
\toprule
design & parameter & bias & emp.\ SD & mean SE & $t_{{15}}$ & wild (M) & wild (W6) & WCR (W6) & WCR (Rad) & JK ($t_{{15}}$) & Pairs ($t$) \\
\midrule
\multicolumn{{12}}{{l}}{{\emph{{(a) Interval performance, M0 ($G = 16$ unequal state clusters; MCSE at $0.95$: ${mcse_tot:.3f}$ at $R={r_tot}$, ${mcse_pairs:.3f}$ for Pairs at $R={r_pairs}$)}}}} \\
{body_a}
\midrule
\multicolumn{{12}}{{l}}{{\emph{{(b) Finite-sample null distributions of over-identification and root-selection statistics}}}} \\
design & statistic & mean & q95 & q99 & max & application & \multicolumn{{5}}{{r}}{{exceed ($p_{{\mathrm{{cal}}}}$)}} \\
\midrule
{body_b}
\midrule
\multicolumn{{12}}{{l}}{{\emph{{(c) Comparator and spatial aggregation performance against known generating truth ($\tau+\kappa={tpk_cal['true']:.3f}$, $\phi={c_cal['beta_true'][2]:.3f}$)}}}} \\
design & method & \multicolumn{{2}}{{l}}{{level (zero share)}} & mean & bias (MCSE) & RMSE & \multicolumn{{5}}{{r}}{{convergence / mechanism note}} \\
\midrule
{body_c}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Calibrated DGP retains the observed BKG VG250 400-county network, 261-week summer holiday schedule, Vink et al.\ (2014) working weights, truth set to primary M0 ($\tau={c_cal['beta_true'][0]:.3f}$, $\kappa={c_cal['beta_true'][1]:.3f}$, $\phi={c_cal['beta_true'][2]:.3f}$), NB2 $k={c_cal['k_nb']:.2f}$ and state $\times$ week log-normal shock SD $\sigma={c_cal['sig_sw']:.3f}$ (\emph{{cal}}), or $k/2$ and $2\sigma$ (\emph{{high-noise}}). All rates assume the maintained working model (Assumption~A7 and parallel trends hold). County M0 achieves {rc_cal['one']}/{r_tot} ({100*rc_cal['share_unique']:.1f}\%, \emph{{cal}}) and {rc_hn['one']}/{r_tot} ({100*rc_hn['share_unique']:.1f}\%, \emph{{high-noise}}) unique interior convergence; sandwich SEs understate $\tau+\kappa$ SD by {100*(1-tpk_cal['se_ratio']):.1f}\% (SE ratio {tpk_cal['se_ratio']:.3f}; {tpk_hn['se_ratio']:.3f} high-noise). Unified delete-one-state jackknife (JK; {tpk_cal['jk_ok']}/{r_tot} converged) is slightly conservative (SE ratio {tpk_cal['jk_se_ratio']:.3f} \emph{{cal}}, {tpk_hn['jk_se_ratio']:.3f} \emph{{high-noise}}), WCR-score ($B=399$) covers at {min(tpk_cal['cov_wcr_webb6'], tpk_cal['cov_wcr_rademacher']):.3f}--{max(tpk_hn['cov_wcr_webb6'], tpk_hn['cov_wcr_rademacher']):.3f}, and nested pairs cluster bootstrap-$t$ ($R={r_pairs}, B={b_pairs}$; $\approx {tpk_cal['pairs_mean_ok']:.0f}/{b_pairs}$ inner convergence) improves on $t_{{15}}$ to {tpk_cal['cov_pairs']:.3f} (\emph{{cal}}) and {tpk_hn['cov_pairs']:.3f} (\emph{{high-noise}}). In Panel~(c), simulated county M2 positive bias under default EpiEstim prior (${cc_cal['M2_county']['bias']:+.3f}$ \emph{{cal}}, ${cc_hn['M2_county']['bias']:+.3f}$ \emph{{high-noise}}) is dominated by prior pseudo-counts under small expected counts (noise-free $Y=\mu$ gives ${cc_cal['M2_county_noisefree']['bias']:+.3f}$; noise/zero realizations add ${cc_cal['M2_county_noise_zero_component']:+.3f}$; moderate prior $a=0.25$ gives ${cc_cal['M2_county_mod']['bias']:+.3f}$, weak $a=0.1$ gives ${cc_cal['M2_county_weak']['bias']:+.3f}$, ultra-low $a=0.01$ gives ${cc_cal['M2_county_ultralow']['bias']:+.3f}$), while M3 estimates the $\log(Y+1)$ TWFE slope rather than $\tau+\kappa$. At NUTS-2, global grid search ($\phi\ge 0$) resolves {cc_cal['M0_nuts2_global_tpk']['n']}/{r_tot} draws with {100*cc_cal['M0_nuts2_global_boundary_share']:.1f}\% (\emph{{cal}}) and {100*cc_hn['M0_nuts2_global_boundary_share']:.1f}\% (\emph{{high-noise}}) selecting $\phi=0$; unconstrained fits do not restrict $\phi\ge 0$, converging on {cc_cal['M0_nuts2_tpk']['n']}/{r_tot} (\emph{{cal}}, of which {cc_cal['M0_nuts2_neg_phi_count']} have $\hat\phi<0$) and {cc_hn['M0_nuts2_tpk']['n']}/{r_tot} (\emph{{high-noise}}, {cc_hn['M0_nuts2_neg_phi_count']} with $\hat\phi<0$), while failed unconstrained draws stalled at finite residuals on the $\phi<0$ side ({fr_cal['tol_neg_phi']}/{f_tot_cal} \emph{{cal}}, {fr_hn['tol_neg_phi']}/{f_tot_hn} \emph{{high-noise}}; {fr_cal['lam_nonpos']}/{f_tot_cal} entered $\Lambda(\phi)\le 0$). Source: \texttt{{experiment/analysis/germany\_calibrated\_sim.json}}.}}
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
        w_str = f"({g['w'][0]:.4f}, {g['w'][1]:.4f}, {g['w'][2]:.4f})"
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
\caption{{Working generation-interval surrogate (Vink et al., 2014 serial interval) and reporting-dispersion sensitivity for the German varicella panel (M0, summer holidays lagged one week).}}
\label{{tab:rev-gi}}
\footnotesize
\setlength{{\tabcolsep}}{{2.8pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{lccrrrr}}
\toprule
interval specification & $(w_1, w_2, w_3)$ & $w_0$ & $\hat\tau+\hat\kappa$ (95\% CI) & $\hat\phi$ (95\% CI) & $J_{{LM,\mathrm{{oi}}}}$ ($p_F$) & $J_{{LM,\mathrm{{inv}}}}$ ($p_F$) \\
\midrule
{chr(10).join(lines)}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Following \citet{{vink2014serial}} (pooled varicella serial interval: mean 14.0 d, SD $\approx 2.5$ d), weights discretise our assumed lognormal density onto weekly lags $s=1,2,3$ under a uniform within-week infection day; $w_0$ is the resulting same-week share (three-decimal main-text rounding $0.141, 0.720, 0.138$ omits 0.001 in rounding). Passing $J_{{LM,\mathrm{{oi}}}}$ only at SD 6.0 d indicates sensitivity to lag dispersion and reporting delays. Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
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
        tL = evL["total_scale"][nm]
        s_dir = fmt_est_ci(cL["est"], *cL["ci"])
        s_tot = f"{fmt_est_ci(tL['est'], *tL['ci'])} [{tL['p_t']:.3f}]"
        s_nolead = fmt_est_ci(evN["total_scale"][nm]["est"], *evN["total_scale"][nm]["ci"]) if nm in evN["total_scale"] else "--"
        role = "Pre-trend lead" if k < 0 else ("Onset week" if k == 0 else ("Holiday weeks" if k <= 5 else "Reopening"))
        lines.append(f"$e = {k:+d}$ & {role} & {s_dir} & {s_tot} & {s_nolead} \\\\")
    for nm, lbl in [("kappa", r"Neighbour control ($\hat\kappa$)"), ("phi", r"Contagion coupling ($\hat\phi$)")]:
        cL = evL["coefs"][nm]
        cN = evN["coefs"][nm]
        lines.append(f"{lbl} & Control / structural & {fmt_est_ci(cL['est'], *cL['ci'])} & -- & {fmt_est_ci(cN['est'], *cN['ci'])} \\\\")
    pw = evL["pretrend_wald"]
    ptot = evL["post_total_wald_e1_e5"]
    rtot = evL["reopen_total_wald_e6_e7"]
    e0t = evL["total_scale"]["e0"]
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Event-time dynamic profile around German summer holiday onset ($e = -2, \dots, 7$): direct bin coefficients $\hat\tau_e$ versus total-policy scale effects $\hat\tau_e + c_e\hat\kappa$ ($t_{{15}}$ 95\% CIs).}}
\label{{tab:rev-event}}
\footnotesize
\setlength{{\tabcolsep}}{{3.5pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{llrrr}}
\toprule
event week & timing role & direct $\hat\tau_e$ (with leads) & total $\hat\tau_e + c_e\hat\kappa$ (with leads) [$p_t$] & total (no leads) \\
\midrule
{chr(10).join(lines)}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Because $\bar D_{{it}}\approx D_{{it}}$ inside states ($c_e \in [0.94, 0.97]$ for $e=0,\dots,5$), $\hat\tau_e$ isolates only within-border contrasts while total policy scale $\hat\tau_e + c_e\hat\kappa$ matches $\tau+\kappa$. Pre-trend joint test ($e_{{-2}}=e_{{-1}}=0$): $\chi^2_2 = {pw['stat']:.2f}$ ($p={pw['p_chi2']:.3f}$, $p_F={pw['p_F']:.3f}$). On the total-policy scale, onset week $e=0$ already drops to ${e0t['est']:.3f}$ (SE {e0t['se']:.3f}, $p_t={e0t['p_t']:.3f}$), followed by $e=1,\dots,5$ ($F_{{5,11}}$ $p_F={ptot['p_F']:.3f}$, $\chi^2_5={ptot['stat']:.2f}$ $p={ptot['p_chi2']:.3f}$) and non-positive reopening weeks $e=6,7$ ($p_F={rtot['p_F']:.3f}$). Immediate onset decline before the 14-day incubation window prevents ruling out holiday-linked reporting drops. Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_14_event_weights.tex").write_text(tex, encoding="utf-8")


def write_table_15(app: dict):
    asc = app["ascertainment"]
    keep_rhos = {-0.20, -0.10, 0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40}
    rows_g = [r for r in asc["grid"] if r["ok"] and round(r["rho"], 2) in keep_rhos]
    lines_a = []
    for r in rows_g:
        lines_a.append(" & ".join([
            f"{r['rho']:+.2f}",
            f"{r['bound_log_shift']:.3f}",
            fmt_est_ci(r["est"]["tpk"], r["lo"], r["hi"]),
            f"{r['est']['phi']:.3f} ({r['se']['phi']:.3f})",
        ]) + r" \\")
    r_hi = asc["rho_star_upper_ci_zero"]
    r_0 = asc["rho_star_tpk_zero"]
    if r_hi is not None:
        f_hi = asc["rho_star_upper_ci_zero_fit"]
        lines_a.append(rf"$\rho^*_{{\mathrm{{CI}}}} = {r_hi:+.3f}$ & {f_hi['bound_log_shift']:.3f} & {fmt_est_ci(f_hi['est']['tpk'], f_hi['lo'], f_hi['hi'])} & {f_hi['est']['phi']:.3f} ({f_hi['se']['phi']:.3f}) \\")
    if r_0 is not None:
        f_0 = asc["rho_star_tpk_zero_fit"]
        lines_a.append(rf"$\rho^*_0 = {r_0:+.3f}$ & {f_0['bound_log_shift']:.3f} & {fmt_est_ci(f_0['est']['tpk'], f_0['lo'], f_0['hi'])} & {f_0['est']['phi']:.3f} ({f_0['se']['phi']:.3f}) \\")

    # Panel (b): Contemporaneous-D0 correction & onset-moment zeroing
    lines_b = []
    for r1, r2 in zip(asc["contemp_grid_lag1"], asc["contemp_grid_lag2"]):
        lines_b.append(" & ".join([
            f"{r1['rho']:.2f}",
            f"{r1['m_onset']:+.4f}",
            fmt_est_ci(r1["est"]["tpk"], r1["lo"], r1["hi"]),
            fmt_est_ci(r2["est"]["tpk"], r2["lo"], r2["hi"]),
        ]) + r" \\")
    f_on1 = asc["fit_at_onset_zero_lag1"]
    f_on2 = asc["fit_at_onset_zero_lag2"]
    if f_on1 and f_on2:
        lines_b.append(" & ".join([
            rf"$\rho^*_{{\mathrm{{on}}}} = {asc['rho_star_onset_zero_lag1']:.3f}$ / ${asc['rho_star_onset_zero_lag2']:.3f}$",
            "0.0000",
            fmt_est_ci(f_on1["est"]["tpk"], f_on1["lo"], f_on1["hi"]),
            fmt_est_ci(f_on2["est"]["tpk"], f_on2["lo"], f_on2["hi"]),
        ]) + r" \\")

    # Panel (c): psi sensitivity to omega
    psi0 = app["primary"]["psi"]
    om_psi = app["primary"]["omega_psi_sensitivity"]
    lines_c = [
        f"$\\omega_i \\times 1.00$ (primary) & {app['primary']['M0']['est']['phi']:.3f} & {psi0['treated']['est']:,.0f} ({psi0['treated']['se']:,.0f}) & {psi0['network']['est']:,.0f} ({psi0['network']['se']:,.0f}) \\\\",
        f"$\\omega_i \\times 0.25$ & {app['sensitivities']['omega_x_0.25']['est']['phi']:.3f} & {om_psi['omega_x_0.25']['treated']['est']:,.0f} ({om_psi['omega_x_0.25']['treated']['se']:,.0f}) & {om_psi['omega_x_0.25']['network']['est']:,.0f} ({om_psi['omega_x_0.25']['network']['se']:,.0f}) \\\\",
        f"$\\omega_i \\times 4.00$ & {app['sensitivities']['omega_x_4']['est']['phi']:.3f} & {om_psi['omega_x_4']['treated']['est']:,.0f} ({om_psi['omega_x_4']['treated']['se']:,.0f}) & {om_psi['omega_x_4']['network']['est']:,.0f} ({om_psi['omega_x_4']['network']['se']:,.0f}) \\\\",
    ]

    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Ascertainment sensitivity (lagged-$D$ tipping points and contemporaneous-$D^{{(0)}}$ onset-moment matching) and sensitivity of $\hat\psi^\dagger$ to background importation $\omega_i$.}}
\label{{tab:rev-ascert}}
\footnotesize
\setlength{{\tabcolsep}}{{3.5pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{rrrr}}
\toprule
\multicolumn{{4}}{{l}}{{\emph{{(a) Lagged-$D$ tipping-point grid: $Y^*_{{it}} = Y_{{it}}/(1-\rho)^{{D_{{it}}}}$}}}} \\
$\rho$ & $|\log(1-\rho)|$ & $\hat\tau+\hat\kappa$ (95\% $t_{{15}}$ CI) & $\hat\phi$ (SE) \\
\midrule
{chr(10).join(lines_a)}
\midrule
\multicolumn{{4}}{{l}}{{\emph{{(b) Contemporaneous-$D^{{(0)}}$ reporting adjustment: $Y^*_{{it}} = Y_{{it}}/(1-\rho)^{{D^{{(0)}}_{{it}}}}$}}}} \\
$\rho$ & onset moment & lag-1 $\hat\tau+\hat\kappa$ (95\% CI) & lag-2 $\hat\tau+\hat\kappa$ (95\% CI) \\
\midrule
{chr(10).join(lines_b)}
\midrule
\multicolumn{{4}}{{l}}{{\emph{{(c) Sensitivity of conditional cases averted $\hat\psi^\dagger$ to background importation $\omega_i$}}}} \\
scaling & $\hat\phi$ & treated $\hat\psi^\dagger$ (SE) & network $\hat\psi^\dagger$ (SE) \\
\midrule
{chr(10).join(lines_c)}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize In Panel (b), matching the holiday-onset moment to zero requires a {100*asc['rho_star_onset_zero_lag1']:.1f}\% contemporaneous reporting reduction, which eliminates the lag-1 effect (${f_on1['est']['tpk']:+.3f}$, SE ${f_on1['se']['tpk']:.3f}$) while leaving an imprecise negative lag-2 estimate (${f_on2['est']['tpk']:+.3f}$, SE ${f_on2['se']['tpk']:.3f}$, covering zero). Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
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
            fmt_ivs(st["S3_inv_F_95"]),
            fmt_ivs(st["S3_oi_chi2_95"]),
            fmt_ivs(st["S3_oi_F_95"]),
            f"({bp[0]:.3f}, {bp[1]:.3f})",
        ]) + r" \\")
    s0 = sw["phi_bar_0.5"]["S1_at_0"]
    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Identification-robust profiled Stock--Wright confidence sets for $\phi$ (including empty $\chi^2$ and uninformative full-grid $F$ sets) and Bonferroni 95\% projection for $\tau+\kappa$.}}
\label{{tab:rev-robust}}
\footnotesize
\setlength{{\tabcolsep}}{{2.4pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{rlllllll}}
\toprule
$\bar\phi$ & $S_1(F)$ & $S_1(\chi^2)$ & $S_{{3,\mathrm{{inv}}}}(\chi^2)$ & $S_{{3,\mathrm{{inv}}}}(F)$ & $S_{{3,\mathrm{{oi}}}}(\chi^2)$ & $S_{{3,\mathrm{{oi}}}}(F)$ & $\tau+\kappa$ (Bonf.\ 95\%) \\
\midrule
{chr(10).join(lines)}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Profiled Stock--Wright statistics evaluated on 304 grid values $\phi\in[0,3]$ ($G=16$ states), presuming strong identification of $(\delta,\tau,\kappa)$ at each tested $\phi$. Notice: $S_{{3,\mathrm{{oi}}}}(\chi^2)$ is empty across all three $\bar\phi$ (reflecting the rejected $J_{{LM,\mathrm{{oi}}}}$ over-identifying test), while $S_{{3,\mathrm{{inv}}}}(F)$, $S_{{3,\mathrm{{oi}}}}(F)$, and the 97.5\% $S_1(F)$ set used in the Bonferroni union all span the entire grid $[0.00,3.00]$. At $\phi=0$, $\bar\phi=0.5$: $S_1={s0['S1']:.2f}$ ($p_F={s0['S1_p_F']:.3f}$). Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_16_robust_sets.tex").write_text(tex, encoding="utf-8")


def write_table_18_placebos(app: dict):
    pl = app["placebos"]
    oh = pl["other_holidays"]
    lines_a = []
    for hkey, hlbl in [("easter", "Easter"), ("autumn", "Autumn"), ("winter", "Winter"), ("pentecost", "Pentecost")]:
        for lg in (0, 1, 2):
            r = oh[hkey][f"lag_{lg}"]
            lines_a.append(" & ".join([
                hlbl if lg == 0 else "",
                str(lg),
                fmt_est_ci(r["est"]["tpk"], *r["ci_t"]["tpk"]),
                f"{r['p_t_tpk'] if 'p_t_tpk' in r else 2*__import__('scipy').stats.t.sf(abs(r['est']['tpk']/r['se']['tpk']), 15):.4f}",
                fmt_est_ci(r["est"]["phi"], *r["ci_t"]["phi"]),
                f"{r['J_LM_oi']['stat']:.2f} ({r['J_LM_oi']['p_F']:.3f})",
            ]) + r" \\")
    ts = pl["time_shifts"]
    lines_b = []
    for sh in (-4, -3, -2, -1, 0, 1, 2, 3, 4):
        r = ts[f"shift_{sh:+d}"]
        pval = 2 * __import__("scipy").stats.t.sf(abs(r["est"]["tpk"] / r["se"]["tpk"]), 15)
        lines_b.append(" & ".join([
            f"{sh:+d} wk (eff lag {1+sh:+d})" if sh >= 0 else f"{sh:+d} wk (lead {-sh} wk)",
            fmt_est_ci(r["est"]["tpk"], *r["ci_t"]["tpk"]),
            f"{pval:.4f}",
            fmt_est_ci(r["est"]["phi"], *r["ci_t"]["phi"]),
            f"{r['J_LM_oi']['stat']:.2f} ({r['J_LM_oi']['p_F']:.3f})",
        ]) + r" \\")
    pm = pl["state_permutations"]
    lines_b_fixed = []
    for sh in (-4, -3, -2, -1, 0, 1, 2, 3, 4):
        r = ts[f"shift_{sh:+d}"]
        pval = 2 * __import__("scipy").stats.t.sf(abs(r["est"]["tpk"] / r["se"]["tpk"]), 15)
        lbl = f"{sh:+d} wk (lag {1+sh})" if 1 + sh >= 0 else f"{sh:+d} wk (lead {-(1+sh)})"
        lines_b_fixed.append(" & ".join([
            lbl,
            rf"\multicolumn{{2}}{{r}}{{{fmt_est_ci(r['est']['tpk'], *r['ci_t']['tpk'])}}}",
            f"{pval:.4f}",
            fmt_est_ci(r["est"]["phi"], *r["ci_t"]["phi"]),
            f"{r['J_LM_oi']['stat']:.2f} ({r['J_LM_oi']['p_F']:.3f})",
        ]) + r" \\")
    lines_c_fixed = []
    for key, lbl in [("lag1", "Summer perm. (lag 1)"), ("lag2", "Summer perm. (lag 2)")]:
        p = pm[key]
        lines_c_fixed.append(f"{lbl} & {p['obs_tpk']:+.3f} & {p['mean']:+.3f} ({p['sd']:.3f}) & [{p['q05']:+.3f}, {p['q95']:+.3f}] & \\multicolumn{{2}}{{r}}{{{p['p_one_sided']:.3f} / {p['p_two_sided']:.3f}}} \\\\")
    body_a = "\n".join(lines_a)
    body_b = "\n".join(lines_b_fixed)
    body_c = "\n".join(lines_c_fixed)

    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Placebo and randomisation diagnostics: non-summer school holidays (lags 0--2), time-shifted summer calendars, and cross-state schedule permutations.}}
\label{{tab:rev-placebo}}
\footnotesize
\setlength{{\tabcolsep}}{{3.0pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{lrrrrr}}
\toprule
\multicolumn{{6}}{{l}}{{\emph{{(a) Non-summer school holiday placebos (M0, lags 0--2)}}}} \\
holiday & lag & $\hat\tau+\hat\kappa$ (95\% $t_{{15}}$ CI) & $p_t$ & $\hat\phi$ (95\% CI) & $J_{{LM,\mathrm{{oi}}}}$ ($p_F$) \\
\midrule
{body_a}
\midrule
\multicolumn{{6}}{{l}}{{\emph{{(b) Time-shifted summer schedules relative to primary lag-1 alignment ($\Delta t \in \{{-4,\dots,+4\}}$)}}}} \\
calendar shift & \multicolumn{{2}}{{r}}{{$\hat\tau+\hat\kappa$ (95\% $t_{{15}}$ CI)}} & $p_t$ & $\hat\phi$ (95\% CI) & $J_{{LM,\mathrm{{oi}}}}$ ($p_F$) \\
\midrule
{body_b}
\midrule
\multicolumn{{6}}{{l}}{{\emph{{(c) Cross-state randomisation test (permuting the 16 states' summer schedules across states, B=200)}}}} \\
specification & observed & perm.\ mean (SD) & [5th, 95th] & \multicolumn{{2}}{{r}}{{$p$ (1-sided / 2-sided)}} \\
\midrule
{body_c}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize All four short holidays (Panel a) exhibit sharp negative lag-0 dips ($-0.20$ to $-0.50$), near-zero lag-1 contrasts, and positive lag-2 rebounds ($+0.07$ to $+0.24$, matching the 14-day serial-interval echo of contemporaneous holiday reporting suppression). Panel (c) permutes whole state schedules ($B=200$, minimum attainable two-sided $p=0.010$; exchangeability is imperfect because BW and BY never rotate, and schedule permutations reflect both contact and reporting changes aligned with holiday calendars). Panel (b) reveals positive pre-window slopes at overlapping calendar shifts $-4$ to $-2$ weeks (leads 3--1) and maximum negative magnitude at $+2$ weeks (lag 3). Source: \texttt{{experiment/analysis/germany\_application\_results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_18_placebos.tex").write_text(tex, encoding="utf-8")


def write_table_envelope(cal: dict):
    res = json.loads((ROOT / "experiment" / "analysis" / "results.json").read_text(encoding="utf-8"))
    bnd = json.loads((ROOT / "verification" / "root_selection" / "full_summary_tsummary.json").read_text(encoding="utf-8"))
    sw = json.loads((ROOT / "verification" / "same_week" / "results.json").read_text(encoding="utf-8"))

    cores = [k for k in res["summary"]["M0"]["convergence"] if k.startswith("core_") and "_spill" in k]
    c_bias = [res["summary"]["M0"]["bias"][f"{k}|tau"]["mean"] for k in cores]
    c_cov = [100.0 * res["summary"]["M0"]["coverage"][f"{k}|tau"]["mean"] for k in cores]

    b0 = [x for x in bnd if x["key"][0] == 0.0][0]
    b_weak_gi = [x["cov_gi_t"] for x in bnd if x["key"][0] in (0.05, 0.1)]
    b_weak_m0 = [x["cov_m0_t"] for x in bnd if x["key"][0] in (0.05, 0.1)]

    sw_tpk = sw["arms"]["same_week"]["tau+kappa"]
    sw_truth = sw["truth"]["tau+kappa"]
    sw_frac = sw_tpk["bias"] / sw_truth

    tpk_c = cal["calibrated"]["params"]["tpk"]
    tpk_h = cal["high_noise"]["params"]["tpk"]
    wcr_lo = min(tpk_c["cov_wcr_webb6"], tpk_c["cov_wcr_rademacher"], tpk_h["cov_wcr_webb6"], tpk_h["cov_wcr_rademacher"])
    wcr_hi = max(tpk_c["cov_wcr_webb6"], tpk_c["cov_wcr_rademacher"], tpk_h["cov_wcr_webb6"], tpk_h["cov_wcr_rademacher"])
    jk_lo = min(tpk_c["cov_jackknife"], tpk_h["cov_jackknife"])
    jk_hi = max(tpk_c["cov_jackknife"], tpk_h["cov_jackknife"])
    jk_cov_str = f"{jk_lo:.3f}" if f"{jk_lo:.3f}" == f"{jk_hi:.3f}" else f"{jk_lo:.3f}--{jk_hi:.3f}"
    jk_se_lo = 100.0 * (min(tpk_c["jk_se_ratio"], tpk_h["jk_se_ratio"]) - 1.0)
    jk_se_hi = 100.0 * (max(tpk_c["jk_se_ratio"], tpk_h["jk_se_ratio"]) - 1.0)
    pr_lo = min(tpk_c["cov_pairs"], tpk_h["cov_pairs"])
    pr_hi = max(tpk_c["cov_pairs"], tpk_h["cov_pairs"])
    t15_lo = min(tpk_c["cov_t15"], tpk_h["cov_t15"])
    t15_hi = max(tpk_c["cov_t15"], tpk_h["cov_t15"])

    x2_tpk = (
        res["summary"]["M0"]["bias"]["X2_reporting_drops_after_adoption|tau"]["mean"]
        + res["summary"]["M0"]["bias"]["X2_reporting_drops_after_adoption|kappa"]["mean"]
    )
    x2_phi = res["summary"]["M0"]["bias"]["X2_reporting_drops_after_adoption|phi"]["mean"]

    r_cal = cal["R"]
    r_pr = cal.get("R_pairs", r_cal)
    mcse_cal = math.sqrt(0.95 * 0.05 / r_cal)
    mcse_pr = math.sqrt(0.95 * 0.05 / r_pr)

    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Operating envelope and diagnostic workflow for weekly QD-GMM across data-generating and observation regimes.}}
\label{{tab:envelope}}
\footnotesize
\setlength{{\tabcolsep}}{{3.4pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{llll}}
\toprule
regime / data structure & estimator & simulation performance & recommended usage \\
\midrule
Multi-week GI, feedback, spillover ($\phi_0>0$) & M0 & bias ${min(c_bias):.3f}$ to ${max(c_bias):.3f}$; cov.\ {min(c_cov):.1f}\%--{max(c_cov):.1f}\% & primary estimator (maintained conds.) \\
Zero/weak contagion ($\phi_0\approx 0$) & M0-GI (post hoc) & cov.\ {b0['cov_gi_z']:.3f}--{b0['cov_gi_t']:.3f} ($\phi_0=0$); {min(b_weak_gi):.2f}--{max(b_weak_gi):.2f} ($\phi_0\in[0.05,0.10]$; vs.\ M0 {min(b_weak_m0):.2f}--{max(b_weak_m0):.2f}) & profile $[0,3]$; post hoc M0-GI; $S_1(\phi)$ \\
Substantial same-week mass ($w_0\gtrsim 0.3$) & sub-weekly Eq. & weekly M0 bias $\approx {sw_frac:.2f}\times$ total (cov.\ {sw_tpk['coverage']:.2f}) & switch to sub-weekly model \\
Moderate/unequal clusters (e.g.\ $G=16$) & M0 + WCR / JK & cal.\ cov.\ {wcr_lo:.3f}--{wcr_hi:.3f} (WCR), {jk_cov_str} (JK, SE $+{jk_se_lo:.0f}$--${jk_se_hi:.0f}$\%); Pairs {pr_lo:.3f}--{pr_hi:.3f} & report WCR-score \& JK; Pairs supp. \\
Policy-linked ascertainment (violates A7) & working M0 + diag. & biases $\tau+\kappa$ ($\approx {x2_tpk:.2f}$) \& $\phi$ ($+{x2_phi:.2f}$) in X2 & onset/placebo diag.; working output only \\
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Summary of operating limits generated dynamically from core ($R=8000$), boundary validation ($R=2000$; Web Appendix~C), same-week ($R=1000$ per arm), and design-calibrated finite-cluster ($R={r_cal}$, MCSE at $0.95$ equals ${mcse_cal:.3f}$; nested Pairs $R={r_pr}$, MCSE at $0.95$ equals ${mcse_pr:.3f}$; vs.\ ${t15_lo:.3f}$--${t15_hi:.3f}$ for $t_{{15}}$, assuming the maintained working model under A7 and $\phi_0={cal['calibration']['beta_true'][2]:.3f}$) experiments. Note that jackknife is slightly conservative (SE ratio $+{jk_se_lo:.0f}$--${jk_se_hi:.0f}$\%) while M0-GI avoids boundary collapse at $\phi_0=0$ (where M0 covers at {b0['cov_m0_t']:.2f}) at a modest weak-contagion coverage cost relative to M0. Diagnostics separate regimes in which M0 serves as the primary estimator under the maintained renewal, timing, and separable-reporting (Assumption~A7) conditions from regimes requiring post hoc boundary heuristics, sub-weekly models, finite-cluster calibration, or exploratory working-model handling.}}
\end{{table}}
"""
    (TAB_DIR / "table_envelope.tex").write_text(tex, encoding="utf-8")


def main():
    app = json.loads(APP_PATH.read_text(encoding="utf-8"))
    cal = json.loads(CAL_PATH.read_text(encoding="utf-8"))
    write_table_envelope(cal)
    write_table_7(app)
    write_table_7b(app)
    write_table_11(app)
    write_table_12(cal)
    write_table_13(app)
    write_table_14(app)
    write_table_15(app)
    write_table_16(app)
    write_table_18_placebos(app)
    print("Generated all 10 German Varicella & envelope LaTeX tables in tables/.")


if __name__ == "__main__":
    main()
