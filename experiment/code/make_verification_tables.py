"""Generate supplementary verification LaTeX tables (Tables 19 and 20) for R9:
  - tables/table_19_nonnested_weakprior.tex : Non-nested daily renewal DGP (R=500),
      same-week w0 grid with matched treated share (R=1000), and weak-prior M2 across
      all 25 core/null/misspecification scenarios (R=8000).
  - tables/table_20_cluster_grid.tex        : Pure-G balanced grid (N=2400, T=12) and
      cluster skewness x T decomposition (N=400, G=16) comparing t_{G-1}, t_{G-2},
      delete-one-cluster jackknife, WCR-Webb6, and WCR-Rademacher (R=250).
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TAB_DIR = ROOT / "tables"


def write_table_19():
    nn = json.loads((ROOT / "verification" / "non_nested_sim" / "results.json").read_text(encoding="utf-8"))
    wk = json.loads((ROOT / "verification" / "weak_prior_m2" / "results.json").read_text(encoding="utf-8"))
    sw = json.loads((ROOT / "verification" / "same_week" / "grid_results.json").read_text(encoding="utf-8"))

    # Panel (a): Non-nested daily Gamma-GI + delay + SEIR depletion (R=500)
    nn_labels = [
        ("NN1_daily_gamma_gi", r"NN1: daily $\Gamma$(12d, 3.5d) + wrong $\hat w,\hat\omega$"),
        ("NN2_delay_2d", r"NN2: + 2-day stochastic reporting delay"),
        ("NN3_delay_4d", r"NN3: + 4-day stochastic reporting delay"),
        ("NN4_seir_depletion", r"NN4: + 2d delay + SEIR depletion + 5\% jitter"),
    ]
    lines_a = []
    for k, lbl in nn_labels:
        sc = nn["scenarios"][k]
        def cell(m):
            r = sc[m]
            return f"{r['bias']:+.3f} ({100*r['coverage']:4.1f}\\%)" if m != "M3" else f"{r['bias']:+.3f} (slope)"
        lines_a.append(" & ".join([
            lbl,
            cell("M0"),
            cell("M0_delay_adj"),
            cell("M1"),
            cell("M2_weak"),
            cell("M2_def"),
            cell("M3"),
            cell("M4"),
        ]) + r" \\")

    # Panel (b): Same-week w0 grid (R=1000, free vs matched treated share)
    lines_b = []
    for k, v in sw["grid"].items():
        f_ = v["free"]
        m_ = v["matched_D"]
        lines_b.append(" & ".join([
            f"{v['w0']:.2f}",
            f"{100*v['mean_reinf_share']:.1f}\\% ({100*v['mean_seeded_inc_share']:.1f}\\%)",
            f"{f_['treated_share']:.3f}",
            f"{f_['tau_bias']:+.3f} ({100*f_['tau_cov']:.1f}\\%)",
            f"{f_['tpk_bias']:+.3f} ({100*f_['tpk_cov']:.1f}\\%)",
            f"{m_['tpk_bias']:+.3f} ({100*m_['tpk_cov']:.1f}\\%)",
            f"{100*f_['j_rej']:.1f}\\%",
            "--",
        ]) + r" \\")

    # Panel (c): Weak-prior M2 across key scenarios (R=8000, identical seeds to Table 1)
    wk_pick = [
        ("1", r"1: $N=500, T=12$, FB, spill"),
        ("2", r"2: $N=500, T=26$, FB, spill"),
        ("4", r"4: $N=2000, T=26$, FB, spill"),
        ("9", r"9: $N=500, T=12$, RD, spill"),
        ("12", r"12: $N=2000, T=26$, RD, spill"),
        ("18", r"18: $N=500, T=12$, null, FB, $k=\infty$"),
        ("21", r"21: X1 wrong generation interval"),
        ("22", r"22: X2 reporting drop after adoption"),
        ("24", r"24: X4 susceptible depletion"),
    ]
    lines_c = []
    for sid, lbl in wk_pick:
        r = wk["summary"][sid]
        lines_c.append(" & ".join([
            lbl,
            f"{r['M2_def']['bias']:+.3f} ({100*r['M2_def']['coverage']:4.1f}\\%)",
            f"{r['M2_weak']['bias']:+.3f} ({100*r['M2_weak']['coverage']:4.1f}\\%)",
            f"{r['M2_ultra']['bias']:+.3f} ({100*r['M2_ultra']['coverage']:4.1f}\\%)",
            f"{r['M2_weak']['empSE']:.3f} / {r['M2_weak']['modelSE']:.3f}",
            f"{r['M2_weak']['rmse']:.3f}",
            f"{r['M2_ultra']['rmse']:.3f}",
            "prior offset" if sid in ("2", "4", "21") else "feedback/noise",
        ]) + r" \\")

    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Non-nested daily renewal simulations ($R=500$), within-week generation grid ($R=1000$), and weak-prior EpiEstim M2 sensitivity across main-protocol scenarios ($R=8000$).}}
\label{{tab:rev-nonnested}}
\footnotesize
\setlength{{\tabcolsep}}{{2.6pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{lrrrrrrr}}
\toprule
\multicolumn{{8}}{{l}}{{\emph{{(a) Non-nested daily Gamma-GI renewal + stochastic reporting delay + SEIR depletion ($N=500, G=50, T=16$, $R={nn['R']}$)}}}} \\
scenario & M0 ($\tau$ bias, cov) & M0 delay-adj & M1 & M2 weak $(0.1,5)$ & M2 def $(1,5)$ & M3 ($\log(Y+1)$) & M4 \\
\midrule
{chr(10).join(lines_a)}
\midrule
\multicolumn{{8}}{{l}}{{\emph{{(b) Same-week generation grid $w_0 \in [0, 0.30]$: free adoption vs.\ treated-share-matched control ($R={sw['R']}$)}}}} \\
$w_0$ & reinf.\ (seeded) & $\bar D$ (free) & M0 $\tau$ (free) & M0 $\tau+\kappa$ (free) & M0 $\tau+\kappa$ (matched $\bar D$) & $J$ rej. & -- \\
\midrule
{chr(10).join(lines_b)}
\midrule
\multicolumn{{8}}{{l}}{{\emph{{(c) Weakening the EpiEstim prior on the $R=8000$ main-protocol seeds (M2 default $(1,5)$ vs.\ weak $(0.1,5)$ and ultra-low $(0.01,100)$)}}}} \\
main scenario & M2 $(1,5)$ & M2 $(0.1,5)$ & M2 $(0.01,100)$ & emp / mod SE $(0.1,5)$ & RMSE $(0.1,5)$ & RMSE $(0.01,100)$ & note \\
\midrule
{chr(10).join(lines_c)}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Panel~(a): data simulated on a daily branching renewal process ($\Gamma$ mean 12d, SD 3.5d) with misspecified working analyst weights $\hat w$ ($\Gamma$ mean 13.5d, SD 4.0d) and $\hat\omega$; M0 delay-adj shifts $D$ by fractional mean delay $\hat\delta/7$. Notice unadjusted M0 bias rises cleanly with delay (highlighting the need for delay alignment), whereas fractional delay alignment restores 85.4\%--99.7\% coverage; M3 targets a distinct outcome-scale slope. Panel~(b): computes runtime same-week reinfection share (\emph{{reinf.}}) and total seeded+reinfection share (\emph{{seeded}}) alongside a treated-share-matched ($\bar D=0.652$) control confirming that same-week bias ($-0.046$ matched vs.\ $-0.057$ free at $w_0=0.30$) stems predominantly from unmodelled within-week nonlinearity rather than shrinking treated share, while Hansen's $J$ stays near nominal size (3.8\%--5.0\%). Panel~(c): removing the positive Gamma(1,5) small-count pseudo-count under higher-count $R=8000$ scenarios exposes severe negative feedback/spillover bias in M2 (e.g.\ Scenario~4 coverage drops from 90.5\% at $(1,5)$ to 0.0\% at $(0.1,5)$ and $(0.01,100)$; X1 drops from 93.8\% to 0.0\%). Sources: \texttt{{verification/non\_nested\_sim/results.json}}, \texttt{{verification/same\_week/grid\_results.json}}, \texttt{{verification/weak\_prior\_m2/results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_19_nonnested_weakprior.tex").write_text(tex, encoding="utf-8")


def write_table_20():
    cg = json.loads((ROOT / "verification" / "root_selection" / "cluster_grid_results.json").read_text(encoding="utf-8"))
    cells = cg["cells"]
    order = [
        ("G10_N2400_T12_bal", r"Grid A: $G=10, N=2400, T=12$, bal.\ (FB)"),
        ("G16_N2400_T12_bal", r"Grid A: $G=16, N=2400, T=12$, bal.\ (FB)"),
        ("G20_N2400_T12_bal", r"Grid A: $G=20, N=2400, T=12$, bal.\ (FB)"),
        ("G30_N2400_T12_bal", r"Grid A: $G=30, N=2400, T=12$, bal.\ (FB)"),
        ("G50_N2400_T12_bal", r"Grid A: $G=50, N=2400, T=12$, bal.\ (FB)"),
        ("G100_N2400_T12_bal", r"Grid A: $G=100, N=2400, T=12$, bal.\ (FB)"),
        ("G16_N400_T12_bal", r"Grid B: $G=16, N=400, T=12$, bal.\ (RD)"),
        ("G16_N400_T60_bal", r"Grid B: $G=16, N=400, T=60$, bal.\ (RD)"),
        ("G16_N400_T12_skew", r"Grid B: $G=16, N=400, T=12$, DE-skew (RD)"),
        ("G16_N400_T60_skew", r"Grid B: $G=16, N=400, T=60$, DE-skew (RD)"),
    ]
    lines = []
    for i, (k, lbl) in enumerate(order):
        c = cells[k]
        sep = "\\midrule\n" if i == 6 else ""
        lines.append(sep + " & ".join([
            lbl,
            f"{c['bias']:+.3f}",
            f"{c['empSD']:.3f}",
            f"{c['meanSE']:.3f} ({c['se_ratio']:.3f})",
            f"{c['cov_tG1']:.3f}",
            f"{c['cov_hwG2']:.3f}",
            f"{c['cov_jk']:.3f} ({c['jk_se_ratio']:.3f})",
            f"{c['cov_wcr_w6']:.3f}",
            f"{c['cov_wcr_rad']:.3f}",
        ]) + r" \\")

    tex = rf"""\begin{{table}}[!t]
\centering
\caption{{Controlled finite-cluster coverage across pure $G \in \{{10,\dots,100\}}$ at fixed $N=2400$ and across German cluster-size skewness $\times T$ at $N=400, G=16$ ($\tau+\kappa$, $R={cg['R']}$).}}
\label{{tab:rev-clustergrid}}
\footnotesize
\setlength{{\tabcolsep}}{{2.6pt}}
\begin{{adjustbox}}{{max width=\linewidth}}
\begin{{tabular}}{{lrrrrrrrr}}
\toprule
design cell & bias & emp.\ SD & mean SE (ratio) & $t_{{G-1}}$ & $t_{{G-2}}$ & JK (SE ratio) & WCR-W6 & WCR-Rad \\
\midrule
{chr(10).join(lines)}
\bottomrule
\end{{tabular}}
\end{{adjustbox}}
\par\smallskip\parbox{{\linewidth}}{{\scriptsize Grid~A isolates pure $G \in \{{10,16,20,30,50,100\}}$ under balanced clusters at fixed $N=2400, T=12$ (feedback adoption): sandwich SE ratios stay near 1 for $G \ge 16$ (0.947--0.984) and $t_{{G-1}}$ avoids under-coverage. Grid~B crosses balanced vs.\ German federal-state cluster sizes ($1$ to $96$ counties out of $400$) against $T \in \{{12, 60\}}$ profiled week effects under staggered random adoption (so late-horizon cases remain active): neither $T=60$ alone (SE ratio 1.012 balanced) nor DE-skewness at $T=12$ (SE ratio 0.934) replicates the German SE shortfall; it is specifically the interaction of severe cluster-size skewness and long-$T$ week profiling that lowers SE ratio to 0.768 at $T=60$, while WCR-score (Webb6: 0.952) and delete-one-cluster jackknife (0.996, SE ratio 1.237) maintain nominal coverage. Note $t_{{G-2}}$ mimics fixed-$G$ $t$-test adjustments \citep{{hwang2021simple}}. Source: \texttt{{verification/root\_selection/cluster\_grid\_results.json}}.}}
\end{{table}}
"""
    (TAB_DIR / "table_20_cluster_grid.tex").write_text(tex, encoding="utf-8")


if __name__ == "__main__":
    write_table_19()
    write_table_20()
    print("Generated tables/table_19_nonnested_weakprior.tex and tables/table_20_cluster_grid.tex.")
