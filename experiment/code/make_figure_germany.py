"""Generate publication-grade figures for the illustrative German varicella
implementation (Main Text Figure 3) and Web Appendix F diagnostics (Web Figure 5).

Reads:  experiment/analysis/germany_application_results.json
Writes: figures/fig_03_application.png (and .pdf)
        figures/fig_web_05_diagnostics.png (and .pdf)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiment" / "code"))
from style_kit import PALETTE, apply_publication_style, save_publication_bundle

RES_PATH = ROOT / "experiment" / "analysis" / "germany_application_results.json"
OUT_MAIN = ROOT / "figures" / "fig_03_application"
OUT_WEB5 = ROOT / "figures" / "fig_web_05_diagnostics"


def build_main_figure(res: dict) -> None:
    apply_publication_style(font_size=8.8)
    plt.rcParams["mathtext.fontset"] = "stix"

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2,
        figsize=(10.8, 4.15),
        gridspec_kw={"width_ratios": [1.38, 1.0]},
    )
    plt.subplots_adjust(wspace=0.28, left=0.175, right=0.985, top=0.89, bottom=0.15)

    # ------------------------------------------------------------------
    # Panel (a): Forest plot + right-aligned Est. (95% CI) table
    # ------------------------------------------------------------------
    m0 = res["primary"]["M0"]
    tc = m0["tcrit"]
    comps = res["comparators"]
    nuts_m0 = res["nuts2"]["M0"]
    nuts_m2 = res["nuts2"]["comparators"]["M2"]
    nuts_m2_ci = [
        nuts_m2["est"]["tau"] - tc * nuts_m2["se"]["tau"],
        nuts_m2["est"]["tau"] + tc * nuts_m2["se"]["tau"],
    ]

    rows = [
        ("QD-GMM M0 ($t_{15}$)", m0["est"]["tpk"], m0["ci_t"]["tpk"][0], m0["ci_t"]["tpk"][1],
         PALETTE["blue_main"], "o", "proposed"),
        ("QD-GMM M0 (WCR Webb6)", m0["est"]["tpk"], res["primary"]["wcr_webb6"]["ci"][0],
         res["primary"]["wcr_webb6"]["ci"][1], PALETTE["blue_secondary"], "P", "proposed"),
        ("QD-GMM M0 (pairs boot.)", m0["est"]["tpk"], res["primary"]["pairs_bootstrap"]["ci_tpk"][0],
         res["primary"]["pairs_bootstrap"]["ci_tpk"][1], PALETTE["blue_secondary"], "D", "proposed"),
        ("M0' (no spillover)", comps["M0_prime"]["est"]["tau"],
         comps["M0_prime"]["est"]["tau"] - tc * comps["M0_prime"]["se"]["tau"],
         comps["M0_prime"]["est"]["tau"] + tc * comps["M0_prime"]["se"]["tau"],
         PALETTE["neutral_dark"], ">", "comp"),
        ("M1 (Poisson FE)", comps["M1"]["est"]["tpk"], comps["M1"]["ci_t_tpk"][0],
         comps["M1"]["ci_t_tpk"][1], PALETTE["red_strong"], "^", "comp"),
        (r"M2 (EpiEstim $\log R$)", comps["M2"]["est"]["tau"], comps["M2"]["ci_t_tau"][0],
         comps["M2"]["ci_t_tau"][1], PALETTE["violet"], "s", "comp"),
        (r"M3 ($\log(Y+1)$ TWFE)", comps["M3"]["est"]["tau"], comps["M3"]["ci_t_tau"][0],
         comps["M3"]["ci_t_tau"][1], PALETTE["green_3"], "X", "comp"),
        ("M4 (Pooled QMLE)", comps["M4"]["est"]["tpk"], comps["M4"]["ci_t_tpk"][0],
         comps["M4"]["ci_t_tpk"][1], PALETTE["neutral_mid"], "<", "comp"),
        ("QD-GMM M0 (NUTS-2)", nuts_m0["est"]["tpk"], nuts_m0["ci_t"]["tpk"][0],
         nuts_m0["ci_t"]["tpk"][1], PALETTE["blue_main"], "v", "nuts"),
        (r"M2 (EpiEstim, NUTS-2)", nuts_m2["est"]["tau"], nuts_m2_ci[0],
         nuts_m2_ci[1], PALETTE["violet"], "s", "nuts"),
    ]

    # Visual y spacing with slight gap between groups
    y_vals = [10.6, 9.6, 8.6, 7.0, 6.0, 5.0, 4.0, 3.0, 1.4, 0.4]

    ax_a.axhspan(7.9, 11.3, color=PALETTE["blue_main"], alpha=0.04, zorder=0)
    ax_a.axhspan(-0.3, 2.1, color=PALETTE["teal"], alpha=0.05, zorder=0)
    ax_a.axhline(7.8, color=PALETTE["neutral_light"], linewidth=0.75, zorder=1)
    ax_a.axhline(2.2, color=PALETTE["neutral_light"], linewidth=0.75, zorder=1)
    ax_a.axvline(0.0, color=PALETTE["neutral_black"], linestyle="--", linewidth=0.85, zorder=1)

    for y, (lbl, est, lo, hi, col, mk, grp) in zip(y_vals, rows):
        mfc = "white" if grp == "nuts" and mk == "s" else col
        ax_a.plot([lo, hi], [y, y], color=col, linewidth=1.55, solid_capstyle="round", zorder=2)
        ax_a.plot(
            est, y, marker=mk, color=col, markerfacecolor=mfc,
            markeredgewidth=1.1, markersize=5.0, zorder=3,
        )
        txt = f"{est:+.3f} ({lo:+.3f}, {hi:+.3f})"
        ax_a.text(0.19, y, txt, va="center", ha="left", fontsize=7.6, color=PALETTE["neutral_black"])

    # Clean arrow showing M2 small-count sign reversal (58.8% zeros -> 2.9% zeros)
    ax_a.annotate(
        "", xy=(nuts_m2["est"]["tau"] + 0.015, 0.4),
        xytext=(comps["M2"]["est"]["tau"], 5.0),
        arrowprops=dict(
            arrowstyle="->,head_width=0.22,head_length=0.28",
            color=PALETTE["violet"], lw=1.0, ls=":",
            connectionstyle="arc3,rad=0.26",
        ),
        zorder=2,
    )
    ax_a.text(
        0.048, 2.85, "Zero share\n58.8% to 2.9%",
        fontsize=7.1, color=PALETTE["violet"], ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.90),
        zorder=4,
    )

    ax_a.text(0.19, 11.45, "Working est. (95% CI)", va="bottom", ha="left", fontsize=7.7, color=PALETTE["neutral_dark"])
    ax_a.text(-0.46, 11.45, "County level (58.8% zeros)", va="bottom", ha="left", fontsize=7.4, color=PALETTE["blue_main"])
    ax_a.text(-0.46, 1.75, "NUTS-2 (2.9% zeros)", va="bottom", ha="left", fontsize=7.2, color=PALETTE["teal"])

    ax_a.set_yticks(y_vals)
    ax_a.set_yticklabels([r[0] for r in rows], fontsize=8.0)
    ax_a.set_ylim(-0.4, 11.9)
    ax_a.set_xlim(-0.48, 0.56)
    ax_a.set_xticks([-0.4, -0.2, 0.0, 0.15])
    ax_a.set_xlabel(r"(a) Working-model log-rate coefficient ($\hat\tau+\hat\kappa$ or $\hat\tau$)")
    ax_a.set_title("Sparse-count estimator comparison", pad=8)

    # ------------------------------------------------------------------
    # Panel (b): Profiled Stock-Wright S1(phi) + zoomed root inset
    # ------------------------------------------------------------------
    sw05 = res["stock_wright"]["phi_bar_0.5"]
    curve = sw05["curve"]
    phis = np.array([c["phi"] for c in curve])
    s1 = np.array([c["S1"] for c in curve])
    crit_f = float(stats.f.ppf(0.95, 1, 15))
    phi_hat = m0["est"]["phi"]
    iv = sw05["sets"]["S1_F_95"][0]

    ax_b.axvspan(iv[0], iv[1], color=PALETTE["blue_main"], alpha=0.08, zorder=0, label=r"95% $F_{1,15}$ set")
    ax_b.axhline(crit_f, color=PALETTE["red_strong"], linestyle="--", linewidth=0.95, zorder=2)
    ax_b.plot(phis, s1, color=PALETTE["blue_main"], linewidth=1.75, zorder=3)
    ax_b.axvline(phi_hat, color=PALETTE["teal"], linestyle="-.", linewidth=1.0, zorder=2)
    ax_b.plot(phi_hat, 0.0, marker="o", color=PALETTE["teal"], markersize=5.4, zorder=4)
    ax_b.plot(0.0, s1[0], marker="o", color=PALETTE["red_strong"], markersize=4.4, zorder=4)

    ax_b.annotate(
        rf"$S_1(0) = {s1[0]:.2f}$",
        xy=(0.0, s1[0]), xytext=(0.44, 5.45),
        color=PALETTE["red_strong"], fontsize=7.8, va="center",
        arrowprops=dict(arrowstyle="->", color=PALETTE["red_strong"], lw=0.75, shrinkA=2, shrinkB=3),
    )
    ax_b.text(
        phi_hat + 0.14, 0.25, rf"$\hat\phi = {phi_hat:.3f}$",
        color=PALETTE["teal"], fontsize=8.3,
        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.88),
        zorder=4,
    )
    ax_b.text(2.94, crit_f + 0.20, rf"$F_{{1,15}}$ 95% crit. ({crit_f:.2f})",
              color=PALETTE["red_strong"], fontsize=7.6, ha="right")
    ax_b.text(1.48, 6.08, rf"95% $S_1(F)$ set: $[{iv[0]:.2f}, {iv[1]:.2f}]$",
              color=PALETTE["blue_main"], fontsize=7.8, ha="center")

    # Zoom inset cleanly lifted above outer bottom x-ticks inside open right space
    ax_in = inset_axes(
        ax_b, width="45%", height="42%", loc="lower right",
        bbox_to_anchor=(0.0, 0.10, 0.98, 0.86), bbox_transform=ax_b.transAxes,
    )
    mask = phis <= 0.62
    ax_in.axvspan(iv[0], 0.62, color=PALETTE["blue_main"], alpha=0.09, zorder=0)
    ax_in.axhline(crit_f, color=PALETTE["red_strong"], linestyle="--", linewidth=0.8, zorder=1)
    ax_in.plot(phis[mask], s1[mask], color=PALETTE["blue_main"], linewidth=1.45, zorder=2)
    ax_in.axvline(phi_hat, color=PALETTE["teal"], linestyle="-.", linewidth=0.85, zorder=2)
    ax_in.plot(phi_hat, 0.0, marker="o", color=PALETTE["teal"], markersize=4.2, zorder=3)
    ax_in.plot(0.0, s1[0], marker="o", color=PALETTE["red_strong"], markersize=3.8, zorder=3)
    ax_in.set_xlim(-0.02, 0.60)
    ax_in.set_ylim(-0.35, 6.5)
    ax_in.set_xticks([0.0, 0.21, 0.33, 0.5])
    ax_in.tick_params(labelsize=7.0)
    ax_in.set_title("Boundary detail", fontsize=7.3, pad=3)
    ax_in.spines["right"].set_visible(True)
    ax_in.spines["top"].set_visible(True)

    ax_b.set_xlim(-0.04, 3.02)
    ax_b.set_ylim(-0.35, 6.8)
    ax_b.set_xlabel(r"(b) Contagion coupling ($\phi$)")
    ax_b.set_ylabel(r"Profiled Stock--Wright $S_1(\phi)$")
    ax_b.set_title(r"Contagion profile & $S_1$ confidence set", pad=8)

    save_publication_bundle(fig, str(OUT_MAIN), dpi=300, check_layout=True)
    plt.close(fig)
    print(f"Saved {OUT_MAIN}.png and {OUT_MAIN}.pdf successfully.")


def build_web_figure_5(res: dict) -> None:
    apply_publication_style(font_size=8.8)
    plt.rcParams["mathtext.fontset"] = "stix"

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.6, 3.95), gridspec_kw={"width_ratios": [1.05, 1.12]})
    plt.subplots_adjust(wspace=0.26, left=0.075, right=0.985, top=0.88, bottom=0.17)

    # ------------------------------------------------------------------
    # Panel (a): Total-scale summer event-time profile
    # ------------------------------------------------------------------
    ts = res["event_time_with_leads"]["total_scale"]
    ks = list(range(-2, 8))
    ests = [ts[f"e{k}"]["est"] for k in ks]
    los = [ts[f"e{k}"]["ci"][0] for k in ks]
    his = [ts[f"e{k}"]["ci"][1] for k in ks]

    ax_a.axhline(0.0, color=PALETTE["neutral_black"], linestyle="--", linewidth=0.85, zorder=1)
    ax_a.axvline(-0.5, color=PALETTE["neutral_mid"], linestyle=":", linewidth=0.9, zorder=1)
    ax_a.axvspan(-0.5, 0.5, color=PALETTE["red_strong"], alpha=0.07, zorder=0)
    ax_a.axvspan(0.5, 5.5, color=PALETTE["blue_main"], alpha=0.06, zorder=0)

    for k, lo, hi in zip(ks, los, his):
        col = PALETTE["neutral_mid"] if k < 0 else (PALETTE["red_strong"] if k == 0 else PALETTE["blue_main"])
        ax_a.plot([k, k], [lo, hi], color=col, linewidth=1.45, solid_capstyle="round", zorder=2)
    ax_a.plot(ks, ests, color=PALETTE["blue_main"], linewidth=1.2, zorder=3)
    for k, e_ in zip(ks, ests):
        col = PALETTE["neutral_mid"] if k < 0 else (PALETTE["red_strong"] if k == 0 else PALETTE["blue_main"])
        ax_a.plot(k, e_, marker="o", color=col, markersize=4.8, zorder=4)

    ax_a.annotate(
        rf"Onset $\hat\tau_0 + c_0\hat\kappa = {ts['e0']['est']:+.3f}$ ($p=0.049$)",
        xy=(0.0, ts["e0"]["est"]),
        xytext=(0.75, 0.16),
        color=PALETTE["red_strong"],
        fontsize=7.6,
        arrowprops=dict(arrowstyle="->", color=PALETTE["red_strong"], lw=0.85),
    )
    ax_a.set_xticks(ks)
    ax_a.set_xlim(-2.5, 7.5)
    ax_a.set_ylim(-0.78, 0.32)
    ax_a.set_xlabel("(a) Weeks relative to summer holiday onset ($e$)")
    ax_a.set_ylabel(r"Total policy scale ($\hat\tau_e + c_e\hat\kappa$)")
    ax_a.set_title("Total-scale summer event-time profile", pad=8)

    # ------------------------------------------------------------------
    # Panel (b): Non-summer school holiday placebos (lags 0-2)
    # ------------------------------------------------------------------
    hols = [
        ("Easter", "easter", PALETTE["blue_main"], "o"),
        ("Autumn", "autumn", PALETTE["teal"], "s"),
        ("Winter", "winter", PALETTE["green_3"], "D"),
        ("Pentecost", "pentecost", PALETTE["violet"], "^"),
    ]
    lags = [0, 1, 2]
    offsets = [-0.18, -0.06, 0.06, 0.18]
    oh = res["placebos"]["other_holidays"]

    ax_b.axhline(0.0, color=PALETTE["neutral_black"], linestyle="--", linewidth=0.85, zorder=1)
    ax_b.axvspan(-0.42, 0.42, color=PALETTE["red_strong"], alpha=0.05, zorder=0)
    ax_b.axvspan(1.58, 2.42, color=PALETTE["teal"], alpha=0.06, zorder=0)

    for (label, key, col, mk), dx in zip(hols, offsets):
        xs = np.array(lags) + dx
        ys = [oh[key][f"lag_{l}"]["est"]["tpk"] for l in lags]
        yerr_lo = [y - oh[key][f"lag_{l}"]["ci_t"]["tpk"][0] for l, y in zip(lags, ys)]
        yerr_hi = [oh[key][f"lag_{l}"]["ci_t"]["tpk"][1] - y for l, y in zip(lags, ys)]
        ax_b.errorbar(
            xs, ys, yerr=[yerr_lo, yerr_hi],
            fmt=mk + "-", color=col, lw=1.15, elinewidth=1.35, markersize=4.8,
            capsize=0, label=label, zorder=3,
        )

    ax_b.text(0.0, 0.30, "Immediate dip\n(lag 0)", ha="center", va="center", fontsize=7.4, color=PALETTE["red_strong"])
    ax_b.text(2.0, -0.52, "Catch-up rebound\n(lag 2)", ha="center", va="center", fontsize=7.4, color=PALETTE["teal"])
    ax_b.set_xticks(lags)
    ax_b.set_xticklabels(["Lag 0\n(holiday week)", "Lag 1\n(primary alignment)", "Lag 2\n(post-break)"])
    ax_b.set_xlim(-0.48, 2.48)
    ax_b.set_ylim(-0.68, 0.49)
    ax_b.set_xlabel("(b) Non-summer school holiday alignment")
    ax_b.set_ylabel(r"Working contrast ($\hat\tau+\hat\kappa$, 95% $t_{15}$ CI)")
    ax_b.set_title("Non-summer holiday placebo cycles", pad=22)
    ax_b.legend(
        frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.005),
        handletextpad=0.35, columnspacing=0.85,
    )

    save_publication_bundle(fig, str(OUT_WEB5), dpi=300, check_layout=True)
    plt.close(fig)
    print(f"Saved {OUT_WEB5}.png and {OUT_WEB5}.pdf successfully.")


def main():
    res = json.loads(RES_PATH.read_text(encoding="utf-8"))
    build_main_figure(res)
    build_web_figure_5(res)


if __name__ == "__main__":
    main()
