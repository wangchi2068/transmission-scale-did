"""Generate publication-grade Figure 3 (figures/fig_03_application.png/.pdf)
for the German county-level Varicella empirical stress test (2015--2019).

Reads:  experiment/analysis/germany_application_results.json
Writes: figures/fig_03_application.png (and .pdf)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiment" / "code"))
from style_kit import PALETTE, apply_publication_style, save_publication_bundle

RES_PATH = ROOT / "experiment" / "analysis" / "germany_application_results.json"
OUT_BASE = ROOT / "figures" / "fig_03_application"


def main():
    res = json.loads(RES_PATH.read_text(encoding="utf-8"))
    apply_publication_style(font_size=8.8)
    plt.rcParams["mathtext.fontset"] = "stix"

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.75), gridspec_kw={"width_ratios": [1.32, 1.15, 1.05]})
    plt.subplots_adjust(wspace=0.35, left=0.16, right=0.985, top=0.88, bottom=0.18)

    # ------------------------------------------------------------------
    # Panel (a): Policy-level effect across estimators, finite-cluster intervals & seasonal/state checks
    # ------------------------------------------------------------------
    ax = axes[0]
    m0 = res["primary"]["M0"]
    tc = m0["tcrit"]
    comps = res["comparators"]
    nuts = res["nuts2"]["M0"]
    sw05 = res["stock_wright"]["phi_bar_0.5"]["sets"]
    reg = res["regional_fe"]

    rows = [
        ("M0 ($t_{15}$)", m0["est"]["tpk"], m0["ci_t"]["tpk"][0], m0["ci_t"]["tpk"][1], PALETTE["blue_main"], "o"),
        ("M0 (WCR Webb6)", m0["est"]["tpk"], res["primary"]["wcr_webb6"]["ci"][0],
         res["primary"]["wcr_webb6"]["ci"][1], PALETTE["blue_secondary"], "P"),
        ("M0 (pairs boot.)", m0["est"]["tpk"], res["primary"]["pairs_bootstrap"]["ci_tpk"][0],
         res["primary"]["pairs_bootstrap"]["ci_tpk"][1], PALETTE["blue_secondary"], "D"),
        ("M0 (SW proj.)", m0["est"]["tpk"], sw05["tpk_bonferroni_95"][0],
         sw05["tpk_bonferroni_95"][1], PALETTE["teal"], "s"),
        (r"M0 (4-reg $\times$ wk)", reg["four_region_x_week_lag1"]["est"]["tpk"],
         reg["four_region_x_week_lag1"]["ci_t"]["tpk"][0],
         reg["four_region_x_week_lag1"]["ci_t"]["tpk"][1], "#C86D2D", "h"),
        ("M0 (Drop BW&BY L1)", reg["drop_bw_by_lag1"]["est"]["tpk"],
         reg["drop_bw_by_lag1"]["ci_t"]["tpk"][0],
         reg["drop_bw_by_lag1"]["ci_t"]["tpk"][1], "#C86D2D", "^"),
        ("M0 (Drop BW&BY L2)", reg["drop_bw_by_lag2"]["est"]["tpk"],
         reg["drop_bw_by_lag2"]["ci_t"]["tpk"][0],
         reg["drop_bw_by_lag2"]["ci_t"]["tpk"][1], "#C86D2D", "v"),
        ("M0' (no spill.)", comps["M0_prime"]["est"]["tau"],
         comps["M0_prime"]["est"]["tau"] - tc * comps["M0_prime"]["se"]["tau"],
         comps["M0_prime"]["est"]["tau"] + tc * comps["M0_prime"]["se"]["tau"], PALETTE["neutral_dark"], ">"),
        ("M1 (Poisson FE)", comps["M1"]["est"]["tpk"], comps["M1"]["ci_t_tpk"][0],
         comps["M1"]["ci_t_tpk"][1], PALETTE["red_strong"], "^"),
        ("M2 (EpiEstim)", comps["M2"]["est"]["tau"], comps["M2"]["ci_t_tau"][0],
         comps["M2"]["ci_t_tau"][1], PALETTE["violet"], "P"),
        ("M3 (log TWFE)", comps["M3"]["est"]["tau"], comps["M3"]["ci_t_tau"][0],
         comps["M3"]["ci_t_tau"][1], PALETTE["green_3"], "X"),
        ("M0 (NUTS-2)", nuts["est"]["tpk"], nuts["ci_t"]["tpk"][0],
         nuts["ci_t"]["tpk"][1], PALETTE["blue_main"], "v"),
    ]
    y_pos = np.arange(len(rows))[::-1]
    ax.axvline(0.0, color=PALETTE["neutral_black"], linestyle="--", linewidth=0.8, zorder=1)
    for y, (lbl, est, lo, hi, col, mk) in zip(y_pos, rows):
        ax.plot([lo, hi], [y, y], color=col, linewidth=1.45, zorder=2)
        ax.plot(est, y, marker=mk, color=col, markersize=4.6, zorder=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.8)
    ax.set_xlabel(r"(a) Policy-level log-rate contrast ($\tau+\kappa$ or $\tau$)")
    ax.set_xlim(-0.56, 0.28)
    ax.set_title("Estimators, intervals & sensitivity", pad=6)

    # ------------------------------------------------------------------
    # Panel (b): Event-time profile ON THE TOTAL-POLICY SCALE (tau_e + c_e * kappa)
    # ------------------------------------------------------------------
    ax = axes[1]
    ts = res["event_time_with_leads"]["total_scale"]
    ks = list(range(-2, 8))
    ests = [ts[f"e{k}"]["est"] for k in ks]
    los = [ts[f"e{k}"]["ci"][0] for k in ks]
    his = [ts[f"e{k}"]["ci"][1] for k in ks]
    ax.axhline(0.0, color=PALETTE["neutral_black"], linestyle="--", linewidth=0.8, zorder=1)
    ax.axvline(-0.5, color=PALETTE["neutral_mid"], linestyle=":", linewidth=0.9, zorder=1)
    ax.axvspan(-0.5, 0.5, color=PALETTE["red_strong"], alpha=0.07, zorder=0)
    ax.axvspan(0.5, 5.5, color=PALETTE["blue_main"], alpha=0.07, zorder=0)
    for k, e_, lo, hi in zip(ks, ests, los, his):
        col = PALETTE["neutral_mid"] if k < 0 else (PALETTE["red_strong"] if k == 0 else PALETTE["blue_main"])
        ax.plot([k, k], [lo, hi], color=col, linewidth=1.35, zorder=2)
    ax.plot(ks, ests, color=PALETTE["blue_main"], linewidth=1.15, zorder=3)
    for k, e_ in zip(ks, ests):
        col = PALETTE["neutral_mid"] if k < 0 else (PALETTE["red_strong"] if k == 0 else PALETTE["blue_main"])
        ax.plot(k, e_, marker="o", color=col, markersize=4.5, zorder=4)
    ax.set_xticks(ks)
    ax.set_xlim(-2.5, 7.5)
    ax.set_ylim(-0.78, 0.32)
    ax.set_xlabel("(b) Weeks from summer holiday onset ($e$)")
    ax.set_ylabel(r"Total policy scale ($\hat\tau_e + c_e\hat\kappa$)")
    ax.set_title("Total-scale event-time profile", pad=6)
    ax.annotate(
        rf"$\hat\tau_0 + c_0\hat\kappa = {ts['e0']['est']:+.3f}$",
        xy=(0.0, ts["e0"]["est"]),
        xytext=(1.1, 0.14),
        color=PALETTE["red_strong"],
        fontsize=7.6,
        arrowprops=dict(arrowstyle="->", color=PALETTE["red_strong"], lw=0.8),
    )

    # ------------------------------------------------------------------
    # Panel (c): Stock-Wright S1(phi) profile and interior root
    # ------------------------------------------------------------------
    ax = axes[2]
    curve = res["stock_wright"]["phi_bar_0.5"]["curve"]
    phis = np.array([c["phi"] for c in curve])
    s1 = np.array([c["S1"] for c in curve])
    crit_f = float(stats.f.ppf(0.95, 1, 15))
    phi_hat = m0["est"]["phi"]
    iv = sw05["S1_F_95"][0]

    ax.axvspan(iv[0], iv[1], color=PALETTE["blue_main"], alpha=0.10, zorder=0)
    ax.axhline(crit_f, color=PALETTE["red_strong"], linestyle="--", linewidth=0.9)
    ax.plot(phis, s1, color=PALETTE["blue_main"], linewidth=1.6, zorder=3)
    ax.axvline(phi_hat, color=PALETTE["teal"], linestyle="-.", linewidth=1.0, zorder=2)
    ax.plot(phi_hat, 0.0, marker="o", color=PALETTE["teal"], markersize=5.0, zorder=4)
    ax.set_xlim(0.0, 3.0)
    ax.set_ylim(-0.35, 6.8)
    ax.set_xlabel(r"(c) Contagion coupling ($\phi$)")
    ax.set_ylabel("Profiled Stock\u2013Wright $S_1(\\phi)$")
    ax.set_title(rf"Interior root ($\hat\phi = {phi_hat:.3f}$) & $S_1$ set", pad=6)
    ax.text(0.60, 0.72, rf"$\hat\phi = {phi_hat:.3f}$", color=PALETTE["teal"], fontsize=8.3)
    ax.text(1.65, crit_f + 0.30, r"$F_{1,15}$ 95% crit.", color=PALETTE["red_strong"], fontsize=8)

    save_publication_bundle(fig, str(OUT_BASE), dpi=300, check_layout=True)
    plt.close(fig)
    print(f"Saved {OUT_BASE}.png and {OUT_BASE}.pdf successfully.")


if __name__ == "__main__":
    main()
