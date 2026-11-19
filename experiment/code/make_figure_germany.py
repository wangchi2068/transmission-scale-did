"""Generate publication-grade Figure 3 (figures/fig_03_application.png/.pdf)
for the German county-level Varicella + staggered summer school holidays application.

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
    apply_publication_style(font_size=9)
    plt.rcParams["mathtext.fontset"] = "stix"

    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.45), gridspec_kw={"width_ratios": [1.22, 1.15, 1.05]})
    plt.subplots_adjust(wspace=0.34, left=0.14, right=0.985, top=0.88, bottom=0.18)

    # ------------------------------------------------------------------
    # Panel (a): Policy-level effect (tau + kappa or tau) across methods
    # ------------------------------------------------------------------
    ax = axes[0]
    m0 = res["primary"]["M0"]
    tc = m0["tcrit"]
    comps = res["comparators"]
    nuts = res["nuts2"]["M0"]
    sw05 = res["stock_wright"]["phi_bar_0.5"]["sets"]

    rows = [
        ("M0 ($t_{15}$)", m0["est"]["tpk"], m0["ci_t"]["tpk"][0], m0["ci_t"]["tpk"][1], PALETTE["blue_main"], "o"),
        ("M0 (pairs boot.)", m0["est"]["tpk"], res["primary"]["pairs_bootstrap"]["ci_tpk"][0],
         res["primary"]["pairs_bootstrap"]["ci_tpk"][1], PALETTE["blue_secondary"], "D"),
        ("M0 (SW proj.)", m0["est"]["tpk"], sw05["tpk_bonferroni_95"][0],
         sw05["tpk_bonferroni_95"][1], PALETTE["teal"], "s"),
        ("M0' (no spill.)", comps["M0_prime"]["est"]["tau"],
         comps["M0_prime"]["est"]["tau"] - tc * comps["M0_prime"]["se"]["tau"],
         comps["M0_prime"]["est"]["tau"] + tc * comps["M0_prime"]["se"]["tau"], PALETTE["neutral_dark"], ">"),
        ("M1 (Poisson FE)", comps["M1"]["est"]["tpk"], comps["M1"]["ci_t_tpk"][0],
         comps["M1"]["ci_t_tpk"][1], PALETTE["red_strong"], "^"),
        ("M2 (EpiEstim)", comps["M2"]["est"]["tau"], comps["M2"]["ci_t_tau"][0],
         comps["M2"]["ci_t_tau"][1], PALETTE["violet"], "P"),
        ("M3 (log TWFE)", comps["M3"]["est"]["tau"], comps["M3"]["ci_t_tau"][0],
         comps["M3"]["ci_t_tau"][1], PALETTE["green_3"], "X"),
        ("M4 (QMLE)", comps["M4"]["est"]["tpk"], comps["M4"]["ci_t_tpk"][0],
         comps["M4"]["ci_t_tpk"][1], PALETTE["neutral_mid"], "<"),
        ("M0 (NUTS-2)", nuts["est"]["tpk"], nuts["ci_t"]["tpk"][0],
         nuts["ci_t"]["tpk"][1], PALETTE["blue_main"], "v"),
    ]
    y_pos = np.arange(len(rows))[::-1]
    ax.axvline(0.0, color=PALETTE["neutral_black"], linestyle="--", linewidth=0.8, zorder=1)
    for y, (lbl, est, lo, hi, col, mk) in zip(y_pos, rows):
        ax.plot([lo, hi], [y, y], color=col, linewidth=1.6, zorder=2)
        ax.plot(est, y, marker=mk, color=col, markersize=5.0, zorder=3)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlabel(r"(a) Policy-level log-rate effect ($\tau+\kappa$ or $\tau$)")
    ax.set_xlim(-0.56, 0.16)
    ax.set_title("Estimators & intervals", pad=6)

    # ------------------------------------------------------------------
    # Panel (b): Event-time profile around Summer holiday onset
    # ------------------------------------------------------------------
    ax = axes[1]
    ev = res["event_time_with_leads"]["coefs"]
    ks = list(range(-2, 8))
    ests = [ev[f"e{k}"]["est"] for k in ks]
    los = [ev[f"e{k}"]["ci"][0] for k in ks]
    his = [ev[f"e{k}"]["ci"][1] for k in ks]
    ax.axhline(0.0, color=PALETTE["neutral_black"], linestyle="--", linewidth=0.8, zorder=1)
    ax.axvline(-0.5, color=PALETTE["neutral_mid"], linestyle=":", linewidth=0.9, zorder=1)
    # Shade post-incubation holiday window (weeks 1..5)
    ax.axvspan(0.5, 5.5, color=PALETTE["blue_main"], alpha=0.07, zorder=0)
    for k, e_, lo, hi in zip(ks, ests, los, his):
        col = PALETTE["neutral_mid"] if k < 0 else (PALETTE["teal"] if k == 0 else PALETTE["blue_main"])
        ax.plot([k, k], [lo, hi], color=col, linewidth=1.4, zorder=2)
    ax.plot(ks, ests, color=PALETTE["blue_main"], linewidth=1.2, zorder=3)
    for k, e_ in zip(ks, ests):
        col = PALETTE["neutral_mid"] if k < 0 else (PALETTE["teal"] if k == 0 else PALETTE["blue_main"])
        ax.plot(k, e_, marker="o", color=col, markersize=4.6, zorder=4)
    ax.set_xticks(ks)
    ax.set_xlim(-2.5, 7.5)
    ax.set_ylim(-0.62, 0.36)
    ax.set_xlabel("(b) Weeks from summer holiday onset ($e$)")
    ax.set_ylabel(r"Direct log-$R$ shift ($\hat\tau_e$)")
    ax.set_title("Dynamic event-time profile", pad=6)

    # ------------------------------------------------------------------
    # Panel (c): Stock-Wright S1(phi) profile and unique interior root
    # ------------------------------------------------------------------
    ax = axes[2]
    curve = res["stock_wright"]["phi_bar_0.5"]["curve"]
    phis = np.array([c["phi"] for c in curve])
    s1 = np.array([c["S1"] for c in curve])
    # F_{1,15} 95% critical value on S1 scale: S1 * (G-1)/(1*(G-1)) = S1 <= F_{0.95}(1, 15)
    crit_f = float(stats.f.ppf(0.95, 1, 15))
    phi_hat = m0["est"]["phi"]
    iv = sw05["S1_F_95"][0]

    ax.axvspan(iv[0], iv[1], color=PALETTE["blue_main"], alpha=0.10, zorder=0)
    ax.axhline(crit_f, color=PALETTE["red_strong"], linestyle="--", linewidth=0.9, label=r"$F_{1,15}$ 95% crit.")
    ax.plot(phis, s1, color=PALETTE["blue_main"], linewidth=1.6, zorder=3)
    ax.axvline(phi_hat, color=PALETTE["teal"], linestyle="-.", linewidth=1.0, zorder=2)
    ax.plot(phi_hat, 0.0, marker="o", color=PALETTE["teal"], markersize=5.2, zorder=4)
    ax.set_xlim(0.0, 3.0)
    ax.set_ylim(-0.35, 6.8)
    ax.set_xlabel(r"(c) Contagion coupling ($\phi$)")
    ax.set_ylabel("Profiled Stock\u2013Wright $S_1(\\phi)$")
    ax.set_title(r"Unique root ($\hat\phi = 0.325$) & $S_1$ set", pad=6)
    ax.text(0.62, 0.72, r"$\hat\phi = 0.325$", color=PALETTE["teal"], fontsize=8.5)
    ax.text(1.70, crit_f + 0.30, r"$F_{1,15}$ 95% crit.", color=PALETTE["red_strong"], fontsize=8)

    save_publication_bundle(fig, str(OUT_BASE), dpi=300, check_layout=True)
    plt.close(fig)
    print(f"Saved {OUT_BASE}.png and {OUT_BASE}.pdf successfully.")


if __name__ == "__main__":
    main()
