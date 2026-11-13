"""Publication figures from results.json (PDF + PNG through style_kit); writes figures/manifest.json.

Every figure is drawn only from results.json, so partial runs work: scenarios without replications are
skipped and the replications behind each panel are stated in the manifest caption. Error bars are
+/- 2 Monte Carlo standard errors (MCSE), never across-replication standard deviations.
Fonts: style_kit requires Times New Roman; fonts/*.ttf (created by fetch_fonts.sh) are registered first.
Run:  python src/make_figures.py [--results results.json] [--outdir figures]
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import matplotlib.ticker
import json
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.dirname(HERE)
sys.path.insert(0, HERE)
for _f in glob.glob(os.path.join(RUN, "fonts", "*.ttf")):
    font_manager.fontManager.addfont(_f)
import scenarios as sc  # noqa: E402
from style_kit import PALETTE, apply_publication_style, save_publication_bundle  # noqa: E402

COL = {"M0": PALETTE["blue_main"], "M0_G": PALETTE["teal"], "M0_ET": PALETTE["blue_secondary"],
       "M1": PALETTE["red_strong"], "M1_BC": PALETTE["red_2"], "M2": PALETTE["violet"], "M3": PALETTE["green_3"],
       "M4": PALETTE["neutral_mid"], "Oracle": PALETTE["neutral_black"], "M0p": PALETTE["neutral_light"],
       "M0_OI": PALETTE["blue_secondary"], "M0_GI": PALETTE["blue_secondary"]}
MARK = {"M0": "o", "M0_G": "D", "M0_ET": "v", "M1": "s", "M1_BC": "^", "M2": "P", "M3": "X", "M4": "<",
        "Oracle": "*", "M0p": ">", "M0_OI": "h", "M0_GI": "d"}
LAB = {"M0": "M0", "M0_G": "M0-G", "M0_ET": "M0-ET", "M1": "M1", "M1_BC": "M1-BC", "M2": "M2", "M3": "M3",
       "M4": "M4", "Oracle": "Oracle", "M0p": "M0'", "M0_OI": "M0-OI", "M0_GI": "M0-GI (post hoc)"}


class Res:
    def __init__(self, res):
        self.res = res
        self.idx = {(r["sid"], r["method"], r["param"], r["target"]): r for r in res.get("metrics_table", [])}
        self.scen = {int(k): v for k, v in res.get("scenarios", {}).items()}

    def g(self, sid, m, par, key, target="model"):
        r = self.idx.get((sid, m, par, target))
        v = None if r is None else r.get(key)
        return np.nan if v is None else float(v)

    def Rtxt(self, sids):
        Rs = [self.scen[s]["R"] for s in sids if s in self.scen]
        if not Rs:
            return None
        return "R = %d" % Rs[0] if min(Rs) == max(Rs) else "R = %d-%d" % (min(Rs), max(Rs))


def _label_note(fig, res):
    lab = (res.res.get("provenance") or {}).get("label")
    if lab:
        ax = fig.axes[-1]           # like style_kit.add_simulated_watermark: inside the last panel
        ax.text(0.99, 0.02, lab, transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="gray",
                alpha=0.8, style="italic")


def _dodge(k, n, width=0.7):
    return (k - (n - 1) / 2) * width / max(n, 1)


def core_labels(scen, sids):
    out = []
    for s in sids:
        p = scen[s].params
        out.append("%s%s\n%d/%d" % ("F" if p.adopt == "feedback" else "R", "S" if p.phi > 0 else "0", p.N, p.T))
    return out


# ------------------------------------------------------------------------------------------------
def fig_core(res, scen, par, meths, fname, ylab):
    sids = [s for s in sorted(scen) if scen[s].block == "core" and s in res.scen]
    if not sids:
        return None
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 4.6), sharex=True, constrained_layout=True)
    x = np.arange(len(sids))
    for k, m in enumerate(meths):
        b = np.array([res.g(s, m, par, "bias") for s in sids])
        bm = np.array([res.g(s, m, par, "bias_mcse") for s in sids])
        c = np.array([res.g(s, m, par, "coverage") for s in sids])
        cm = np.array([res.g(s, m, par, "coverage_mcse") for s in sids])
        if not np.isfinite(b).any():
            continue
        xx = x + _dodge(k, len(meths))
        axes[0].errorbar(xx, b, yerr=2 * bm, fmt=MARK[m], color=COL[m], ms=3.5, lw=0.8, capsize=0, label=LAB[m])
        axes[1].errorbar(xx, c, yerr=2 * cm, fmt=MARK[m], color=COL[m], ms=3.5, lw=0.8, capsize=0)
    axes[0].axhline(0, color=PALETTE["neutral_dark"], lw=0.6)
    axes[1].axhspan(0.93, 0.97, color=PALETTE["green_1"], zorder=0, lw=0)
    axes[1].axhline(0.95, color=PALETTE["neutral_dark"], lw=0.6)
    axes[0].set_ylabel("Bias of %s" % ylab)
    axes[1].set_ylabel("95% CI coverage")
    axes[1].set_ylim(-0.02, 1.02)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(core_labels(scen, sids), fontsize=7)
    axes[0].set_xlabel("(a) Scenario")
    axes[1].set_xlabel("(b) Scenario (F/R = feedback/random adoption, S/0 = spillover/none; N/T)")
    axes[0].legend(ncol=len(meths), loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=7)
    _label_note(fig, res)
    return fig, sids


def fig_null(res, scen):
    sids = [s for s in sorted(scen) if scen[s].block == "null" and s in res.scen]
    if not sids:
        return None
    meths = ["M0", "M0_G", "M1", "M1_BC", "M2", "M3", "M4", "Oracle"]
    fig, ax = plt.subplots(figsize=(5.0, 2.8), constrained_layout=True)
    x = np.arange(len(sids))
    for k, m in enumerate(meths):
        r = np.array([res.g(s, m, "tau", "reject0") for s in sids])
        rm = np.array([res.g(s, m, "tau", "reject0_mcse") for s in sids])
        if not np.isfinite(r).any():
            continue
        ax.errorbar(x + _dodge(k, len(meths)), r, yerr=2 * rm, fmt=MARK[m], color=COL[m], ms=3.5, lw=0.8,
                    label=LAB[m])
    ax.axhline(0.05, color=PALETTE["neutral_dark"], lw=0.6)
    ax.set_yscale("symlog", linthresh=0.1)
    ax.set_ylim(0, 1.05)
    lab = []
    for s in sids:
        p = scen[s].params
        lab.append("%s, k=%s" % ("feedback" if p.adopt == "feedback" else "random",
                                 "10" if math.isfinite(p.k) else "inf"))
    ax.set_xticks(x)
    ax.set_xticklabels(lab)
    ax.set_ylabel("Rejection rate of H0: tau = 0")
    ax.set_xlabel("Null scenario (N = 500, T = 12)")
    ax.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=7)
    _label_note(fig, res)
    return fig, sids


def fig_roots(res, scen):
    sids = [s for s in sorted(scen) if s in res.scen]
    if not sids:
        return None
    keys = [("M0G_multiple_roots", "more than one interior root", PALETTE["blue_main"]),
            ("M0G_at_boundary", "boundary selected (phi = 0)", PALETTE["red_strong"]),
            ("M0G_differs_from_M0", "M0-G differs from M0", PALETTE["teal"])]
    fig, ax = plt.subplots(figsize=(7.0, 2.6), constrained_layout=True)
    x = np.arange(len(sids))
    for k, (key, lab, c) in enumerate(keys):
        v = np.array([np.nan if res.scen[s].get(key, {}).get("p") is None else res.scen[s][key]["p"] for s in sids])
        e = np.array([np.nan if res.scen[s].get(key, {}).get("mcse") is None else res.scen[s][key]["mcse"]
                      for s in sids])
        ax.bar(x + _dodge(k, 3, 0.8), v, width=0.8 / 3, color=c, label=lab, yerr=2 * e,
               error_kw=dict(lw=0.6, ecolor=PALETTE["neutral_dark"]))
    ax.set_xticks(x)
    ax.set_xticklabels([str(s) for s in sids], fontsize=7)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Share of replications")
    ax.set_xlabel("Scenario id (1-16 core; 5-8, 13-16 without spillover; 17-20 null; 21-25 misspecified)")
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=7)
    _label_note(fig, res)
    return fig, sids


def fig_misspec(res, scen):
    sids = [s for s in sorted(scen) if scen[s].block == "misspec" and s in res.scen]
    if not sids:
        return None
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), constrained_layout=True,
                             gridspec_kw=dict(width_ratios=[1.3, 1]))
    meths = ["M0", "M0_G", "M1", "M2", "M3", "Oracle"]
    x = np.arange(len(sids))
    for k, m in enumerate(meths):
        b = np.array([res.g(s, m, "tau", "bias") for s in sids])
        bm = np.array([res.g(s, m, "tau", "bias_mcse") for s in sids])
        if np.isfinite(b).any():
            axes[0].errorbar(x + _dodge(k, len(meths)), b, yerr=2 * bm, fmt=MARK[m], color=COL[m], ms=3.5, lw=0.8,
                             label=LAB[m])
    axes[0].axhline(0, color=PALETTE["neutral_dark"], lw=0.6)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([scen[s].name.split("_")[0] for s in sids])
    axes[0].set_ylabel("Bias of tau")
    axes[0].set_xlabel("(a) Misspecification scenario")
    axes[0].legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=7)
    ax = axes[1]
    if 25 in res.scen:
        dyn = scen[25].params.dyn
        e = np.arange(7)
        truth = np.array([dyn[min(k, len(dyn) - 1)] for k in e])
        est = np.array([res.g(25, "M0_ET", "tau_e%d" % k, "mean_est") for k in e])
        se = np.array([res.g(25, "M0_ET", "tau_e%d" % k, "bias_mcse") for k in e])
        ax.plot(e, truth, color=PALETTE["neutral_dark"], lw=1.0, label="truth")
        ax.errorbar(e, est, yerr=2 * se, fmt="v", color=COL["M0_ET"], ms=3.5, lw=0.8, label="M0-ET mean")
        tau_c = res.g(25, "M0", "tau", "mean_est")
        if np.isfinite(tau_c):
            ax.axhline(tau_c, color=COL["M0"], lw=0.8, ls="--", label="M0 constant tau")
        ax.set_xticks(e)
        ax.set_xticklabels([str(k) for k in e[:-1]] + ["6+"])
        ax.legend(ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=7)
    ax.set_ylabel("Event-time effect")
    ax.set_xlabel("(b) Weeks since adoption (X5)")
    _label_note(fig, res)
    return fig, sids


def fig_application(res):
    app = res.res.get("application")
    if not app or not app.get("table"):
        return None
    rows = [r for r in app["table"] if r.get("tau") is not None and r.get("figure", True)]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 0.35 * len(rows) + 1.4), constrained_layout=True,
                             gridspec_kw=dict(width_ratios=[1.4, 1]))
    ax = axes[0]
    y = np.arange(len(rows))[::-1]
    # x-range from the estimates and the intervals that are not more than 4x the median width; wider
    # intervals are clipped at the axis limit and marked with an open arrow head
    widths = [r["tau_hi"] - r["tau_lo"] for r in rows if r.get("tau_lo") is not None and r.get("tau_hi") is not None]
    wmed = float(np.median(widths)) if widths else 0.0
    xs = [r["tau"] for r in rows]
    for r in rows:
        if r.get("tau_lo") is not None and r.get("tau_hi") is not None and r["tau_hi"] - r["tau_lo"] <= 4 * wmed:
            xs += [r["tau_lo"], r["tau_hi"]]
    pad = 0.08 * (max(xs) - min(xs) + 1e-9)
    xl, xh = min(min(xs), 0.0) - pad, max(max(xs), 0.0) + pad
    for yy, r in zip(y, rows):
        m = r["method"]
        lo, hi = r.get("tau_lo"), r.get("tau_hi")
        c = COL.get(m, PALETTE["neutral_mid"])
        if lo is not None and hi is not None:
            ax.plot([max(lo, xl), min(hi, xh)], [yy, yy], color=c, lw=1.0)
            if lo < xl:
                ax.plot(xl, yy, "<", color=c, ms=4, mfc="none")
            if hi > xh:
                ax.plot(xh, yy, ">", color=c, ms=4, mfc="none")
        ax.plot(r["tau"], yy, MARK.get(m, "o"), color=c, ms=4)
    ax.set_xlim(xl, xh)
    ax.axvline(0, color=PALETTE["neutral_dark"], lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(["%s: %s" % (LAB.get(r["method"], r["method"]), r["analysis"].replace("_", " "))
                        for r in rows], fontsize=7)
    ax.set_xlabel("(a) tau (log rate ratio), 95% CI")
    ax = axes[1]
    prof = app.get("profile") or {}
    if prof.get("phi"):
        ph = np.asarray(prof["phi"], float)
        Q = np.asarray([np.nan if v is None else v for v in prof.get("Q", [])], float)
        F = np.asarray([np.nan if v is None else v for v in prof.get("F", [])], float)
        if len(F) == len(ph):
            ax.plot(ph, F, "o-", color=COL["M0"], ms=3, lw=1.0, label="profiled phi-moment")
            ax.axhline(0, color=PALETTE["neutral_dark"], lw=0.6)
        hs, ls = ax.get_legend_handles_labels()
        if len(Q) == len(ph) and np.isfinite(Q).any():
            ax2 = ax.twinx()
            ax2.plot(ph, Q, "s--", color=COL["M0_G"], ms=3, lw=0.8, label="selection criterion Q")
            ax2.set_yscale("log")
            _fmt = matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % v)
            ax2.yaxis.set_major_formatter(_fmt)
            ax2.yaxis.set_minor_formatter(_fmt)
            ax2.set_ylabel("Q")
            ax2.grid(False)
            h2, l2 = ax2.get_legend_handles_labels()
            hs, ls = hs + h2, ls + l2
        ax.legend(hs, ls, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=7, ncol=2)
    ax.set_ylabel("Scaled moment")
    ax.set_xlabel("(b) phi")
    return fig, []


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(RUN, "results.json"))
    ap.add_argument("--outdir", default=os.path.join(RUN, "figures"))
    ap.add_argument("--only", default=None, help="comma-separated figure ids (used by figures/make_fig_*.py); "
                                                  "the manifest entries of the other figures are kept")
    a = ap.parse_args(argv)
    with open(a.results, "rb") as fh:
        res_sha = hashlib.sha256(fh.read()).hexdigest()
    with open(a.results) as fh:
        res = Res(json.load(fh))
    scen = sc.build_scenarios()
    apply_publication_style()
    os.makedirs(a.outdir, exist_ok=True)
    _gi = [m["R"] for m in res.res.get("metrics_table", []) if m.get("method") == "M0_GI" and m.get("param") == "tau"
           and m.get("target") == "model" and m.get("sid", 99) <= 16 and m.get("R")]
    gitxt = ("M0-GI R = %d" % _gi[0] if min(_gi) == max(_gi) else "M0-GI R = %d-%d" % (min(_gi), max(_gi))) \
        if _gi else "M0-GI not available"
    specs = [
        ("fig_01_core_tau", lambda: fig_core(res, scen, "tau", ["M0", "M0_G", "M0_GI", "M1", "M1_BC", "M2", "M3",
                                                                    "Oracle"],
                                             "fig_01", "tau"),
         "Bias (a) and coverage of 95%% t_{G-1} intervals (b) for the policy effect tau in the 16 core scenarios; "
         "error bars +/- 2 MCSE; shaded band 0.93-0.97. %s per scenario; M0-GI is the post-hoc root-selection "
         "variant (not pre-specified) evaluated on the same data sets over the replications its supplementary run "
         "reached ({GI} per scenario; Table tab:sim-m0gi).", "simulation"),
        ("fig_02_core_phi", lambda: fig_core(res, scen, "phi", ["M0", "M0_G", "M0_GI", "M1", "M4", "Oracle"], "fig_02",
                                             "phi"),
         "Bias (a) and coverage (b) for the spillover-scale parameter phi in the core scenarios (truth 0.4 with "
         "spillover, 0 without; M0-G and M0-GI coverage over interior selections). Error bars +/- 2 MCSE. %s per "
         "scenario; M0-GI = post-hoc variant ({GI} per scenario).",
         "simulation"),
        ("fig_03_core_kappa", lambda: fig_core(res, scen, "kappa", ["M0", "M0_G", "M1", "M4", "Oracle"], "fig_03",
                                               "kappa"),
         "Bias (a) and coverage (b) for the neighbour-policy effect kappa in the core scenarios. Error bars +/- 2 "
         "MCSE. %s per scenario.", "simulation"),
        ("fig_04_null_typeI", lambda: fig_null(res, scen),
         "Rejection rate of H0: tau = 0 at nominal 5%% (horizontal line) in the null block; error bars +/- 2 MCSE. "
         "%s per scenario.", "simulation"),
        ("fig_05_roots", lambda: fig_roots(res, scen),
         "Multiple roots and boundary solutions of the profiled QD-GMM equations by scenario (M0-G grid search); "
         "error bars +/- 2 MCSE. %s per scenario.", "simulation"),
        ("fig_06_misspec", lambda: fig_misspec(res, scen),
         "(a) Bias of tau under misspecification X1-X5 (+/- 2 MCSE). (b) X5 dynamic effects: mean M0-ET event-time "
         "estimates (+/- 2 MCSE) against the generating effects and the constant-tau M0 estimate. %s per scenario.",
         "simulation"),
        ("fig_07_application", lambda: fig_application(res),
         "Application to US county COVID-19 cases and state mask orders: (a) tau estimates with 95%% intervals "
         "across methods and pre-specified sensitivity analyses (open arrow heads: interval clipped at the axis limit); "
         "(b) main specification: profiled phi-moment over the phi grid (estimating-equation roots are its zero "
         "crossings) and the M0-G selection criterion Q evaluated along the profile (M0-G compares Q only at the "
         "roots and at the boundary phi = 0).%s", "application"),
    ]
    only = set(a.only.split(",")) if a.only else None
    man_path = os.path.join(a.outdir, "manifest.json")
    manifest = []
    if only and os.path.exists(man_path):
        with open(man_path) as fh:
            manifest = [m for m in json.load(fh) if m.get("id") not in only]
    for fid, fn, cap, section in specs:
        if only and fid not in only:
            continue
        out = fn()
        if out is None:
            print("[figures] %s skipped (no data yet)" % fid)
            continue
        fig, sids = out
        save_publication_bundle(fig, os.path.join(a.outdir, fid))
        plt.close(fig)
        rt = res.Rtxt(sids) if sids else ""
        def _sha(ext):
            pth = os.path.join(a.outdir, fid + ext)
            if not os.path.exists(pth):
                return None
            with open(pth, "rb") as fh:
                return hashlib.sha256(fh.read()).hexdigest()
        manifest.append(dict(id=fid, path_pdf=fid + ".pdf", path_png=fid + ".png",
                             sha256_pdf=_sha(".pdf"), sha256_png=_sha(".png"),
                             caption=(cap % (rt or "")).replace("{GI}", gitxt), section=section, generator="src/make_figures.py",
                             script="make_fig_%s.py" % fid[4:6], data_sources=["results.json"],
                             results_sha256=res_sha, scenarios=sids,
                             label=(res.res.get("provenance") or {}).get("label")))
        print("[figures] wrote %s" % fid)
    manifest.sort(key=lambda m: m["id"])
    with open(man_path, "w") as fh:
        json.dump(manifest, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
