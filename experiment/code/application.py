"""Application: state public mask orders and county-level COVID-19 transmission, US, June-November 2020
(explore/method_brief.md section 9; data sources and processing in data_contract.md).

Input  : data/processed/panel.npz (built by src/data_prep.py from the raw public files).
Output : application/application_results.json (read by aggregate.py -> results.json["application"],
         make_tables.py table_7 and make_figures.py fig_07).

Pre-specified analyses
  main          W = row-normalised 2010 Census county adjacency; generation-interval (GI) weights from the
                Ferretti et al. (2020) Weibull generation time discretised to weeks; omega_i = pop_i / 1e5;
                phi_bar = 0.5; clusters = states (G = 49).  Methods: M0, M0-G (with phi profile), M0-OI
                (Hansen J), M0', M0-ET (event time), M1, M1-BC, M2, M3, M4; cases averted psi for M0 / M0-G.
  sensitivity   GI from the Nishiura et al. (2020) serial interval; GI grid w = (a, 1-a), a in {0.6,...,1.0};
                commuting-flow W (ACS 2011-2015); policy lagged one week; omega x 0.25 and x 4;
                phi_bar in {0.25, 1}; excluding Kansas and Texas (county-level opt-outs/exemptions not coded
                in the CDC file); conventional lead placebo (invalid under feedback, Proposition 2(c);
                reported for comparison only).
Every sensitivity analysis is fitted with M0 and M0-G.
Usage: python src/application.py [--out application/application_results.json] [--quick]
"""
from __future__ import annotations

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    __import__("os").environ.setdefault(_v, "1")

import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np
import scipy.sparse as sps
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comparators as cmp  # noqa: E402
import dgp  # noqa: E402
import psi as psimod  # noqa: E402
import qdgmm  # noqa: E402

RUN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PANEL = os.path.join(RUN, "data", "processed", "panel.npz")
VERSION = "1.0.0"
S = dgp.S_LAGS
STATE_KS, STATE_TX = 20, 48

# Published distributions (parameters as reported in the cited papers)
GI_SOURCES = {
    "ferretti": dict(dist="weibull", shape=2.826, scale=5.665,
                     what="generation time, Weibull (Ferretti et al. 2020, Science 368:eabb6936, Table 1)"),
    "nishiura": dict(dist="lognormal", mean=4.7, sd=2.9,
                     what="serial interval, lognormal, mean 4.7 d, SD 2.9 d (Nishiura, Linton & Akhmetzhanov "
                          "2020, Int J Infect Dis 93:284-286)"),
}
GI_GRID = (0.6, 0.7, 0.8, 0.9, 1.0)
PSI_HORIZONS = (4, 8, 13)


# ------------------------------------------------------------------------------------------------
# Generation interval: daily continuous distribution -> weekly lag weights
# ------------------------------------------------------------------------------------------------
def _cdf(src: dict):
    if src["dist"] == "weibull":
        return stats.weibull_min(src["shape"], scale=src["scale"]).cdf
    if src["dist"] == "lognormal":
        s2 = np.log(1.0 + (src["sd"] / src["mean"]) ** 2)
        mu = np.log(src["mean"]) - s2 / 2
        return stats.lognorm(np.sqrt(s2), scale=np.exp(mu)).cdf
    raise ValueError(src["dist"])


def gi_weekly(src: dict, S: int = S, n_u: int = 7000) -> dict:
    """P(lag = s weeks) when the infector's case falls uniformly within its week and the delay X has the
    published distribution: p_s = E_u[F(7(s+1) - u) - F(7s - u)], u ~ U[0, 7).  Same-week transmission
    (s = 0) cannot enter a lagged renewal term; weights for s = 1..S are renormalised to sum to one
    (methods.tex: sum_s w_s = 1) and the excluded mass is reported."""
    F = _cdf(src)
    u = (np.arange(n_u) + 0.5) / n_u * 7.0
    p = np.array([np.mean(F(np.clip(7 * (s + 1) - u, 0, None)) - F(np.clip(7 * s - u, 0, None)))
                  for s in range(S + 1)])
    w = p[1:] / p[1:].sum()
    return dict(w=w.tolist(), p_by_week_lag=p.tolist(), mass_same_week=float(p[0]),
                mass_beyond_S=float(1.0 - p.sum()), source=src)


# ------------------------------------------------------------------------------------------------
# Data
# ------------------------------------------------------------------------------------------------
def row_normalise(A: sps.spmatrix) -> sps.csr_matrix:
    A = sps.csr_matrix(A, dtype=float)
    A.setdiag(0.0)
    A.eliminate_zeros()
    deg = np.asarray(A.sum(1)).ravel()
    inv = np.where(deg > 0, 1.0 / np.where(deg > 0, deg, 1.0), 0.0)
    return (sps.diags(inv) @ A).tocsr()


def load_data(path: str = PANEL) -> dict:
    z = np.load(path, allow_pickle=False)
    d = {k: z[k] for k in z.files}
    N = d["Y"].shape[1]
    e = d["edges"].astype(int)
    A = sps.coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(N, N))
    d["W_adj"] = row_normalise(A + A.T)
    c = d["comm"]
    if len(c):
        C = sps.coo_matrix((c[:, 2], (c[:, 0].astype(int), c[:, 1].astype(int))), shape=(N, N)).tocsr()
        d["W_comm"] = row_normalise(C + C.T)      # flows in both directions between i and j
    else:
        d["W_comm"] = None
    states = np.unique(d["state"])
    d["cl"] = np.searchsorted(states, d["state"])
    d["states"] = states
    return d


def verify(d: dict) -> dict:
    """Key-variable checks (fail loudly on violation)."""
    Y, D, pop = d["Y"], d["D"], d["pop"]
    T_all, N = Y.shape
    chk = dict(shape_ok=Y.shape == D.shape == d["Dshare"].shape and pop.shape == (N,),
               Y_nonneg_integer=bool(np.all(Y >= 0) and np.all(Y == np.round(Y))),
               D_binary=bool(np.isin(D, (0.0, 1.0)).all()),
               pop_positive=bool(np.all(pop > 0)), weeks=int(T_all), units=int(N),
               states=int(len(d["states"])), no_isolated_units=bool(np.all(np.asarray(d["W_adj"].sum(1)) > 0)))
    bad = [k for k, v in chk.items() if v is False]
    if bad:
        raise RuntimeError("data checks failed: %s" % bad)
    return chk


def make_panel(d, w, W, D_all=None, omega_mult=1.0, units=None, drop_presample=False):
    D_all = d["D"] if D_all is None else D_all
    Yall, cl, pop = d["Y"], d["cl"], d["pop"]
    if drop_presample:      # M0-ET2: remove units already under the policy in the pre-sample weeks
        keep = ~D_all[:S].any(0)
        units = keep if units is None else units[keep[units]]
    if units is not None:
        W = W[:, units]     # subnetwork adjacency restricted to kept units; row-normalised after subsetting
        Yall, D_all, pop = Yall[:, units], D_all[:, units], pop[units]
        W = row_normalise(W[units])
        cl = np.unique(cl[units], return_inverse=True)[1]
    omega = omega_mult * pop / 1e5
    panel = dgp.build_panel(Yall, D_all, W, omega, cl, w=w, S=S)
    return qdgmm.attach_presample(panel, D_all[:S]), Yall   # Dpre: pre-sample policy path for spec_m0_et2


# ------------------------------------------------------------------------------------------------
# Fitting helpers
# ------------------------------------------------------------------------------------------------
def _f(x):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def ci_row(analysis, method, res, tc, names=("tau", "kappa", "phi"), figure=False, extra=None, panel=None, spec=None):
    r = dict(analysis=analysis, method=method, ok=bool(res.get("ok")), figure=figure)
    est, se = res.get("est", {}), res.get("se", {})
    for n in names:
        e, s = _f(est.get(n)), _f(se.get(n))
        r[n], r[n + "_se"] = e, s
        r[n + "_lo"] = None if e is None or s is None else e - tc * s
        r[n + "_hi"] = None if e is None or s is None else e + tc * s
    if panel is not None and spec is not None and res.get("ok") and {"tau", "kappa"} <= set(est):
        # total effect tau+kappa from the joint covariance (Reviewer M3: report the policy-level target)
        try:
            mod = qdgmm.QDModel(panel, spec)
            inf = mod.inference(np.array([est[n] for n in spec.names]))
            if inf is not None:
                i, j = list(spec.names).index("tau"), list(spec.names).index("kappa")
                V = inf["Vb"]
                tot = est["tau"] + est["kappa"]
                s_tot = float(np.sqrt(V[i, i] + V[j, j] + 2 * V[i, j]))
                r.update(tau_plus_kappa=float(tot), tau_plus_kappa_se=s_tot,
                         tau_plus_kappa_lo=float(tot - tc * s_tot), tau_plus_kappa_hi=float(tot + tc * s_tot))
        except Exception:  # noqa: BLE001 - the column stays empty if the covariance is unavailable
            pass
    if extra:
        r.update(extra)
    return r


def slim(res):
    keep = ("ok", "est", "se", "reason", "resid", "start_used", "cond", "iters", "J", "J_df", "J_p", "n_roots",
            "roots", "at_boundary", "chosen", "Q", "profile")
    out = {k: res[k] for k in keep if k in res}
    return json.loads(json.dumps(out, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))


def timed(fn, *a, **k):
    t0 = time.perf_counter()
    try:
        out = fn(*a, **k)
    except Exception as e:  # noqa: BLE001 - recorded, never fatal
        out = dict(ok=False, est={}, se={}, reason="exception: %r" % (e,))
    out["time_sec"] = time.perf_counter() - t0
    return out


def phi_profile(panel, phi_bar=0.5, grid=qdgmm.PHI_GRID):
    """Profiled phi-moment F(phi) (scaled) and the selection criterion Q at each grid value (M0-G ingredients)."""
    sp = qdgmm.spec_m0(panel, phi_bar)
    mod = qdgmm.QDModel(panel, sp)
    mod_oi = qdgmm.QDModel(panel, qdgmm.spec_m0_oi(panel, phi_bar))
    x = np.zeros(mod.K)
    out = dict(phi=[], F=[], Q=[], tau=[], kappa=[], ok=[])
    for ph in grid:
        b, ok, Fp = qdgmm._solve_given_phi(mod, ph, x)
        if ok:
            x = b
        out["phi"].append(float(ph)); out["ok"].append(bool(ok))
        out["F"].append(_f(Fp) if ok else None)
        out["Q"].append(_f(qdgmm.selection_criterion(mod_oi, np.r_[b, ph])) if ok else None)
        out["tau"].append(_f(b[0]) if ok else None); out["kappa"].append(_f(b[1]) if ok else None)
    return out


def psi_block(Yall, panel, res, label):
    if not res.get("ok") or res.get("beta") is None:
        return dict(method=label, ok=False, reason="estimator failed")
    beta = np.asarray(res["beta"], float)

    def idx(b):
        return b[0] * panel["D"] + b[1] * panel["Db"], b[2]
    out = psimod.psi_with_se(Yall, panel, beta, res.get("Vb"), idx)
    # cumulative psi over the first h window weeks (Theorem 4 with P = {(i, t): t < h}); the full-window
    # psi compounds the per-week multiplier through renewal and contagion, so shorter horizons are reported too
    Vb = res.get("Vb")
    T = panel["Y"].shape[0]
    hs = [h for h in PSI_HORIZONS if h < T] + [T]

    def cum(b):
        ix, ph = idx(b)
        nu = psimod.nu_recursion(Yall, ix, panel["W"], panel["omega"], panel["w"], ph, S)
        dw = (nu - Yall[S:]).sum(1)
        return np.array([dw[:h].sum() for h in hs])
    base = cum(beta)
    se = np.full(len(hs), np.nan)
    if Vb is not None and np.all(np.isfinite(Vb)):
        hstep = 1e-5
        grad = np.column_stack([(cum(beta + hstep * e) - cum(beta - hstep * e)) / (2 * hstep)
                                for e in np.eye(len(beta))])
        se = np.sqrt(np.clip(np.einsum("ij,jk,ik->i", grad, Vb, grad), 0, None))
    out["by_horizon"] = [dict(weeks=int(h), est=_f(v), se=_f(s)) for h, v, s in zip(hs, base, se)]
    out["method"] = label
    out["observed_cases_treated"] = float(Yall[S:][panel["D"] == 1].sum())
    out["observed_cases_all"] = float(Yall[S:].sum())
    return out


def event_time_block(res, spet, tc):
    """Event-time bins (tau_e*), equal-weight average over the single-week bins (open-ended last bin excluded)
    with its t_{G-1} interval, and the cluster-robust Wald test of H0: all event-time effects = 0
    (F form, (G-1) denominator; chi-square form reported alongside)."""
    G = int(spet.G)   # clusters of the panel the spec was built on (the M0-ET2 panel drops some states)
    names = [n for n in spet.names if n.startswith("tau_e")]
    out = dict(ok=res.get("ok"), bins=names, G=G,
               est={n: _f(res.get("est", {}).get(n)) for n in spet.names},
               se={n: _f(res.get("se", {}).get(n)) for n in spet.names},
               n_cells={n: float(spet.lin[n].sum()) for n in names})
    if not (res.get("ok") and res.get("Vb") is not None):
        return out
    ix = [spet.names.index(n) for n in names]
    Vb = np.asarray(res["Vb"])
    bvec = np.array([res["est"][n] for n in names])
    a_w = np.zeros(len(names))
    a_w[:-1] = 1.0 / (len(names) - 1)
    Ve = Vb[np.ix_(ix, ix)]
    avg, avg_se = float(a_w @ bvec), float(np.sqrt(a_w @ Ve @ a_w))
    W = float(bvec @ np.linalg.solve(Ve, bvec))
    q = len(names)
    Fst = W / q * (G - q) / (G - 1)
    out["average_e0_to_e%d" % (len(names) - 2)] = dict(
        est=avg, se=avg_se, lo=avg - tc * avg_se, hi=avg + tc * avg_se,
        weights="equal over single-week bins, open-ended last bin excluded")
    out["joint_wald"] = dict(W=W, q=q, F=float(Fst), df1=q, df2=int(G - q),
                             p_F=float(stats.f.sf(Fst, q, G - q)), p_chi2=float(stats.chi2.sf(W, q)))
    return out


def unit_counts(d, panel):
    """Sample composition of the application panel by treatment timing (counts of units)."""
    D, ev = panel["D"], panel["ev"]
    pre = qdgmm.pre_exposed(panel)
    ever = D.any(0)
    Dpre = np.asarray(panel["Dpre"])
    return dict(
        n_units_full_sample=int(D.shape[1]),
        n_units_m0et2=int((~pre).sum()), n_states=int(d["states"].size),
        n_never_treated_window=int((~ever).sum()),
        n_ever_treated_window=int(ever.sum()),
        n_presample_exposed=int(pre.sum()),
        n_presample_exposed_ever_window=int((pre & ever).sum()),
        n_presample_exposed_not_window=int((pre & ~ever).sum()),
        n_switch_window_week0_not_presample=int((~pre & (ev == 0)).sum()),
        n_switch_window_week_ge1=int((ev >= 1).sum()),
        presample_weeks=int(Dpre.shape[0]),
        week0_treated_window=int(D[0].sum()))


def lead_indicator(D):
    """1{policy switches on next week} (the week before adoption); last window week has no lead -> 0."""
    L = np.zeros_like(D)
    L[:-1] = np.clip(D[1:] - D[:-1], 0, None)
    return L


# ------------------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", default=PANEL)
    ap.add_argument("--out", default=os.path.join(RUN, "application", "application_results.json"))
    ap.add_argument("--quick", action="store_true", help="main analysis only (smoke test)")
    a = ap.parse_args(argv)
    t_start = time.time()
    d = load_data(a.panel)
    checks = verify(d)
    N = d["Y"].shape[1]
    G = len(d["states"])
    tc = qdgmm.tcrit(G)
    gis = {k: gi_weekly(v) for k, v in GI_SOURCES.items()}
    w_main = gis["ferretti"]["w"]
    print("[app] %d units, %d states, GI weights ferretti=%s nishiura=%s" %
          (N, G, np.round(w_main, 4), np.round(gis["nishiura"]["w"], 4)), flush=True)

    panel, Yall = make_panel(d, w_main, d["W_adj"])
    analyses, table, psi, notes = {}, [], {}, []

    # ---------------- main analysis: all methods ----------------
    main_res = {}
    m0 = timed(qdgmm.fit, panel, qdgmm.spec_m0(panel))
    main_res["M0"] = m0
    mg = timed(qdgmm.fit_global, panel)
    main_res["M0_G"] = mg
    # M0-GI: post-hoc root-selection variant (run-sim-application, experiment_design.md deviation 7)
    mgi = timed(qdgmm.fit_global, panel, selection="inv")
    main_res["M0_GI"] = mgi
    main_res["M0_OI"] = timed(qdgmm.fit_two_step, panel, qdgmm.spec_m0_oi(panel), m0) if m0.get("ok") else \
        dict(ok=False, est={}, se={}, reason="M0 failed")
    main_res["M0p"] = timed(qdgmm.fit, panel, qdgmm.spec_m0prime(panel))
    spet = qdgmm.spec_m0_et(panel)
    met = timed(qdgmm.fit, panel, spet)
    # M0-ET2 (revision): pre-sample-exposed units (Dall[:S] = 1, ev = 0 only inside the window) form their
    # own bin tau_pre instead of being pooled into tau_e0; spec_m0_et and event_time are kept unchanged.
    panel2, Yall2 = make_panel(d, w_main, d["W_adj"], drop_presample=True)
    spet2 = qdgmm.spec_m0_et2(panel2)
    met2 = timed(qdgmm.fit, panel2, spet2)
    m1 = timed(cmp.poisson_fe, panel)
    main_res["M1"] = m1
    main_res["M1_BC"] = timed(cmp.poisson_fe_bc, panel, m1)
    main_res["M2"] = timed(cmp.m2_logr_did, panel)
    main_res["M3"] = timed(cmp.m3_logy_twfe, panel)
    main_res["M4"] = timed(cmp.poisson_qmle, panel)
    for k, v in main_res.items():
        print("[app] main %-6s ok=%s tau=%s (%.1fs)" % (k, v.get("ok"), _f(v.get("est", {}).get("tau")),
                                                         v.get("time_sec", 0)), flush=True)
    analyses["main"] = {k: slim(v) for k, v in main_res.items()}
    analyses["main"]["M0_ET"] = slim(met)
    analyses["main"]["M0_ET2"] = slim(met2)
    for k in ("M0", "M0_G", "M0_GI", "M0_OI", "M0p", "M1", "M1_BC", "M2", "M3", "M4"):
        names = ("tau",) if k in ("M2", "M3") else ("tau", "kappa", "phi")
        if k == "M0p":
            names = ("tau",)
        extra = {"at_boundary": main_res[k].get("at_boundary")} if k in ("M0_G", "M0_GI") else None
        panel_k, spec_k = (panel, qdgmm.spec_m0(panel)) if k in ("M0", "M0_G", "M0_GI") else (None, None)
        table.append(ci_row("main", k, main_res[k], tc, names, figure=True, extra=extra,
                            panel=panel_k, spec=spec_k))
    J = {k: _f(main_res["M0_OI"].get(k)) for k in ("J", "J_df", "J_p")}
    profile = phi_profile(panel)
    profile.update(roots=mg.get("roots"), chosen=mg.get("chosen"), at_boundary=mg.get("at_boundary"))
    event_time = event_time_block(met, spet, tc)
    event_time2 = event_time_block(met2, spet2, qdgmm.tcrit(int(panel2["G"])))
    event_time2["spec"] = ("pre-exposed units (Dall[:S] = 1) removed from the sample; tau_e0 = switches on in "
                           "window week 0; bins tau_e1..tau_e6 and kappa as in M0-ET")
    event_time2["sample_counts"] = unit_counts(d, panel)
    print("[app] main M0_ET2 ok=%s fields=%s (%.1fs)" % (met2.get("ok"), spet2.names, met2.get("time_sec", 0)),
          flush=True)

    for lab, res in (("M0", m0), ("M0_G", mg), ("M0_GI", mgi)):
        psi[lab] = json.loads(json.dumps(timed(psi_block, Yall, panel, res, lab), default=float))
    desc = dgp.panel_descriptives(panel)

    # ---------------- sensitivity analyses: M0 and M0-G ----------------
    sens = []
    if not a.quick:
        sens.append(("gi_nishiura", dict(w=gis["nishiura"]["w"]), True))
        for g in GI_GRID:
            sens.append(("gi_w1_%.1f" % g, dict(w=[g, 1.0 - g]), False))
        if d["W_comm"] is not None:
            sens.append(("W_commuting", dict(W=d["W_comm"]), True))
        else:
            notes.append("commuting W unavailable (ACS file missing); sensitivity skipped")
        D_lag = np.vstack([d["D"][:1], d["D"][:-1]])
        sens.append(("policy_lag1", dict(D_all=D_lag), True))
        sens.append(("omega_x0.25", dict(omega_mult=0.25), False))
        sens.append(("omega_x4", dict(omega_mult=4.0), False))
        sens.append(("phibar_0.25", dict(phi_bar=0.25), False))
        sens.append(("phibar_1", dict(phi_bar=1.0), False))
        keep = np.where(~np.isin(d["state"], (STATE_KS, STATE_TX)))[0]
        sens.append(("excl_KS_TX", dict(units=keep), True))
    for name, kw, fig in sens:
        phi_bar = kw.pop("phi_bar", 0.5)
        pnl, Ya = make_panel(d, kw.pop("w", w_main), kw.pop("W", d["W_adj"]), **kw)
        r0 = timed(qdgmm.fit, pnl, qdgmm.spec_m0(pnl, phi_bar))
        rg = timed(qdgmm.fit_global, pnl, phi_bar)
        rgi = timed(qdgmm.fit_global, pnl, phi_bar, selection="inv")
        analyses[name] = {"M0": slim(r0), "M0_G": slim(rg), "M0_GI": slim(rgi), "n_units": int(pnl["Y"].shape[1]),
                          "G": int(pnl["G"])}
        tcs = qdgmm.tcrit(int(pnl["G"]))
        tab_spec = qdgmm.spec_m0(pnl, phi_bar)
        table.append(ci_row(name, "M0", r0, tcs, figure=fig, panel=pnl, spec=tab_spec))
        table.append(ci_row(name, "M0_G", rg, tcs, figure=False, extra={"at_boundary": rg.get("at_boundary")},
                            panel=pnl, spec=tab_spec))
        table.append(ci_row(name, "M0_GI", rgi, tcs, figure=False, extra={"at_boundary": rgi.get("at_boundary")},
                            panel=pnl, spec=tab_spec))
        print("[app] %-12s M0 tau=%s phi=%s | M0-G tau=%s phi=%s (%s) | M0-GI tau=%s (%s)" % (
            name, _f(r0.get("est", {}).get("tau")), _f(r0.get("est", {}).get("phi")),
            _f(rg.get("est", {}).get("tau")), _f(rg.get("est", {}).get("phi")), rg.get("chosen"),
            _f(rgi.get("est", {}).get("tau")), rgi.get("chosen")), flush=True)
    lead = None
    if not a.quick:
        sp = qdgmm.QDSpec(lin={"tau": panel["D"], "kappa": panel["Db"], "lead": lead_indicator(panel["D"])},
                          phi_free=True, phi_bar=0.5)
        rl = timed(qdgmm.fit, panel, sp)
        analyses["lead_placebo"] = {"M0": slim(rl)}
        e, s = _f(rl.get("est", {}).get("lead")), _f(rl.get("se", {}).get("lead"))
        lead = dict(est=e, se=s, t=None if e is None or s is None else e / s,
                    note="conventional lead placebo; not a valid test under feedback (Proposition 2(c))")
        table.append(ci_row("lead_placebo", "M0", rl, tc, figure=False))
        # the row must report the lead coefficient itself, not tau (Reviewer m3, round 2)
        table[-1]["lead"] = e
        table[-1]["lead_se"] = s
        table[-1]["lead_lo"] = None if e is None or s is None else e - tc * s
        table[-1]["lead_hi"] = None if e is None or s is None else e + tc * s

    out = dict(
        task="epi-novel-method", version=VERSION, created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        runtime_sec=round(time.time() - t_start, 1), mode="quick" if a.quick else "full",
        n_units=N, G=G, T=int(panel["Y"].shape[0]), S=S, first_week=str(d["week_start"][S]),
        last_week=str(d["week_start"][-1]), ci="95%% t_{G-1}, G = %d" % G, tcrit=tc,
        data_checks=checks, descriptives=desc, gi=gis, omega="pop/1e5", phi_bar=0.5,
        table=table, analyses=analyses, profile=profile, J=J, psi=psi, event_time=event_time,
        event_time2=event_time2, lead_placebo=lead,
        kansas_subset="not run: the CDC file codes state orders only; Kansas county opt-outs are not coded",
        notes=notes, git_sha=_git_sha())
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(_clean(out), fh, indent=1)
    os.replace(tmp, a.out)
    print("[app] wrote %s in %.1f s" % (a.out, time.time() - t_start), flush=True)
    return 0


def _clean(x):
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return _clean(x.tolist())
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return _f(x)
    return x


def _git_sha():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=RUN, capture_output=True, text=True,
                              timeout=10).stdout.strip() or None
    except Exception:  # noqa: BLE001
        return None


if __name__ == "__main__":
    sys.exit(main())
