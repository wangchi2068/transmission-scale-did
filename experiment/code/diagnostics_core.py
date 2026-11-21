"""Exact diagnostic checks for core simulation edge behaviors (Minors 5, 6, 7, 11):
1. Scenario 5 (N=500, T=12, FB, no-spill): distribution of M0 tau SE (median, IQR, p95, p99, max)
   explaining relerr_SE mean=5.089 / MCSE=16523 from near-singular boundary Jacobians.
2. Scenarios 14 & 16: unweighted Oracle Poisson QMLE cluster mass concentration (max cluster share / HHI)
   vs M0 weighted 1/(Lambda+1) scaling explaining Oracle under-coverage (80.5%, 80.2%).
3. Scenarios 2 & 4: exact M1-BC half-panel failure share due to zero within-unit variation in D under feedback.
4. Scenario 22 (X2): exact psi^dagger vs psi_hat medians and absolute bias distribution explaining +1004.5 rel. bias.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiment" / "code"))

import comparators as cmp  # noqa: E402
import dgp  # noqa: E402
import psi as psi_mod  # noqa: E402
import qdgmm  # noqa: E402
import scenarios as sc  # noqa: E402

OUT_PATH = ROOT / "experiment" / "analysis" / "core_diagnostics.json"


def _run_scen_5(scen_map: dict, R: int = 500) -> dict:
    s5 = scen_map[5]
    geo = dgp.make_geography(s5.params.N)
    ses, taus, phis, min_svals = [], [], [], []
    n_fail = 0
    for r in range(R):
        seed = sc.seed_for(5, r)
        sim = dgp.simulate(s5.params, geo, [seed, 0], [seed, 1])
        pan = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl, w=s5.w_est)
        res = qdgmm.fit(pan, qdgmm.spec_m0(pan))
        if not res["ok"]:
            n_fail += 1
            continue
        ses.append(float(res["se"]["tau"]))
        taus.append(float(res["est"]["tau"]))
        phis.append(float(res["est"]["phi"]))
        # Check minimum singular value of profiled Jacobian if available via QDModel
        mod = res["model"]
        mm = mod.moments(res["beta"]) if hasattr(mod, "moments") else None
        if mm is not None:
            min_svals.append(float(np.linalg.svd(mm["J"], compute_uv=False)[-1]))
    arr_se = np.asarray(ses)
    arr_t = np.asarray(taus)
    emp_sd = float(arr_t.std(ddof=1))
    emp_iqr_sd = float((np.quantile(arr_t, 0.75) - np.quantile(arr_t, 0.25)) / 1.349)
    top5 = sorted([(float(se), float(ph), float(ta)) for se, ph, ta in zip(ses, phis, taus)], reverse=True)[:5]
    return dict(
        R=R, n_ok=len(ses), n_fail=n_fail,
        emp_sd=emp_sd, emp_iqr_sd=emp_iqr_sd,
        se_mean=float(arr_se.mean()), se_median=float(np.median(arr_se)),
        se_iqr=float(np.quantile(arr_se, 0.75) - np.quantile(arr_se, 0.25)),
        se_p95=float(np.quantile(arr_se, 0.95)), se_p99=float(np.quantile(arr_se, 0.99)),
        se_max=float(arr_se.max()),
        share_se_gt_0_5=float((arr_se > 0.5).mean()),
        relerr_median_se=float(np.median(arr_se) / emp_sd - 1.0),
        top5_outliers=[dict(se=a, phi=b, tau=c) for a, b, c in top5],
    )


def _run_oracle_mass(scen_map: dict, sid: int, R: int = 100) -> dict:
    s_obj = scen_map[sid]
    geo = dgp.make_geography(s_obj.params.N)
    cl = geo.cl
    G = geo.G
    max_raw, hhi_raw, max_m0, hhi_m0, se_ratios = [], [], [], [], []
    cov_cnt = 0
    ses_or, taus_or = [], []
    for r in range(R):
        seed = sc.seed_for(sid, r)
        sim = dgp.simulate(s_obj.params, geo, [seed, 0], [seed, 1])
        pan = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl, w=s_obj.w_est)
        # Cluster mass in unweighted Poisson score ~ sum_t Y_{it}, vs M0 weighted score ~ sum_t Y_{it} / (L + omega + 1)
        Y_g = np.bincount(cl, pan["Y"].sum(0), minlength=G)
        sh_raw = Y_g / max(float(Y_g.sum()), 1.0)
        W_m0 = pan["Y"] / (pan["L"] + pan["omega"][None, :] + 1.0)
        Wm_g = np.bincount(cl, W_m0.sum(0), minlength=G)
        sh_m0 = Wm_g / max(float(Wm_g.sum()), 1.0)
        max_raw.append(float(sh_raw.max()))
        hhi_raw.append(float(np.sum(sh_raw ** 2)))
        max_m0.append(float(sh_m0.max()))
        hhi_m0.append(float(np.sum(sh_m0 ** 2)))
        ora = cmp.poisson_qmle(pan, offset=sim["alpha"])
        if ora["ok"]:
            taus_or.append(float(ora["est"]["tau"]))
            ses_or.append(float(ora["se"]["tau"]))
    arr_t = np.asarray(taus_or)
    arr_se = np.asarray(ses_or)
    emp_sd = float(arr_t.std(ddof=1))
    return dict(
        sid=sid, R=R,
        raw_max_cluster_share_mean=float(np.mean(max_raw)),
        raw_max_cluster_share_max=float(np.max(max_raw)),
        raw_hhi_mean=float(np.mean(hhi_raw)),
        raw_effective_G_mean=float(np.mean([1.0 / h for h in hhi_raw])),
        m0_max_cluster_share_mean=float(np.mean(max_m0)),
        m0_max_cluster_share_max=float(np.max(max_m0)),
        m0_effective_G_mean=float(np.mean([1.0 / h for h in hhi_m0])),
        oracle_emp_sd=emp_sd,
        oracle_mean_se=float(arr_se.mean()),
        oracle_se_ratio=float(arr_se.mean() / emp_sd),
    )


def _run_x2_psi(scen_map: dict, R: int = 200) -> dict:
    s22 = scen_map[22]
    geo = dgp.make_geography(s22.params.N)
    psi_dag, psi_hat, abs_errs, rel_errs = [], [], [], []
    for r in range(R):
        seed = sc.seed_for(22, r)
        sim = dgp.simulate(s22.params, geo, [seed, 0], [seed, 1])
        pan = dgp.build_panel(sim["Yall"], sim["Dall"], geo.W, sim["omega"], geo.cl, w=s22.w_est)
        res = qdgmm.fit(pan, qdgmm.spec_m0(pan))
        if not res["ok"]:
            continue
        idx0 = s22.params.tau * pan["D"] + s22.params.kappa * pan["Db"]
        p_dag = psi_mod.psi_values(
            sim["Yall"], idx0, geo.W, sim["omega"], s22.w_est, s22.params.phi, pan["D"]
        )["treated"]
        idxh = res["est"]["tau"] * pan["D"] + res["est"]["kappa"] * pan["Db"]
        p_hat = psi_mod.psi_values(
            sim["Yall"], idxh, geo.W, sim["omega"], s22.w_est, res["est"]["phi"], pan["D"]
        )["treated"]
        psi_dag.append(float(p_dag))
        psi_hat.append(float(p_hat))
        abs_errs.append(float(p_hat - p_dag))
        rel_errs.append(float((p_hat - p_dag) / max(abs(p_dag), 1.0)))
    ad, ah, ae, re = map(np.asarray, (psi_dag, psi_hat, abs_errs, rel_errs))
    return dict(
        R=R, n_ok=len(ae),
        psi_dagger_median=float(np.median(ad)), psi_dagger_mean=float(ad.mean()),
        psi_hat_median=float(np.median(ah)), psi_hat_mean=float(ah.mean()),
        abs_bias_mean=float(ae.mean()), abs_bias_median=float(np.median(ae)),
        rel_bias_median=float(np.median(re)), rel_bias_mean=float(re.mean()),
    )


def main():
    cal = sc.load_calibration()
    scen_map = sc.build_scenarios(cal)
    res_json = json.loads((ROOT / "experiment" / "analysis" / "results.json").read_text(encoding="utf-8"))
    m1bc_fails = {
        str(sid): dict(
            conv=res_json["scenarios"][str(sid)]["convergence"]["M1_BC"],
            reasons=res_json["scenarios"][str(sid)]["failure_reasons"]["M1_BC"],
        )
        for sid in (1, 2, 3, 4, 21, 22, 23, 24, 25)
    }
    out = dict(
        scen5_se_outliers=_run_scen_5(scen_map, R=400),
        oracle_cluster_mass={
            "scen14": _run_oracle_mass(scen_map, 14, R=80),
            "scen16": _run_oracle_mass(scen_map, 16, R=80),
        },
        m1bc_failures=m1bc_fails,
        x2_psi_abs_bias=_run_x2_psi(scen_map, R=100),
    )
    OUT_PATH.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
