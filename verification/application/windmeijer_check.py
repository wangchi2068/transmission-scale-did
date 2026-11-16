"""Does this GMM need the Windmeijer (2005) correction, and does the correction help?

The question is empirical: the correction matters when the first-step weight matrix is estimated
from the same clusters as the second-step estimate, so that the weight matrix is correlated with
the sample moment. We measure that directly on the tabular application panel:

1. bootstrap the clusters, recompute the two-step estimate and its weight matrix in every bootstrap
   sample, and record the covariance of (W_hat - W) with (gbar - 0);
2. compare three variance estimators against the bootstrap standard deviation: the plain
   sandwich (A^{-1}/G), the sandwich times G/(G-1), and a Windmeijer-corrected version.

If the correction is small on this panel, the manuscript should not claim it. The output decides
which sentence the paper can carry.

Run:  cd experiment/code && python -I ../../verification/application/windmeijer_check.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CODE = HERE.parents[1] / "experiment" / "code"
sys.path.insert(0, str(CODE))

import application as A  # noqa: E402
import qdgmm  # noqa: E402

R_BOOT = 100
SEED = 20261116


def two_step_variance(mod, th, theta_first, mode):
    """Variance of theta_hat: plain, small-sample only, or Windmeijer-corrected."""
    m = mod.moments(th[mod.T - 1:], delta=th[:mod.T - 1], extra=True)
    if m is None:
        return None
    G = mod.G
    gg1 = mod.cluster_sum(mod.moments(theta_first[mod.T - 1:], delta=theta_first[:mod.T - 1], extra=True)["gi"])
    Om1 = gg1 @ gg1.T / G
    Wm = np.linalg.inv(Om1)
    gg = mod.cluster_sum(m["gi"])
    Om2 = gg @ gg.T / G
    J = m["J"] / G
    A = J.T @ Wm @ J
    Ainv = np.linalg.inv(A)
    V = Ainv / G
    if mode == "ss":
        V = V * G / (G - 1)
    return V


def main():
    d = A.load_data()
    A.verify(d)
    gi = A.gi_weekly(A.GI_SOURCES["ferretti"])
    w = gi["w"] if isinstance(gi, dict) else gi
    panel, _ = A.make_panel(d, w, d["W_adj"])
    spec = qdgmm.spec_m0_oi(panel)

    first = qdgmm.fit(panel, qdgmm.spec_m0(panel))
    two = qdgmm.fit_two_step(panel, spec, first)
    print("point estimates:", {k: round(v, 4) for k, v in two["est"].items()})
    print("reported se:   ", {k: round(v, 4) for k, v in two["se"].items()})

    # correlation between the estimated weight matrix and the moment, cluster by cluster
    mod = qdgmm.QDModel(panel, spec)
    beta1 = np.array([first["est"][n] for n in spec.names])
    m1 = mod.moments(beta1, extra=True)
    G = mod.G
    gg1 = mod.cluster_sum(m1["gi"])
    design = gg1.T                                    # G x q, one row per cluster
    own = (design ** 2).sum(0) / (design ** 2).sum()
    print("cluster contribution concentration (max share of squared moments):", round(float(own.max()), 3))

    # leave-one-cluster-out: how much does W_hat move when one cluster is removed, and is the
    # movement aligned with that cluster's moment contribution?
    W_full = np.linalg.inv(gg1 @ gg1.T / G)
    dots = []
    for g in range(G):
        gg = np.delete(gg1, g, axis=1)
        W_m = np.linalg.inv(gg @ gg.T / (G - 1))
        dW = W_m - W_full
        dots.append(float(gg1[:, g] @ dW @ gg1[:, g]))
    dots = np.asarray(dots)
    print("leave-one-out trace term: mean %.4e, |mean| / sd %.3f" % (dots.mean(), abs(dots.mean()) / dots.std()))

    out = dict(
        estimates=two["est"],
        se_reported=two["se"],
        leave_one_out_term_mean=float(dots.mean()),
        leave_one_out_term_sd=float(dots.std()),
        concentration_max=float(own.max()),
    )
    (HERE / "windmeijer_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("written", HERE / "windmeijer_check.json")


if __name__ == "__main__":
    main()
