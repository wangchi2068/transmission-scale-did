"""Delete-one-state jackknife and wild cluster bootstrap for the application estimate.

Why the jackknife standard error is larger than the sandwich here
---------------------------------------------------------------
Deleting one state is not a small perturbation: Texas contributes 254 of the 3091 counties, and
removing it also removes its policy from every neighbour's pressure term. The delete-one jackknife
therefore measures sensitivity to the composition of the cluster set, not the sampling variability
that the clustered sandwich estimates. Both are reported, and the text says which notion each one
describes: the sandwich is the reference for the reported interval, the jackknife indicates how much
the estimate moves when a large cluster is dropped, and the wild bootstrap-t is applied to the
sandwich so that its reference distribution accounts for the number of clusters.

Run: cd experiment/code && python -I ../../verification/application/wild_bootstrap.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CODE = HERE.parents[1] / "experiment" / "code"
sys.path.insert(0, str(CODE))

import application as A  # noqa: E402
import qdgmm  # noqa: E402


def delete_one(d, w, names, est):
    """Refit with each state removed; returns (deviations, sizes) both indexed by state."""
    theta = np.array([est[n] for n in names])
    states = np.unique(d["state"])
    dev = np.full((len(names), states.size), np.nan)
    sizes = np.array([(d["state"] == s).sum() for s in states])
    for k, s in enumerate(states):
        keep = np.where(d["state"] != s)[0]
        pnl, _ = A.make_panel(d, w, d["W_adj"], units=keep)
        f = qdgmm.fit(pnl, qdgmm.spec_m0(pnl), starts=[theta])
        if f.get("ok"):
            dev[:, k] = np.array([f["est"][n] for n in names]) - theta
    return dev, sizes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=999)
    a = ap.parse_args()
    d = A.load_data()
    A.verify(d)
    gi = A.gi_weekly(A.GI_SOURCES["ferretti"])
    w = gi["w"] if isinstance(gi, dict) else gi
    panel, _ = A.make_panel(d, w, d["W_adj"])
    spec = qdgmm.spec_m0(panel)
    names = list(spec.names)
    m0 = qdgmm.fit(panel, spec)
    Vb = qdgmm.QDModel(panel, spec).inference(np.array([m0["est"][n] for n in names]))["Vb"]
    se_s = np.sqrt(np.diag(Vb))
    est = np.array([m0["est"][n] for n in names])

    dev, sizes = delete_one(d, w, names, m0["est"])
    ok = np.all(np.isfinite(dev), axis=0)
    X, sz = dev[:, ok], sizes[ok]
    G = X.shape[1]
    se_j = np.sqrt((G - 1) / G * ((X - X.mean(1, keepdims=True)) ** 2).sum(1))

    # wild bootstrap-t on the sandwich: Rademacher weights applied to the cluster deviations,
    # studentised by the weighted sandwich standard error
    rng = np.random.default_rng(20261124)
    t = np.zeros((a.B, len(names)))
    for b in range(a.B):
        wt = rng.choice([-1.0, 1.0], size=G)
        d_b = (X * wt[None, :]).sum(1) / G
        om_b = ((X * wt[None, :]) ** 2).sum(1) / (G * (G - 1))
        t[b] = d_b / np.sqrt(np.maximum(om_b, 1e-300))
    q = np.quantile(np.abs(t), [0.9, 0.95, 0.975], axis=0)

    res = dict(names=names, est=est.tolist(), se_sandwich=se_s.tolist(), se_jack=se_j.tolist(),
               q90=q[0].tolist(), q95=q[1].tolist(), q975=q[2].tolist(),
               lo=(est - q[1] * se_s).tolist(), hi=(est + q[1] * se_s).tolist(),
               largest_state=int(sz.max()), n_states=int(G), B=a.B)
    print("states", G, "largest", int(sz.max()))
    print("sandwich se", se_s.round(4), "| jackknife se", se_j.round(4))
    for i, n in enumerate(names):
        print("  %-6s est %8.4f  CI(t48) [%8.4f,%8.4f]  wild q95 %5.3f  CI [%8.4f,%8.4f]"
              % (n, est[i], est[i] - 2.0106 * se_s[i], est[i] + 2.0106 * se_s[i], q[1][i],
                 est[i] - q[1][i] * se_s[i], est[i] + q[1][i] * se_s[i]))

    tot = float(est[0] + est[1])
    se_t = float(np.sqrt(Vb[0, 0] + Vb[1, 1] + 2 * Vb[0, 1]))
    xd = X[0] + X[1]
    rng = np.random.default_rng(20261125)
    tt = np.zeros(a.B)
    for b in range(a.B):
        wt = rng.choice([-1.0, 1.0], size=G)
        db = (xd * wt).sum() / G
        omb = ((xd * wt) ** 2).sum() / (G * (G - 1))
        tt[b] = db / np.sqrt(max(omb, 1e-300))
    qt = float(np.quantile(np.abs(tt), 0.95))
    res["tau_plus_kappa"] = dict(est=tot, se=se_t, q95=qt, lo=tot - qt * se_t, hi=tot + qt * se_t,
                                 ci_t48=[tot - 2.0106 * se_t, tot + 2.0106 * se_t])
    print("  tau+k  est %8.4f  CI(t48) [%8.4f,%8.4f]  wild q95 %5.3f  CI [%8.4f,%8.4f]"
          % (tot, tot - 2.0106 * se_t, tot + 2.0106 * se_t, qt, tot - qt * se_t, tot + qt * se_t))
    (HERE / "wild_bootstrap.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print("written")


if __name__ == "__main__":
    main()
