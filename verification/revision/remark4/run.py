"""Review Comment 1 ("a point in the authors' favour") / Remark 4: adjacent-cluster score dependence in the
application, computed with the paper's own code and data (port of the reviewer's check to the revision tree).
At the M0 estimate, cluster score sums g_c (q x G):
  Omega      = sum_c g_c g_c' / G                       (cluster-robust middle matrix used in the paper)
  Omega_adj  = Omega + sum_{c ~ c'} g_c g_c'' / G        (adds cross-products of adjacent states, both orders)
Reports the adjacent-pair score correlation per parameter moment and the SEs of tau, kappa, phi, tau+kappa.
Deterministic.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import env  # noqa
import numpy as np
import application as A
import qdgmm

HERE = Path(__file__).resolve().parent


def main():
    d = A.load_data()
    w = A.gi_weekly(A.GI_SOURCES["ferretti"])["w"]
    panel, _ = A.make_panel(d, w, d["W_adj"])
    spec = qdgmm.spec_m0(panel)
    m0 = qdgmm.fit(panel, spec)
    mod = qdgmm.QDModel(panel, spec)
    inf = mod.inference(m0["beta"])
    G = mod.G
    gg = inf["gg"]
    e = d["edges"].astype(int)
    cl = d["cl"]
    Adj = np.zeros((G, G), bool)
    a, b = cl[e[:, 0]], cl[e[:, 1]]
    Adj[a, b] = True
    Adj[b, a] = True
    np.fill_diagonal(Adj, False)
    Om = gg @ gg.T / G
    cross = sum(np.outer(gg[:, i], gg[:, j]) for i, j in zip(*np.nonzero(Adj))) / G
    Tm1 = mod.T - 1
    names = m0["names"]
    Gi = np.linalg.inv(inf["Gam"])

    def se(Omega):
        Vb = (Gi @ Omega @ Gi.T / G * G / (G - 1))[Tm1:, Tm1:]
        out = {n: float(np.sqrt(Vb[k, k])) if Vb[k, k] > 0 else float("nan") for k, n in enumerate(names)}
        v = Vb[0, 0] + Vb[1, 1] + 2 * Vb[0, 1]
        out["tau+kappa"] = float(np.sqrt(v)) if v > 0 else float("nan")
        return out
    corr = {}
    for k, n in enumerate(names):
        r = gg[Tm1 + k]
        num = sum(r[i] * r[j] for i, j in zip(*np.nonzero(Adj))) / Adj.sum()
        corr[n] = float(num / (r @ r / G))
    res = dict(G=G, adjacent_state_pairs=int(Adj.sum() // 2), mean_neighbours=float(Adj.sum(1).mean()),
               est={n: float(v) for n, v in zip(names, m0["beta"])}, se_cluster=se(Om), se_with_adjacent=se(Om + cross),
               adj_pair_score_corr=corr, git=env.git_sha())
    (HERE / "results.json").write_text(json.dumps(env.clean(res), indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
