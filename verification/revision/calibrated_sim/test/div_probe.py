import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common")); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import env  # noqa
import numpy as np, calib
from qdx import groups_for_units
cal = calib.calibrate()
p = cal["panel"]; Y = p["Y"]; G = int(p["G"]); cl = p["cl"]
lam = p["L"] + cal["phi"] * p["Nb"] + p["omega"]
mu = cal["ea"][None, :] * cal["eg"][:, None] * np.exp(cal["tau"] * p["D"] + cal["kappa"] * p["Db"]) * lam
Yg = np.stack([np.bincount(cl, Y[t], G) for t in range(Y.shape[0])]); mg = np.stack([np.bincount(cl, mu[t], G) for t in range(Y.shape[0])])
div, reg = groups_for_units(cal["d"]["state"])
sdiv = np.array([div[cl == g][0] for g in range(G)])
lr = np.log(np.maximum(Yg, 0.5) / mg)
for lab, grp in (("division", sdiv), ("region", np.array([reg[cl == g][0] for g in range(G)]))):
    R = grp.max() + 1
    means = np.stack([lr[:, grp == r].mean(1) for r in range(R)], 1)   # T x R
    within = np.concatenate([(lr[:, grp == r] - means[:, [r]]).ravel() for r in range(R)])
    nr = np.array([(grp == r).sum() for r in range(R)])
    var_w = within.var() * G / (G - R)
    var_between = means.var() - np.mean(var_w / nr)
    print(lab, "within sd", np.sqrt(var_w), "between sd", np.sqrt(max(var_between, 0)))
# correlation of lr among adjacent states
import scipy.sparse as sps
W = cal["W"]; A_ = (W > 0).astype(float)
S = np.zeros((G, G))
coo = sps.coo_matrix(A_)
for i, j in zip(coo.row, coo.col):
    if cl[i] != cl[j]: S[cl[i], cl[j]] = 1
pairs = np.argwhere(np.triu(S) > 0)
c = [np.corrcoef(lr[:, a], lr[:, b])[0, 1] for a, b in pairs]
print("adjacent-state corr of lr over weeks: mean", np.mean(c), "n pairs", len(pairs))
nonadj = np.argwhere(np.triu(1 - S - np.eye(G)) > 0)
c2 = [np.corrcoef(lr[:, a], lr[:, b])[0, 1] for a, b in nonadj]
print("non-adjacent mean", np.mean(c2))
