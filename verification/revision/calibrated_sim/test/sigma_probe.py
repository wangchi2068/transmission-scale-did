import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common")); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import env  # noqa
import numpy as np, qdgmm, calib, application as A
cal = calib.calibrate()
p = cal["panel"]; Y = p["Y"]
lam = p["L"] + cal["phi"] * p["Nb"] + p["omega"]
mu = cal["ea"][None, :] * cal["eg"][:, None] * np.exp(cal["tau"] * p["D"] + cal["kappa"] * p["Db"]) * lam
G = int(p["G"]); cl = p["cl"]
Yg = np.stack([np.bincount(cl, Y[t], G) for t in range(Y.shape[0])]); mg = np.stack([np.bincount(cl, mu[t], G) for t in range(Y.shape[0])])
lr = np.log(Yg / mg)
# NB noise contribution approx: var(log sum) ~ (sum mu + sum mu^2/k)/ (sum mu)^2
mu2 = np.stack([np.bincount(cl, mu[t] ** 2, G) for t in range(Y.shape[0])])
noise = (mg + mu2 / cal["k"]) / mg ** 2
print("var log ratio", lr.var(), "mean noise", noise.mean(), "sigma est", np.sqrt(max(lr.var() - noise.mean(), 0)))
print("sd by week", np.round(lr.std(0)[:5], 3), "acf1", np.corrcoef(lr[1:].ravel(), lr[:-1].ravel())[0, 1])
ok = Yg > 0
print("zero state-weeks", (~ok).sum(), "of", ok.size)
x = lr[ok]
print("var", x.var(), "sigma", np.sqrt(max(x.var() - noise[ok].mean(), 0)))
q = np.quantile(x, [0.05, 0.25, 0.5, 0.75, 0.95]); print("quantiles", np.round(q, 3))
iqr_sigma = (q[3] - q[1]) / 1.349; print("robust sigma", iqr_sigma)
both = ok[1:] & ok[:-1]
print("acf1", np.corrcoef(lr[1:][both], lr[:-1][both])[0, 1])
