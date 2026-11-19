"""Check that QDX (national week effects, common phi) reproduces the paper's M0 estimates, sandwich SEs,
J (M0-OI, Hansen two-step) and Q_G / Q_GI; and report the score-form J_LM next to Hansen J."""
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import env  # noqa
import numpy as np
import application as A, qdgmm
from qdx import QDX, oi_extras, inv_extras

t0 = time.time()
d = A.load_data(); w = A.gi_weekly(A.GI_SOURCES["ferretti"])["w"]
panel, _ = A.make_panel(d, w, d["W_adj"])
sp = qdgmm.spec_m0(panel)
m0 = qdgmm.fit(panel, sp)
beta = np.array([m0["est"][n] for n in sp.names])
ex, exn = oi_extras(panel)
X = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, extras=ex, extra_names=exn)
x, ok, res = X.fit(starts=[beta])
inf = X.inference(x)
se_x = np.sqrt(np.diag(inf["Vb"]))
oi = qdgmm.fit_two_step(panel, qdgmm.spec_m0_oi(panel), m0)
jl = X.j_lm(x); jlc = X.j_lm(x, centered=True)
qg = X.q_full(x)
exi, exin = inv_extras(panel)
XI = QDX(panel, {"tau": panel["D"], "kappa": panel["Db"]}, extras=exi, extra_names=exin)
qgi = XI.q_full(x); jli = XI.j_lm(x)
gG = qdgmm.fit_global(panel); gGI = qdgmm.fit_global(panel, selection="inv")
out = dict(
    m0_paper_code=dict(est=m0["est"], se=m0["se"]),
    qdx=dict(est=dict(zip(X.names, x)), se=dict(zip(X.names, se_x)), ok=ok),
    max_abs_est_diff=float(np.max(np.abs(x - beta))),
    max_rel_se_diff=float(np.max(np.abs(se_x / np.array([m0["se"][n] for n in sp.names]) - 1))),
    hansen_J_two_step=oi["J"], J_LM_uncentered=jl["J_LM"], J_LM_centered=jlc["J_LM"], J_LM_p=jl["p"],
    Q_G_paper_code=gG["Q"], Q_G_qdx=qg, Q_GI_paper_code=gGI["Q"], Q_GI_qdx=qgi, J_LM_inv=jli["J_LM"],
    seconds=time.time() - t0)
out = env.clean(out)
Path(__file__).with_name("results.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
