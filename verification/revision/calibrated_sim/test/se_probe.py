import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common")); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import env  # noqa
import numpy as np, qdgmm, calib
from qdx import tsum
for sig in (float(a) for a in sys.argv[1:]):
    cal = calib.calibrated(sig); print("rho", cal["shock_rho"], cal["shock_adj_corr"], cal["adj_corr_data"], cal["nonadj_corr_data"])
    rows = []
    for r in range(6):
        pnl, _ = calib.simulate(cal, 777 + r)
        f = qdgmm.fit(pnl, qdgmm.spec_m0(pnl), starts=[[cal["tau"], cal["kappa"], cal["phi"]]])
        if f["ok"]:
            rows.append([f["se"]["tau"], f["se"]["kappa"], f["se"]["phi"], tsum(f["Vb"], 0, 1), f["est"]["tau"] + f["est"]["kappa"]])
    print(sig, "est_sigma", round(cal["sigma_sw_est"], 3), "median se tau,kappa,phi,tpk; tpk", np.round(np.median(rows, 0), 3), "recal", np.round(cal["recal_hist"], 3), flush=True)
