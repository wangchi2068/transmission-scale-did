"""应用部分新增：州级政策下的 tau 与 kappa 共线性，以及总效应 tau+kappa。

背景：处理变量是州级口罩令，分析单元是县。腹地县的邻居全在同州，故 D_it = Db_it，
tau 与 kappa 只靠跨州县识别。

用法：python tau_kappa.py <repo>/experiment/code
"""
import json
import os
import sys
from pathlib import Path

if len(sys.argv) > 1:
    sys.path.insert(0, sys.argv[1])

import numpy as np

import application as A
import qdgmm

HERE = Path(__file__).parent


def main():
    d = A.load_data()
    A.verify(d)
    gi = A.gi_weekly(A.GI_SOURCES['ferretti']) if hasattr(A, 'GI_SOURCES') else {'w': (0.85, 0.15)}
    w = gi['w'] if isinstance(gi, dict) else gi
    panel, _ = A.make_panel(d, w, d['W_adj'])
    sp = qdgmm.spec_m0(panel, 0.5)
    fit = qdgmm.fit(panel, sp)
    mod = qdgmm.QDModel(panel, sp)
    beta = np.array([fit['est'][n] for n in ('tau', 'kappa', 'phi')])
    inf = mod.inference(beta)
    Vb = inf['Vb']
    tau, kap = fit['est']['tau'], fit['est']['kappa']
    total = tau + kap
    se_tot = float(np.sqrt(Vb[0, 0] + Vb[1, 1] + 2 * Vb[0, 1]))
    corr = float(Vb[0, 1] / np.sqrt(Vb[0, 0] * Vb[1, 1]))
    tc = qdgmm.tcrit(mod.G)
    D, Db = panel['D'], panel['Db']
    out = dict(tau=tau, kappa=kap, phi=fit['est']['phi'],
               se_tau=float(np.sqrt(Vb[0, 0])), se_kappa=float(np.sqrt(Vb[1, 1])),
               corr_tau_kappa=corr, tau_plus_kappa=float(total), se_total=se_tot,
               ci_total=[float(total - tc * se_tot), float(total + tc * se_tot)],
               t_total=float(total / se_tot), df=int(mod.G - 1),
               share_D_equals_Db=float(np.mean(np.abs(D - Db) < 1e-12)))
    print(json.dumps(out, ensure_ascii=False, indent=1))
    (HERE / 'tau_kappa.json').write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding='utf-8')


if __name__ == '__main__':
    main()
