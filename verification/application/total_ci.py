"""对应用中每个方法/敏感性分析，用联合协方差计算 tau+kappa 的区间。

重跑每个拟合以取回 Vb（application.py 未保存协方差矩阵）。
用法：python total_ci.py <repo>/experiment/code
"""
import json
import sys

sys.path.insert(0, sys.argv[1])

import numpy as np

import application as A
import qdgmm

OUT = []


def record(label, panel, spec):
    fit = qdgmm.fit(panel, spec)
    names = list(spec.names)
    if not fit.get('ok') or 'tau' not in names or 'kappa' not in names:
        return
    mod = qdgmm.QDModel(panel, spec)
    beta = np.array([fit['est'][n] for n in names])
    inf = mod.inference(beta)
    if inf is None:
        return
    Vb = inf['Vb']
    i, j = names.index('tau'), names.index('kappa')
    tot = fit['est']['tau'] + fit['est']['kappa']
    se = float(np.sqrt(Vb[i, i] + Vb[j, j] + 2 * Vb[i, j]))
    tc = qdgmm.tcrit(mod.G)
    OUT.append(dict(label=label, tau=fit['est']['tau'], kappa=fit['est']['kappa'],
                    total=float(tot), se=se, lo=float(tot - tc * se), hi=float(tot + tc * se)))


def main():
    d = A.load_data()
    A.verify(d)
    gi = A.gi_weekly(A.GI_SOURCES['ferretti'])
    w = gi['w'] if isinstance(gi, dict) else gi
    panel, _ = A.make_panel(d, w, d['W_adj'])
    record('main', panel, qdgmm.spec_m0(panel, 0.5))
    print(json.dumps(OUT, ensure_ascii=False, indent=1))


if __name__ == '__main__':
    main()
