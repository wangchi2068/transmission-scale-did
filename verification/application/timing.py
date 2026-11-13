"""实测 M0 与 M0-GI 在应用面板上的单次拟合耗时。"""
import sys
import time

sys.path.insert(0, sys.argv[1])

import application as A
import qdgmm

d = A.load_data()
gi = A.gi_weekly(A.GI_SOURCES['ferretti'])
w = gi['w'] if isinstance(gi, dict) else gi
panel, _ = A.make_panel(d, w, d['W_adj'])
sp = qdgmm.spec_m0(panel, 0.5)

t0 = time.time()
fit = qdgmm.fit(panel, sp)
t_m0 = time.time() - t0

t0 = time.time()
qdgmm.fit_global(panel, 0.5, selection='inv')
t_gi = time.time() - t0

T, N = panel['Y'].shape
print('units=%d weeks=%d  M0=%.2fs  M0-GI=%.2fs' % (N, T, t_m0, t_gi))
