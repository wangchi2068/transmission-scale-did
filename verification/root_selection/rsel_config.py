"""根选择规则验证：新设计配置。

目标：在现有 25 个情景之外，检验 M0-GI（事后设计的根选择规则）在
  (a) φ0 靠近零但不等于零（0.05 / 0.10 / 0.15），
  (b) 不同参考值 φ̄（0.25 / 0.5 / 1.0），
  (c) 不同簇数 G（20 / 100），
下的表现。必须使用与开发阶段完全独立的种子。

用法：python rsel_run.py --R 200 --workers 8
"""
from __future__ import annotations

import json
from pathlib import Path

# 开发阶段的种子基是 20261112；这里用独立的基
BASE_SEED = 20270101
PSI_SEED_OFFSET = 950_000_000
ALPHA_SCALE = 0.25
TAU0 = -0.2231435513142097  # log(0.8)
GBAR = {  # 从 analysis/results.json 的工况参数读出（T, adopt）-> gbar
    (12, "feedback"): 0.1016,
    (26, "feedback"): 0.2118,
    (12, "random"): 0.0056,
    (26, "random"): 0.0231,
}

# 新设计的格子：(phi0, phi_bar, G_scale, adopt, N, T)
# G_scale 为 1.0 表示沿用默认的 50 个簇；0.4 -> 20，2.0 -> 100
CELLS = [
    (0.05, 0.50, 1.0, "feedback", 500, 12),
    (0.10, 0.50, 1.0, "feedback", 500, 12),
    (0.15, 0.50, 1.0, "feedback", 500, 12),
    (0.05, 0.25, 1.0, "feedback", 500, 12),
    (0.05, 1.00, 1.0, "feedback", 500, 12),
    (0.10, 0.50, 0.4, "feedback", 2000, 12),
    (0.10, 0.50, 2.0, "feedback", 2000, 12),
    (0.10, 0.50, 1.0, "random", 500, 12),
    (0.00, 0.50, 1.0, "feedback", 500, 12),   # 与既有场景同设定，仅换种子：内部对照
]


def build(calib_gbar=None):
    from dgp import DGPParams
    out = {}
    for i, (phi, phibar, gs, adopt, N, T) in enumerate(CELLS):
        sid = 101 + i
        key = (T, adopt)
        out[sid] = dict(sid=sid, phi0=phi, phi_bar=phibar, g_scale=gs, adopt=adopt, N=N, T=T,
                        params=DGPParams(N=N, T=T, adopt=adopt, tau=TAU0, kappa=-0.1 * (phi > 0),
                                         phi=phi, k=10.0, gbar=GBAR[key],
                                         alpha_state_sd=0.2 * ALPHA_SCALE,
                                         alpha_unit_sd=0.2 * ALPHA_SCALE,
                                         alpha_pi=0.3 * ALPHA_SCALE))
    return out


def seed_for(sid, rep):
    return BASE_SEED + sid * 100_000 + rep
