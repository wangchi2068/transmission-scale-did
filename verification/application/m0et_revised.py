"""M0-ET 事件时间分箱修正：窗口前已暴露单位（Dall[:2] = 1，窗口内 ev = 0）单独处理。

原 spec_m0_et 把两类 ev = 0 的单位并入 tau_e0："窗口首周开始" (136) 与 "窗口前已暴露" (435)。
修正版 spec_m0_et2 把窗口前已暴露单位从估计样本中剔除，tau_e0 只包含首周开始的 136 个单位，
tau_e1..e6 与 kappa 的定义不变。旧的 spec_m0_et / event_time 保留，用于对照。

Dall[:2] 只通过窗口内的 L、Nb、Db 进入估计（Lag L 用 Yall 的两行 pre-sample Y；Db = W@D 用窗口内的 D）。
被剔除单位作为邻居仍出现在 W 中，两种替换方式：
  rowwise     被剔除单位的 Y 列为 0：贡献 W@L 时 L = 0，不贡献政策；等于把 D_all 的 pre-sample 两行清零。
  columnwise  被剔除单位整列保留：只作为邻居贡献 W@L 与 Db = W@D。
两种方式都用 panel.npz 的 Dall / Yall 重建，不改变窗口内任何取值。

用法（在仓库根目录）：
    python -I verification/application/m0et_revised.py experiment/code

输出 verification/application/m0et_revised.json。各版本都走 application.py 的路径
（make_panel + qdgmm.fit + event_time_block + unit_counts）。
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CODE = sys.argv[1] if len(sys.argv) > 1 else str(REPO / 'experiment' / 'code')
sys.path.insert(0, CODE)

import application as A  # noqa: E402
import qdgmm  # noqa: E402

OUT = HERE / 'm0et_revised.json'
COMMAND = "python -I verification/application/m0et_revised.py experiment/code"


def columnwise_panel(d, keep, w):
    """Estimation panel for units `keep`, with the dropped units kept as neighbours only (L, Nb, Db)."""
    Wf = A.row_normalise(d['W_adj'])
    S = A.S
    T = d['Y'].shape[0] - S
    L = sum(w[s - 1] * d['Y'][S - s:S - s + T] for s in range(1, len(w) + 1))     # full sample
    Db = np.asarray((Wf @ d['D'][S:].T).T)
    Nb = np.asarray((Wf @ L.T).T)
    pnl = dict(A.make_panel(d, w, d['W_adj'], units=keep)[0])
    pnl['L'], pnl['Nb'], pnl['Db'] = L[:, keep], Nb[:, keep], Db[:, keep]
    return pnl


def blocks(ev, spec):
    return dict(ok=ev.get('ok'), bins=ev.get('bins'), fields=spec.names,
                est={k: ev['est'].get(k) for k in spec.names},
                se={k: ev['se'].get(k) for k in spec.names},
                n_cells={k: ev['n_cells'].get(k) for k in ev.get('bins', [])},
                joint_wald=ev.get('joint_wald'),
                average={k: v for k, v in ev.items() if k.startswith('average_')} or None)


def run(panel, name, G):
    spec = qdgmm.spec_m0_et(panel) if name == 'M0_ET_stock' else qdgmm.spec_m0_et2(panel)
    tc = qdgmm.tcrit(G)
    res = A.timed(qdgmm.fit, panel, spec)
    ev = A.event_time_block(res, spec, tc)
    print('[m0et_revised] %-26s ok=%s in %.1fs' % (name, res.get('ok'), res.get('time_sec', 0)), flush=True)
    return dict(table=blocks(ev, spec), ok=bool(res.get('ok')), cond=res.get('cond'),
                time_sec=res.get('time_sec'), n_units=int(panel['Y'].shape[1]), G=int(panel['G']),
                tcrit=float(tc))


def main():
    t0 = time.time()
    d = A.load_data(A.PANEL)
    A.verify(d)
    w = A.gi_weekly(A.GI_SOURCES['ferretti'])['w']
    G = len(d['states'])
    S = A.S
    pre = d['D'][:S].any(0)
    keep = np.where(~pre)[0]
    panel_full, _ = A.make_panel(d, w, d['W_adj'])
    panel_row, _ = A.make_panel(d, w, d['W_adj'], drop_presample=True)       # rowwise: W columns subset
    panel_col = columnwise_panel(d, keep, w)                                 # columnwise: full W columns
    # check: zeroing the pre-sample rows of D_all reproduces the rowwise panel
    panel_zero, _ = A.make_panel(d, w, d['W_adj'], D_all=d['D'] * (~pre)[None])
    for k in ('Y', 'D', 'Db', 'L', 'Nb', 'Ypre', 'Dpre'):
        panel_zero[k] = panel_zero[k][:, keep]
    panel_zero['omega'] = panel_zero['omega'][keep]
    panel_zero['ev'] = panel_zero['ev'][keep]
    panel_zero['cl'] = np.unique(panel_zero['cl'][keep], return_inverse=True)[1]
    panel_zero['G'] = int(panel_zero['cl'].max()) + 1
    identical = bool(np.allclose(panel_zero['L'], panel_row['L'])
                     and np.allclose(panel_zero['Db'], panel_row['Db'])
                     and np.allclose(panel_zero['Nb'], panel_row['Nb']))

    runs = dict(
        old_full_sample=run(panel_full, 'M0_ET_stock', G),
        new_rowwise=run(panel_row, 'M0_ET2', panel_row['G']),
        new_columnwise=run(panel_col, 'M0_ET2', panel_col['G']),
        check_Dpre_zeroed_rowwise=run(panel_zero, 'M0_ET2', panel_zero['G']),
        old_bins_within_kept_panel=run(panel_row, 'M0_ET_stock', panel_row['G']),
    )

    D = panel_full['D']
    switch = (D[1:] == 1) & (D[:-1] == 0)
    counts = A.unit_counts(d, panel_full)
    counts.update(n_units_m0et2=int(len(keep)), n_states_m0et2=int(panel_row['G']),
                  switch_units_by_week={int(t + 1): int(switch[t].sum())
                                        for t in range(switch.shape[0]) if switch[t].any()})
    states_drop = sorted({int(x) for x in d['state'][pre]} - {int(x) for x in d['state'][keep]})
    counts['states_fully_dropped'] = states_drop
    counts['states_partially_present'] = sorted({int(x) for x in d['state'][pre]} - set(states_drop))

    def cmp(a, b):
        o, n = runs[a]['table'], runs[b]['table']
        return dict(F_old=o['joint_wald']['F'], F_new=n['joint_wald']['F'],
                    p_F_old=o['joint_wald']['p_F'], p_F_new=n['joint_wald']['p_F'],
                    W_old=o['joint_wald']['W'], W_new=n['joint_wald']['W'],
                    df_old=[o['joint_wald']['df1'], o['joint_wald']['df2']],
                    df_new=[n['joint_wald']['df1'], n['joint_wald']['df2']],
                    p_chi2_old=o['joint_wald']['p_chi2'], p_chi2_new=n['joint_wald']['p_chi2'],
                    tau_e0_cells_old=o['n_cells'].get('tau_e0'), tau_e0_cells_new=n['n_cells'].get('tau_e0'),
                    tau_e0_est_old=o['est'].get('tau_e0'), tau_e0_est_new=n['est'].get('tau_e0'),
                    tau_e0_se_old=o['se'].get('tau_e0'), tau_e0_se_new=n['se'].get('tau_e0'),
                    other_bins_est_old={k: o['est'].get(k) for k in o['bins'] if k != 'tau_e0'},
                    other_bins_est_new={k: n['est'].get(k) for k in n['bins'] if k != 'tau_e0'},
                    other_bins_delta={k: n['est'][k] - o['est'][k] for k in o['bins'] if k != 'tau_e0'},
                    kappa_phi_old={k: o['est'].get(k) for k in ('kappa', 'phi')},
                    kappa_phi_new={k: n['est'].get(k) for k in ('kappa', 'phi')})

    out = dict(
        created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        command=COMMAND, runtime_sec=round(time.time() - t0, 1),
        data=dict(panel='experiment/data/processed/panel.npz', n_units=int(D.shape[1]), G=int(G),
                  T=int(panel_full['Y'].shape[0]), S=int(S)),
        spec=dict(
            old='spec_m0_et: tau_e0 = {D_it = 1, ev_i = 0} pools 435 pre-exposed and 136 first-week switchers',
            new='spec_m0_et2: pre-exposed units (Dall[:2] = 1) removed; tau_e0 = {D_it = 1, ev_i = 0} among '
                'the remaining 136 units; tau_e1..e6 and kappa unchanged',
            pre_exposed_path='all 435 pre-exposed units are treated in every window week, so a tau_pre bin inside '
                             'the full panel has an identically zero Jacobian column (collinear with the '
                             'quasi-difference intercept) and M0-ET does not identify; removal is the '
                             'identified option',
        ),
        sample_counts=counts,
        data_replacement=dict(
            rowwise='removed units are neighbours with L = 0 (their pre-sample Y is outside the kept columns) and '
                    'with their pre-sample policy dropped from the neighbour sums: Db here is A.make_panel\'s '
                    'W_subset @ D_subset (W restricted to the kept columns), i.e. W_kept,kept @ D_kept. '
                    'Zeroing the pre-sample rows of D_all keeps W_kept,kept @ D_pre instead and differs '
                    '(check_Dpre_zeroed_rowwise)',
            columnwise='removed units are neighbours only: contribute W@L (L = 0 for them) and W@D, not their own Y',
            Dpre_zeroed_equals_rowwise=identical,
        ),
        runs=runs,
        comparison_rowwise=cmp('old_full_sample', 'new_rowwise'),
        comparison_within_kept_sample=cmp('old_bins_within_kept_panel', 'new_rowwise'),
        comparison_rowwise_vs_columnwise=cmp('old_full_sample', 'new_columnwise'),
    )
    OUT.write_text(json.dumps(A._clean(out), ensure_ascii=False, indent=1), encoding='utf-8')
    print('[m0et_revised] wrote %s in %.1f s' % (OUT, time.time() - t0))
    return 0


if __name__ == '__main__':
    sys.exit(main())
