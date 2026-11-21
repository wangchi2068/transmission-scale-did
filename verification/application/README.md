# verification/application —— 应用部分的辅助检验

> **Provenance note (E8)**: `wild_bootstrap.py` (`score_wild_bootstrap`, `pairs_moment_bootstrap`) and `STROBE_checklist.md` are active components for the German varicella implementation (`experiment/code/application.py`, $N=400, G=16, T=261$). The earlier scripts `tau_kappa.py`, `total_ci.py`, `timing.py`, `m0et_revised.py` / `m0et_revised.json` below document the **retired US COVID-19 mask-order prototype** ($N=3{,}091, G=49$) and are retained for historical audit only.

| 文件 | 内容 |
|---|---|
| `tau_kappa.py` / `tau_kappa.json` | 州级政策的 tau 与 kappa 共线性、tau+kappa 合并效应 |
| `total_ci.py`、`timing.py` | tau+kappa 区间、单次拟合耗时 |
| `m0et_revised.py` / `m0et_revised.json` | **M0-ET 事件时间分箱修正**（窗口前已暴露单位） |
| `STROBE_checklist.md` | 写作侧 STROBE 观察性研究核对清单 |

## M0-ET 修正（m0et_revised.json）

缺陷：`panel["ev"]` 由 `dgp.py` 的 `first = np.where(D.any(0), D.argmax(0), -1)` 在窗口内计算，窗口前已暴露单位（`Dall[:2] = 1`）的 ev 也是 0。原 `spec_m0_et` 的 tau_e0 因而把 435 个窗口前已暴露单位与 136 个窗口首周开始的单位混在一起，"in the week of adoption" 的解释不成立。435 个单位在窗口全部 26 周都是 D = 1。

修正：`qdgmm.spec_m0_et2(panel)`，配合 `application.make_panel(..., drop_presample=True)`。窗口前已暴露单位从估计样本中剔除（3,091 → 2,656 个单位，49 → 37 个州；有 12 个州只由这类县组成，另有 1 个州（州编号 35）只有 4 个县被剔除、其余保留）。tau_e0 只含窗口首周开始的 136 个单位，tau_e1..tau_e6 与 kappa 的分箱定义不变。
这类单位在窗口内始终处理，若在完整面板中单设 tau_pre 分箱，其分箱列在准差分后恒为零（与周效应截距共线），M0-ET 无法识别，所以只能剔除而不能另设分箱。剔除方式见 json 的 `data_replacement`：`rowwise` 为默认（邻居关系限于保留单位，按 `W_kept,kept @ D_kept` 计算 Db），`columnwise` 保留被剔除单位作为邻居（Db = W @ D，用全样本）。两种结果接近，见 `comparison_rowwise_vs_columnwise`。

旧 `spec_m0_et` 与 `application.py` 中的 `event_time` 未改动；新结果写入 `event_time2`，并在 `analyses.main.M0_ET2` 中保存完整拟合对象。

运行（仓库根目录，约 1 秒）：

```bash
python -I verification/application/m0et_revised.py experiment/code
```

`application.py` 的 quick 路径同时拟合两个版本（约 10–30 秒）：

```bash
cd experiment
python -I code/application.py --quick --out application/application_results.json
```

关键数字（panel.npz，G = 49，T = 26，S = 2，Ferretti GI）：

| 项目 | 旧 M0-ET（全样本） | 新 M0-ET2（剔除 pre-exposed） |
|---|---|---|
| tau_e0 系数（SE） | -0.232 (0.269) | -0.343 (0.341) |
| tau_e0 观测单元-周数 | 2,291（571 个单位） | 1,856（136 个单位） |
| tau_e1..tau_e6 | -0.416, -0.411, -0.503, -0.334, -0.577, -0.557 | -0.619, -0.672, -0.823, -0.680, -0.924, -0.923 |
| kappa / phi | -0.507 (0.315) / 0.771 (0.086) | -0.579 (0.406) / 0.757 (0.094) |
| 联合 Wald W（q = 7） | 21.135 | 17.546 |
| 联合 F（p） | 2.642（df 7, 42；p = 0.0234） | 2.089（df 7, 30；p = 0.0759） |
| chi-square p | 0.0036 | 0.0142 |
| e0–e5 平均（SE） | -0.412 (0.369) | -0.677 (0.521) |

说明：新样本只剩 37 个州聚类，F 的分母自由度由 42 变为 30。在同一 2,656 个单位的样本内，旧分箱（tau_e0 含 435 个窗口前已暴露单位）与 M0-ET2 的每个分箱本来就相同：这些单位在窗口内不切换状态，对分箱系数没有贡献。因此系数变化完全来自样本剔除，旧 tau_e0 −0.232 与 435 个常数处理单位的窗口前期混在一起不能解释为首周效应。所有数字见 `m0et_revised.json` 的 `runs`、`comparison_rowwise`、`comparison_within_kept_sample`。
