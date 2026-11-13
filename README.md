# 传播尺度双重差分（QD-GMM）：定稿与补充验证

日期：2026-11-13。本文件夹收录论文定稿、LaTeX 源码，以及本轮新增的两组验证。

## 目录

```
QD-GMM_transmission_DiD_2026-11-13/
├── paper.pdf                  定稿（65 页，含 Web Appendix A–C）
├── paper.html                 网页阅读版
├── main.tex, sections/, tables/, figures/
├── bibliography.bib, biom.cls, biom.bst, endrotfloat.sty, compile.sh
├── experiment/                模拟、比较估计量与应用的全部代码与数据
│   ├── code/                  dgp、qdgmm、comparators、psi、simulate、application …
│   ├── data/                  raw（JHU、CDC 口罩令、Census 邻接与通勤、ACS）与 processed 面板
│   ├── analysis/              results.json、metrics_*.csv、实验报告与设计
│   └── logs/                  运行日志
└── verification/
    ├── root_selection/        根选择规则在独立设计上的验证（18000 次重复）
    └── application/           τ+κ 共线性、联合区间与计算耗时
```

实验复现：`cd experiment/code` 后按 `README.md` 运行；应用部分可用
`python application.py --quick` 做三秒级烟测（应与论文一致：τ = 0.040、κ = −0.506、φ = 0.772）。

编译：`cd` 到本目录后运行 `bash compile.sh`（xelatex + bibtex）。依赖 XITS Math 字体
（用于 Assumption A4′ 的 𝒞 符号），脚本已在 main.tex 中通过 `\setmathfont{XITS Math}[range={\mathcal,\mathscr}]` 指定。

## 本轮修改

**1. 正文与附录结构调整**

- 独立验证（Table 3）移到 Web Appendix C，正文只保留一段精简结论并引用 Web Table 7；
- “Departures from the protocol”压缩到三句，删去重复次数从 1000 提到 8000 的记述；
- 实证表格 Table 4 不设 τ+κ 合计列，总效应只在正文与摘要中报告主分析的联合区间。

**2. 字符与引用修复**

- 讨论段 `\tau` 曾被写成制表符，导致 PDF 显示 “au”；已修复，现为 τ 与 κ；
- Assumption A4′ 的 𝒞 通过 XITS Math 渲染；
- 两处内部文件名（`experiment_design.md` / `data_contract.md`）已从表注中移除，改为指向 Web Appendix；
- 表格中的代码变量名（`core_N500_T12_fb_spill` 等）全部替换为 `$N=500$, $T=12$, FB, spill` 这类排版形式。

## 新增证据

### 根选择规则的独立验证（verification/root_selection）

9 个新情景，种子基 `20270101`（主研习为 `20261112`），每情景 2000 次重复。

| φ₀ | M0 覆盖率 | M0-GI 覆盖率 | M0-GI 选边界比例 |
|---|---|---|---|
| 0 | 0.361 | 0.888 | 64% |
| 0.05 | 0.957–0.967 | 0.912–0.929 | 56%–79% |
| 0.10 | 0.930–0.961 | 0.890–0.940 | 1%–30% |
| 0.15 | 0.944 | 0.946 | 0% |

结论：规则在 φ₀=0 时有效，在 φ₀=0.05 时损失覆盖率，φ₀≥0.10 时与 M0 无差别。因此只作为多根时的补救。

复现：`python rsel_run.py --R 2000 --workers 8 --code <原仓库>/experiment/code`

### τ 与 κ 的共线性（verification/application）

- 腹地县（所有邻居同州）占 61.9%；`D` 与邻居政策均值完全相等的县周占 83.3%；
- τ 与 κ 的估计相关系数 +0.238；
- **τ+κ = −0.466**，联合协方差 95% 区间 **−1.091 到 0.159**（t = −1.50，48 自由度）；该总和是匹配州级政策的模型条件目标，重参数化本身不增加信息。
- 应用工作模型被过度识别检验拒绝（$J=17.873$, $p=0.0001$）；估计区间为常规州聚类结果。49 个州集群中约 24 个有县在窗口内切换，州规模不均，模拟未针对这种处理集群结构校准覆盖率。
- 同周传播会改变周频模型的条件均值，且可能改变准差分协方差；现有预设模拟不含同周传播，不能从中量化相应偏差或方差影响。
- 计算耗时：3091 县 × 26 周，M0 拟合 0.05 秒，M0-GI 全局搜索 0.73 秒。

复现：`total_ci.py`（联合区间）、`tau_kappa.py`（共线性结构）、`timing.py`（耗时），均需指向原实验仓库的 `experiment/code`。

## 复现与依赖

- 论文：`bash compile.sh`（xelatex + bibtex），需 XITS Math 字体（A4′ 的 𝒞）。
- 实验：Python 3.11，依赖见 `experiment/code/requirements.txt`。全量模拟约 3.5 小时（16 核）；
  应用部分约 3 秒。
- `experiment/calibration.json` 原包未附，已从 `analysis/results.json` 复原
  （gbar：T12_feedback 0.1016、T26_feedback 0.2118、T12_random 0.0056、T26_random 0.0231；
  depletion_pop：34391.9）。原目录也已补上同一文件。
- 原始数据快照在 `experiment/data/raw/`，校验和见 `experiment/data/checksums.tsv` 与 `data_contract.md`。

## 投稿前仍需处理

作者、单位、通讯、基金与利益冲突声明仍是模板占位（Apevon Science、页码 000 0000），需替换。
