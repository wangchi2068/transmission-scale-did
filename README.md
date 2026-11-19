# 传播尺度双重差分（QD-GMM）

论文、LaTeX 源码、实验代码与验证材料。按 *Biometrics* Regular Article 规范拆分为两份独立 PDF：
- 英文主稿 `main.pdf`：首页为标题/摘要（不计入限页，≤ 225 词），含正文、Back Matter、References 及 `endfloat` 排出的图表（未强行压缩，保留完整限定条件与独立的 Data Availability Statement）。
- 网络补充材料 `supporting_information.pdf`：Web Appendices A--F（含 Web Tables 1--15 与 Web Figures 1--4）及附带文献，与主稿通过 `xr` 双向交叉引用。
- 中文阅读稿 `paper_zh/main.pdf`（已与最新英文稿逐节同步）。

## 目录

```
QD-GMM_transmission_DiD_2026-11-13/
├── main.pdf, main.tex, supporting_information.pdf, supporting_information.tex
├── sections/, tables/, figures/
├── bibliography.bib, biom.cls, biom.bst, endrotfloat.sty, compile.sh
├── paper_zh/                  中文阅读稿（结构与英文逐节对应）
├── experiment/
│   ├── code/                  dgp、qdgmm、comparators、psi、simulate、data_prep、application、calibrated_sim_germany、make_figure_germany、make_tables_germany …
│   ├── theory/                numeric_checks.py（Web Appendix B 的 8 项数值核验）
│   ├── data/                  raw/germany（RKI SurvStat 水痘、OpenHolidays 16 州学校假期、德国 400 县 AGS 边界）与 processed 面板
│   ├── analysis/              results.json、germany_application_results.json、germany_calibrated_sim.json …
│   └── logs/                  运行日志
├── verification/
│   ├── root_selection/        根选择规则在 9 个独立设计上的验证（18000 次重复）
│   ├── same_week/             同周传播的匹配模拟（w0 = 0 对 w0 = 0.3，每臂 1000 次）
│   ├── revision/common/       qdx.py（分组周效应、满秩得分形式 J_LM / S3、全局求根）
│   ├── revision/population_roots/ 零传染设计总体多根验证
│   └── application/           wild_bootstrap.py（得分野 bootstrap-t 与非线性成对集群矩 bootstrap-t）
└── scripts/                   构建与质量检查脚本
```

## 编译与复现

主稿与网络补充材料交叉编译（`xelatex` + `bibtex`，依赖 XITS Math 与 Times New Roman 字体）：
```bash
xelatex -interaction=nonstopmode main.tex
xelatex -interaction=nonstopmode supporting_information.tex
bibtex supporting_information
xelatex -interaction=nonstopmode supporting_information.tex
xelatex -interaction=nonstopmode supporting_information.tex
xelatex -interaction=nonstopmode main.tex
bibtex main
xelatex -interaction=nonstopmode main.tex
xelatex -interaction=nonstopmode main.tex
xelatex -interaction=nonstopmode supporting_information.tex
```
中文稿：`cd paper_zh && xelatex main && bibtex main && xelatex main && xelatex main`。

德国水痘应用一键复现（Python 3.11，见 `experiment/code/README.md`）：
```bash
python experiment/code/data_prep.py                 # 校验 SHA-256 并构建 processed/panel.npz
python experiment/code/application.py               # 生成 germany_application_results.json
python experiment/code/calibrated_sim_germany.py    # 400 县 x 261 周校准模拟 (R=500 x 2)
python experiment/code/make_figure_germany.py       # 生成 Figure 3 (fig_03_application.png/.pdf)
python experiment/code/make_tables_germany.py       # 生成 Table 3 及 Web Appendix F 全部附表
```

## 主要结果

- **正溢出核心模拟（8 个情景，$R=8000$）**：M0 对 $\tau$ 偏倚 $-0.002$ 至 $0.001$，95\% 区间覆盖率 94.5\%--95.4\%，同时存在反馈实施、过度离散与小计数；常用比较量 M1--M4 与 M0$'$ 均存在偏倚。
- **无传染边界情景与同周压力测试**：零传染下总体矩存在第二根，M0-G/M0-GI 配合剖面 Stock--Wright $S$ 集控制选根风险；同周传播匹配实验（$w_0=0.3$）将周度 QD-GMM 的有效适用域明确界定在多周代际间隔病原体。
- **实证应用（德国 400 县、16 联邦州 2015--2019 年周度水痘与交错暑假，115,129 例）**：
  - 水痘 14 天代际间隔满足周度聚合（$w_0 = 2.6\times 10^{-6}$），且 16 州每年 5--6 周交错轮换暑假提供非吸收式重复开关暴露。
  - QD-GMM 定位唯一内部传染耦合根 $\hat\phi = 0.325$（95\% CI $0.206$--$0.445$），估计政策层面传播率下降 21.2\%（$\hat\tau+\hat\kappa = -0.238$，SE 0.058，95\% $t_{15}$ 区间 $[-0.362, -0.114]$），与欧洲历史水痘基准（22\%--31\%）一致，对应受处理周避免 2,440 例（SE 670）、全网络避免 2,995 例（SE 900）。
  - 全部有限簇与识别稳健区间一致不含零：Mammen 得分野 bootstrap-$t$ $[-0.357, -0.120]$、Webb 六点 $[-0.363, -0.113]$、非线性成对集群矩 bootstrap-$t$（999/999 收敛）$[-0.383, -0.094]$、删一州刀切 $[-0.374, -0.103]$（16 次删一范围 $[-0.268, -0.202]$）、Stock--Wright Bonferroni 95\% 投影 $[-0.492, -0.081]$；$S_1(F)$ 95\% 集 $[0.21, 2.33]$ 单一连通且不含零，$S_{3,\mathrm{inv}}(\chi^2)$ 与 $S_{3,\mathrm{oi}}(F)$ 均非空。
  - 事件时间剖面显示平缓前趋势（$p=0.741$）与零起始周效应（$\hat\tau_0 = -0.008$，排除即时暑期漏报），随后在潜伏期后（$e_1$--$e_5$）持续下降。
  - 真实数据直接验证定理 3：两步 EpiEstim 对数 $R$ 回归（M2）在 58.8\% 县级零计数下反号为 $+0.064$，而在 38 个 NUTS-2 区域聚合（零计数 2.9\%）下立即反转为 $-0.102$；$\log(Y+1)$ TWFE（M3）严重衰减至 $-0.102$。

## 投稿前仍需处理

作者、单位、通讯、致谢细节与资助信息仍是占位（`Apevon Science`、页码 `000 0000`、`[Funding statement to be completed before submission.]`），请在正式投递 *Biometrics* ScholarOne 前补齐真实作者与基金信息。
