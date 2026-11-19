# 传播尺度双重差分（QD-GMM）

论文、LaTeX 源码、实验代码与验证材料。按 *Biometrics* 结构拆分并保留完整边界条件披露：
- 英文主稿 `main.pdf`（共 38 页；按 *Biometrics* 计页规则第 1--5 节 + Back Matter + References 为 31 页，未强行压缩至 25 页）：首页为标题/摘要（225 词，≤ 225），含第 1--5 节、Back Matter、References 及 `endfloat` 排出的图表，完整保留方案披露（Protocol & Departures）与局限/边界条件。
- 网络补充材料 `supporting_information.pdf`（共 45 页）：Web Appendices A--F（含 Web Tables 1--16 与 Web Figures 1--4）及附带文献，与主稿通过 `xr` 双向交叉引用。
- 中文阅读稿 `paper_zh/main.pdf`（共 26 页，已与最新英文稿逐节同步）。

## 目录

```
QD-GMM_transmission_DiD_2026-11-13/
├── .gitattributes             锁定 experiment/data/raw/** 为 binary，防止跨平台 CRLF 改变 SHA-256
├── main.pdf, main.tex, supporting_information.pdf, supporting_information.tex
├── sections/, tables/, figures/
├── bibliography.bib, biom.cls, biom.bst, endrotfloat.sty, compile.sh
├── paper_zh/                  中文阅读稿（结构与英文逐节对应）
├── experiment/
│   ├── code/                  dgp、qdgmm、comparators、psi、simulate、data_prep、application、calibrated_sim_germany、make_figure_germany、make_tables_germany …
│   ├── theory/                numeric_checks.py（Web Appendix B 的 8 项数值核验）
│   ├── data/                  fetch_germany.py（RKI SurvStat@RKI 2.0 / OpenHolidays / AGS GeoJSON 下载与校验脚本）、checksums.tsv、data_contract.md、raw/germany 与 processed 面板
│   ├── analysis/              results.json、germany_application_results.json、germany_calibrated_sim.json …
│   └── logs/                  运行日志
├── verification/
│   ├── root_selection/        根选择规则在 9 个独立设计上的验证（18000 次重复）
│   ├── same_week/             同周传播的匹配模拟（w0 = 0 对 w0 = 0.3，每臂 1000 次）
│   ├── revision/common/       qdx.py（分组周效应、满秩得分形式 J_LM / S3、全局求根）
│   ├── revision/population_roots/ 零传染设计总体多根验证
│   └── application/           wild_bootstrap.py（得分野 bootstrap-t 与非线性成对集群矩 bootstrap-t）
├── response/                  各轮审稿意见回复（r1r2、r3、r4、r5）
└── scripts/                   scan_control_chars.py、qa_checks.py、section_pages.py
```

## 编译与复现

主稿、网络补充材料与中文伴随稿编译（`xelatex` + `bibtex`）：
```bash
bash compile.sh
# 随后运行质量与分节计页核查（compile.sh 已自动调用）：
python scripts/scan_control_chars.py
python scripts/qa_checks.py
python scripts/section_pages.py
```

德国水痘应用一键复现（Python 3.11，见 `experiment/code/README.md` 与 `experiment/data/data_contract.md`）：
```bash
python experiment/data/fetch_germany.py --check     # 校验 LF 原始文件与 SHA-256（3/3 OK）
python experiment/code/data_prep.py                 # 校验 SHA-256（0 个插补县）并生成 processed/panel.npz
python experiment/code/application.py               # 生成 germany_application_results.json
python experiment/code/calibrated_sim_germany.py    # 400 县 x 261 周校准模拟 (R=500 x 2，含 WCR-score)
python experiment/code/make_figure_germany.py       # 生成 Figure 3 (fig_03_application.png/.pdf)
python experiment/code/make_tables_germany.py       # 生成 Table 3 及 Web Appendix F 全部附表 (含 Web Table 16)
```

## 主要结果

- **正溢出核心模拟（8 个情景，$R=8000$）**：M0 对 $\tau$ 偏倚 $-0.002$ 至 $0.001$，95\% 区间覆盖率 94.5\%--95.4\%，同时存在反馈实施、过度离散与小计数；常用比较量 M1--M4 与 M0$'$ 均存在偏倚。
- **无传染边界情景与同周压力测试**：零传染下总体矩存在第二根，M0-G/M0-GI（第 3 节明确披露为事后规则，参见 *Departures from the protocol*）配合剖面 Stock--Wright $S$ 集揭示边界与多根风险；同周传播匹配实验（$w_0=0.3$）将周度 QD-GMM 的有效适用域界定在多周代际间隔病原体。
- **实证压力测试（德国 400 县、16 联邦州 2015--2019 年周度水痘与交错暑假，115,129 例，事后分析）**：
  - 真实监测数据验证定理 3：两步 EpiEstim 对数 $R$（M2）在 58.8\% 县级零计数下反号为 $+0.064$，而在 38 个 NUTS-2 区域聚合（零计数 2.9\%）下立即反转为 $-0.102$；$\log(Y+1)$ TWFE（M3）严重衰减至 $-0.102$。
  - 基准工作模型 M0 给出政策层面负向对比 $\hat\tau+\hat\kappa = -0.239$（SE $0.0584$，$t_{15}$ CI $[-0.364, -0.115]$；内部耦合根 $\hat\phi = 0.326$），有限簇校准下 $t_{15}$ 覆盖率为 0.902（SE 比值 0.855），由限制得分野 bootstrap WCR-score 恢复至 0.934--0.946（观测 CI $[-0.419, -0.139]$ Webb6，$[-0.429, -0.119]$ Rademacher）；邻近州交叉项协方差（Remark 4）呈非正定（$\min\mathrm{eig} = -0.020$）。
  - 同时，内置诊断明确阻止将其归因为无混杂的因果减毒：滞后工具过度识别检验显著拒绝（$J_{LM,\mathrm{oi}}=9.12$，校准 $p=0.004$），$S_{3,\mathrm{oi}}$ 微弱偏向边界（`M0-G` 落于 $\phi=0$），$S_{3,\mathrm{oi}}(\chi^2)$ 全为空集且 $F_{3,14}$ 集中部分铺满 $[0,3]$；总政策尺度事件时间 $\hat\tau_0+c_0\hat\kappa$ 在放假起始周即跌至 $-0.208$（$p=0.049$）；非暑期短期假期（复活节、秋假、冬假、圣灵降临节）均现滞后 0 负降与滞后 2 正回补，同期 $D^{(0)}$ 21.5\% 报告率修正可完全消除滞后 1 效应（$+0.007$）；剔除不轮换的巴符与巴伐利亚后滞后 1 CI 含零（$-0.211, [-0.431, 0.008]$）、滞后 2 为 $+0.033$，四大区 $\times$ 周效应后亦含零（$-0.169, [-0.376, 0.038]$）。

## 投稿前仍需处理

作者、单位、通讯、致谢细节与资助信息仍是占位（`Apevon Science`、页码 `000 0000`、`[Funding statement to be completed before submission.]`），请在正式投递前补齐真实作者与基金信息。
