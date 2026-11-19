# 传播尺度双重差分（QD-GMM）

论文、LaTeX 源码、实验代码与验证材料。按 *Biometrics* 结构拆分，以方法理论、25 组预设模拟、适用包络界定与稀疏真实面板说明性实现为主线，完整德国水痘敏感性与诊断归档于网络支持材料：
- 英文主稿 `main.pdf`（共 35 页；*Biometrics* 正文+致谢+参考文献计 27 页，含 `endfloat` 排出的 Figure 1--3 与 Table 1--4）：首页为标题/摘要（201 词，≤ 225），含完整方案与偏离披露（Protocol & Departures）、适用包络指南表（Table 3）以及德国稀疏水痘面板说明性实现（§4，Table 4）。
- 网络补充材料 `supporting_information.pdf`（共 46 页）：Web Appendices A--F（含 Web Tables 1--16 与 Web Figures 1--5）及附带文献，与主稿通过 `xr` 双向交叉引用。
- 中文阅读稿 `paper_zh/main.pdf`（共 24 页，已与最新英文稿逐节同步）。

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

德国水痘面板与图表复现（Python 3.11，见 `experiment/code/README.md` 与 `experiment/data/data_contract.md`）：
```bash
python experiment/data/fetch_germany.py --check     # 校验 LF 原始文件与 SHA-256（3/3 OK）
python experiment/code/data_prep.py                 # 校验 SHA-256（0 个插补县）并生成 processed/panel.npz
python experiment/code/application.py               # 生成 germany_application_results.json
python experiment/code/calibrated_sim_germany.py    # 400 县 x 261 周校准模拟 (R=500 x 2，含 WCR-score)
python experiment/code/make_figure_germany.py       # 生成 Figure 3 (fig_03_application) 与 Web Figure 5 (fig_web_05_diagnostics)
python experiment/code/make_tables_germany.py       # 生成 Table 4 及 Web Appendix F 全部附表 (含 Web Table 16)
```

## 主要结果

- **正溢出核心模拟（8 个情景，$R=8000$）**：M0 对 $\tau$ 偏倚 $-0.002$ 至 $0.001$，95\% 区间覆盖率 94.5\%--95.4\%，同时存在反馈实施、过度离散与小计数；常用比较量 M1--M4 与 M0$'$ 均存在偏倚。
- **适用边界与方法指南（Table 3 Operating Envelope）**：零传染下总体矩存在第二根，M0-G/M0-GI（第 3 节 *Departures from the protocol* 显式披露为修订事后规则）配合剖面 Stock--Wright $S_1$ 集揭示多根风险；同周传播匹配实验（$w_0=0.3$）将周度 QD-GMM 的有效适用域界定在多周代际间隔病原体；非暑期假期安慰剂直接检出对假设 A7 可分离报告约束（$c_{it}=c_i c_t$）的偏离。
- **稀疏监测计数说明性实现（德国 400 县、16 联邦州 2015--2019 年水痘面板，115,129 例，事后分析）**：
  - **计算与小计数比较量**：无需剔除零计数县周（58.8\% 为零），M0 稳定收敛至唯一内部耦合根 $\hat\phi=0.326$（95\% CI $[0.205,0.446]$），工作模型对比 $\hat\tau+\hat\kappa = -0.239$（$[-0.364,-0.115]$）；真实监测数据直观展示了定理 3 所刻画的小计数符号反转现象——两步 EpiEstim（M2）在县级反号为 $+0.064$，而在 38 个 NUTS-2 区域聚合（零计数 2.9\%）下立即反转为 $-0.102$（$[-0.150,-0.055]$），$\log(Y+1)$ TWFE（M3）衰减至 $-0.102$。
  - **有限簇设计校准**：在按该 400 县 $\times$ 16 州结构校准的模拟（$R=500\times 2$）中，常规 $t_{15}$ 对 $\tau+\kappa$ 覆盖率为 0.902（SE 比值 0.855），而限制得分野 bootstrap（WCR-score）显著改善并恢复覆盖率至 0.934--0.946（观测 CI $[-0.419,-0.139]$ Webb6，$[-0.429,-0.119]$ Rademacher）。
  - **完整透明诊断（Web Appendix E--F）**：鉴于公开数据按报告周（*Meldewoche*）记录，正文不作暑假因果政策归因；Web Appendices E--F 完整汇报了邻近州非正定协方差（$\min\mathrm{eig}=-0.020$）、$J_{LM,\mathrm{oi}}=9.12$（校准 $p=0.004$）、$S_{3,\mathrm{oi}}(\chi^2)$ 空集、总尺度起始周 $-0.208$（Web Figure 5a）、非暑期假期先降后弹安慰剂（Web Figure 5b）、21.5\% 起始周归零校准敏感性、BW/BY 联合剔除与四大区 $\times$ 周控制结果。

## 投稿前仍需处理

作者、单位、通讯、致谢细节与资助信息仍是占位（`Apevon Science`、页码 `000 0000`、`[Funding statement to be completed before submission.]`），请在正式投递前补齐真实作者与基金信息。
