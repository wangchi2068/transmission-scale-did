# 传播尺度双重差分（QD-GMM）

论文、LaTeX 源码、实验代码与验证材料。当前仓库为未按期刊 25 页正文限制压缩的完整工作稿（正式投稿版将根据目标期刊要求另行排版），以方法理论、25 组预设模拟、适用包络界定与稀疏真实面板说明性实现为主线，完整德国水痘敏感性与诊断透明归档于网络补充材料：
- 英文主稿 `main.pdf`：完整工作稿（含摘要/关键词、Sections 1--5 + Back Matter + References 及 `endfloat` 排出的 Figure 1--3 与 Table 1--4），含研究方案与事后层级区分（Protocol Audit）、适用包络指南表（Table 3）以及德国稀疏水痘面板说明性测试平台（§4，Table 4）。
- 网络补充材料 `supporting_information.pdf`：Web Appendices A--F（含 Web Tables 1--18 与 Web Figures 1--5）及附带文献，与主稿通过 `xr` 双向交叉引用。
- 中文阅读稿 `paper_zh/main.pdf`：已与最新英文稿逐节同步。

## 目录

```
QD-GMM_transmission_DiD_2026-11-13/
├── .gitattributes             锁定 experiment/data/raw/** 为 binary，防止跨平台 CRLF 改变 SHA-256
├── main.pdf, main.tex, supporting_information.pdf, supporting_information.tex
├── sections/, tables/, figures/
├── bibliography.bib, biom.cls, biom.bst, endrotfloat.sty, compile.sh
├── paper_zh/                  中文阅读稿（结构与英文逐节对应）
├── experiment/
│   ├── code/                  dgp、qdgmm、comparators、psi、simulate、data_prep、application、calibrated_sim_germany、diagnostics_core、make_tables、make_figure_germany、make_tables_germany、make_verification_tables …
│   ├── theory/                numeric_checks.py 与 web_appendix_b_mc.py（Web Appendix B 数值与 MC 核验）
│   ├── data/                  fetch_germany.py（RKI SurvStat@RKI 2.0 / OpenHolidays ODbL-1.0 / BKG VG250 dl-de/by-2-0 下载与校验脚本）、checksums.tsv、data_contract.md、raw/germany 与 processed 面板
│   ├── analysis/              results.json、core_diagnostics.json、germany_application_results.json、germany_calibrated_sim.json …
│   └── logs/                  运行日志
├── verification/
│   ├── root_selection/        根选择规则 9 组验证（18000 次）+ cluster_grid.py（纯 G=10..100 与偏度 × T 网格）
│   ├── same_week/             同周传播匹配模拟（w0=0 对 0.3）+ grid_w0.py（w0∈[0,0.30] 分级与处理比例匹配网格）
│   ├── non_nested_sim/        run.py（非同构日度 Gamma 分支 + 随机报告延迟 + SEIR 耗竭模拟）
│   ├── weak_prior_m2/         run.py（在主协议 R=8000 同种子上重新评估弱先验与超弱先验 EpiEstim M2）
│   ├── revision/              common/qdx.py、population_roots/（N=20k×40 + N=50k×15 共 540 次零传染总体第二根核查）、README_PROVENANCE.md
│   └── application/           wild_bootstrap.py（得分野 bootstrap-t 与非线性成对集群矩 bootstrap-t）
├── response/                  各轮审稿意见回复（r1r2 至 r9）
└── scripts/                   scan_control_chars.py、qa_checks.py、section_pages.py、rebuild_manifest.py
```

## 编译与复现

主稿、网络补充材料与中文伴随稿编译（`xelatex` + `bibtex`；图表脚本 `experiment/code/style_kit.py` 严格校验并要求系统安装 **Times New Roman** 常规字体，缺字体时会直接报错而不静默回退）：
```bash
bash compile.sh
# 随后运行质量与分节计页核查（compile.sh 已自动调用）：
python scripts/scan_control_chars.py
python scripts/qa_checks.py
python scripts/section_pages.py
```

德国水痘面板与图表复现（Python 3.11，见 `experiment/code/README.md` 与 `experiment/data/data_contract.md`；存档 SurvStat 数据提取日期 **2026-11-19**，官方 BKG VG250 县界按 `dl-de/by-2-0` 与 OpenHolidays 按 `ODbL-1.0` 发布，易变运行日期写入 `fetch_metadata.json` 使规范化重下载与存档 SHA-256 一致；`calibrated_sim_germany.py` 默认最多启用 18 个工作进程）：
```bash
python experiment/data/fetch_germany.py --check            # 校验 2026-11-19 提取的 LF 原始文件与 SHA-256（3/3 OK）
python experiment/data/fetch_germany.py --verify-download  # 从公共接口在线重下载并规范化比对存档 SHA-256（3/3 OK）
python experiment/code/data_prep.py                        # 校验 SHA-256（0 个插补县）并生成 processed/panel.npz
python experiment/code/application.py                      # 生成 germany_application_results.json（含精确 WCR brentq、A4 置换检验、w/延迟与 B=2000 置换）
python experiment/code/calibrated_sim_germany.py           # 400 县 x 261 周校准模拟 (R=500 x 2，含 WCR-score、统一刀切、全量嵌套 pairs R=500/B=199 与 M2 先验敏感性/NUTS-2 全局检查)
python experiment/code/make_figure_germany.py              # 生成 Figure 3 (fig_03_application) 与 Web Figure 5 (fig_web_05_diagnostics)
python experiment/code/make_tables.py                      # 动态生成 8 张核心模拟表
python experiment/code/make_tables_germany.py              # 动态生成 Table 3--4 及 Web Appendix F 全部附表（共 10 张表）
python experiment/code/make_verification_tables.py         # 动态生成 Web Tables 7--8（非同构/同周网格/弱先验 M2 与纯 G/偏度 x T）
python scripts/rebuild_manifest.py                         # 以 LF SHA-256 注册全部 22 张活动表格至 tables/manifest.json
```

## 主要结果

- **正溢出核心模拟与边界界定（$R=8000$ + 非同构/分级验证）**：在工作更新模型、时间对齐与可分离报告成立下，M0 在全部 8 个正溢出核心情景与 4 个原假设情景中对 $\tau$ 偏倚 $-0.002$ 至 $0.001$（95\% 覆盖率 94.5\%--95.4\%），$\hat\phi,\hat\kappa$ 绝对偏倚至多 $0.004$（情景 9 为 $0.0035$）；反而在 $\phi_0=0$ 边界（全部 540 次大 $N$ 核查均见随 $\bar\phi$ 移动的第二根）、实施后报告偏移（X2）或动态效应（X5）下预设估计量失效。非同构日度 Gamma 分支（$R=500$）、弱先验 M2 重新评估（暴露负反馈/溢出偏倚）与分级 $w_0\in[0,0.30]$ 及处理比例匹配试验精确定位方法退化阈值。
- **适用包络与使用指南（Table 3 Operating Envelope）**：零传染边界下 M0-G/M0-GI 属无保证的启发式选根且 $\phi_0\in[0.05,0.10]$ 存在覆盖代价（推荐使用 $S_1(\phi)$ 识别稳健集为主）；周度聚合在 $w_0\le 0.10$ 稳定而在 $w_0\ge 0.20$ 明显恶化；$\omega_i\times[0.1,10]$ 与 GI/延迟卷积敏性检定规模/延迟依赖；非暑期假期安慰剂直接检出对假设 A7 可分离报告约束（$c_{it}=c_i c_t$）的偏离。
- **稀疏监测计数说明性测试平台（德国 BKG VG250 400 县、16 联邦州 2015--2019 年水痘面板，115,129 例，事后分析）**：
  - **稀疏面板可计算性与结构稳健诊断**：在 58.8\% 县周为零计数的面板上无需删零即可收敛至内部根 $\hat\phi=0.327$（$S_1(F)=[0.21,2.45]$），但 $F$ 校准联合/不变性集铺满 $[0,3]$ 且 $S_{3,\mathrm{oi}}(\chi^2)$ 为空；$\omega_i\times[0.1,10]$ 剖面使 $\hat\phi$ 跨 $0.783$ 至 $-0.118$、$\hat\psi^\dagger$ 变化逾 7 倍。
  - **两步法小计数符号反转（实测对比 + 已知真值校准模拟）**：真实面板中 M2（EpiEstim）在县级为 $+0.064$、NUTS-2 聚合（零计数 2.9\%）为 $-0.102$；在**按德国设计校准的维持工作模型已知真值模拟（Web Table 13(c)）**中，县级 M2 偏倚达 $+0.227$（主要由 Gamma(1,5) 先验伪计数主导），但实证 M2 相对校准均值跨所有先验高出 $+1.93$ 至 $+4.33$ 个 SD；NUTS-2 下全局 $\phi\ge 0$ 搜索 500/500 收敛而无约束拟合常见 $\hat\phi<0$。
  - **有限簇机制解析与 A4 正式检验**：在**德国校准模拟**（$R=500\times 2$）中常规 $t_{15}$ 覆盖率 0.898/0.888（SE 比值 0.850）；纯 $G$ 与偏度 $\times T$ 网格（Web Table 8）证实仅当严重簇规模偏度遇上长 $T$ 周剖面时 SE 比值才降至 $0.768$（$T=60$）。WCR-score（精确 `brentq` $<10^{-4}$ 端点 $[-0.422,-0.131]$ Webb6、$[-0.438,-0.119]$ Rademacher）、删一州刀切与嵌套 pairs 均改善校准覆盖（0.928--0.956）；29 对相邻州正式 A4 置换检验（$B=2000$）不拒绝（$Z_{A4}=-0.295, p=0.767$），A4$'$ 下 PSD 与四大区 SE 为 $0.0703$ 与 $0.0617$。
  - **工作模型数值输出与完整诊断归档**：基线对比 $\hat\tau+\hat\kappa = -0.240$ 仅作非因果工作输出报告，因数据按报告周记录、日历前置呈 $+0.247$--$+0.328$ 且 $B=2000$ 假期安慰剂与置换（$p=0.0075$）检出假期间违背 A7；完整诊断透明归档于网络补充材料（Web Appendices D--F）。

## 投稿前仍需处理

作者、单位、通讯、致谢细节与资助信息仍是占位（`Apevon Science`、页码 `000 0000`、`[Funding statement to be completed before submission.]`），请在正式投递 ScholarOne 前补齐真实作者、基金信息并按期刊要求确定最终版式。
