# 传播尺度双重差分（QD-GMM）

论文、LaTeX 源码、实验代码与验证材料。按 *Biometrics* Regular Article 规范拆分为两份独立 PDF：
- 英文主稿 `main.pdf`（共 31 页）：首页为标题/摘要（不计入限页），正文、Back Matter 与 References 为第 1–24 印刷页（严格 ≤ 25 页上限），末 6 页为 `endfloat` 排出的图表（不计入限页）。
- 网络补充材料 `supporting_information.pdf`（共 30 页）：Web Appendices A–E 与附带文献，与主稿通过 `xr` 双向交叉引用。
- 中文阅读稿 `paper_zh/main.pdf`（22 页）。

## 目录

```
QD-GMM_transmission_DiD_2026-11-13/
├── main.pdf, main.tex, supporting_information.pdf, supporting_information.tex
├── sections/, tables/, figures/
├── bibliography.bib, biom.cls, biom.bst, endrotfloat.sty, compile.sh
├── paper_zh/                  中文阅读稿（结构与英文逐节对应，含最新文献核验）
├── experiment/
│   ├── code/                  dgp、qdgmm、comparators、psi、simulate、application、make_tables …
│   ├── theory/                numeric_checks.py（Web Appendix B 的 8 项数值核验）
│   ├── data/                  raw（JHU、CDC 口罩令、Census 邻接与通勤、ACS）与 processed 面板
│   ├── analysis/              results.json、metrics.csv、实验报告与设计文档
│   └── logs/                  运行日志
├── verification/
│   ├── root_selection/        根选择规则在 9 个独立设计上的验证（18000 次重复）、边界推断诊断
│   ├── application/           τ+κ 联合区间、M0-ET 分箱修正、STROBE 清单、稳健推断方案、少簇检查（删除一州刀切 + 野聚类 bootstrap-t）
│   └── same_week/             同周传播的匹配模拟（w0 = 0 对 w0 = 0.3，每臂 1000 次）
└── scripts/                   构建与质量检查脚本（见下）
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
也可直接执行 `bash compile.sh`（或 `bash compile.sh --clean` 编译并清理中间临时文件），依赖 XITS Math 与 Times New Roman 字体、`xelatex`/`bibtex` 及 `pypdf`。

中文稿：`cd paper_zh && xelatex main && bibtex main && xelatex main && xelatex main`。

实验：见 `experiment/code/README.md`，依赖锁在 `experiment/code/requirements.txt`（Python 3.11）。应用部分可用
`python -I application.py --quick` 做烟测，应与论文一致：τ = 0.040、κ = −0.506、φ = 0.772、τ+κ = −0.466。

## 主要结果

- 正溢出核心情景（8 个）：M0 对 τ 的偏倚 −0.002 到 0.001，95% 区间覆盖率 94.5%–95.4%，反馈实施、过度离散与小计数同时在场。
- 无传染边界情景：局部求解约一半重复选到非相合根；剖面规则把它修回。Oracle 在 T=26 随机实施下覆盖率为 80.2%/80.5%，系无尺度权重 $s_{it}$ 导致的单个州集群极端集中（情景 16 单簇占全样本病例 95 分位达 58.5%）所致；使用 $s_{it}$ 的 M0$'$ 与 M0-GI 保持约 95% 覆盖率。
- 同周传播匹配实验（w0 = 0 对 0.3，每臂 1000 次）：单元内周内再感染（非线性周内份额 14.5%，A4 跨集群独立性仍成立而条件均值 A1 误设）使 τ 偏倚 −0.033、覆盖率 0.890；τ+κ 偏倚 −0.053、覆盖率 0.783；J 检验 4.9% 对 4.8%，无功效。
- 应用：政策层面总效应 τ+κ = −0.466（95% 区间 −1.091 到 0.159），来自联合协方差；工作模型被过度识别检验拒绝（J = 17.873，p = 0.0001）。
- 少簇检查（网络附录 E）：删除一州刀切标准误比夹心大 13%–30%；Mammen 得分野聚类 bootstrap-t 临界值 1.72–1.94（τ+κ 区间 −0.999 到 0.067）；纳入雅可比簇间变异的非线性成对集群矩 bootstrap-t 临界值 1.59（τ）、2.22（κ）、2.44（φ）、2.78（τ+κ），政策层面合计区间 −1.330 到 0.397 仍含零。

## 目录内工具

| 脚本 | 作用 |
|---|---|
| `experiment/theory/numeric_checks.py` | Web Appendix B 的 S1、S3、S4、S5、S7、S8 数值核验，8/8 通过 |
| `verification/root_selection/build_rsel_table.py` | 由 `full.jsonl` 重算 Web Table 7（t 临界值、Multi 列） |
| `experiment/code/make_tables.py` | 由 `analysis/results.json` 生成正文与应用表 |
| `scripts/regen_app_table.py` | 由应用重跑结果重建 Table 3 |
| `scripts/split_app_table.py` | 主分析表与敏感性表（Web Table 8）分别生成 |
| `scripts/scan_control_chars.py` | 扫描源文件里被 shell 转义破坏的控制字符 |
| `scripts/qa_checks.py` | running head 长度、摘要词数、表注指向 |
| `scripts/section_pages.py` | 报告各节在 PDF 中的起始页 |

## 版本说明

`5eaa6a7` 之后的这一版按 2026-11-16 的同行评审意见做了三类改动：

**错误修正**
- Theorem 3(a) 后的衰减表述改写：上界随 T 增长，不推出衰减速率，也未证明下界。
- M0-ET 剔除窗口前已暴露的 435 个单位（剩 2656 个单元、37 个州），联合检验 F = 2.089（7, 30），p = 0.076，卡方形式 p = 0.014。
- 描述统计区分 1768（进入处理的次数）与 1720（改变状态的单位）。
- 参数空间上界改用 `\bar\Phi`，与工具参考值 `\bar\phi` 分开。

**新增证据**
- 同周传播匹配实验（`verification/same_week/`），写入正文 3.2.1 节。
- 边界推断诊断（`verification/root_selection/boundary_inference.md`），Oracle 重尾诊断写入正文 3.2 节。
- Web Table 7 按 t_{G−1} 重算并补 Multi 列；Table 3 以 τ+κ 为首列。

**结构**
- 方法节四小节（伪真值、异质与动态效应、避免病例数、周内传播与 A4）移入 Web Appendix A.9–A.12，正文保留一段概述。
- 模拟、应用、讨论、引言的篇幅与语域按 Biometrics 同类论文重写，正文由压缩前的 40 页降到 26 页。
- 新增 Web Appendix D（数据来源与清洗）；补利益冲突、资助、报告规范声明与 STROBE 清单。
- 旧的 `paper.pdf`、`paper.html` 已删除，`main.pdf` 是唯一编译产物。

2026-11-17 又按第三轮意见及 *Biometrics* 冲刺要求完成全面优化：

- **少簇与数值修正**：新增 Web Appendix E（应用少簇推断：删除一州刀切、Mammen/Webb6 得分野聚类 bootstrap-t 与 999/999 收敛的非线性成对集群矩 bootstrap-t，见 `verification/application/wild_bootstrap.py`）；修正 Theorem 3(b) 反馈条件均值偏移与 Lemma 5 逐行中值标记、Proposition 3 有界/无界支撑表述及 `numeric_checks.py` S1/S3/S5/S7/S8（8/8 通过）；明确 Scenario 16 未加权 Oracle 覆盖率下降系单簇集中度所致（中位 21.5%，95 分位 58.5%），同周实验明确区分非线性再感染份额（14.5%）与合计份额（34%）。
- **结构拆包与合规（≤ 25 页）**：将 Web Appendix A–E 拆为独立 `supporting_information.tex`（30 页），主稿启用 `endfloat` 将图表移至文献之后并通过 `\externaldocument[][nocite]` 实现 `xr` 双向交叉引用；精简正文 1–5 节与 Back Matter 后，计入 *Biometrics* 限页的正文+声明+References 为 24 印刷页（严格 ≤ 25 页，零 `Overfull \vbox/\hbox`、零未定义引用）。
- **参考文献核验**：为全部 7 篇 2021–2026 年 arXiv 预印本补齐已核验的 `eprint` 与 `journal` 条目，`bibtex` 零警告。

## 投稿前仍需处理

作者、单位、通讯、致谢细节与资助信息仍是占位（`Apevon Science`、页码 `000 0000`、`[Funding statement to be completed before submission.]`），请在正式投递 *Biometrics* ScholarOne 前补齐真实作者与基金信息。
