# -*- coding: utf-8 -*-
"""Round-3 revision: bring the repository README in line with the built
documents. Only facts are corrected: page counts, appendix letters, the
verification tree, and the version notes for the third-round changes."""

import io

P = "README.md"
s = io.open(P, encoding="utf-8").read()

old_head = """论文、LaTeX 源码、实验代码与验证材料。按 *Biometrics* Regular Article 规范拆分为两份独立 PDF：
- 英文主稿 `main.pdf`（共 31 页）：首页为标题/摘要（不计入限页），正文、Back Matter 与 References 为第 1–24 印刷页（严格 ≤ 25 页上限），末 6 页为 `endfloat` 排出的图表（不计入限页）。
- 网络补充材料 `supporting_information.pdf`（共 30 页）：Web Appendices A–E 与附带文献，与主稿通过 `xr` 双向交叉引用。
- 中文阅读稿 `paper_zh/main.pdf`（22 页）。"""

new_head = """论文、LaTeX 源码、实验代码与验证材料。按 *Biometrics* Regular Article 规范拆分为两份独立 PDF：
- 英文主稿 `main.pdf`（共 32 页）：首页为标题/摘要（不计入限页），正文、Back Matter 与 References 为第 1–19 印刷页，加上 6 页 `endfloat` 排出的图表（不计入限页），计入限页的共 25 页（上限 25，恰好用满）。
- 网络补充材料 `supporting_information.pdf`（共 43 页）：Web Appendices A–F 与附带文献，与主稿通过 `xr` 双向交叉引用。
- 中文阅读稿 `paper_zh/main.pdf`（26 页）。
- 回复信 `response/`：`response-to-reviewers-r1r2.md`（第一、二轮）与 `response-to-reviewers-r3.md`（第三轮小修）。"""

assert old_head in s
s = s.replace(old_head, new_head)

old_tree = """├── verification/
│   ├── root_selection/        根选择规则在 9 个独立设计上的验证（18000 次重复）、边界推断诊断
│   ├── application/           τ+κ 联合区间、M0-ET 分箱修正、STROBE 清单、稳健推断方案、少簇检查（删除一州刀切 + 野聚类 bootstrap-t）
│   └── same_week/             同周传播的匹配模拟（w0 = 0 对 w0 = 0.3，每臂 1000 次）
└── scripts/                   构建与质量检查脚本（见下）"""

new_tree = """├── verification/
│   ├── revision/              修订轮新增的应用诊断（区域效应、校准模拟、代际重组、稳健集、事件权重、报告率、总体根），各子目录的 results.json 是正文数字的出处
│   ├── root_selection/        根选择规则在 9 个独立设计上的验证（18000 次重复）
│   ├── application/           τ+κ 联合区间、M0-ET 分箱修正、STROBE 清单、少簇检查（删除一州刀切 + 野聚类 bootstrap-t）
│   └── same_week/             同周传播的匹配模拟（w0 = 0 对 w0 = 0.3，每臂 1000 次）
├── response/                  两封逐点回复信
└── scripts/                   构建与质量检查脚本（见下）"""

assert old_tree in s
s = s.replace(old_tree, new_tree)

io.open(P, "w", encoding="utf-8", newline="").write(s)
print("head and tree updated")