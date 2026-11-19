# -*- coding: utf-8 -*-
"""Round-3 revision: rewrite the README version notes so they describe the
current build and record this round's changes."""

import io

P = "README.md"
s = io.open(P, encoding="utf-8").read()

start = "2026-11-17 又按第三轮意见及 *Biometrics* 冲刺要求完成全面优化："
end_marker = "\n## 投稿前仍需处理"
i = s.index(start)
j = s.index(end_marker)

new = """2026-11-18 按第二轮评审（大修 6/10）的意见补做了应用诊断，2026-11-19 再按第三轮评审（小修 7/10）改写了措辞。合并后的当前状态：

- **应用改为对周度工作模型的诊断。** 加入区域 × 周、分区 × 周效应与区域特定 φ 后，$J_{LM}$ 仍在 16.9–17.8，校准 95 分位 6.38–6.66；$J=17.873$ 与 $Q_G=40.64$ 超过全部 600 次校准抽样，$Q_{GI}=27.62$ 被超过 8 次。τ+κ = −0.466（95% 区间 −1.091 到 0.159）因此报为被拒绝模型的伪真值，不作因果解读。
- **新增 Web Appendix F 与网络表 10–17**：区域效应、按应用校准的模拟（并给出 $J$、$Q$ 的参照分布）、代际分布同周质量重组、Stock–Wright 识别稳健集、事件时间隐含权重、报告率临界点。各表注标明对应的 `verification/revision/*/results.json`。
- **少簇推断成文。** Web Appendix E 报告删除一州刀切（标准误比夹心大 14%–30%）与两种 bootstrap-t；正文 §4 给出成对集群矩区间 $[-1.330, 0.397]$。
- **第三轮只改措辞，不动数字。** 讨论的校准覆盖率句限定为"工作模型成立时"，并写出同周情形 τ+κ 覆盖率 0.783；正文补出 τ+κ 投影所依赖的"给定 φ 强识别"条件（φ=0 处条件标准误 2.13）；事件时间写出隐含权重非凸（2.109 至 −0.471）并把 ω 加权对比标为线性化核对；正文补出集群规模 1–254 县与 τ̂ 随 φ̄ 的移动（0.074 到 0.018）；摘要 "regional week effects" 改 "division × week effects"。
- **排版与一致性。** 0.292 改为 0.291（二次舍入）；表 3 注的均值改为 $e_0$–$e_5$ 的均值；网络附录 C 印出 C.1／C.2 编号；A.12 去掉代码变量名；网络表 11 的 "stress design" 改名 "high-noise design"，避免与同周的 "stress regime" 混淆；主稿与 SI 均零 Overfull、零未定义引用。
- **篇幅。** 为容纳新增披露，压掉了约 380 词冗余（致谢与数据可得性声明合段，利益冲突与资助合用一个标题），未删任何数字、条件或披露。计入限页仍恰好 25 页。

"""

s = s[:i] + new + s[j:]
io.open(P, "w", encoding="utf-8", newline="").write(s)
print("version notes rewritten")