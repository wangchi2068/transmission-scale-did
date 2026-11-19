# -*- coding: utf-8 -*-
"""Round-3 revision: bring the Chinese reading copy in line with the English
manuscript after the second review (M1-M5 and the minor items)."""

import io

P = "paper_zh/main.tex"
s = io.open(P, encoding="utf-8").read()

pairs = []

# M2 / abstract scope: regional -> division
pairs.append((
 "在美国州级口罩令与县级 COVID-19 病例的应用中，这些诊断拒绝了周度工作模型：过度识别统计量高于全部校准抽样，加入区域周效应后依然如此。",
 "在美国州级口罩令与县级 COVID-19 病例的应用中，这些诊断拒绝了周度工作模型：过度识别统计量高于全部校准抽样，加入分区 $\\times$ 周效应后依然如此；加入过度识别矩后，耦合参数的稳健集为空。"))

# M2: robust sets do not depend on the selected root
pairs.append((
 "设计还附带能够显示自身工作模型何时失效的诊断：剖面传染方程揭示多根，过度识别检验对照校准参照分布解读，识别稳健集无需选根。",
 "设计还附带能够显示自身工作模型何时失效的诊断：剖面传染方程揭示多根，过度识别检验对照校准参照分布解读，$\\phi$ 的识别稳健集不依赖所选根。"))

# M2: same-week regime comes from the generation-interval arithmetic
pairs.append((
 "在美国州级口罩令的应用中，这些诊断在全国周效应与区域周效应下都拒绝了周度工作模型，所报告的政策层面总效应因此是周度拟合的概括量，不作因果解读；诊断同时指出，同周传播是周度模型无法处理的情形。",
 "在美国州级口罩令的应用中，这些诊断在全国周效应与区域周效应下都拒绝了周度工作模型，所报告的政策层面总效应因此是周度拟合的概括量，不作因果解读。代际权重把约 0.30 的传播落在与感染者同一周内，正处于周度模型无法处理、$J$ 也无法察觉聚合误差的区间（第~\\ref{ss:sameweek} 节）；判定周度规格不可用的是这一算术，而非某项诊断检验。这里需要说明一处限制：$\\phi$ 的稳健集含一个边界分量和一个内部分量，加入过度识别矩后集合为空，因此在本应用中稳健集确认了拒绝，却未把选根问题替换为可用的区间。"))

# M5: unequal cluster sizes in the main text
pairs.append((
 "49 个集群中只有 24 个包含切换县。",
 "49 个集群中只有 24 个包含切换县，且集群规模相差悬殊：最小 1 个县，最大 254 个县，最大的州在每个矩条件中都占很大权重。"))

# M4: non-convex implicit weights, linearisation check
pairs.append((
 "合并的 M0 在该样本上给出 $\\hat\\tau=0.056$（标准误 0.181）；减去隐含加权的分箱（网络附录~\\ref{ss:hetero}）为 0.079（标准误 0.240，$p=0.742$），减去等权分箱为 0.733（0.457，$p=0.117$；网络表 14）。这一差异与动态加权相符，但无法与抽样误差区分。",
 "合并的 M0 在该样本上给出 $\\hat\\tau=0.056$（标准误 0.181）。合并估计目标的隐含权重在实施分箱上为 2.109，其后各分箱为负，最低到 $e6+$ 的 $-0.471$，因此合并估计不是各分箱效应的凸组合。减去隐含加权分箱后为 0.079（标准误 0.240，$p=0.742$），这是对线性化的核对而非两个估计目标的比较，按构造两边在一阶上相等；减去等权分箱后为 0.733（0.457，$p=0.117$），这一对比才提供信息（网络表 14）。这一差异与动态加权相符，但无法与抽样误差区分。"))

# M1: discussion coverage sentence
pairs.append((
 "集群数是第三个条件：定理~\\ref{thm:asy} 是大 $G$ 结果；在 49 个州、其中 24 个有切换的情形下，按应用校准的模拟支持 $t_{G-1}$ 夹心区间（$\\tau+\\kappa$ 覆盖率 0.935 到 0.950），而受处理集群较少时推断可能失效\\citep{ferman2019inference,canay2021wild}。",
 "集群数是第三个条件：定理~\\ref{thm:asy} 是大 $G$ 结果；在 49 个州、其中 24 个有切换的情形下，若工作模型成立，按应用校准的模拟支持 $t_{G-1}$ 夹心区间，$\\tau+\\kappa$ 的覆盖率在基准设计中为 0.935，把州 $\\times$ 周冲击调大到与应用的标淮误相当后为 0.950；受处理集群较少时推断一般可能失效\\citep{ferman2019inference,canay2021wild}。这一校准不覆盖本应用所处的同周情形，此时同一区间对 $\\tau+\\kappa$ 的覆盖率降至 0.783，$J$ 对聚合没有检验力（第~\\ref{ss:sameweek} 节）。"))

# minor 1: SE 0.292 -> 0.291
pairs.append(("由 0.311 变为 0.292", "由 0.311 变为 0.291"))

# M3 + phi_bar movement: robust-set paragraph
pairs.append((
 "加入过度识别矩后集合为空，与 $J$ 一致。$\\tau+\\kappa$ 的 Bonferroni 投影为 $[-3.55,6.29]$，内部分量为 $[-1.19,0.26]$（网络表 13）。内部根与 $\\hat\\phi$ 依赖工具参考值 $\\bar\\phi$。",
 "加入过度识别矩后集合为空，与 $J$ 一致。$\\tau+\\kappa$ 的 Bonferroni 投影为 $[-3.55,6.29]$，内部分量为 $[-1.19,0.26]$（网络表 13）。投影把 $\\phi$ 的 97.5\\% 集合与给定 $\\phi$ 时 $\\tau+\\kappa$ 的 97.5\\% Wald 区间合并，其第二步要求 $(\tau,\kappa)$ 在该 $\\phi$ 处强识别。边界分量上这一前提并不可靠：$\\phi=0$ 处拟合给出 $\\hat\\tau=-2.00$、$\\hat\\kappa=3.37$，$\\tau+\\kappa$ 的条件标准误为 2.13。因此本文把投影报告为在该条件下算出的界，而不是具备覆盖保证的区间，并推荐不依赖它的 $\\phi$ 的 $S$ 集。内部根与 $\\hat\\phi$ 依赖工具参考值 $\\bar\\phi$；$\\hat\\tau$ 也随之移动，从 $\\bar\\phi=0.25$ 时的 0.074 到 $\\bar\\phi=1$ 时的 0.018，而 $\\hat\\tau+\\hat\\kappa$ 保持在 $-0.476$ 到 $-0.466$ 之间。"))

# high-noise rename in Chinese application paragraph
pairs.append(("在标准误与应用相当的压力设计中为 0.950",
              "在把冲击调大到与应用标准误相当的高噪声设计中为 0.950"))
pairs.append(("在压力设计中（$R=300$）为 0.950、0.930 与 0.923",
              "在高噪声设计中（$R=300$）为 0.950、0.930 与 0.923"))
pairs.append(("压力设计调高冲击标准差，使 $\\hat\\tau+\\hat\\kappa$ 的标准误接近应用中的 0.311。",
              "高噪声设计调高冲击标准差，使 $\\hat\\tau+\\hat\\kappa$ 的标准误接近应用中的 0.311。"))

# C.1 numbering and code variable in the Chinese appendix
pairs.append(("论文表 7", "网络表 7"))

n_missing = 0
for a, b in pairs:
    if a not in s:
        print("MISSING: %r" % a[:60])
        n_missing += 1
        continue
    s = s.replace(a, b)

io.open(P, "w", encoding="utf-8", newline="").write(s)
print("missing=%d" % n_missing)