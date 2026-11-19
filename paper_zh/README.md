# 中文稿（英文稿的阅读版）

`main.pdf`：中文全文（ctexart，25 页，页数随系统中文字体略有浮动），结构与修订后的英文稿 `../main.pdf` 逐节对应，数字与英文稿及结果文件一致。供阅读，不作为投稿版本。

## 与英文稿的对应关系

| 英文稿 | 中文稿 |
|---|---|
| Abstract | 摘要（七句，问题、缺口、做法与主要结果在前） |
| 1 Introduction | 1 引言 |
| 2 Methods（6 小节，含命题 1 报告率） | 2 方法（6 小节） |
| 3 Simulation Study（设计、结果、3.2.1 周内传播） | 3 模拟研究（同结构） |
| 4 Application: Diagnosing a Weekly Transmission Model for State Mask Orders | 4 应用：诊断州级口罩令的周度传播模型 |
| 5 Discussion | 5 讨论 |
| Back matter（致谢、数据可得性、支持材料、利益冲突、资助） | 同名声明 |
| Web Appendix A.9–A.12 | 网络附录 A.9–A.12 |
| Web Appendix D（数据来源与清洗） | 网络附录 D |
| Web Appendix E（少簇推断，含校准覆盖率） | 网络附录 E |
| Web Appendix F（应用诊断 F.1–F.6） | 网络附录 F |

定理证明与正则条件（Web Appendix A.1–A.8）、数值核对（B）、补充模拟表图（C）未翻译，正文引用处写作“网络附录 A”“网络附录 B”“网络附录 C”。

## 排版说明

- 图 1–3 直接使用英文稿的图片文件（`../figures/`，图内文字为英文），图题译为中文，编号与英文稿相同。
- 中文稿不排版表格；正文提到表时直接写英文稿编号（表 1–3、网络表 1–15），表格见英文稿与网络支持材料。
- 数学字体与英文稿相同，用 Latin Modern Math（unicode-math）。
- 定理、命题、引理、推论、注的编号与英文稿一致（命题 1 报告率、命题 2 识别来源、命题 3 异质效应；注 2–4 为动态效应、过程噪声与耗竭、可检验性）。
- 术语：准差分、鞅差阵列、聚类夹心方差、避免病例数、伪真值、剖面矩、过度识别、政策层面总效应、识别稳健集、压力区间。

编译：`xelatex main` → `bibtex main` → `xelatex main` ×2，依赖 ctex、unicode-math 与系统中文字体；图片从 `../figures/` 读取。`bibliography.bib` 与英文稿相同。
