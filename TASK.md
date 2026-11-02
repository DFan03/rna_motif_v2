# TASK (verbatim brief from the user, kept so the README can answer it point by point)

## 研究背景
病毒 RNA 包含多种可能参与 RNA 稳定性、翻译调控和分子相互作用的局部序列及结构模式。
现有公共数据库提供了大量病毒核酸序列和注释，但绝大部分局部片段没有直接的功能实验标签。
请设计一个计算工具，从公开病毒 RNA 相关数据中识别、统计和排序潜在的 RNA 调控 motif，并说明如何评价候选结果的可靠性。

## 任务1：数据获取与质量控制
选择一个公开病毒序列数据库，获取少量适合分析的 RNA 序列记录（或使用离线 FASTA/元数据）。
- Q1.1 解析序列及其元数据。
- Q1.2 识别序列类型、来源和已有注释。
- Q1.3 处理重复序列、非标准碱基、长度异常和缺失数据。
- Q1.4 记录序列来源、物种/分类信息、数据库版本。
- Q1.5 解释如何避免由于同源病毒被重复采样而产生统计偏差。

## 任务2：motif 发现算法
实现候选序列 motif 识别方法（k-mer、PWM、统计富集、序列聚类等均可）。
- Q2.1 找到候选局部序列模式。
- Q2.2 统计 motif 的出现频率和覆盖范围。
- Q2.3 与合理背景进行比较。
- Q2.4 输出候选 motif 的排名和统计证据。
- Q2.5 不要求发现真实功能元件，但必须说明候选 motif 与已验证功能元件的区别。

## Overall: also "说明如何评价候选结果的可靠性" (how to judge reliability of candidates).
## User constraints: <1 h total, local M2 Pro CPU, results + workflow in GitHub README; principle papers: Jin 2025, Razumova 2025, Bernett 2024.
