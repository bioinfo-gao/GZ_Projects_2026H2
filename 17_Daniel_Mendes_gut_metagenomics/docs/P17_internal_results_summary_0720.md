# P17 内部项目结果总结文档 — Daniel_Mendes_gut_metagenomics

**创建日期：2026-07-20**（内部文档，创建日期永久不变；后续更新记于文末"更新记录"）
**用途**：全项目结果的一站式内部总览（三条分析路线 + 交付情况 + 技术坑点），供内部复盘/后续项目参考。客户可见内容以两份 `custom_research_report_*` 为准，本文档不对外发送。

---

## 1. 项目基本信息

| 项 | 内容 |
| :--- | :---: |
| Project ID | QTE_26_06_25_001_Daniel_Mendes |
| Species | *Mus musculus*（host genome GRCm39） |
| Tissue/Material | Stool（fecal pellets） |
| 研究设计 | **HFD**（High-Fat Diet，高脂饮食）小鼠，case-control：**AL**（ad libitum，自由采食）vs **IF**（intermittent fasting，间歇禁食），每臂 n=5，共 10 样本 |
| 测序平台 | Illumina NovaSeq X Plus，PE150，shotgun metagenomics |
| 数据量 | 单样本 17.2–26.5 M read pairs（均值 22.8M），**合计 227.5 M read pairs（~68 Gbp）**。原 plan 用 gzip 文件大小反推 ~13M/样本，实测几乎翻倍，已在 plan 文档订正回填 |
| Host content | 2.2%–18.2%（均值 ~6.8%），符合粪便样本预期，微生物信号充足 |
| Plan 文档 | `docs/P17_research_plan_and_review_0717.md`（living doc，创建日 0717 永不改名） |

## 2. 研究设计：三条正交互补路线

本项目是 [[reference_shotgun_meta_three_routes]] skill 家族的实战 worked example，三条路线全部跑完并交付：

| Phase | 路线 | Pipeline | 交付日期 | 交付目录 |
| :--- | :--- | :--- | :---: | :--- |
| Phase 1 | 组成 taxonomy + diversity | nf-core/taxprofiler 2.0.1（Kraken2+Bracken / MetaPhlAn 4） | 分析完成 2026-07-18 | `custom_research_report_20260722/`（原 `_20260718/`，0722 因发出前 review 措辞改动而改名，见更新记录） |
| Phase 1b | 功能 function | HUMAnN 3.9（mag_biobakery env） | 分析完成 2026-07-19（并入同目录 `function/`） | 同上 |
| Phase 2 | 基因组 MAG | nf-core/mag 5.4.2 | 2026-07-20 | `custom_research_report_20260720/` |

编排原则：Phase1 → Phase1b 复用去宿主 reads，构成标准分析完整交付；Phase2 是增值项，在数据深度评估（~4 Gbp/样本，不足以单样本组装，改为按 diet arm group co-assembly）确认后才开跑，且未与 HUMAnN 满载线程叠跑。

## 3. Phase 1 结果摘要：组成 + 多样性

- **群落整体结构 AL vs IF 无统计学可分辨差异**：Bray-Curtis PERMANOVA R²=0.149, p=0.25（999 permutations）；between-group distance（0.295）≈ within-group distance（0.298），ratio=0.99。PCoA 前两轴解释 83.7%（PCo1 66.9% + PCo2 16.8%）。
- ***Akkermansia muciniphila* 是绝对优势菌种**，均值占比 AL 48.7% / IF 64.4%，IF 组呈上升趋势但不显著（Wilcoxon p=0.42）——方向与间歇禁食文献一致。
- **60 个 core species 差异丰度检验，无一通过 FDR<0.05**（最小 adjusted p=0.60）。非显著趋势：↑IF 为 *Akkermansia*、*Paramuribaculum intestinale*；↑AL 为 *Lactobacillus johnsonii*、*Lactococcus* spp.。
- **Alpha diversity**：AL 三项指数（Observed/Shannon/Simpson）均略高于 IF，但均不显著（p=0.31–0.42），趋势方向一致。
- **两工具交叉验证**：Bracken vs MetaPhlAn 中等一致（Spearman ρ=0.64），Bracken 报告丰度系统性偏高，主导菌种归属上有分歧（MetaPhlAn 把大量信号分给未命名 SGB），但均认同群落 *Akkermansia*-rich。
- **样本层面异质性是本研究真正的限制因素**（而非"n=5 太小"这么简单）：组内距离与组间距离相当甚至更大；层级聚类不按 diet 分组，反而 AL_7_06_12 与 IF_6_05_22（距离 0.08）等跨组样本对最相似。
- **两个明确离群样本**：HFD_AL_4_02_25（组内均距 0.479，全队列最偏离）、HFD_IF_4_03_11（0.317）。经检验二者偏离与测序深度（r=−0.11）、host 比例（r=+0.07）均不相关，是生物学真实差异而非技术假象；剔除离群样本后 between/within ratio 反而从 0.99 降到 0.88（更不可分），说明零结果是 robust 的，不是被噪声样本掩盖的假阴性。
- **⚑ 关键待澄清点（已在客户报告中作为 flag 提出）**：Aitchison/CLR PCA（Fig 10）显示 PC1（59% 方差）按 replicate ID 前缀（`4_…` vs `6_…`/`7_…`）分离，而非按 diet 分离——疑似隐藏的 cage/litter/batch/collection-timepoint 效应。已请客户澄清这些编号是否代表跨臂匹配的 cage/litter/timepoint；若是，可免费改跑 paired/blocked 设计，可能改变显著性结论。**此问题截至今日（0720）客户尚未回复，是项目唯一悬而未决事项。**

## 4. Phase 1b 结果摘要：功能 pathway（HUMAnN）

- HUMAnN 重建 **1,243 个 MetaCyc pathway features**，聚焦 128 个 well-represented community pathways。
- **功能层面同样无统计学可分辨差异**：pathway-level Bray-Curtis PERMANOVA R²=0.155, p=0.22（PCo1=59.2%）——diet 仍解释约 16% 的功能变异，但该 5 vs 5 设计不足以把它 resolve 出显著性（underpowered，非"零效应"结论，见 [[feedback_underpowered_not_null_report_framing]]）。
- **无 pathway 通过 FDR 校正**（最小 adjusted p=0.76）。方向性趋势：↑IF 为氨基酸/核苷酸/peptidoglycan 生物合成通路（生物合成能力偏强）；↑AL 为 glycolysis IV + 一条精氨酸合成通路。与 Phase 1 分类学趋势方向自洽。
- Top-25 高丰度 pathway（核苷酸/氨基酸/辅因子合成 + 中心碳代谢）在两组间共享核心代谢谱，量级相当。

## 5. Phase 2 结果摘要：基因组重建（MAG）

- **Pipeline 全绿**：3,328/3,328 Nextflow tasks 完成，0 failed。
- **Co-assembly 策略**：单样本深度（~4 Gbp）不足以支持个体组装，按 diet arm 分组 co-assembly（group-AL 汇总 5 样本 ~117M pairs，group-IF ~111M pairs）——因此 Phase 2 的比较是 AL-vs-IF group 级别，不是 animal-by-animal。
- **组装统计**：group-AL 297,332 contigs / 553.9 Mb / N50=7,366 bp；group-IF 302,440 contigs / 518.7 Mb / N50=5,593 bp（AL 略大略连续，与其略高的 pooled input reads 一致）。
- **Binning**：三 binner（MetaBAT2/MaxBin2/SemiBin2）+ DASTool 去冗余，**184 个非冗余 final MAG**（AL 95 / IF 89）——⚠ 注意原始 per-binner bins 是 1060 行（跨 binner 重复计数），**184 才是正确的"基因组数"**，见下文技术坑记录。
- **MIMAG 质量分级**（CheckM2 post-hoc 打分，不采信 DASTool 内部分数）：**High 79 / Medium 89 / Low 16**（AL: 44 High/43 Medium/8 Low；IF: 35 High/46 Medium/8 Low）。GTDB 分类率 167/184（91%）。两组质量分布无实质差异。
- **分类组成**：GTDB-classified MAG 中 **Bacillota 主导两组**（AL 84% 73/87，IF 80% 64/80），Bacteroidota/Actinomycetota 为稳定少数门——门水平无 diet 效应，与 Phase 1 群落层面 null 结果吻合。
- **Akkermansia muciniphila 表面矛盾，已解释**：仅在 AL 组装出 genome（100% complete, Medium quality, 5.18% contamination），IF 组未能干净分箱——尽管 Phase 1 读段丰度 IF（64.4%）高于 AL（48.7%）。检查 AL bin 的 per-sample coverage 发现跨 5 只小鼠有 8 倍差异（31×–238×），个体/株系异质性使该菌的 assembly graph 碎片化，这在 IF 组更严重导致分箱失败。**这不是矛盾而是信息量**：read-based abundance 与 genome assembly 是两种不同灵敏度，两者互补而非互斥（此点已写入客户报告 §6.3 作为诚实讨论，未掩盖）。
- **交叉验证 Phase 1 发现**：属水平最常被回收的 genera 包括 *Oscillibacter*、*Paramuribaculum*（均是 Phase 1 flag 出的 diet-associated trend taxa）——基因组层面独立佐证了群落层面的观察。

## 6. 跨路线整合结论

三条路线（组成/功能/基因组）在"AL vs IF 无统计学可分辨的全局差异"这一点上**完全一致**——不是巧合，是同一套生物学限制（组内异质性 + 可能的隐藏 batch 结构）在三个层面上的一致体现：

| 层面 | 效应量估计 (R²) | 显著性 | 结论 |
| :--- | :---: | :---: | :--- |
| 群落组成（taxonomy PERMANOVA） | 0.149 | p=0.25 | Underpowered，非零效应 |
| 功能通路（HUMAnN PERMANOVA） | 0.155 | p=0.22 | Underpowered，非零效应 |
| 基因组分类组成（MAG 门水平） | 定性一致 | — | 无 phylum-level shift |

三层 R² 都落在 0.15 附近，方向一致（*Akkermansia*/生物合成趋势 IF 偏高），**这是一个真实但当前样本量下无法统计学确证的 candidate effect**，不是简单的"没有差异"。

## 7. 局限性、统计效力与后续建议（已写入客户报告 §7.1）

限制因素排序：**组内异质性（个体/cage 效应）> 样本量本身**。建议的三个 lever（已提供给客户）：
1. 扩大样本量至 ~10–12 只/臂（最可靠的提升 power 方式）
2. 强化/标准化 fasting 干预强度，拉大组间对比
3. cage/litter-matched 配对 + 固定采样时间点，压缩组内方差（AL 组组内方差 0.340 高于 IF 的 0.256，优先收紧 AL）

**⚑ 唯一悬而未决事项**：客户尚未回复 replicate ID（`4_…`/`6_…`/`7_…`）是否代表 cage/litter/timepoint 匹配。若确认匹配，可免费用现有样本改跑 paired/blocked 设计，可能改变显著性结论——**下次与客户联系时应主动跟进此问题**。

## 8. 交付物清单

| 交付轮次 | 目录 | 内容 |
| :--- | :--- | :--- |
| 第一轮（分析完成 2026-07-18，报告 0722 发出前 review 改名） | `custom_research_report_20260722/`（原 `_20260718/`） | 报告 + `qc/`（multiqc）+ `taxonomy/`（组成图表+丰度表）+ `diversity/`（alpha/beta/差异丰度/离群样本分析）+ `function/`（HUMAnN pathway 结果） |
| 第二轮（2026-07-20） | `custom_research_report_20260720/` | 报告 + `qc/`（MAG pipeline multiqc）+ `assembly_binning/`（final_MAG_catalog.tsv + 质量分布图 + 分类组成图） |

**两轮均尚未实际发给客户**（截至 0722，用户明确表示"未全面审查完毕不发出"）。第一轮因发出前 review 时改进了 §7.1/Key Findings 的措辞（"ask"→"suggest" + 补充配对分析为何能还原真实处理效果的解释），按"发出前 review 原地改+改名到今天"规则改名为 `_20260722/`；第二轮尚未审阅，暂未改动/改名。两轮目录一旦实际发出，即视为冻结的历史记录，之后如需新分析则另开新日期目录，不再改写。

## 9. 关键技术坑与经验教训（内部记录，不对客户展示）

1. **HUMAnN 3.9 ↔ MetaPhlAn 4.2.4 三连坑**（详见 [[reference_shotgun_meta_three_routes]]）：① version-check 解析崩溃需 PATH 前置 shim；② db 代际不匹配需用 vJun23（非新装的 vJan25）；③ `--bowtie2out` 需 shim 翻译为 `--mapout`。且必须整段 `conda activate mag_biobakery`，不能逐条 `conda run`，否则子进程 PATH 丢失重新踩坑①。黄金参考脚本：`scripts/7_run_humann.sh` + `scripts/mpa_shim/`。
2. **DASTool final bin 与原始 per-binner bin 命名不同、数量差 5.8 倍的陷阱**：CheckM2/BUSCO/GTDB-Tk 这轮跑在原始 1060 行 per-binner bins（跨 MetaBAT2/MaxBin2/SemiBin2 有重复计数同一 genome）上，不是跑在 DASTool 去冗余后的 184 个 final bins 上；二者命名规则不同（`XxxRefined-group-XX.NNN[_sub].fa` vs 原始命名，SemiBin2 还多一层 `_`→`.` 转换）。写了 `scripts/15_mag_final_bin_catalog.py` 做名字映射把 184 个最终 bin 对回原始 QC/分类记录（100% 映射成功）——**任何人复用这批 MAG 输出都必须先踩这个坑**，否则会误把 1060（含重复）当作"基因组数"上报。
3. **数据量估算陷阱**：plan 阶段用 gzip 文件大小反推 read pairs 数（~13M/样本），实测（22.8M/样本均值）几乎翻倍，已订正回填 plan 文档 — 提醒今后优先用 `seqkit stats` 或实测而非文件大小反推。
4. **co-assembly 分组决策**：单样本深度 ~4 Gbp 低于个体组装门槛，改为按 diet arm（而非 cohort-wide 或 per-sample）co-assembly，是深度/对比度之间的标准折中——此决策在 Phase2 启动前已确认，避免了组装完成后才发现无法回答 AL/IF 对比的返工风险。

## 10. 项目状态

**三路线（组成/功能/MAG）全部交付完毕，项目 17 分析与交付流程已实质完结。** 唯一开放项：等待客户回复 replicate ID 是否代表 cage/litter/timepoint 匹配（§7）；若回复确认匹配，需追加一轮免费 paired/blocked 重分析并生成新一轮带日期的交付目录（不覆盖现有两轮）。

大文件已按 [[reference_git_skill]] 加入 `.gitignore`（GenomeBinning/QC 原始 BUSCO/CheckM2/DASTool 产物 13G+ 排除，只留 bin_summary.tsv）；两轮 commit 均已 push 成功（3c11996 + 7199920），仓库与 origin 完全同步。

---

## 更新记录

- 2026-07-20 — 文档创建，汇总 Phase1（组成+多样性）、Phase1b（功能 HUMAnN）、Phase2（MAG）全部结果。
- 2026-07-21 — §1 补全 HFD 缩写全称（High-Fat Diet）。
- 2026-07-22 — 第一轮报告发出前 review：改进 §7.1/Key Findings 措辞（"ask the client"→更礼貌的
  "suggest"，并给实验学家补充"为什么配对分析能还原真实处理效果"的解释）；因未实际发出，原地改+
  文件夹/报告文件改名 `_20260718/…_0718.md` → `_20260722/…_0722.md`（发出前 review 规则）。同步修正
  0720 报告里指向旧路径的交叉引用。§2/§8 表格路径同步更新。
