# Ellen 人源基因敲入小鼠 WGS — 分析方案 (Analysis Plan)

- **Project**: 13_Ellen_knockin_wgs
- **Plan Date**: 2026-07-05
- **Prepared by**: Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
- **Data**: `/home/gao/Dropbox/Ellen/` — 6 samples, WGS PE150, NovaSeq X Plus
- **参考方法来源**: `4_wgs_human_immu/` (nf-core/sarek 已装好并测试通过) 提供比对/QC/资源管理的基础做法；本项目的整合位点/拷贝数/SV 分析为在其之上的定制流程。

---

## 1. 样本清单 (来自 Sample_Sheet_for_interal_use.xlsx)

**样本命名约定（客户确认）**：样本名下划线前的部分 = 敲入的人源基因标识；下划线后为个体/founder 编号。

| 编号 | 样本名 | FASTQ 前缀 (S 号) | Sample Sheet 标注 | R1/R2 大小 | 角色 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| 1 | CD1A_B125 | S83 | Human knock in | 26G / 26G | **CD1A 敲入** |
| 2 | RAGH_153 | S89 | Mouse (普通) | 17G / 17G | **对照 (transgene 阴性)** |
| 3 | RAGH_273 | S90 | Mouse (普通) | 20G / 20G | **对照 (transgene 阴性)** |
| 4 | MTTH_284 | S84 | Human knock in | 17G / 17G | **MTTH 敲入 (founder 1)** |
| 5 | MTTH_412 | S85 | Human knock in | 17G / 17G | **MTTH 敲入 (founder 2)** |
| 6 | MTTH_524 | S86 | Human knock in | 18G / 18G | **MTTH 敲入 (founder 3)** |

\* 预估覆盖度 ~20–35×/样本（PE150，2.7 Gb 小鼠基因组），QC 后核实。对整合位点与拷贝数分析足够。

**两个敲入人源基因 = CD1A 与 MTTH；RAGH 两样本为对照。基因身份核查结果：**
- **CD1A** ✅ 已确认是真实人类蛋白编码基因（ENSG00000158477，GENCODE v45，免疫相关，呈递脂类抗原）。
  **利好**：小鼠仅保留 CD1d，不含 group-1 CD1（CD1a/b/c），因此人源 CD1A **几乎无小鼠同源交叉比对**，整合信号会很干净。可直接从本机 GRCh38 提取 CD1A 序列作 bootstrap。
- **MTTH** ⚠️ **不是标准人类基因符号**（不在 GENCODE v45 中），可能是客户内部构建体/品系简称。**需向 Ellen 确认 MTTH 对应的确切基因/序列**，否则无法为该组建混合参考。
- **RAGH** 非基因符号（人类有 RAG1/RAG2）。Sample Sheet 标为普通 "Mouse"，判定为对照。仍建议确认其遗传背景（是否 RAG 相关免疫缺陷背景），但只要 transgene 阴性即可作阴性对照与拷贝数基线。

---

## 2. 对朋友策略的评估：**核心思路完全合理，是行业标准做法**

朋友提出的三步——(1) 混合参考基因组，(2) 用嵌合读段（split reads + discordant pairs）定位整合位点，(3) 深度比值估拷贝数——**在方法学上正确且是转基因整合定位的标准范式**。逐条确认：

| 朋友的观点 | 评价 | 补充/修正 |
| :--- | :---: | :--- |
| 构建混合参考 (mm39 + 外源序列) | ✅ 正确且必需 | mm39=GRCm39，本机已就绪。**外源序列必须是客户提供的真实构建体全序列**（见下方 Blocker） |
| Split reads 定位碱基级断点 | ✅ 正确 | 软剪切读段跨越 mouse↔transgene 结合部，给出单碱基精度 breakpoint |
| Discordant pairs (R1 小鼠 / R2 人源) | ✅ 正确 | 与 split reads 互为佐证，提升灵敏度与置信度 |
| 深度比值估拷贝数 | ✅ 原理正确 | 需注意：**人源基因若与小鼠有同源区，会交叉比对污染深度**；需用唯一比对 (MAPQ 过滤) + 用野生型对照区域做基线 |
| 评估是否破坏小鼠基因/引起 SV | ✅ 正确且是重要交付物 | 断点注释到 GENCODE vM35 基因模型；结合点局部做 SV/CNV 检测 |

**朋友没有讲到、但决定成败的三个点：**
1. **必须拿到外源构建体的精确序列**（见 Blocker）。没有它，整合定位退化为"未知插入检测"，难度剧增、精度骤降。
2. **人源基因身份**：两个基因是什么？若与小鼠有直系同源（如很多免疫基因），需对同源小鼠区做屏蔽或强制唯一比对，否则整合信号与内源信号混淆。
3. **对照的意义**：野生型对照有两大作用——(a) 确认 transgene contig 在野生型中零覆盖（阴性对照，证明信号特异）；(b) 提供单拷贝背景深度用于拷贝数归一化。这要求先厘清第 1 节的样本角色歧义。

---

## 🚧 执行前的头号 Blocker（必须先要到）

**外源构建体序列文件** —— 来自客户 Ellen 的：
- 两个人类基因的**精确序列**（cDNA 还是基因组片段？含哪些外显子/内含子？）
- **载体/质粒骨架**全序列：启动子、增强子、polyA、筛选标记 (如 Neo/Puro)、LoxP/FRT、接头等
- 最好是 **GenBank (.gb) 或注释过的 FASTA / 质粒图谱 (SnapGene .dna)**

只要这个到位，其余全部可自动化推进。

**当前进展（基于样本名=基因名的确认）**：
- **CD1A 组**：可立即用本机 GRCh38 提取 CD1A 基因序列作 bootstrap 混合参考，先跑通流程。
- **MTTH 组**：MTTH 非标准基因符号，**必须先向客户确认它对应什么**，否则无法建该组混合参考。
- 无论哪组，用 GRCh38 提取的仅是**人源基因本体**，**会漏掉载体骨架（启动子/polyA/标记等）处的整合信号**，精度打折——真实质粒图谱仍应尽量索取。

---

## 3. 分析流程 (Pipeline)

> 遵循 CLAUDE.md：拆分为独立顺序步骤（不用三通管道），每步各自 tmux，参数按本机资源上限（≤28 物理核 / 56 线程）。

### Step 0 — QC
- `fastp` 每样本单独跑（`-w 8`），产出 HTML/JSON；`MultiQC` 汇总。
- 确认真实覆盖度、接头、重复率。

### Step 1 — 构建混合参考 (Hybrid Reference)
- `cat GRCm39.primary_assembly.genome.fa  construct.fa  >  GRCm39_plus_construct.fa`
- 外源序列作为独立 contig（如 `>TG_construct`）。
- `bwa-mem2 index` + `samtools faidx` + `samtools dict`。
- 本机已有 GRCm39 (GENCODE M35) 及其 bwa-mem2 索引，只需对**加了 contig 的新参考重新建索引**。

### Step 2 — 比对 + 排序 + 去重
- Step 2a `bwa-mem2 mem -t 20 ... | samtools sort -@ 8 -m 8G`（2 通管道，每样本单独）
- Step 2b `samtools markdup`（或 sambamba/GATK MarkDuplicates）
- Step 2c `samtools index -@ 8`
- **同时只跑 2 个样本**（沿用项目 4 经验，防内存紧张）。

### Step 3 — 整合位点检测 (核心交付物 1)
两条互补路线，结果交叉验证：
- **(A) 定向嵌合读段法**：`samtools view TG_construct` 抓所有落在外源 contig 上的读段及其配偶 → 看配偶/软剪切段落回小鼠基因组何处 → 聚类成整合位点，取 split-read 断点为碱基级坐标。
- **(B) 通用 SV caller**：`manta` 或 `delly`（BND/translocation 模式）在混合参考上直接调用，输出跨 contig 断点。
- 汇总每个整合位点：小鼠侧染色体坐标、外源侧坐标、支持读段数（split + discordant）、方向。

### Step 4 — 拷贝数估算 (核心交付物 2)
- `mosdepth` 计算：外源 contig 平均深度、全基因组常染色体中位深度、若干已知单拷贝小鼠区域深度。
- 拷贝数 ≈ (外源 contig 唯一比对深度) / (单拷贝小鼠基线深度)。
- **用野生型对照样本**校验基线；对人源/小鼠同源区做 MAPQ≥某阈值的唯一比对过滤，避免交叉比对虚高。

### Step 5 — 基因破坏 / 局部 SV 评估 (核心交付物 3)
- 整合断点坐标注释到 `gencode.vM35.annotation.gtf` → 是否落在某小鼠基因外显子/内含子/调控区。
- 结合点局部检查缺失/重复/倒位（Step 3B 的 SV 输出）——整合常伴随靶位点小片段缺失。

### Step 6 (可选) — 全基因组种系变异
- 若客户还想要常规 SNV/indel，可套用项目 4 的 **nf-core/sarek** 流程（germline，`--genome` 换成 GRCm39/自定义）。**注意：sarek 面向人类 iGenomes，跑小鼠需自定义 reference，且它本身不做整合定位——整合分析仍走上面的定制流程。** 属加分项，非核心目标。

---

## 4. 本机工具与参考现状

| 组件 | 状态 |
| :--- | :---: |
| GRCm39 (GENCODE M35) FASTA + bwa-mem2 索引 | ✅ 已就绪 `/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/` |
| GENCODE vM35 GTF (基因注释) | ✅ 已就绪 |
| GRCh38 (提取人源基因备用) | ✅ 已就绪 |
| bwa-mem2 / samtools / bcftools / blastn / minimap2 / seqkit / fastp | ✅ regular_bioinfo |
| **mosdepth / manta / delly** (深度+SV) | ❌ **需补装**（conda 一条命令即可） |
| nf-core/sarek 3.8.1 (可选 germline) | ✅ 项目 4 已装并验证 |

---

## 5. 执行前需客户/你确认的清单

1. **MTTH 到底是什么基因/构建体**（CD1A 已确认，MTTH 不在人类基因库中）—— 现头号 Blocker。
2. **外源构建体载体序列**（GenBank/质粒图谱/FASTA，含骨架）—— CD1A 组可先用 GRCh38 基因本体 bootstrap，但有骨架就能定位更准。
3. **RAGH 对照的遗传背景确认**（是否 RAG 相关；只要 transgene 阴性即可用作阴性对照/拷贝数基线）。
4. 是否需要 Step 6 的常规 germline SNV/indel（加分项）。
5. 交付物形式：整合位点表 + 拷贝数表 + 断点注释 + IGV 截图/可视化 + 英文报告。

---

*方案待你审阅确认。拿到外源序列后即可搭 Step 1–5 的定制流程并先补装 mosdepth/manta/delly。*
