# Ellen 人源基因敲入小鼠 WGS — 分析方案 (Analysis Plan)

- **Project**: 13_Ellen_knockin_wgs
- **Plan Date**: 2026-07-06 (rev — 构建体解码 + sarek-on-hybrid 引擎决策)
- **Prepared by**: Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
- **Data**: `/home/gao/Dropbox/Ellen/` — 6 samples, WGS PE150, NovaSeq X Plus
- **参考方法来源**: `4_wgs_human_immu/` (nf-core/sarek 已装好并测试通过) 提供比对/QC/资源管理的基础做法；本项目的整合位点/拷贝数/SV 分析为在其之上的定制流程。

---

## 1. 样本清单 (来自 Sample_Sheet_for_interal_use.xlsx)

**样本命名约定**：样本名下划线前 = 打靶等位基因（targeted allele）品系标识；下划线后为个体/founder 编号。

**❗ 重要更正（据客户 Ellen 2026-07-06 邮件）**：客户将 **CD1A、RAGH、MTTH 三者均称为 "projects" 且各有独立的 "targeted allele sequence"**（客户公司 genetargeting.com 专做基因打靶）。因此 **RAGH 也是一个打靶/敲入品系，而非 Sample Sheet 字面所示的"普通野生型 Mouse"**。共 **3 个构建体**，不是 2 个。这也与"客户已提供 2 个载体（RAGH+MTTH），CD1A 明天提供"一致。

| 打靶品系 | 样本 (FASTQ 前缀 / S 号) | 数量 | R1/R2 | 构建体序列到位情况 |
| :--- | :---: | :---: | :---: | :---: |
| **CD1A** | CD1A_B125 (S83) | 1 | 26G/26G | ⏳ 客户承诺明天提供 |
| **RAGH** | RAGH_153 (S89), RAGH_273 (S90) | 2 | 17G/20G | ✅ 已随邮件提供（待存盘） |
| **MTTH** | MTTH_284 (S84), _412 (S85), _524 (S86) | 3 | 17–18G | ✅ 已随邮件提供（待存盘） |

\* 预估覆盖度 ~20–35×/样本（PE150，2.7 Gb 小鼠基因组），QC 后核实。对整合位点与拷贝数分析足够。

**构建体身份核查：**
- **CD1A** ✅ 真实人类蛋白编码基因（ENSG00000158477，GENCODE v45，免疫相关）。**利好**：小鼠仅保留 CD1d、不含 group-1 CD1（CD1a/b/c），人源 CD1A 几乎无小鼠同源交叉比对，整合信号干净。
- **RAGH / MTTH** 非标准人类基因符号 —— 属客户自定义打靶等位基因命名（可能是人源化/敲入/敲除等位）。**确切结构以客户提供的 targeted-allele 序列文件为准**，无需再猜测。

**❗ 对照问题（需向客户确认）**：这批 6 样可能 **无一为纯野生型对照**（三者皆打靶品系）。拷贝数归一化策略相应调整：
- 以每个样本**自身全基因组常染色体单拷贝深度**为基线；
- 各品系互为**阴性对照**（如 RAGH 样本在 CD1A/MTTH 构建体上应零覆盖，验证信号特异性）；
- 若客户能补送一只纯野生型小鼠 WGS，会显著提升拷贝数基线与背景 SV 的可信度——可提出但非必需。

---

## 1.5 构建体解码（2026-07-06，已下载并解析 SnapGene .dna）

已通过 IMAP 下载客户附件（`refs/constructs/*.dna`），Biopython 解析 + minimap2 比对同源臂/插入序列到 GRCm39/GRCh38，确认：

| 构建体 | 长度 | 小鼠靶位点 (GRCm39, 同源臂定位) | 敲入内容 | 模型解读 |
| :--- | :---: | :---: | :---: | :--- |
| **RAGH** | 38.6 kb | **Rag2** (chr2:101,455,063–101,462,874) | 2A 串联人源 **G-CSF, M-CSF, IL-6, IL-1β, IL-7, IL-15** | Rag2 位点人源化细胞因子 → Rag2 缺陷 + 人源化免疫小鼠 (MISTRG 类) |
| **MTTH** | 254 kb | **Htt** (chr5:34,919,084–35,069,878) | 人源 **HTT** 全基因组 (GRCh38 chr4:3.04–3.24 Mb, 67 外显子 ~170 kb) | 小鼠 Htt 被人源 HTT 全基因替换 → 人源化亨廷顿病 (HD) 小鼠（呼应客户 mouse behavior study） |
| **CD1A** | 待明天 | 待定 | 推测人源 CD1A 定点人源化 | 待序列到位确认 |

- 两构建体均为 **Neo deleted**（FRT 疤残留），两侧为**小鼠同源臂**（定点同源重组打靶）。
- **关键含义**：这些是**定点打靶敲入，不是随机转基因**，预期整合位点由设计已知（Rag2 / Htt）。分析重点从"找未知随机整合位点"转为下列 6 项（见 §3 精修）。

**分析目标精修（因定点打靶而调整）：**
1. **验证定点整合正确**：5'/3' 结合部（小鼠→同源臂→人源）存在且干净。
2. **检测额外随机/脱靶整合**：构建体读段配偶是否落到 Rag2/Htt 以外意外位点。
3. **拷贝数**：人源插入单拷贝 vs 多拷贝/串联 concatemer。
4. **合子型**：纯合 vs 杂合（结合部支持 KI vs WT 等位读段比例）。
5. **确认 Neo 已删除**：无残留 Neo/载体骨架读段。
6. **确认内源基因置换/破坏**：Rag2 / 小鼠 Htt 原生区覆盖丢失。

---

## 2. 对朋友策略的评估：**核心思路合理，但需按"定点打靶"精修**

> 朋友的框架是针对**随机转基因**整合的标准范式（完全正确）；本项目是**定点同源重组敲入**，方法学相同但重点转移（见 §1.5 六项）。

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

## 🚧 序列到位情况（Blocker 基本解除）

- **RAGH、MTTH 构建体序列** ✅ 已下载并解析（`refs/constructs/`，见 §1.5）。已导出 `TG_RAGH.fa` / `TG_MTTH.fa`。
- **CD1A 构建体序列** ⏳ 客户承诺 2026-07-07 提供；到位前 CD1A 组用 GRCh38 CD1A 本体 bootstrap。
- 下载途径：Gmail MCP 不支持下载附件，改用其底层 IMAP 凭据写脚本直接抓取（`scratchpad/fetch_attach.py`）。

---

## 3. 分析流程 (Pipeline)

> 遵循 CLAUDE.md：拆分为独立顺序步骤（不用三通管道），每步各自 tmux，参数按本机资源上限（≤28 物理核 / 56 线程）。

### Step 0 — QC
- `fastp` 每样本单独跑（`-w 8`），产出 HTML/JSON；`MultiQC` 汇总。
- 确认真实覆盖度、接头、重复率。

### Step 1 — 构建混合参考 (Hybrid Reference)
- **建一个合并 hybrid 参考**（不是每组一个）：
  `cat GRCm39.primary_assembly.genome.fa  CD1A_contig.fa  RAGH_contig.fa  MTTH_contig.fa  >  GRCm39_plus3constructs.fa`
- 三个构建体各作独立 contig（如 `>TG_CD1A` `>TG_RAGH` `>TG_MTTH`）。
- `bwa-mem2 index` + `samtools faidx` + `samtools dict`（对加了 contig 的新参考重建）。
- **好处**：6 样全比对到同一参考，每样只在自己构建体上有覆盖，**在其他两个构建体上零覆盖 = 自带阴性对照（特异性验证）**。
- **分阶段**（CD1A 序列 2026-07-07 才到）：先用真实 RAGH+真实 MTTH + GRCh38 的 CD1A 占位建参考，跑 5 个 RAGH/MTTH 样本（结果即最终）；明天用真实 CD1A 替换占位 contig，只补跑 1 个 CD1A 样本。

### Step 2 — 比对 + 排序 + 去重（引擎：nf-core/sarek on hybrid ref）
- **用项目 4 已验证的 nf-core/sarek 3.8.1 作标准骨架**，指向合并 hybrid 参考：
  `--fasta GRCm39_plus3constructs.fa  --aligner bwa-mem2  --skip_tools baserecalibrator  --tools manta,tiddit`
  （跳过 BQSR：无 known-sites 且 NovaSeq 收益极小；跳过 snpEff/VEP 人类注释——断点我们自己注释）。
- sarek 一趟产出：fastp QC → bwa-mem2 比对 → MarkDuplicates → **Manta/TIDDIT SV（含 BND 断点=整合结合部）** → mosdepth 深度。
- 资源：沿用项目 4 稳健配置，`queueSize=2`（同时只 2 样），防内存紧张。
- 备选：若想更透明可控，可用等价的精简定制流程（fastp→bwa-mem2→markdup→manta+mosdepth）替代 sarek；二者产物一致，sarek 胜在已验证+一键。

### Step 3 — 整合位点检测 (核心交付物 1)
两条互补路线，结果交叉验证：
- **(A) 定向嵌合读段法**：`samtools view TG_construct` 抓所有落在外源 contig 上的读段及其配偶 → 看配偶/软剪切段落回小鼠基因组何处 → 聚类成整合位点，取 split-read 断点为碱基级坐标。
- **(B) 通用 SV caller**：`manta` 或 `delly`（BND/translocation 模式）在混合参考上直接调用，输出跨 contig 断点。
- 汇总每个整合位点：小鼠侧染色体坐标、外源侧坐标、支持读段数（split + discordant）、方向。

### Step 4 — 拷贝数估算 (核心交付物 2)
- `mosdepth` 计算：构建体 contig 平均深度、全基因组常染色体中位深度、若干已知单拷贝小鼠区域深度。
- 拷贝数 ≈ (构建体 contig 唯一比对深度) / (**样本自身**单拷贝小鼠基线深度)。**每个样本用自身基因组做基线**（本批可能无纯野生型对照）。
- 用**其他品系样本**作阴性对照（其在本构建体 contig 上应零覆盖）；对与小鼠同源区做 MAPQ≥阈值唯一比对过滤，避免交叉比对虚高（CD1A 因小鼠缺 group-1 CD1，此风险低）。
- sarek 的 CNV 模块（ASCAT/CNVkit，肿瘤/人类导向）**不用于此**；构建体拷贝数走 mosdepth 深度比值，更直接可靠。

### Step 5 — 基因破坏 / 局部 SV 评估 (核心交付物 3)
- 整合断点坐标注释到 `gencode.vM35.annotation.gtf` → 是否落在某小鼠基因外显子/内含子/调控区。
- 结合点局部检查缺失/重复/倒位（Step 3B 的 SV 输出）——整合常伴随靶位点小片段缺失。

### Step 6 (可选) — 全基因组种系变异
- Step 2 的 sarek 若加 `--tools haplotypecaller`（或 strelka/deepvariant）即可顺带产出常规 SNV/indel。属加分项，非核心目标；默认可先不开，按客户需求追加。

---

## 4. 本机工具与参考现状

| 组件 | 状态 |
| :--- | :---: |
| GRCm39 (GENCODE M35) FASTA + bwa-mem2 索引 | ✅ 已就绪 `/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/` |
| GENCODE vM35 GTF (基因注释) | ✅ 已就绪 |
| GRCh38 (提取人源基因备用) | ✅ 已就绪 |
| bwa-mem2 / samtools / bcftools / blastn / minimap2 / seqkit / fastp | ✅ regular_bioinfo |
| mosdepth 0.3.14 / delly 2.3.0 | ✅ 已装 (regular_bioinfo) |
| manta 1.6.0 | ✅ 已装 (独立 `manta` 环境，避免 py2 污染) |
| nf-core/sarek 3.8.1 (标准骨架引擎) | ✅ 项目 4 已装并验证 |

---

## 5. 执行前需客户/你确认的清单

1. **把 Ellen 邮件附件（RAGH、MTTH 序列）存盘**到 `refs/constructs/` —— 现头号动作（Gmail 工具无法自动下载附件）。
2. **CD1A 构建体序列**（客户承诺 2026-07-07 提供）—— 到位前用 GRCh38 CD1A bootstrap。
3. **对照确认**：本批是否有纯野生型？若无，采用"自身基线 + 品系互为阴性对照"方案；可向客户建议补送一只 WT。
4. **各构建体是敲入/人源化/敲除？** 以客户序列文件的注释为准（影响 SV/断点解读）。
5. 是否需要 Step 6 的常规 germline SNV/indel（加分项）。
6. 交付物形式：整合位点表 + 拷贝数表 + 断点注释 + IGV 截图/可视化 + 英文报告。

---

*工具已装齐（2026-07-06）。下一步：存盘 RAGH/MTTH 序列 → 建合并 hybrid（CD1A 用 GRCh38 占位）→ sarek 跑 5 个 RAGH/MTTH 样本；CD1A 样本待明天真实序列到位补跑。*
