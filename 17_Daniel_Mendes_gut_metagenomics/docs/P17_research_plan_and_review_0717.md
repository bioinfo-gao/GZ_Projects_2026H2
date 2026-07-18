# P17 — Daniel Mendes 小鼠肠道 shotgun metagenomics · Research Plan & Review

- **Project ID**: 17 (`17_Daniel_Mendes_gut_metagenomics`)
- **Plan Date（创建）**: 2026-07-17（文件名 MMDD=0717 为创建日，永不改名）
- **Client**: Daniel Mendes ｜ 送样单 `QTE_26_06_25_001_Daniel_Mendes`
- **Species / 材料**: *Mus musculus* (GRCm39) ｜ **Stool（粪便）** ｜ shotgun metagenomics
- **平台**: NovaSeq X Plus, PE150
- **状态**: ⏸ **执行已暂停，等待客户/用户 review 本文档后批准再跑**（脚本已就位、已 dry-check，未启动任何计算）

## 更新记录 / change-log
- 2026-07-17 — 初稿：完成数据勘察、路线抉择、资源规划、脚本落地；一度启动 Phase 1 预热（GTDB 下载/建索引），随即按用户要求**全部停止并清理**，改为先出完整思考文档待批。

---

## 0. 本文档目的

把"该走 assembly 依赖路线还是 mapping 为主路线"这个抉择背后的**完整思考过程、备选方案、取舍依据、风险与待定问题**摊开，供 review。批准后再执行第 5 节的脚本。下面第 2 节是核心决策论证，第 8 节是需要你拍板的 open questions。

---

## 1. 样品与数据量（实测 on-disk，2026-07-17）

数据源：`/home/gao/Dropbox/QTE_26_06_25_001_Daniel_Mendes/HFD_*`（10 个样品目录，PE150 双端）。
研究设计：**case-control**——HFD 背景下 **AL（ad libitum 自由采食）vs IF（intermittent fasting 间歇禁食）**，两臂各 5 个生物学重复。

| Sample | Group | R1 (GiB) | R2 (GiB) |
| :--- | :---: | :---: | :---: |
| HFD_AL_4_02_25 | AL | 1.5 | 1.5 |
| HFD_AL_4_03_11 | AL | 1.5 | 1.5 |
| HFD_AL_6_05_12 | AL | 1.6 | 1.5 |
| HFD_AL_6_05_22 | AL | 1.7 | 1.7 |
| HFD_AL_7_06_12 | AL | 1.2 | 1.2 |
| HFD_IF_4_02_25 | IF | 1.3 | 1.4 |
| HFD_IF_4_03_11 | IF | 1.1 | 1.1 |
| HFD_IF_6_05_12 | IF | 1.5 | 1.5 |
| HFD_IF_6_05_22 | IF | 1.6 | 1.6 |
| HFD_IF_7_06_12 | IF | 1.5 | 1.5 |

- **数据集总量**：10 samples，gzip FASTQ 合计 **≈ 28 GiB**（`du -shc HFD_*` → 28G）。
- 每样品 PE150 ≈ 1.1–1.7 GiB/端 → **粗估 ~13M read pairs、~4 Gbp/样品**（gzip 反推，非精确；准确值由 fastp 运行时确认后回填）。
- 送样单 `Analysis Requirement` 字段：**"Shot gun meta std analysis"**（标准 shotgun 宏基因组分析）。这一点对范围界定很关键，见第 8 节。

## 2. 核心抉择：assembly-free 优先，assembly-based MAG 作二期（论证）

### 2.1 先问"客户到底要回答什么"
本项目是干预对照（IF vs AL）下的**肠道菌群比较**。这类研究的主问题几乎总是三件事：
1. **组成**——两臂菌群物种构成有无差异，哪些菌 IF 富集/耗竭；
2. **多样性**——alpha（组内丰富度/均匀度）、beta（组间群落结构分离）；
3. **功能**——代谢通路/基因家族丰度差异（禁食最受关注的正是代谢层面）。

这三件事，**assembly-free 的 reference-based profiling 直接、稳健地回答**，不需要先拼基因组。

### 2.2 两条路线各自能/不能回答什么

| 维度 | Assembly-free（Kraken2/Bracken + MetaPhlAn + HUMAnN） | Assembly-based（MAG：组装+分箱） |
| :--- | :--- | :--- |
| 回答 | 有哪些**已知**物种、丰度、功能通路差异 | 能否拼出近完整**基因组草图**、它们是谁、编码哪些基因 |
| 与主问题匹配 | **直接命中**组成/多样性/功能三问 | 补充——基因组层面证据、发现新/未培养菌 |
| 对深度的要求 | 低即可（~4 Gbp/样品足够定量常见菌） | 偏高；per-sample ~4 Gbp 对高质量 MAG 偏薄 |
| 参考依赖 | 依赖库；小鼠肠道 reference-rich，覆盖好 | 不依赖库，能抓库里没有的菌 |
| 速度/成本 | 轻、快、当天出主结果 | 重（GTDB-Tk r226 ~102GB、组装/分箱耗时耗内存） |
| 交付确定性 | 高（成熟标准流程） | 中（MAG 数量/质量随样本复杂度波动） |

### 2.3 深度这一条为什么决定"MAG 别当入口"
经验阈值：per-sample **≥10 Gbp** 才较稳地拿到中低丰度菌的高质量 MAG；本数据 ~4 Gbp/样品明显偏薄。若逐样本组装+分箱，多半只能拿到少数高丰度菌的 bin，且完整度/污染度参差。**对策不是放弃 MAG，而是 group co-assembly**：把同臂 5 个样本合并（AL、IF 各 ~20 Gbp），提升低丰度菌 contig 连续性与跨样本 co-abundance 分箱信号——但这天然是"二期增值"定位，不是回答主问题的入口。

### 2.4 "std analysis" 的范围含义
送样单写的是 **std analysis**。业界"标准 shotgun 宏基因组分析"通常 = **taxonomy + diversity + functional profiling（assembly-free）**；**MAG 重建属于 advanced/增值**，一般单独报价。因此：
- **Phase 1（assembly-free）本身就完整交付了客户下单的 std analysis。**
- **Phase 2（MAG）是超出 std 的加做**——我把它排进来，一是有科学价值（基因组级证据），二是正好吃满机器（你要求"充分利用"），但它要额外拉 102GB 库 + 长时算。**是否真的要做 MAG，请见第 8 节 open question，由你/客户拍板。**

### 2.5 结论
**Phase 1 先 assembly-free 出主交付**（直接回答 AL vs IF 的组成/多样性/功能），**Phase 2 视批准再 assembly-based MAG**（group co-assembly 重建基因组、GTDB 分类、跨臂丰度）。两者互补、非二选一。这也符合 `/taxnom` 与 `/tax-resemb-mag` 两个 skill 的选型边界（先看组成→需要基因组再上组装）。

## 3. 备选方案与为何不采（Alternatives considered）

| 方案 | 评价 | 采纳? |
| :--- | :--- | :---: |
| **A. 只做 MAG（assembly-first 为入口）** | 深度偏薄、库重、耗时，且回答不了"相对丰度/多样性/功能"主问题；本末倒置 | ✗ |
| **B. 只做 assembly-free、不做 MAG** | 完整满足 std analysis；但放弃基因组级证据，也没吃满机器 | 作为**保底**（若客户不要 MAG） |
| **C. assembly-free 优先 + MAG group co-assembly 二期**（本plan） | 主问题稳交付 + 基因组增值 + 机器吃满 | ✓ |
| **D. MAG 用 per-sample metaSPAdes** | 深度偏薄下逐样本 SPAdes 质量有限，且 SPAdes 内存易撞 125GB 墙 | ✗（改 co-assembly + MEGAHIT） |
| **E. Kraken2 换更大库（Standard-full/nt）** | 本机已有 Standard-8GB，够小鼠肠道；大库要 100+GB 内存与下载，性价比低 | ✗（保 8GB，必要时再议） |
| **F. 单一 classifier** | 单工具易受库偏差影响；双工具（Kraken2+MetaPhlAn）交叉验证更稳 | ✗（坚持双工具） |

## 4. 分析设计与每个决策的理由

### Phase 1 — Taxonomic + Functional profiling（assembly-free）
- **Pipeline**: nf-core/taxprofiler 2.0.1（`taxprofiler` env 只调度，process 走 singularity 容器）。理由：已在 `8_taxprofiler_setup` 全 classifier 验证通过。
- **QC + 去宿主**: fastp；Bowtie2 vs mouse GRCm39。粪便是宿主源样本，**必须去宿主**去掉小鼠 read 再分类，否则丰度被宿主污染。→ 预建共享 Bowtie2 索引复用。
- **Classifier 双工具交叉验证**:
  - **Kraken2 + Bracken**（Standard-8GB，属/种级丰度；Bracken `-r 150` 匹配 PE150 的 kmer_distrib）。
  - **MetaPhlAn**（CHOCOPhlAn SGB vJan25，marker-gene 种级）。两法原理不同，一致性高则结论稳，冲突处报告需标注。
- **合并/可视化**: taxpasta 跨样品标准化合并 + Krona。
- **功能通路**: taxprofiler **不含**功能模块 → 用 `mag_biobakery` env 原生 **HUMAnN 3.9**（ChocoPhlAn 16G + UniRef 34G）出 pathway/gene-family 丰度。禁食研究最关心代谢功能，这步不能省。

### Phase 2 — MAG recovery（assembly-based，待批）
- **Pipeline**: nf-core/mag 5.4.2（绝对路径 nextflow ≥25.04.2 + JAVA_HOME，避 PATH 版本遮蔽坑）。
- **组装**: `--coassemble_group`（AL、IF 各一个 co-assembly，~20 Gbp/组）+ **MEGAHIT**，`--skip_spades`。理由：125GB 内存是本机唯一硬墙，co-assembly 级 metaSPAdes 常要 ~200GB 会 OOM；MEGAHIT 内存友好、co-assembly 首选。
- **分箱**: MetaBAT2 + MaxBin2 + SemiBin2 → **DAS Tool** 精炼；关掉较慢的 CONCOCT/COMEBin/MetaBinner（小项目性价比低）。
- **质控**: BUSCO（本地 bacteria_odb10）+ CheckM2（需预下载库）；高污染(>10%)bin 视为潜在嵌合，不当干净基因组交付。
- **分类**: GTDB-Tk r226，**split-tree（非 full_tree）** 省内存；需预下载 ~102GB 库。
- **注释**: Prodigal/Prokka 基因预测。

### 下游统计（两 Phase 共用）
- 相对丰度堆叠图（属/种 Top N + Other）。
- **alpha 多样性**（Shannon/Simpson/richness）AL vs IF 箱线图 + 组间检验。
- **beta 多样性**（Bray-Curtis + PCoA），量化组内 vs 组间距离比（排序图**必须配数值解读**，不只放图）。
- **差异丰度**（Maaslin2 / ALDEx2 / Wilcoxon）找 IF 富集/耗竭 taxa 与 pathway。
- **跨工具交叉验证**：Kraken2/Bracken vs MetaPhlAn Top 物种一致性（Spearman/并排表）。

## 5. 执行编排（充分利用机器 · 守 28 核/56 线程 · 125GB 内存硬墙）

批准后按此并行预热 + 串行主算：

| 阶段 | tmux | 资源 | 时机 |
| :--- | :---: | :---: | :--- |
| 预下载 GTDB-Tk + CheckM2 库 | `db17` | 网络/IO | 与 Phase 1 并行（仅 MAG 获批才需要） |
| 预建 mouse Bowtie2 去宿主索引 | `idx17` | 20 threads | 立即，存共享参考目录复用 |
| Phase 1 taxprofiler | `tax17` | 24 核/96GB | 索引就绪即启动 |
| Phase 1b HUMAnN | 后续 | 分批 | taxprofiler 完成后 |
| Phase 2 MAG（待批） | `mag17` | 24 核/110GB | Phase 1 完成后 |

- taxprofiler 与 MAG **不同时抢核**（MAG 排在其后）；预热（下载/建索引）与 Phase 1 属轻量重叠，短时可接受。
- 全程 tmux + Phase-1 前 3 分钟早失败监控 + Phase-2 常驻看门狗。

## 6. 脚本清单（`scripts/`，编号即执行顺序 · 已 `bash -n` 通过 · 20/20 fastq 路径校验通过）

| 脚本 | 作用 |
| :--- | :--- |
| `0_produce_inputs.py` | 生成 taxprofiler samplesheet + databases.csv + mag samplesheet（已跑，10 样品无误） |
| `1_predownload_dbs.sh` | 预下载 GTDB-Tk r226 + CheckM2 库（MAG 依赖） |
| `2_prebuild_host_index.sh` | 预建 mouse GRCm39 Bowtie2 去宿主索引 |
| `3_run_taxprofiler.sh` | Phase 1 taxprofiler（self-relaunch tmux + resume） |
| `3b_wait_then_run_taxprofiler.sh` | 等索引就绪再启动 taxprofiler |
| `4_run_mag.sh` | Phase 2 nf-core/mag（co-assembly + 分箱 + 质控 + GTDB-Tk） |
| `5_monitor.sh` | 各 tmux 状态 + 输出目录 + 负载快照 |

## 7. 数据库与环境（就位情况）

| 资源 | 状态 | 路径 |
| :--- | :---: | :--- |
| taxprofiler env / nf-core/taxprofiler 2.0.1 | ✅ | `taxprofiler` conda env |
| Kraken2 Standard-8GB | ✅ | `/Work_bio/references/Metagenomics/kraken2/k2_standard_08gb_20260226/` |
| MetaPhlAn CHOCOPhlAn SGB vJan25 | ✅ | `/Work_bio/references/Metagenomics/metaphlan/` |
| HUMAnN ChocoPhlAn + UniRef | ✅ | `/Work_bio/references/Metagenomics/humann/` |
| mouse GRCm39 host fasta | ✅ | `/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/` |
| mouse Bowtie2 去宿主索引 | ❌ 待建 | `.../mouse_gencode_M35/bowtie2_index/`（脚本 2） |
| mag_biobakery env / nf-core/mag 5.4.2 | ✅ | `mag_biobakery` conda env（nextflow 26.04.4） |
| BUSCO bacteria_odb10 | ✅ | `/Work_bio/references/Metagenomics/busco/` |
| GTDB-Tk r226 (~102GB) | ❌ 待下载 | `/Work_bio/references/Metagenomics/gtdbtk/release226/`（仅 MAG 获批才拉） |
| CheckM2 DIAMOND DB | ❌ 待下载 | `/Work_bio/references/Metagenomics/checkm2/`（仅 MAG 获批才拉） |

## 8. 待你拍板的问题（Open questions）

1. **MAG 二期做不做？** 送样单只写了 std analysis；assembly-free（Phase 1）已完整交付 std。Phase 2 MAG 要额外拉 **~102GB GTDB 库**并长时算，属超出 std 的增值。选项：
   - (a) 只做 Phase 1（std，保底、最快）；
   - (b) Phase 1 + Phase 2 MAG（本 plan 默认，机器吃满、基因组级证据）。
2. **实验设计是独立两组还是配对？** 两臂重复编号（4_02_25 / 4_03_11 / 6_05_12 / 6_05_22 / 7_06_12）在 AL 与 IF **一一对应**。若这是 litter/cage/批次配对（而非同一只鼠——鼠不可能既 AL 又 IF），差异分析应把它当 **blocking factor / paired 设计**（统计功效更高）。请确认这些编号的含义。默认按独立两组分析，若确认配对则改用配对模型。
3. **功能通路 HUMAnN 是否要做？** 我默认做（禁食研究重代谢，几乎必需）。若客户只要物种组成可省，能省不少算力。
4. **参考库版本**：Kraken2 用本机 Standard-8GB（够小鼠肠道常见菌）。若客户要更高灵敏度可换更大库，但要额外下载+内存。默认 8GB。

## 9. 交付（endpoint = delivery folder）
- `custom_research_report_YYYYMMDD/`：英文 report（签名 Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics）+ 分类/功能/多样性/(MAG) 结果表 + 图。
- 子目录：`qc/`、`taxonomy/`、`function/`、`diversity/`、`mag/`。
- Report header 注明 Species=*Mus musculus* (GRCm39) + Tissue=Stool。
