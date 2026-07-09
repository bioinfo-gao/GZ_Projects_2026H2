# nf-core/taxprofiler 2.0.1 — 完整操作步骤 · Nextflow 内部调用流程 · 参数详解

**记录日期**：2026-07-08
**核对方式**：直接阅读本机已 `nextflow pull` 的 pipeline 源码
（`~/.nextflow/assets/nf-core/taxprofiler`，对应 release `2.0.1`），逐个子工作流/`conf/modules.config`
核对参数拼接逻辑，不是转述官方文档大意。

**关联文档**：
- 安装 + test profile 验证记录（怎么装好、踩过什么坑）：`taxprofiler_安装与测试记录_0701.md`
- 英文教程（安装步骤对应版）：`README_taxprofiler_tutorial.md`
- 可复用操作 skill（信息收集/目录结构/脚本模板/交付前自查）：`/taxnom`（`~/.claude/commands/taxnom.md`）

**本文档定位**：上面两篇文档记录的是"怎么装好、怎么跑通官方 test profile"；这篇聚焦 pipeline **内部**
到底做了什么——从一条 `nextflow run` 命令发出开始，到调度器实际按什么顺序调用哪些 process、每个
process 收到的命令行参数具体是什么、为什么这样设计，逐层拆解，供以后设计真实项目参数时有依据可查，
不用每次都重新翻源码。

---

## 目录

1. 我该操作的步骤（外部视角 checklist）
2. Nextflow 内部调用流程（按真实 DAG 顺序逐步展开）
3. 参数详解与选择依据（分类别表格）
4. `8_taxprofiler_setup` test run 实际用的参数 vs 真实项目该怎么设
5. 快速参考：命令行模板

---

## 1. 我该操作的步骤

真实项目从零到出结果，我这边需要做的事情，按顺序：

| 步骤 | 做什么 | 产出/确认点 |
|:---|:---|:---:|
| ① 收集样品信息 | 确认短读长/长读长、样品来源（是否需要去宿主）、目标 classifier | 见 `/taxnom` 第一步 |
| ② 写 `samplesheet.csv` | 固定列 `sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta`，短读长填 fastq_1/2，长读长填 fasta，不能混填 | 每样品至少一行，路径存在 |
| ③ 写 `databases.csv` | 固定列 `tool,db_name,db_params,db_type,db_path`；本地已有 Kraken2 Standard-8GB / MetaPhlAn CHOCOPhlAn SGB 库直接可用 | 见本文档第 3 节 db_params 语法 |
| ④ 决定要开哪些 `--perform_*`/`--run_*` 开关 | 每个都默认关闭，必须显式打开；`--run_<tool>` 要和 databases.csv 里的 `tool` 行对应，否则该工具静默不跑 | 见本文档第 3 节参数表 |
| ⑤ 若要去宿主，先确认/预建宿主索引 | 首次用某宿主可以让 pipeline 现场建 Bowtie2/minimap2 索引（慢，一次性）；反复用的宿主建议预建一次长期复用 | `--shortread_hostremoval_index`/`--longread_hostremoval_index` |
| ⑥ tmux 内启动（`2_run_taxprofiler.sh`，self-relaunch + `-resume`） | `nextflow run nf-core/taxprofiler -r 2.0.1 -profile singularity -c local_resources.config ...` | tmux session 已建，日志开始写入 |
| ⑦ 两阶段监控 | 前 3 分钟 30 秒轮询查错；健康后 10–20 分钟轮询；>15 分钟无新日志查 `.nextflow.log` 是否反复打印 `No more task to compute` | `Succeeded/Cached/Failed` 统计 |
| ⑧ 完成后核对 `output_results/multiqc/multiqc_report.html` | 各样品 QC/去宿主/分类是否有异常 | 见 `/taxnom` 交付前自查 |
| ⑨ 下游丰度/多样性分析 + 报告 | 读 `output_results/taxpasta/` 合并表 | `4_downstream_analysis.R` / `5_generate_report.R` |

---

## 2. Nextflow 内部调用流程

以下按 `workflows/taxprofiler.nf` 里真实的执行顺序（不是我猜的，是源码 `main:` 块里 include 的调用顺序）逐段展开。每段说明：**触发条件**（对应哪个 `--param`）、**实际调用的 process/工具**、**命令行参数怎么拼出来的**、**为什么**。

### Stage 0 — 参数校验（nf-schema 插件）

启动时先用 `nextflow_schema.json` 校验 `samplesheet.csv`、`databases.csv` 和所有 CLI 参数。**关键行为**：schema 里不认识的参数名会被**静默丢弃**——不报错、不出现在打印的参数汇总里（这是 `8_taxprofiler_setup` 项目实测踩过的坑，如 `--max_cpus`/`--remove_host` 这类旧模板/臆测参数名）。任何新参数上线前建议先 `nextflow run ... -preview` 核对它确实出现在参数汇总里。

### Stage 1 — 输入分流（`workflows/taxprofiler.nf:87-102`）

`samplesheet.csv` 每一行按 `instrument_platform` + 是否填了 `fasta` 列，分流成 5 路 channel：
`fastq`（短读长）、`nanopore`、`pacbio`、`fasta_short`、`fasta_long`。**这一步决定了后面走短读长处理链还是长读长处理链**——FASTA 输入（不管长短读长哪种平台标注）会直接跳过 QC/host removal 一路到 profiling，因为组装后的 contig 本身已经是"干净"序列。

### Stage 2 — 数据库准备（`databases.csv` → UNTAR）

- 每行数据库先按 `db_type` 拆分（未填默认 `short;long`，即同时供短读长和长读长使用）。
- **只有对应 `--run_<tool>` 确实打开的数据库才会进入下一步**（`params["run_${db_meta.tool}"]` 过滤），未启用工具的数据库行即使写在 CSV 里也不会被解压/使用——不用手动删掉不用的行。
- `db_path` 以 `.tar.gz` 结尾的才走 `UNTAR`；本地已有的 Kraken2 tar 包、或直接是目录的 MetaPhlAn 库，处理方式不同但结果一致。

### Stage 3 — 原始 reads QC（FastQC / Falco）

`--skip_preprocessing_qc`（默认不跳过）+ `--preprocessing_qc_tool`（默认 `fastqc`，可选更快的 `falco`）。**在任何清洗之前**先跑一次，用于 MultiQC 里的 raw vs processed 对比。

### Stage 4 — 短读长预处理 `SHORTREAD_PREPROCESSING`（`--perform_shortread_qc`）

工具由 `--shortread_qc_tool` 选（默认 `fastp`，备选 `adapterremoval`）。以 fastp 为例，实际拼给容器的命令行（`conf/modules.config` FASTP_PAIRED 块）：

```
--length_required {shortread_qc_minlength}          # 默认 15
[--detect_adapter_for_pe  |  --adapter_sequence.. --adapter_sequence_r2..]  # 未指定 adapter 时 PE 默认自动探测
[--include_unmerged]                                 # 仅 shortread_qc_mergepairs 开且 includeunmerged 也开时
[--disable_adapter_trimming]                         # 仅 skipadaptertrim 开时
[--low_complexity_filter --complexity_threshold 30]  # 仅 complexityfilter_tool=fastp 时（非默认工具）
[--dedup]                                             # 仅 shortread_qc_dedup 开时
```

**为什么用 fastp 而不是 adapterremoval**：fastp 是 pipeline 默认值，速度快、能一步完成 adapter trim + 质量过滤 + 可选的低复杂度过滤 + 可选去重，产出的 JSON/HTML 也是 MultiQC 直接识别的标准格式；adapterremoval 是备选（更贴近 aDNA/古菌群落分析习惯，有独立的 `--collapse` 序列合并逻辑）。

清洗完成后再跑一次 FASTQC/Falco（`FASTQC_PROCESSED`），供 MultiQC 里 before/after 对比。

### Stage 5 — 长读长预处理 `LONGREAD_PREPROCESSING`（`--perform_longread_qc`）

- **Adapter 去除**：porechop 或 porechop_abi（`--longread_qc_predictadapters` 开启时用 porechop_abi 的 ab-initio 模式直接从 reads 里推断 adapter，不依赖已知库——适合用了非标准/自定义 barcode 的项目）。
- **长度/质量过滤**：filtlong 或 nanoq（默认 nanoq）。
  - `filtlong`: `--min_length {minlength=1000} --keep_percent {keeppercent=90} --target_bases {targetbases=5e8}`
  - `nanoq`: `--min-len {minlength=1000} --min-qual {minquality=7}`
- `--longread_qc_skipadaptertrim`/`--longread_qc_skipqualityfilter` 两个开关组合决定走"只 trim"/"只 filter"/"两步都做"哪条子路径。

**为什么 minlength 默认 1000bp**：短于 1000bp 的 Nanopore reads 信息量往往不足以可靠定位到种/株级别，是官方给的经验阈值，一般项目不需要改；若做病毒宏基因组（基因组本身就短）可能需要调低。

### Stage 6 — 测序冗余度估计 `NONPAREIL`（`--perform_shortread_redundancyestimation`）

仅短读长。**跟物种分类本身无关**，是回答"这批测序深度是否已经把群落测饱和、还要不要加测"的独立 QC 指标（Nonpareil curve）。默认关闭，项目讨论"是否测够深"时才开。

### Stage 7 — 复杂度过滤 `SHORTREAD_COMPLEXITYFILTERING`（`--perform_shortread_complexityfilter`）

**互斥设计**：若 `--shortread_complexityfilter_tool=fastp`（即用 fastp 自带的低复杂度过滤），这一步会被整个跳过——因为 Stage 4 的 fastp 调用里已经顺带做了，避免同一件事做两遍。只有工具选 `bbduk`（默认）或 `prinseqplusplus` 时才作为独立 process 跑：

- `bbduk`: `entropy=0.3 entropywindow=50 entropymask=f`（默认**丢弃**低复杂度 reads，而非仅打 mask 标记）
- `prinseqplusplus`: entropy 或 dust 两种打分模式（`--shortread_complexityfilter_prinseqplusplus_mode`）

**为什么要做这一步**：低复杂度序列（poly-A/poly-G、简单重复）会被 Kraken2 等 k-mer 工具误判成"能唯一比对"从而产生假阳性物种命中，宏基因组分析里是标准前处理步骤，不建议跳过。

### Stage 8 — 短读长去宿主 `SHORTREAD_HOSTREMOVAL`（`--perform_shortread_hostremoval`）

**前提**：必须同时给 `--hostremoval_reference <fasta>`，只开开关是 no-op（`8_taxprofiler_setup` 已验证的坑）。

1. 若未给 `--shortread_hostremoval_index`，现场 `BOWTIE2_BUILD` 建索引（这一步每个项目都要重新做一次，耗时且吃内存——反复用同一宿主的项目建议预建一次长期复用）。
2. `BOWTIE2_ALIGN`：比对到宿主参考，同时输出 (a) 未比对上的 reads（FASTQ，供下游分类用）(b) 全部比对结果 BAM。
3. `SAMTOOLS_INDEX` + `SAMTOOLS_STATS`：统计比对率——**这个比对率就是"去宿主比例"的直接依据**，宿主源样品该比例应该在合理区间（因样品类型而异），交付前必须核对。

### Stage 9 — 长读长去宿主 `LONGREAD_HOSTREMOVAL`（`--perform_longread_hostremoval`）

同理，用 minimap2 代替 Bowtie2，`--longread_hostremoval_index` 对应预建的 `.mmi` 文件。

### Stage 10 — Run merging（`--perform_runmerging`，可选）

同一 `sample` 下若有多个 `run_accession`（比如一个样品分了两条 lane 测），在 **QC + 去宿主之后**才合并成一份 reads。**关键点是顺序**：先各自独立 QC/去宿主，再合并——这样每个 run 各自的 QC 质量仍可单独追溯（不会因为提前合并而混在一起看不出是哪个 run 的问题）。若两个 run 的 single/paired-end 类型不一致又没开 `shortread_qc_mergepairs`，会各自独立跑完 profiling，文件名加 `_se`/`_pe` 后缀区分，不强行合并。

### Stage 11 — Profiling 主体（`subworkflows/local/profiling`，核心）

这是唯一把"清洗好的 reads"和"databases.csv 里每一行"做笛卡尔组合、真正调用各分类算法的子工作流。三个关键设计点：

**(a) 读长类型必须匹配数据库的 `db_type`。** 例如某个 Kraken2 库标了 `db_type=short`，长读长样品即使 `--run_kraken2` 开着也不会被送进这个库跑（`db_type` 和 `meta.type` 的 join 逻辑，`profiling/main.nf:47-69`）。

**(b) 几乎所有分类器的实际算法参数，直接等于 `databases.csv` 那一行的 `db_params` 字段**（源码里遍布 `ext.args = "${meta.db_params}"` 这种写法：Kraken2/Centrifuge/MetaPhlAn/KrakenUniq 等全是这个模式）。**这是本 pipeline 参数设计的核心思想**：Nextflow CLI 层面的 `--run_<tool>` 只是"要不要跑这个工具"的开关，**真正的算法调参责任下放到 `databases.csv` 的 `db_params` 列**（比如 MetaPhlAn 要不要用 `--mpa3` 老版本数据库格式、Kraken2 要不要 `--quick` 模式、DIAMOND 要不要 `--long-reads`，全部写在这一列，不是 nextflow 参数）。

**(c) Bracken / KMCP / sylph 三个工具是两步流程**，`db_params` 用分号拆成两段，源码精确 `.split(";")` 取用：
- Bracken：`<kraken2的参数>;<bracken自己的参数(如 -r 150)>`——第一段在 Stage 11 里先跑 Kraken2，第二段在 Kraken2 报告出来后单独喂给 `BRACKEN_BRACKEN`。
- KMCP：`<search参数>;<profile参数>`。
- sylph：`<profile参数>;<taxprof参数>`。

**(d) 各分类器有硬编码的输入类型保护**（`profiling/main.nf` 里散落的 `log.warn` + `filter`），不需要我自己记住哪些组合会出问题——pipeline 自己会跳过并警告：
- Centrifuge / sylph / mOTUs：不接受 FASTA 输入，遇到会跳过该样品并 `log.warn`。
- Bracken / KMCP / Ganon：长读长（Nanopore/PacBio）样品会被跳过（未在长读长上验证/灵敏度低）。
- DIAMOND：不接受双端，双端样品只取 R1 跑，并打 warning。
- Melon：只接受长读长，短读长样品直接跳过。

Kraken2+Bracken 联动细节：Bracken 需要先有 Kraken2 的 `.report.txt`，所以 `run_bracken` 单独开着已经隐含跑一次 Kraken2（即使 `--run_kraken2` 没开，Bracken 也能跑通），两者同时开不冲突、也不会重复计算——**这也是 `--run_bracken` 官方 schema 描述"自动带出 Kraken2 前置步骤"的真实含义**。

### Stage 12 — Krona 可视化 `VISUALIZATION_KRONA`（`--run_krona`）

只有部分工具的输出格式支持转 Krona 交互图：**Kraken2（含 Bracken 重新分配后的结果）、Centrifuge、Kaiju** 走 `KRAKENTOOLS_KREPORT2KRONA`/`KAIJU_KAIJU2KRONA` + `KRONA_KTIMPORTTEXT`，`--run_krona` 打开即可直接出图；**MALT** 额外需要同时给 `--krona_taxonomy_directory <NCBI taxonomy 目录>`（走 `MEGAN_RMA2INFO_KRONA` + `KRONA_KTIMPORTTAXONOMY` 这条独立支线），只开 `--run_krona` 而不给这个目录时 MALT 不会出 Krona 图。其余分类器（MetaPhlAn、Bracken 原始输出、KMCP、Ganon、sylph、Melon、MetaCache、mOTUs、KrakenUniq）没有对应 Krona 转换逻辑，即使 `--run_krona` 打开也不会为它们生成图。

### Stage 13 — 标准化合并 `STANDARDISATION_PROFILES`（`--run_profile_standardisation`）

用 [taxpasta](https://taxpasta.readthedocs.io/) 把每个分类器格式不一的原始输出，转成统一的"taxonomy ID + 每样品 read 计数"表格：
- `TAXPASTA_STANDARDISE`：单个数据库单独一份标准化表。
- `TAXPASTA_MERGE`：同一工具跨样品合并成一张宽表（跨工具**不**合并，Kraken2 一份表、MetaPhlAn 一份表，是本 skill 交付前要求"跨工具交叉验证"的直接原因——两个工具的结果天然是分开的两张表，不会自动被这条 pipeline 合成一个"共识"结果，需要下游脚本自己做比对）。
- `--taxpasta_add_name/rank/lineage/idlineage/ranklineage`：默认全关，只输出 taxonomy ID + 计数；打开这些需要额外给 `--taxpasta_taxonomy_dir <NCBI taxdump 目录>`，会让表格直接带上物种名/分类阶元/完整谱系，省去下游脚本自己再查 taxonomy ID 的步骤（推荐真实项目打开，方便下游直接按物种名画图）。

### Stage 14 — MultiQC 汇总 + 版本/流程信息落盘

汇总所有前面步骤产出的 QC 文件（FastQC/fastp/Bowtie2/各分类器统计等）成一份 `multiqc_report.html`；同时收集所有 process 实际用到的软件版本写入 `pipeline_info/`，供报告 Methods 部分核对具体版本号。

---

## 3. 参数详解与选择依据

### 3.1 短读长 QC（`shortread_qc_*`）

| 参数 | 默认值 | 说明 / 为什么 |
|:---|:---:|:---|
| `--perform_shortread_qc` | 关 | **必须显式打开**才会跑 fastp/adapterremoval |
| `--shortread_qc_tool` | `fastp` | 见 Stage 4；一般项目不需要改成 adapterremoval |
| `--shortread_qc_minlength` | `15` | 过短的 read 分类学信噪比差；15bp 是官方偏保守的下限，一般不用改 |
| `--shortread_qc_mergepairs` | 关 | overlap 区域长的建库（如短片段 aDNA）可开，普通 shotgun 一般不需要 |
| `--shortread_qc_dedup` | 关 | PCR 重复去除；宏基因组丰度定量一般**不需要**开（不是变异检测，PCR 重复对相对丰度估计影响有限，开了反而可能损失低丰度物种的信号） |
| `--shortread_qc_skipadaptertrim` | 关 | 只有明确知道数据已经去过接头才跳过 |

### 3.2 短读长复杂度过滤（`shortread_complexityfilter_*`）

| 参数 | 默认值 | 说明 |
|:---|:---:|:---|
| `--perform_shortread_complexityfilter` | 关 | 建议真实项目打开（见 Stage 7 理由：防低复杂度序列产生 k-mer 假阳性） |
| `--shortread_complexityfilter_tool` | `bbduk` | 默认工具，速度快 |
| `--shortread_complexityfilter_entropy` | `0.3` | 熵值越低过滤越严格；0.3 是 bbduk 官方推荐的宏基因组默认阈值 |

### 3.3 长读长 QC（`longread_qc_*`）

| 参数 | 默认值 | 说明 |
|:---|:---:|:---|
| `--perform_longread_qc` | 关 | Nanopore/PacBio 项目必须打开 |
| `--longread_qc_qualityfilter_minlength` | `1000` | 见 Stage 5；短片段/病毒项目可能要调低 |
| `--longread_qc_qualityfilter_minquality` | `7` | nanoq 专属，Q7 约等于 80% 碱基准确率，Nanopore 常规过滤阈值 |

### 3.4 去宿主（`*hostremoval*`）

| 参数 | 默认值 | 说明 |
|:---|:---:|:---|
| `--perform_shortread_hostremoval` / `--perform_longread_hostremoval` | 关 | 宿主源样品（组织/粪便/拭子）必须打开；环境样品（水/土）通常不需要 |
| `--hostremoval_reference` | 无 | **两个开关的必需搭档**，只给开关不给参考基因组是 no-op |
| `--shortread_hostremoval_index` / `--longread_hostremoval_index` | 无 | 给了就跳过现场建索引，直接用预建索引——反复用同一宿主的项目强烈建议预建复用 |
| `--save_hostremoval_unmapped` / `--save_hostremoval_bam` | 关 | 需要留存去宿主后的 reads/比对 BAM 供复查时打开 |

### 3.5 分类器开关与 `databases.csv` 的 `db_params` 语法

| 参数 | 默认值 | 说明 |
|:---|:---:|:---|
| `--run_<tool>`（14 个，见下） | **全部关** | 必须显式打开；且要和 `databases.csv` 里的 `tool` 行对应，否则静默不跑 |
| `db_params`（databases.csv 第 3 列） | 空 | **真正的算法参数写在这里**，不是 nextflow CLI 参数（见 Stage 11(b)），如 Kraken2 的 `--quick`、MetaPhlAn 的 `--mpa3` |
| `db_params` 分号语法 | — | Bracken `<kraken2参数>;<bracken参数>`；KMCP `<search>;<profile>`；sylph `<profile>;<taxprof>` |
| `db_type`（databases.csv 第 4 列） | `short;long` | 限定这个库只给短读长/长读长/两者用 |

**14 个 classifier 短读长/长读长兼容性**（决定 `--run_<tool>` 该不该开）：

| 工具 | 短读长 | 长读长 | 备注 |
|:---|:---:|:---:|:---|
| Kraken2 | ✅ | ✅ | 通用，k-mer 法，速度快 |
| Bracken | ✅ | ❌ | 依赖 Kraken2 report 二次分配丰度；长读长自动跳过 |
| MetaPhlAn | ✅ | ❌ | marker gene 法，种级精度高 |
| Centrifuge | ✅ | ❌ | 不接受 FASTA |
| Kaiju | ✅ | ✅ | 蛋白比对法，对新颖/发散序列更敏感 |
| DIAMOND | ✅ | ❌ | 不接受双端，只用 R1 |
| KrakenUniq | ✅ | ✅ | k-mer 唯一性校验，低丰度判读更保守可靠 |
| MALT | ✅ | ❌ | 不接受双端 |
| mOTUs (v3) | ✅ | ✅* | 长读长需先 `motus prep_long`（>2kb） |
| ganon | ✅ | ❌(未评估) | |
| KMCP | ✅ | ❌(灵敏度低) | |
| sylph | ✅ | ✅ | 不接受 FASTA |
| Melon | ❌ | ✅ | 专为长读长设计 |
| MetaCache | ✅ | ✅ | |

### 3.6 taxpasta 标准化

| 参数 | 默认值 | 说明 |
|:---|:---:|:---|
| `--run_profile_standardisation` | 关 | 建议真实项目打开，否则下游拿不到跨样品合并表 |
| `--taxpasta_add_name/rank/lineage` | 全关 | 建议打开并配 `--taxpasta_taxonomy_dir`，省去下游脚本自查 taxonomy ID |
| `--run_krona` | 关 | 建议开，交付时给客户一个可交互点开的物种组成图 |

### 3.7 资源限制

| 机制 | 是否对 taxprofiler 有效 | 说明 |
|:---|:---:|:---|
| `process.resourceLimits`（`local_resources.config`） | ✅ 有效 | 唯一有效的资源控制方式 |
| `--max_cpus` / `--max_memory` | ❌ 无效（静默忽略） | taxprofiler 2.0.1 是新版 nf-core 模板，已弃用旧模板的 `check_max()` 机制 |

---

## 4. `8_taxprofiler_setup` test run 实际用的参数 vs 真实项目该怎么设

### 4.1 test profile（`conf/test.config`）实际参数一览

test profile 的设计目的是**验证代码路径全覆盖**（几乎每个工具都开一遍），不是"推荐的真实项目默认值"：

```groovy
perform_shortread_qc = true
perform_longread_qc = true
perform_shortread_redundancyestimation = true
shortread_qc_mergepairs = true
perform_shortread_complexityfilter = true
perform_shortread_hostremoval = true
perform_longread_hostremoval = true
hostremoval_reference = <极小的测试参考基因组>
perform_runmerging = true
run_kaiju = true
run_kraken2 = true
run_bracken = true
run_malt = false
run_metaphlan = true
run_centrifuge = true
run_diamond = true
run_krakenuniq = true
run_motus = false
run_ganon = true
run_krona = true
run_kmcp = true
run_sylph = true
run_melon = true
run_metacache = true
run_profile_standardisation = true
```

**这解释了为什么 `8_taxprofiler_setup` 的 test run 要拉 ~40 个 Singularity 镜像、且触发过一次 dataflow 停滞**（见 `taxprofiler_安装与测试记录_0701.md` 第 7.1 节）——14 个分类器全开，DAG 分支数远超真实项目会用到的规模。

### 4.2 真实项目推荐参数（对应 `/taxnom` skill 默认方案）

```groovy
perform_shortread_qc = true
perform_shortread_complexityfilter = true
perform_shortread_hostremoval = true          // 仅宿主源样品；环境样品设 false 并去掉 hostremoval_reference
hostremoval_reference = '<宿主基因组fasta>'
run_kraken2 = true
run_bracken = true
run_metaphlan = true
run_krona = true
run_profile_standardisation = true
taxpasta_add_name = true
taxpasta_add_rank = true
```

**理由**：
- Kraken2+Bracken（k-mer 法，速度快、覆盖广，本地已有 Standard-8GB 库）与 MetaPhlAn（marker gene 法，种级精度更高、假阳性率更低）方法学互补，两者头部物种排名一致时结果更可信；只跑一个工具没有交叉验证。
- 其余 12 个分类器默认不开——多数是本地暂无现成数据库、或者是为特定场景优化（Kaiju/DIAMOND 长于发散/新颖序列，KrakenUniq 长于低丰度判读，Melon/MetaCache 专为长读长），客户明确有对应需求（如"怀疑有新颖病毒"、"需要长读长数据支持"）才按需加开，不作为标准配置——多开一个工具就多一份数据库准备成本和运行时间，且默认组合已经能覆盖大多数常规群落物种组成问题。

---

## 5. 快速参考：命令行模板

### Test profile（仅用于环境验证，已在 `8_taxprofiler_setup` 跑通）

```bash
nextflow run nf-core/taxprofiler -r 2.0.1 -profile test,singularity \
  -c configs/local_resources.config \
  --databases test_run/testdata/database_v2.1_taxprofiler2.0.1.csv \
  --outdir test_run/outdir -work-dir test_run/work -resume
```

### 真实项目（短读长、宿主源样品、Kraken2+Bracken+MetaPhlAn 双验证）

```bash
nextflow run nf-core/taxprofiler -r 2.0.1 -profile singularity \
  -c local_resources.config \
  --input samplesheet.csv \
  --databases databases.csv \
  --perform_shortread_qc \
  --perform_shortread_complexityfilter \
  --perform_shortread_hostremoval --hostremoval_reference <host_genome.fasta> \
  --run_kraken2 --run_bracken --run_metaphlan \
  --run_krona --run_profile_standardisation \
  --taxpasta_add_name --taxpasta_add_rank \
  --outdir output_results -work-dir work -resume
```

`databases.csv` 对应最小可用模板（本地库路径）：
```csv
tool,db_name,db_params,db_type,db_path
kraken2,k2standard8gb,,short;long,/Work_bio/references/Metagenomics/kraken2/k2_standard_08gb_20260226/k2_standard_08_GB_20260226.tar.gz
bracken,k2standard8gb,;-r 150,short,/Work_bio/references/Metagenomics/kraken2/k2_standard_08gb_20260226/k2_standard_08_GB_20260226.tar.gz
metaphlan,mpa_vJan25,,short,/Work_bio/references/Metagenomics/metaphlan/
```
（`-r 150` 按实际测序读长改；本地无 kmer_distrib 对应读长时 Bracken 会用最接近的一档，建库时已覆盖 50/75/100/150/200/250/300bp 七档。）
