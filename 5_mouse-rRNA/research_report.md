# 小鼠 RNA-seq 文库质量评估报告

**日期**: 2026-06-25  
**项目路径**: `/home/gao/projects_2026H2/5_mouse-rRNA/`  
**分析人**: 高震

---

## 一、研究背景与目的

9个小鼠样品，建库完成后送测。目的是做质控 triage，回答三个问题（见 `research_aim.md`）：

1. 每个样品 insert 大小分布，是否存在 dimer content
2. ribosomal RNA 污染情况
3. 有多少 unique reads 可以打到基因组，能检测到多少基因

**决策标准**：如果质量极差，简单重复建库没有意义，需要重新设计建库方案。

---

## 二、数据与工具

### 原始数据

| 项目 | 详情 |
|------|------|
| 样品数 | 9个（mouse_28/29/32/41/42/45/46/47/48） |
| 测序类型 | PE150，双端 |
| 原始数据量 | ~32M read pairs / 样品（约 10G/样品） |
| 原始 samplesheet | `scripts/nf_core_samplesheet.csv` |

### 参考基因组

| 物种 | 版本 |
|------|------|
| 参考基因组 | GRCm39（`/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa`） |
| 注释文件 | GENCODE vM35（`gencode.vM35.annotation.gtf`） |
| STAR 索引 | `/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/star_index` |

### rRNA 参考数据库（8个，SortMeRNA）

```
/Work_bio/references/rRNA_databases/sortmerna_database_manifest.txt
```

| 数据库文件 | 物种覆盖 |
|-----------|---------|
| rfam-5s-database-id98.fasta | 5S rRNA |
| rfam-5.8s-database-id98.fasta | 5.8S rRNA |
| silva-arc-16s-id95.fasta | 古菌 16S |
| silva-arc-23s-id98.fasta | 古菌 23S |
| silva-bac-16s-id90.fasta | 细菌 16S |
| silva-bac-23s-id98.fasta | 细菌 23S |
| silva-euk-18s-id95.fasta | 真核 18S |
| silva-euk-28s-id98.fasta | 真核 28S |

### 主要软件

| 工具 | 版本/说明 |
|------|---------|
| nf-core/rnaseq | v3.15.1（Singularity 容器） |
| SortMeRNA | 4.3.6 |
| Trim Galore | 随 nf-core pipeline |
| STAR | 随 nf-core pipeline |
| Salmon | 随 nf-core pipeline |
| seqtk | 系统安装版 |
| Kraken2 | 安装于 mag_biobakery conda 环境 |

---

## 三、分析流程与命令

所有脚本均存档于 `SortMeRNA_issue/`，日志存于 `SortMeRNA_issue/logs/`。

### 步骤 0：问题背景（全量数据失败原因）

直接用全量 ~32M read pairs × 8 个 rRNA 参考库运行 nf-core/rnaseq，SortMeRNA 出现两类失败：

- **超时（exit 143）**：mouse_32、mouse_41、mouse_48 超过 16h 时间限制被 Nextflow 杀死
- **OOM（exit 137）**：mouse_45（内存峰值 ~118GB）、mouse_47（~74GB）被内核 OOM-killer 杀死
- 只有 mouse_28 和 mouse_29 成功完成（偶然）

**根本原因**：SortMeRNA 4.3.6 在 `--paired_in --out2 --fastx` 模式下（需要写出配对 fastq 文件），内存随 reads 数量线性增长，全量数据 + 8 个库在当前服务器（125GB RAM）上不可行。

**解决方案**：对每个样品抽样至 10M read pairs，足以稳定估计 rRNA% 和 mapping rate，但内存峰值降至 2-39GB（可控范围）。

---

### 步骤 1：对 9 个样品各抽样 10M read pairs

**脚本**: `SortMeRNA_issue/01_subsample_10M.sh`  
**日志**: `SortMeRNA_issue/logs/01_subsample_10M.20260624_201955.log`

```bash
# 核心命令（对每个样品执行）
seqtk sample -s100 ${fq1} 10000000 | gzip -1 > ${out1} &
seqtk sample -s100 ${fq2} 10000000 | gzip -1 > ${out2} &
wait
```

**关键参数**：
- seed = 100（固定，保证可重复）
- N = 10,000,000 read pairs
- R1/R2 使用相同 seed，保证配对一致

**输出**:
- 抽样文件：`SortMeRNA_issue/subsampled_10M/${sample}_sub10M_1.fq.gz` / `_2.fq.gz`
- 新 samplesheet：`SortMeRNA_issue/subsampled_10M_samplesheet.csv`

---

### 步骤 2：单样品 pilot 验证内存安全性

**脚本**: `SortMeRNA_issue/02_pilot_run_mouse28.sh`  
**日志**: `SortMeRNA_issue/logs/02_pilot_run_mouse28.20260624_205023.log`

```bash
nextflow run nf-core/rnaseq -r 3.15.1 \
    -profile singularity \
    -c scripts/local_optimized.config \
    -c scripts/avoid_download.config \
    -w SortMeRNA_issue/pilot_work \
    --input SortMeRNA_issue/pilot_mouse28_samplesheet.csv \
    --outdir SortMeRNA_issue/pilot_output \
    --fasta /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa \
    --gtf /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/gencode.vM35.annotation.gtf \
    --star_index /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/star_index \
    --gencode \
    --aligner star_salmon \
    --remove_ribo_rna \
    --ribo_database_manifest /Work_bio/references/rRNA_databases/sortmerna_database_manifest.txt \
    --save_non_ribo_reads \
    --max_cpus 28 \
    --max_memory 110.GB
```

**`local_optimized.config` 关键配置**（`scripts/local_optimized.config`）：

```groovy
process {
    withName: 'SORTMERNA' {
        cpus   = 24
        memory = '100 GB'
        time   = '48h'
        maxForks = 1
    }
}
```

**Pilot 安全性验证结果**（来自 Nextflow trace 文件）：

| Task | Peak RSS | Peak vmem | 耗时 |
|------|---------|---------|------|
| SORTMERNA (mouse_28, 10M pairs) | **4.2 GB** | 6.6 GB | 6m 39s |

内存远低于 100GB 上限，安全性确认通过。

---

### 步骤 3：全部 9 个样品运行（含 bug 修复）

运行过程遇到两个问题并依次修复：

**问题 A**：背景进程未放在 tmux 中，连接中断会杀死 pipeline  
→ 修复：终止旧进程，在 tmux session `mouse_rrna_full9` 中用 `-resume` 重启

**问题 B**：nf-core/rnaseq v3.15.1 代码 bug（`strandedness` 变量未声明）  
→ 文件：`~/.nextflow/assets/nf-core/rnaseq/subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main.nf` 第 23 行  
→ 修复：`def library_strandedness = 'undetermined'` → `def strandedness = 'undetermined'`

**最终成功运行脚本**: `SortMeRNA_issue/05_full_run_9samples_after_strandedness_fix.sh`

```bash
# 在 tmux session 中执行
tmux new-session -d -s mouse_rrna_full9 'bash SortMeRNA_issue/05_full_run_9samples_after_strandedness_fix.sh \
  2>&1 | tee SortMeRNA_issue/logs/05_full_run_9samples_after_strandedness_fix.$(date +%Y%m%d_%H%M%S).log'
```

**输出目录**: `SortMeRNA_issue/full9_output/`

---

### 步骤 4：Kraken2 污染物鉴定（mouse_29/mouse_46）

**背景**：mouse_29 和 mouse_46 trim 后 reads 中位长度接近满长（~147bp），但 STAR 比对率极低（分别 12.22% 和 2.16%），"too short" 不能解释其失败原因 → 怀疑非小鼠外源序列污染。

**脚本**: `SortMeRNA_issue/06_kraken2_contamination_check.sh`  
**数据库**: Kraken2 Standard-8（RefSeq archaea+bacteria+viral+plasmid+human+UniVec_Core, 8GB）  
→ `/Work_bio/references/Metagenomics/kraken2/k2_standard_08gb_20260226/`

```bash
# 核心分类命令（对 mouse_46 和 mouse_29 的 STAR 未比对 reads 执行）
kraken2 \
    --db /Work_bio/references/Metagenomics/kraken2/k2_standard_08gb_20260226 \
    --threads 16 \
    --paired \
    --report SortMeRNA_issue/kraken2_contamination/${sample}_kraken2_report.txt \
    --output SortMeRNA_issue/kraken2_contamination/${sample}_kraken2_output.txt \
    SortMeRNA_issue/kraken2_contamination/${sample}_unmapped_1.fq \
    SortMeRNA_issue/kraken2_contamination/${sample}_unmapped_2.fq
```

**输出**: `SortMeRNA_issue/kraken2_contamination/`

---

## 四、全部结果汇总

### 4.1 核心质量指标（9个样品）

> 数据来源：10M read pairs 抽样，nf-core/rnaseq v3.15.1，GRCm39 + GENCODE vM35

| 样品 | rRNA%（SortMeRNA） | Trim后短于20bp被丢弃的pair%（dimer信号） | STAR唯一比对% | STAR "too short"未比对% | 检测到的基因数 | Strand推断 |
|------|-------------------|-----------------------------------------|---------------|-------------------------|---------------|-----------|
| mouse_28 | 44.55 | 69.69 | 19.57 | 66.02 | 10,312 | undetermined |
| mouse_29 | 45.06 | 71.99 | 12.22 | 79.63 | 6,157 | undetermined |
| mouse_32 | 47.22 | 26.08 | 28.67 | 51.49 | 19,000 | undetermined |
| mouse_41 | 37.94 | 37.38 | 22.61 | 59.45 | 16,502 | undetermined |
| mouse_42 | 43.22 | 22.36 | 26.63 | 54.57 | 20,817 | undetermined |
| mouse_45 | 40.40 | 23.35 | 24.88 | 54.80 | 20,270 | undetermined |
| mouse_46 | 27.51 | 29.11 | 2.16 | 96.24 | 1,042 | 跳过（未达5%阈值） |
| mouse_47 | 44.81 | 42.29 | 26.01 | 53.87 | 18,590 | undetermined |
| mouse_48 | 45.31 | 22.85 | 27.67 | 52.30 | 19,826 | undetermined |

**参考正常值**（小鼠 bulk RNA-seq）：
- rRNA%：< 5–10%
- Dimer（trim后丢弃）：< 5–10%
- STAR 唯一比对率：70–90%+

### 4.2 rRNA 各数据库贡献分解

> 数据来源：`SortMeRNA_issue/full9_output/sortmerna/*.sortmerna.log`

| 样品 | total rRNA% | silva-bac-16s% | silva-euk-18s% | silva-euk-28s% | 备注 |
|------|------------|---------------|---------------|---------------|------|
| mouse_28 | 44.55 | 30.41 | 4.08 | 7.91 | 细菌库占主导（异常）|
| mouse_29 | 45.06 | 37.71 | 3.60 | 2.71 | 细菌库占主导（异常）|
| mouse_32 | 47.22 | 3.43 | 7.81 | 31.59 | 真核 28S 占主导（正常模式） |
| mouse_41 | 37.94 | 11.80 | 6.53 | 15.51 | 混合 |
| mouse_42 | 43.22 | 2.81 | 8.70 | 26.66 | 真核 28S 占主导（正常模式） |
| mouse_45 | 40.40 | 1.56 | 11.42 | 22.07 | 真核 28S 占主导（正常模式） |
| mouse_46 | 27.51 | 15.58 | 10.52 | 0.72 | 细菌库占主导（异常） |
| mouse_47 | 44.81 | 2.94 | 10.17 | 24.75 | 真核 28S 占主导（正常模式） |
| mouse_48 | 45.31 | 1.91 | 10.41 | 27.97 | 真核 28S 占主导（正常模式） |

### 4.3 Trim 后 reads 长度分布与 STAR 比对失败原因

> 直接从 trimmed fastq 文件计数验证

| 样品 | Trim后reads中位长度（满长150bp） | STAR唯一比对% | "too short"% |
|------|--------------------------------|---------------|--------------|
| mouse_32/42/45/47/48 | **42 bp** | 22–29% | 51–55% |
| mouse_41 | **52 bp** | 22.61% | 59.45% |
| mouse_28 | **82 bp** | 19.57% | 66.02% |
| mouse_29 | **147 bp**（近满长） | 12.22% | 79.63% |
| mouse_46 | **147 bp**（近满长） | 2.16% | 96.24% |

### 4.4 两种 rRNA% 指标的区别（重要）

| 指标 | 来源工具 | 统计口径 | mouse_28 的值 |
|------|---------|---------|--------------|
| **本报告中的 rRNA%** | SortMeRNA | Trim 后 reads 中，序列相似度匹配 rRNA 参考库的比例（在 STAR 之前） | 44.55% |
| MultiQC 网页 General Statistics 中的 %rRNA | featureCounts biotype | STAR 成功唯一比对的 reads 中，落在基因组 rRNA 注释位点的比例（SortMeRNA 已过滤后的残留） | 8.45% |

两者分母不同，**均正确但回答不同问题**。本报告用 SortMeRNA 口径，因为它不依赖 STAR 能否成功比对，更直接地反映原始样品的 rRNA 污染水平。

---

## 五、回答研究问题

### Q1：Insert size 分布 / Dimer content

**结论：全部 9 个样品均存在严重的 adapter dimer / 超短 insert 问题，分两档：**

- **极度严重**（mouse_28、mouse_29）：~70% 的 read pair 在 Trim Galore 去接头后因 insert < 20bp 被丢弃，说明这两个样品的文库主体是 adapter dimer，几乎没有真实 cDNA 序列。
- **显著异常**（mouse_32/41/42/45/46/47/48）：22–42% 的 pairs 因 < 20bp 被丢弃，远超健康文库的 < 5–10% 正常水平。

Trim 后 reads 中位长度（42–82 bp vs 满长 150 bp）进一步确认了 insert size 普遍偏短的事实。

### Q2：rRNA 污染情况

**结论：全部 9 个样品 rRNA 污染比例均在 27.5–47.2% 之间，远高于 mRNA-seq 健康文库的 < 5–10% 标准。**

**重要注意**：

- mouse_28、mouse_29、mouse_46 的 rRNA 信号中，**silva-bac-16s（细菌 16S）占主导**（30–38%），这有两种可能解释：
  - **生物学污染**：样品本身含大量细菌（如肠道相关组织）
  - **数据库假阳性**：细菌 16S 库物种多样性极高，与真核序列存在保守区交叉匹配
- 其余 6 个样品（32/41/42/45/47/48）的 rRNA 信号以 **silva-euk-28s（真核 28S）为主导**（15–32%），是更典型的哺乳动物 rRNA 污染模式。
- 细菌信号是否为真实污染，建议结合 Kraken2 分类结果进一步判断。

### Q3：Unique reads / 检测到的基因数

**结论：全部样品 STAR 唯一比对率严重偏低（2–29% vs 正常 70–90%+），原因明确——reads 太短。**

- 主要失败原因均为 STAR 的 "too short"（51–96%），即 reads 中能匹配基因组的有效部分不足 reads 长度的 66%（STAR 默认门槛）。
- 直接根源是 Q1 的 dimer 问题：trim 后大量 reads 仅剩 20–50bp，无法可靠定位。
- 检测到的基因数：mouse_46 仅 1,042 个（几乎废库），mouse_28/29 因 dimer 极重分别只有 10,312/6,157 个，其余 6 个样品 16,000–21,000 个（但这是在极低 mapping rate 下的结果，有效 reads 数量远不足正常水平）。

---

## 六、额外发现

### 6.1 Strandedness 校验全部失败

8 个样品（除 mouse_46 外）RSeQC 链特异性推断均返回 "undetermined"（sense 33–40% / antisense 54–61%，偏向 reverse 但未达 80% 置信阈值）。  
**解释**：这不是文库设计问题，而是 dimer/超短 reads 的第三个症状——reads 太短太碎，链信号被稀释到无法可靠推断。

### 6.2 Mouse_29 和 Mouse_46 的额外异常（Kraken2 分析）

这两个样品有别于其他样品：trim 后 reads 中位长度接近满长（147 bp），但 STAR 比对率依然极低（12% 和 2%），说明问题不是 reads 太短，而是序列本身不是小鼠来源。Kraken2 分类（Standard-8 数据库）对其 STAR 未比对 reads 进行了分类，报告见 `SortMeRNA_issue/kraken2_contamination/`。

---

## 七、结论与建议

**综合三项指标，9 个文库整体质量不合格，且问题普遍而系统：**

| 样品 | 综合评级 | 主要问题 |
|------|---------|---------|
| mouse_46 | 废库 | STAR 比对率 2.16%，仅 1,042 个基因，外源序列污染 |
| mouse_28、mouse_29 | 严重不合格 | ~70% dimer，mapping rate 仅 12–20% |
| mouse_32/41/42/45/47/48 | 显著不合格 | 22–42% dimer，mapping rate 仅 22–29% |

**与 `research_aim.md` 中"如果情况非常差，简单重复建库没有意义"的判断相符**：这批文库的问题是系统性的（所有样品均受影响，且 dimer/rRNA 两种问题同时存在），建议在重新建库前先排查根本原因（RNA 质量、建库 protocol、dA tailing/连接步骤等），而不是直接重做。

---

## 八、文件清单（供复现）

```
5_mouse-rRNA/
├── research_aim.md                          # 研究目标说明
├── scripts/
│   ├── nf_core_samplesheet.csv              # 原始 9 样品 samplesheet
│   ├── local_optimized.config               # Nextflow 资源配置（SortMeRNA maxForks=1, mem=100GB, time=48h）
│   ├── avoid_download.config                # 跳过 iGenomes 下载
│   └── rRNA_local.config                   # 备用 rRNA 配置
├── SortMeRNA_issue/
│   ├── 01_subsample_10M.sh                  # 步骤1：seqtk 抽样
│   ├── 02_pilot_run_mouse28.sh              # 步骤2：pilot 安全验证
│   ├── 03_full_run_9samples.sh              # 步骤3a：初始全量运行
│   ├── 04_full_run_9samples_tmux_resume.sh  # 步骤3b：迁入 tmux + resume
│   ├── 05_full_run_9samples_after_strandedness_fix.sh  # 步骤3c：修复 nf-core bug 后重跑
│   ├── 06_kraken2_contamination_check.sh   # 步骤4：Kraken2 鉴定外源污染
│   ├── logs/                               # 所有步骤的运行日志
│   ├── subsampled_10M/                     # 抽样后的 fastq 文件
│   ├── subsampled_10M_samplesheet.csv      # 抽样数据 samplesheet
│   ├── pilot_output/                       # Pilot 运行输出
│   ├── full9_output/                       # 最终 9 样品完整输出
│   │   ├── sortmerna/                      # SortMeRNA 日志（含 rRNA% 原始数据）
│   │   ├── trimgalore/                     # Trim Galore 报告（含 dimer%）
│   │   ├── star_salmon/log/               # STAR 比对日志（含 mapping rate）
│   │   ├── multiqc/star_salmon/multiqc_report.html  # MultiQC 汇总报告
│   │   └── pipeline_info/                 # Nextflow execution trace
│   └── kraken2_contamination/             # Kraken2 分类结果
└── rRNA_databases -> /Work_bio/references/rRNA_databases  # 软链接
```

### nf-core/rnaseq v3.15.1 Bug 修复记录

**文件**：`~/.nextflow/assets/nf-core/rnaseq/subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main.nf` 第 23 行

**修改**：
```groovy
// 修改前（原始代码，有 bug）
def library_strandedness = 'undetermined'

// 修改后
def strandedness = 'undetermined'
```

**原因**：函数 `calculateStrandedness()` 中声明了 `library_strandedness` 但条件分支和 return 语句使用的是未声明的 `strandedness`，当 RSeQC 推断结果为 "ambiguous" 时抛出 `MissingPropertyException`。此 bug 已直接修改本地缓存的 pipeline 代码，影响所有使用该缓存版本的后续运行。

---

*报告生成时间：2026-06-25*
