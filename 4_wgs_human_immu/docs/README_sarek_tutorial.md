# nf-core/sarek 安装、调试与使用教程

> 本文档记录了在本机环境中安装、验证 nf-core/sarek（nf-core 官方人类 WGS/WES 变异检测标准流程）的全过程，
> 并给出后续用真实人类 WGS 样本运行该流程的操作指南。所有关键操作的命令与时间点见同目录下
> `../logs/key_operations_log.md`。

## 1. 流程选择：为什么是 nf-core/sarek

nf-core 官方生态里，**人类 WGS/WES germline + somatic 变异检测**的标准流程就是
[nf-core/sarek](https://nf-co.re/sarek)（区别于之前 1_opossum_YuFan 项目用的 nf-core/rnaseq，那是转录组定量流程）。
sarek 覆盖：

- QC：FastQC / fastp
- 比对：BWA-MEM2 / DragMap / Sentieon BWA
- BAM 处理：去重(MarkDuplicates/biobambam2)、BQSR
- 变异检测：
  - Germline SNV/Indel: GATK HaplotypeCaller、DeepVariant、Sentieon DNAscope
  - Somatic SNV/Indel: GATK Mutect2、Strelka2
  - CNV: CNVkit、ControlFREEC
  - SV: Manta、TIDDIT
- 注释：snpEff / VEP
- 报告：MultiQC

本机已安装：nextflow 25.10.4、docker、singularity/apptainer，均可正常调用（已用 `docker ps` 验证免 sudo 权限）。

## 2. 安装步骤（已完成）

```bash
# nf-core CLI 工具（用于浏览/下载流程、生成 samplesheet 等）
/home/gao/.conda/envs/regular_bioinfo/bin/pip install nf-core
# 验证
nf-core --version        # 4.0.2
nf-core pipelines list | grep sarek   # 3.8.1，4 个月前发布
```

nextflow 本身已存在于 `regular_bioinfo` conda 环境（`/home/gao/.conda/envs/regular_bioinfo/bin/nextflow`），无需重装。
sarek 流程代码本身不需要手动 clone —— `nextflow run nf-core/sarek -r 3.8.1` 会自动从 GitHub 拉取并缓存到
`~/.nextflow/assets/nf-core/sarek`。

## 3. 容器与缓存约定（沿用既有项目习惯）

参考 `1_opossum_YuFan/script_opposum/2_nextflow.sh` 的既有做法：

- 容器引擎用 **singularity**（`-profile test,singularity`），而不是 docker，原因是 HPC/多用户场景下
  singularity 不需要 root daemon，且已有缓存目录约定。
- 全局环境变量已经配置好（写在 shell 启动文件里，整机生效）：
  - `SINGULARITY_CACHEDIR=/Work_bio/gao/configs/.singularity`
  - `APPTAINER_CACHEDIR=/Work_bio/singularity_cache`
  - `NXF_SINGULARITY_PULL_TIMEOUT=30m`
- 额外显式导出 `NXF_SINGULARITY_CACHEDIR=/Work_bio/gao/configs/.singularity`（Nextflow singularity profile
  专用的缓存变量，与上面的 APPTAINER_CACHEDIR 不是同一个，需要分别设置），避免重复下载容器、占用 `/home` 空间。
- `NXF_OPTS="-Xms512m -Xmx2g"` 限制 Nextflow 自身 JVM 内存，避免把它跑挂。

## 4. 资源限制配置

`configs/local_resources.config`：按本机 64 核 / 125G 内存，分配 48 核 + 96G 给 nextflow local executor，
按 process label（process_low/medium/high/high_memory）分级限制单进程资源，预留给系统和其他并发任务。

## 5. 验证运行：官方 test profile

```bash
bash scripts/1_run_sarek_test_profile.sh
```

用的是 sarek 内置的 `-profile test`：流程会自动从 nf-core/test-datasets 下载一份极小的 chr21 测试 fastq +
参考序列，全程不需要真实样本或本地参考基因组，是验证"安装+容器+执行器"是否打通的标准做法。
运行在 tmux 会话 `sarek_test` 中后台执行，日志同步写入 `logs/sarek_test_run.log`。

结果与运行结论见 `logs/key_operations_log.md`。

## 6. 后续真实人类 WGS 样本怎么跑

### 6.1 本机已有的人类参考资源（/Work_bio/references/Homo_sapiens/GRCh38/）

| 用途 | 路径 |
| --- | --- |
| 参考基因组 fasta | `human_gencode_v45/GRCh38.primary_assembly.genome.fa`（已有 `.fai` 及 bwa-mem2 索引） |
| GTF 注释 | `human_gencode_v45/gencode.v45.annotation.gtf` |
| dbSNP | `gatk_bundle/dbsnp_146.hg38.vcf.gz` (+`.tbi`) |
| Mills/1000G indel | `gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` (+`.tbi`) |

注意：sarek 默认推荐用 `--genome GATK.GRCh38`（自动用 AWS iGenomes 资源，含更完整的 known-sites 集合），
若要复用本机已有资源而不重新下载，需要手动指定下列参数组合（而不是 `--genome`）：

```bash
--fasta /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/GRCh38.primary_assembly.genome.fa \
--fasta_fai /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/GRCh38.primary_assembly.genome.fa.fai \
--bwa /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45 \
--dbsnp /Work_bio/references/Homo_sapiens/GRCh38/gatk_bundle/dbsnp_146.hg38.vcf.gz \
--known_indels /Work_bio/references/Homo_sapiens/GRCh38/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```

注意 dbsnp_146 版本较旧（sarek/iGenomes 默认是更新的 dbSNP 版本），如需更严谨的 BQSR/过滤建议升级，
但用于先跑通流程没有问题。

### 6.2 Samplesheet 格式

sarek 要求 CSV，最简单的 germline 单样本双端 fastq 示例：

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
SUBJ01,XX,0,SUBJ01_N,L001,/path/to/SUBJ01_N_R1.fastq.gz,/path/to/SUBJ01_N_R2.fastq.gz
```

- `status`: 0=normal, 1=tumor（做 somatic 配对分析时同一 `patient` 下放一行 normal 一行 tumor）
- 可用 `nf-core pipelines download/launch` 或手写脚本（参考 1_opossum_YuFan 里
  `1_产生_nf-core_Samplesheet.py` 的思路）按目录批量生成

### 6.3 起始步骤与变异检测工具选择

- `--step mapping`（默认，从 fastq 开始）；如果已有 BAM 可以 `--step variant_calling` 跳过比对
- `--tools haplotypecaller,deepvariant,strelka,mutect2,manta,cnvkit,snpeff,vep`（按需组合，逗号分隔）
- 全外显子用 `--wes` + `--intervals` 指定 panel/exome bed；全基因组不用 `--wes`

### 6.4 免疫相关（human_immu 项目命名提示）的补充说明

sarek 本身是通用 germline/somatic 变异检测流程，**不包含 HLA 分型**。如果这个项目的"immu"目标包括：

- **HLA 基因型分型**：sarek 跑出来的 BAM 在 HLA 区域（chr6 MHC 区）比对质量通常不够，需要额外用专门工具
  （如 nf-core/hlatyping、HLA-LA、OptiType 等）单独跑，而不是指望 sarek 给出 HLA 分型结果。
- **免疫基因组体细胞突变/免疫相关基因 panel 注释**：sarek 标准 germline/somatic calling + VEP/snpEff
  注释即可覆盖，跑完后再按免疫相关基因列表（如 IMGT/HLA、KIR、TCR/BCR 基因）做下游筛选过滤。
- **TCR/BCR 免疫组库（repertoire）分析**：这是完全不同的分析类型（需要 RNA-seq 或专门的 repertoire-seq
  数据 + MiXCR/Immcantation 等工具），sarek（DNA-seq 变异检测）解决不了这个问题。

**这部分需要你确认具体分析目标**，我才能给出针对性的真实样本运行命令（见下方待确认问题）。
