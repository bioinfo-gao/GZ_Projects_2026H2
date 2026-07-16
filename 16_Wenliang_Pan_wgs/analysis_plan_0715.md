# 项目 16 — Wenliang Pan · Human Germline WGS · 分析计划

- **Plan Date（创建，immutable）**: 2026-07-15
- **Client**: Wenliang Pan
- **Analyst**: Zhen Gao, PhD — Athenomics
- **模式**: `/wgs` Mode A —— standard germline WGS（nf-core/sarek）
- **原始数据**: `/home/gao/Dropbox/Quote_06202601_Wenliang_Pan/`

### 更新记录 / change-log
- 2026-07-15 — 初稿：项目创建、sample table、模式判定、sarek 参数、annotation/rare-variant/HLA/report 步骤。
- 2026-07-15 — 按 user directive 将本内部计划改为中文撰写（technical terms 保留 English 原文）；内容不变。

---

## 1. Sample information 与 input data volume

Canonical sample table：[`sample_info.tsv`](sample_info.tsv)（single source of truth）。
文件大小由 `ls -l` 于 **2026-07-15** 从磁盘实测，source path 见上。

| sample | client | species | seq | machine | R1 (GiB) | R2 (GiB) | per-sample (GiB) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Sample_A | Wenliang Pan | Human | WGS PE150 | NovaSeq X Plus | 17.58 | 17.67 | 35.26 |
| Sample_B | Wenliang Pan | Human | WGS PE150 | NovaSeq X Plus | 21.56 | 21.59 | 43.15 |

- **数据集总量：2 samples，gzip FASTQ 合计 78.40 GiB（`du -shc` = 79 G）。**
- 单一 flowcell/lane `CKDL260012882-1A_23JCJ2LT3_L4`。估计 depth ~30–40×（待 alignment 后由 mosdepth 确认）。
- **Provenance 核查（已做）**：两个样本文件的 MD5 与交付 manifest `MD5.txt` 中的
  `01.RawData/Sample_A|Sample_B/` 一致；manifest 里大量 `HFD_*/ISRIB_*/Tg_*` 条目是**同一 lane 上多路复用的
  其他客户样本**，不属于本项目。无 mislabeling，本项目只有 Sample_A / Sample_B 两个样本。
- **Species sanity check 已排入流程**：该 lane 以看似 mouse metabolic 的样本为主，客户虽标注 Human，
  仍需在信任 human-specific annotation（ClinVar/gnomAD/HLA）之前，先确认 **GRCh38 primary mapping rate
  （预期 >95%）**（来自 sarek / samtools flagstat）。

## 2. Objectives（来自客户邮件）

Standard germline WGS：(1) raw-data QC，(2) alignment 到 human reference，(3) SNV + indel calling，
(4) SV + CNV analysis，(5) 对 ClinVar/gnomAD/其他 DB 的 variant annotation，(6) rare + potentially
functional variant 识别，(7) HLA typing（if feasible），(8) 汇总 biologically significant variants 的 report。

## 3. Pipeline 与 step map

| # | step | tool | script |
| :--- | :--- | :---: | :--- |
| 0 | 下载 gnomAD-AF + ClinVar DB | wget + bcftools | `scripts/0_prep_annotation_dbs.sh` |
| 1 | 生成 sarek samplesheet | python | `scripts/1_make_samplesheet.py` → `scripts/samplesheet.csv` |
| 2 | QC + align + dedup + SNV/indel + SV + CNV + VEP | **nf-core/sarek 3.8.1** | `scripts/2_run_sarek.sh` |
| 3 | 运行状态检查 | tail | `scripts/3_monitor.sh` |
| 4 | 给 calls 叠加 gnomAD AF + ClinVar | bcftools annotate | `scripts/4_annotate_gnomad_clinvar.sh` |
| 5 | rare + functional variant 优先级筛选 | python | `scripts/5_rare_functional_filter.py` |
| 6 | HLA typing（if feasible） | T1K | `scripts/6_hla_typing.sh` |
| 7 | client report | R/py + markdown | （交付时撰写） |

**sarek 一次运行即覆盖 objectives 的第 1–5 项：**
- **QC**：FastQC + fastp adapter/quality trimming（`--trim_fastq`）+ MultiQC 汇总。
- **Alignment**：bwa-mem2 → GATK MarkDuplicates → CRAM，比对到 **GATK.GRCh38**（analysis-set，即
  ClinVar/gnomAD/VEP 全部所依据的 reference）。
- **SNV/indel**：GATK **HaplotypeCaller**（germline 标准，per-sample GVCF → genotyped VCF）。
- **SV**：**Manta** + **TIDDIT**（两个 orthogonal caller；若 Manta 卡住则 TIDDIT 兜底）。
- **CNV**：**CNVkit**（germline，flat reference —— 无 matched normal）。
- **Functional annotation**：**VEP**（consequence、gene、SIFT、PolyPhen、已知 dbSNP ID）；cache 由
  `--download_cache` 自动下载。

## 4. Analysis rationale 与 decision criteria

- **`--genome GATK.GRCh38`**（而非本地 GENCODE fasta）：analysis-set 的 contig naming + 自带 known-sites
  与 intervals，正是每个 clinical DB（ClinVar、gnomAD、VEP）所构建的坐标系；避免 contig-name mismatch。
  iGenomes bundle 首次运行自动获取。
- **`--skip_tools baserecalibrator`**：NovaSeq X Plus 的 base quality 已良好校准，BQSR 收益极小且耗时约
  +20%（proj4 已确认）；本项目也无 project-specific known-sites 需求。
- **`--aligner bwa-mem2`**：比 bwa 快；RAM 已覆盖（60 GB/task）。
- **Annotation 拆分（sarek 内 VEP，gnomAD+ClinVar 放到 post-step 4）**：sarek 内 VEP 给 consequence/impact；
  gnomAD **AF-only**（GATK resource，~3 GB）+ ClinVar（~0.2 GB）由 `bcftools annotate` 叠加。用 AF-only
  gnomAD（相对完整 sites VCF >1 TB）在磁盘几乎无成本，同时正好提供 rarity filtering 需要的 allele-frequency 字段。
- **Rare + functional 定义（step 5）：**
  - **RARE** = gnomAD_AF < 0.001 **或** absent from gnomAD。
  - **FUNCTIONAL** = VEP HIGH/MODERATE consequence（LoF：stop-gain / frameshift / splice-donor·acceptor /
    start-loss；missense；inframe indel）**或** ClinVar Pathogenic/Likely_pathogenic。
  - **Flagged** = RARE **且** FUNCTIONAL。ClinVar P/LP 无论 frequency 一律 surface（clinically actionable）。
  - 仅评估 FILTER=PASS 的 variant。
- **HLA "if feasible"**：从 CRAM 提取 MHC-region（chr6:28–34 Mb）+ unmapped reads，交给 T1K 做 class I + II
  分型 —— 标准 WGS depth 可行。标注为 feasible，最终以 depth 确认为准。

## 5. Resource allocation

- `scripts/local_resources.config`：`queueSize=2`，`cpus=16`/task，`BWAMEM2 memory=60.GB`，
  防御性 `CNNSCOREVARIANTS cpus=9`。CLI `--max_memory 120.GB --max_cpus 56`。
- 2 samples × 16 = **32 threads sustained**（≤ 56 cap ✓）。最坏 RAM 2×60 = 120 GB / 125 GB（+65 GB swap）——
  即 proj4 已验证的 envelope。
- Disk：`/home` free 696 G，承载 `work/` + `output_results/` + VEP cache（2 samples 的 WGS work dir 约
  200–400 GB —— 充足）。共享 annotation DB 放在 `/Work_bio/.../annotation/`。
- 当前（2026-07-15）**独占**运行（无 co-running pipeline）。

## 6. 运行顺序

```bash
conda run -n regular_bioinfo python scripts/1_make_samplesheet.py   # samplesheet.csv
bash scripts/0_prep_annotation_dbs.sh    # （tmux，并行）sarek 运行期间下载 DB
bash scripts/2_run_sarek.sh              # （自动进 tmux 'pan_wgs'，失败自动 -resume）
# sarek 完成后：
bash scripts/4_annotate_gnomad_clinvar.sh
conda run -n regular_bioinfo python scripts/5_rare_functional_filter.py
bash scripts/6_hla_typing.sh
```

## 7. 时间线估计

2 samples、queueSize=2，约 ~10–20 h（bwa-mem2 align + MarkDuplicates ~3–4 h/sample 为主；HaplotypeCaller
scatter-gather；Manta 为不确定项 —— TIDDIT 兜底）。随后 annotation + filtering <1 h；HLA ~0.5 h。

## 8. 已知局限 / 待客户确认项

- **SNV/indel caller** = 仅 HaplotypeCaller（GATK 标准）。如客户需 orthogonal 验证可加 DeepVariant
  （CPU-only → 慢很多）。
- **无 matched normal 的 CNV**：germline CNVkit 用 flat reference；大 / mosaic CNV 可靠，小 focal event
  较弱；无 panel-of-normals。
- **gnomAD** = AF-only resource（v2-derived，genome+exome joint AF）。对 rarity 足够；不区分 subpopulation。
  如需完整 gnomAD v4 per-population AF 可另行添加（视 disk 而定）。
- **HLA** 标注为 feasible，depth 确认后定案。T1K 给 2-field typing；若需更高 resolution / 验证，可改用
  HLA-LA（graph-based，更重）。
- **Report 侧重**：待与客户确认 —— 通用 variant catalogue，还是聚焦某 gene/phenotype（邮件仅泛泛提到
  "variants of potential biological significance"）。

## 9. 交付前 self-audit checklist（mandatory）

- MultiQC：per-sample mapping rate（>95% human）、duplication、coverage（mosdepth）—— 任一 outlier ⇒ 排查。
- 信任 human annotation 前确认 GRCh38 mapping rate（species sanity check，§1）。
- SNV 全基因组 Ti/Tv ratio ~2.0–2.1（SNV calls 的 sanity）。
- 在 IGV 抽查一个 ClinVar P/LP 命中；确认 gnomAD_AF annotation 确实已填充。
- 两个样本间数值 sanity —— variant count 差异悬殊 ⇒ 写入 report 前先排查。

---
*Reference 实现：`4_wgs_human_immu/`（Mode A sarek）、`/wgs` skill。交付签名：
Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics。*
