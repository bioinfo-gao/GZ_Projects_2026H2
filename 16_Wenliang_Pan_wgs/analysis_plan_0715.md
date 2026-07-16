# 项目 16 — Wenliang Pan · Human Germline WGS · 分析计划

- **Plan Date（创建，immutable）**: 2026-07-15
- **Client**: Wenliang Pan
- **Analyst**: Zhen Gao, PhD — Athenomics
- **Species**: Homo sapiens (GRCh38 / GATK.GRCh38)
- **Tissue/Cell**: not specified by client
- **模式**: `/wgs` Mode A —— standard germline WGS（nf-core/sarek）
- **原始数据**: `/home/gao/Dropbox/Quote_06202601_Wenliang_Pan/`

### 更新记录 / change-log

- 2026-07-15 — 初稿：项目创建、sample table、模式判定、sarek 参数、annotation/rare-variant/HLA/report 步骤。
- 2026-07-15 — 按 user directive 将本内部计划改为中文撰写（technical terms 保留 English 原文）；内容不变。
- 2026-07-15 — 按 user directive 在 header 补入 **Species** + **Tissue/Cell**（报告/计划头部必填项）。

---

## 1. Sample information 与 input data volume

Canonical sample table：[`sample_info.tsv`](sample_info.tsv)（single source of truth）。
文件大小由 `ls -l` 于 **2026-07-15** 从磁盘实测，source path 见上。

| sample   |    client    | species |    seq    |    machine    | R1 (GiB) | R2 (GiB) | per-sample (GiB) |
| :------- | :----------: | :-----: | :-------: | :------------: | :------: | :------: | :--------------: |
| Sample_A | Wenliang Pan |  Human  | WGS PE150 | NovaSeq X Plus |  17.58  |  17.67  |      35.26      |
| Sample_B | Wenliang Pan |  Human  | WGS PE150 | NovaSeq X Plus |  21.56  |  21.59  |      43.15      |

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

| # | step                                            |             tool             | script                                                           |
| :- | :---------------------------------------------- | :---------------------------: | :--------------------------------------------------------------- |
| 0 | 下载 gnomAD-AF + ClinVar DB                     |        wget + bcftools        | `scripts/0_prep_annotation_dbs.sh`                             |
| 1 | 生成 sarek samplesheet                          |            python            | `scripts/1_make_samplesheet.py` → `scripts/samplesheet.csv` |
| 2 | QC + align + dedup + SNV/indel + SV + CNV + VEP | **nf-core/sarek 3.8.1** | `scripts/2_run_sarek.sh`                                       |
| 3 | 运行状态检查                                    |             tail             | `scripts/3_monitor.sh`                                         |
| 4 | 给 calls 叠加 gnomAD AF + ClinVar               |       bcftools annotate       | `scripts/4_annotate_gnomad_clinvar.sh`                         |
| 5 | rare + functional variant 优先级筛选            |            python            | `scripts/5_rare_functional_filter.py`                          |
| 6 | HLA typing（if feasible）                       |              T1K              | `scripts/6_hla_typing.sh`                                      |
| 7 | client report                                   |        R/py + markdown        | （交付时撰写）                                                   |

**sarek 一次运行即覆盖 objectives 的第 1–5 项：**

- **QC**：FastQC + fastp adapter/quality trimming（`--trim_fastq`）+ MultiQC 汇总。
- **Alignment**：bwa-mem2 → GATK MarkDuplicates → CRAM，比对到 **GATK.GRCh38**（analysis-set，即
  ClinVar/gnomAD/VEP 全部所依据的 reference）。
- **SNV/indel**：GATK **HaplotypeCaller**（germline 标准，per-sample GVCF → genotyped VCF）。
- **SV**：**Manta** + **TIDDIT**（两个 orthogonal caller；若 Manta 卡住则 TIDDIT 兜底）。
- **CNV**：**CNVkit**（germline，flat reference —— 无 matched normal）。
- **Functional annotation**：**VEP**（consequence、gene、SIFT、PolyPhen、已知 dbSNP ID）；cache 由
  `--download_cache` 自动下载。

===== 07-15 ZG

**HLA 分型** （人类白细胞抗原分型）是 =检测人体细胞表面标志物的基因检测技术，用于评估组织相容性或疾病风险 。

在器官/骨髓移植中，配型越一致，排斥反应风险越低；它还常用于辅助诊断自身免疫疾病（如强直性脊柱炎检测 HLA-B27）及肿瘤免疫治疗评估。 [[1](https://www.uptodate.com/contents/zh-Hans/human-leukocyte-antigens-hla-a-roadmap), [2](https://baike.baidu.com/item/HLA%E5%88%86%E5%9E%8B/5556320), [3](https://www.mskcc.org/zh-hans/pdf/cancer-care/patient-education/human-leukocyte-antigen-typing?mode=large), [4](https://www.wehelpinc.com/news/detail/id/16.html), [5](https://www.nhc.gov.cn/wjw/s9492/202109/732a604482c64166a93c0e6f966aeab6/files/1739781623585_41303.pdf), [6](https://www.apollohospitals.com/zh-CN/diagnostics-investigations/hla-b27-test), [7](https://pmc.ncbi.nlm.nih.gov/articles/PMC6000529/)]

核心要点

* **检测位点** **：主要分为I类（HLA-A, B, C）和II类（HLA-DR, DQ, DP）基因。**
* **样本类型** **：通常通过抽血或口腔黏膜拭子采集样本。**
* **分型技术** **：现代分型以二代测序（NGS）和SBT（基于测序分型）等DNA分子生物学方法为主，准确度极高。** [[1](https://www.immuquad.com/hla-typing/), [2](https://www.catb.org.cn/%E3%80%90%E5%8E%9F%E5%88%9B%E3%80%91%E5%85%B3%E4%BA%8Ehla%E5%88%86%E5%9E%8B%E7%9A%84%E4%B8%B4%E5%BA%8A%E5%BA%94%E7%94%A8%EF%BC%8C%E8%AE%A9%E6%88%91%E4%BB%AC%E7%9C%8B%E7%9C%8B%E6%9C%80%E6%96%B0/), [3](https://www.abace-biology.com/hla-typing.htm), [4](https://www.mskcc.org/zh-hans/pdf/cancer-care/patient-education/human-leukocyte-antigen-typing?mode=large), [5](https://www.nhc.gov.cn/wjw/s9492/202109/732a604482c64166a93c0e6f966aeab6/files/1739781623585_41303.pdf)]

如果你想了解更多，请告诉我你做HLA分型的具体原因，例如：

* **是为了** **寻找造血干细胞/器官捐献者** **？**
* **是出于** **某种疾病（如风湿免疫类疾病）的辅助排查** **？**
* **或是为了进行** **肿瘤免疫治疗评估** **？**

==================

**这个实验要做HLA分形，是不是基本可以判断样品来自富集白细胞？ 或者某一个白细胞株系？**

不能这样推断。HLA 分型能做，跟样品是不是白细胞、是不是某个白细胞株系，没有必然关系。

**核心原因：HLA 是 germline 基因，在体内每个有核细胞里都一模一样。**

* HLA/MHC 位点在 chr6p21.3，属于基因组的固有部分。你做 DNA 层面的 HLA typing（WES/WGS），本质是在读这个人 germline 的 HLA 等位基因序列——血液、颊拭子、唾液、肿瘤、成纤维细胞、任何 cell line，测出来的 genotype 都是同一套（如 `HLA-A*02:01`）。它告诉你的是 **这个个体的基因型，不是组织来源** 。
* "Leukocyte"（白细胞抗原）只是 **历史命名** ：当年最早在白细胞表面发现这些抗原（因为白细胞好从血里分离、且 MHC 表达高）。但  **MHC class I 几乎在所有有核细胞表面都表达** ，编码基因更是每个细胞都带。所以"能做 HLA 分型"跟"富集了白细胞"是两码事。

**一个关键区分——DNA vs RNA：**

* 你这个项目是  **WGS（DNA）** ，HLA typing 对任何组织都可行，推不出组织类型。
* 如果是  **RNA-seq** ，情况才有一点信息量：HLA class II（HLA-DR/DQ/DP）主要在 APC（B cell、DC、macrophage）里高表达，若某样本 class II 表达特别高，可以*提示*免疫细胞成分多。但那是看 **表达量** ，不是分型本身，且不适用于你的 DNA 数据。

**那么，从这批 WGS 里能不能反推"是不是（克隆性）白细胞/淋巴细胞"？可以，但靠的不是 HLA，而是另外两类 DNA signature：**

| 想判断                                | 真正的 DNA 证据                                                                                                                                                                     | 在本 pipeline 里                                                                           |
| ------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| 是否 lymphoid（T/B 细胞）来源         | **clonal V(D)J recombination** ：TCR（TRA/TRB）、IG（IGH/IGK/IGL）位点的克隆性重排/缺失。单克隆淋巴细胞株系会出现**特定**的重排断点；多克隆全血则重排多样、平均掉看不出 | 需额外分析（如`MiXCR`/`TRUST4` 直接从 reads 抓，或看 TCR/IGH 位点的 SV/coverage 断点） |
| 是否 cell line（而非 primary tissue） | **aneuploidy / 大片段 CNV / LOH** ：传代细胞系常有特征性非整倍体和拷贝数异常                                                                                                  | 已在跑的**CNVkit + TIDDIT** 结果里能看出来                                           |

所以结论： **"要做 HLA 分型"这个事实，既不能判断样品富集了白细胞，也不能判断是某个白细胞株系** ，因而也没法用它去填 header 里那个 `Tissue/Cell: not specified`。

如果你想让数据自己"说"组织/来源，我可以在 sarek 跑完后加一步轻量检查：(1) 扫 TCR/IGH 位点看有没有 clonal V(D)J 重排（判 lymphoid / 单克隆株系），(2) 看 CNVkit 的 aneuploidy/LOH 画像（判 primary vs cell line）。要的话我就把它写成 `scripts/8_origin_check.sh`。

===================

T1K 是一款在生物信息学界非常强大且相对较新的工具，专门用于从高通量测序数据中进行 

**HLA（人类白细胞抗原） 和 KIR（杀伤细胞免疫球蛋白样受体） 的基因分型。**

它由开发了著名免疫组库分析工具 TRUST4 的团队 **（Bo Li 实验室等）**

于 2023 年发表（Nature Communications）。



作为生信工程师，你可以从以下几个核心维度来理解 T1K：

1. 为什么需要 T1K？（HLA 分型的痛点）
   HLA 和 KIR 区域是人类基因组中多态性最高（高度突变）且同源性极高的区域。
   传统比对的失败： 如果你用 BWA 或 STAR 把 Reads 比对到标准参考基因组（如 hg38），由于 HLA 有成千上万种等位基因（Alleles），Reads 会发生严重的多重比对（Multi-mapping），导致无法确定样本到底属于哪种 HLA 亚型。
   T1K 的作用： 它是专门为了解决这个“超级多态性”问题而设计的，能够精准给出样本的 HLA 分型（如 HLA-A*02:01）。
2. T1K 是如何工作的？（核心算法）
   T1K 放弃了传统的“线性比对”思路，采用了一套非常适合高度多态性序列的算法：
   K-mer 提取： 它首先利用 IMGT/HLA 和 IPD-KIR 数据库，提取所有已知等位基因的特异性 K-mers。
   序列图 (Sequence Forest)： 它将那些高度相似的等位基因序列构建成一个图/森林结构。
   EM 算法 (期望最大化)： 当把样本的 Reads 映射到这个图中后，T1K 使用 EM 算法来计算每一种可能等位基因的表达丰度 (Abundance)。丰度最高、最合理的那个，就被判定为真实的基因型。
3. T1K 的主要优势（你的应用场景）
   结合你周报中提到的工作，T1K 对你可能有极大的实用价值：
   极其通用的数据兼容性：
   它是目前少数能完美处理 Bulk RNA-seq、单细胞 RNA-seq (scRNA-seq) 以及 WES/WGS 数据的全能工具。（注：它最初是主打 RNA-seq 分型的，因为 RNA-seq 数据中 HLA 表达量极高，分型非常准。）
   计算资源消耗极低：
   相比于老牌的 OptiType 或 POLYSOLVER 需要复杂的整数线性规划或漫长的比对，T1K 速度极快，内存占用很小。
   同时支持 KIR 分型：
   很多工具只能做 HLA，而 T1K 还能做 KIR 分型，这对于肿瘤免疫学（如 NK 细胞研究）非常重要。
4. 如何在你的 Pipeline 中集成 T1K？
   如果你想在你的 RNA-seq 流程（比如你给 Yue Liu 或 Qiuchen Li 做的项目，如果是人类样本的话）中加入 HLA 分型作为一个“增值服务”，T1K 是极其理想的插件：

===== Above 07-15 ZG

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
- 

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
