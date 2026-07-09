
这是一份关于鸟枪法宏基因组学（Shotgun Metagenomics）的详细教程。根据您的要求，本教程分为三个部分：宏基因组学总论、依赖组装的方法学（基于 nf-core/mag 5.4.2）以及不依赖组装的方法学（基于 nf-core/taxprofiler 2.0.1）。

### 一、 鸟枪法宏基因组学总论


鸟枪法宏基因组学（Shotgun Metagenomics）是一种对环境、临床或动植物样本中所有微生物（包括细菌、古菌、真菌和病毒）的基因组 DNA 进行无差别、随机打断并进行高通量测序的技术。
与传统的 16S rRNA 等扩增子测序不同，鸟枪法不仅能回答 “样本里有什么物种？(Who is there?)”，还能回答 “它们能做什么？(What can they do?)”，因为它测序的是微生物的全部功能基因。

在数据分析策略上，宏基因组学主要分为两大流派：

依赖组装的方法 (Assembly-based): 也就是常说的“从头组装（De novo assembly）”。这种方法试图将测序得到的短序列（reads）像拼图一样拼接成长片段（contigs），甚至恢复出完整的单菌基因组（MAGs）。它非常适合发现未知新物种和进行深度的功能及代谢通路分析，但对计算资源要求极高，且低丰度的物种很难被完整拼接。

不依赖组装的方法 (Assembly-free / Taxonomic Profiling): 直接将测序得到的 reads 与已知的参考数据库（如 NCBI RefSeq）进行比对，以确定样本中的物种组成及其相对丰度。它速度快、节省计算资源，适合大批量样本的快速分类鉴定，但高度依赖现有数据库，难以鉴定数据库中没有的未知物种（即“微生物暗物质”）。

### 二、 依赖组装的方法学：基于 nf-core/mag (v5.4.2)


nf-core/mag 是一个用于宏基因组组装、分箱 (binning) 和注释的最佳实践生物信息学分析流程。它可以从复杂的群落中提取出高质量的宏基因组组装基因组（MAGs）。



#### 核心分析步骤 (v5.4.2 架构)



##### 数据质控与预处理 (QC & Pre-processing):

该流程同时支持二代短读长 (Short-reads) 和三代长读长 (Long-reads，如 Nanopore)。

自动调用 fastp, AdapterRemoval 或 trimmomatic 进行接头去除和低质量序列修剪（长读长使用 Porechop）。



##### 组装 (Assembly):

利用 MEGAHIT 或 SPAdes (支持 metaSPAdes 模式) 将散乱的 reads 拼接成连续的重叠群 (Contigs)。

利用 Quast 评估组装的连续性（如 N50指标），短读长还可使用 ALE 评估质量。

古DNA特性： 支持使用 PyDamage 验证古DNA的损伤模式。

基因预测 (Gene Prediction):

使用 Prodigal 在组装好的序列上预测蛋白质编码基因。



##### 分箱 (Binning) - 核心步骤:

“分箱”是利用序列的核苷酸组成频率 (四聚体频率) 和测序深度覆盖度，将属于同一物种的 contigs 聚类到一起的过程。

流程支持多种主流软件：MetaBAT2, MaxBin2, CONCOCT, COMEBin, MetaBinner, 和 SemiBin2。可以开启共组装丰度计算（co-abundance），利用多样本间的丰度协方差提高分箱准确率。

分箱质控与去冗余 (Bin QC & Refinement):

使用 Busco, CheckM 或 CheckM2 评估提取到的基因组的完整度 (Completeness) 和污染度 (Contamination)。

可开启 DAS Tool 将上述多种分箱软件的结果进行合并择优，获得更高质量的 bins。


##### 物种分类与下游注释 (Taxonomy & Annotation):

使用 GTDB-Tk 和/或 CAT 将提取到的高质量基因组比对到 GTDB 数据库，赋予其界门纲目科属种的分类学地位。

支持识别特定序列：使用 geNomad 挖掘潜在的病毒序列，使用 Tiara 识别真核生物序列。

2. 如何运行
   你需要准备一个 samplesheet.csv 来告诉流程你的数据在哪里：

Bash

##### 运行示例命令

nextflow run nf-core/mag -profile docker 
  --input samplesheet.csv 
  --outdir ./mag_results 
  --run_spades 
  --binning_map_mode group # 启用基于群组的丰度映射，优化分箱效果

### 三、 不依赖组装的方法学：基于 nf-core/taxprofiler (v2.0.1)


nf-core/taxprofiler 是一个高度并行的多重物种分类和丰度谱分析流程，专为直接从 reads 层面进行物种定量而设计。



#### 核心分析步骤 (v2.0.1 架构)



##### * 数据质控与宿主去除 (QC & Host-read removal):

除了常规的接头和低复杂度过滤外，该流程着重强调宿主序列去除。在临床（如人粪便、血液）或宿主相关环境样本中，大量 reads 来自宿主。

流程利用 BowTie2 或 Minimap2 将数据比对到宿主参考基因组并剔除，防止假阳性干扰下游分析。

宏基因组覆盖度评估 (Coverage Estimation):

使用 Nonpareil 工具评估当前测序深度是否已经覆盖了样本中的大部分微生物多样性。



##### 高度并行的物种分类与定量 (Taxonomic Classification & Profiling):

这是该流程的最大亮点。你可以同时运行多个算法并引入多个数据库，进行交叉验证。v2.0.1 支持高达 10 余种业界顶级的 Profiler：

基于 K-mer: Kraken2, KrakenUniq, Centrifuge, MetaCache, sylph （速度极快，适合超大数据库）。

基于标记基因 (Marker-gene): MetaPhlAn, mOTUs （特异性极高，丰度估计极其准确，但只评估已知物种的核心基因）。

基于氨基酸/DNA比对: DIAMOND, MALT, Kaiju （在远缘未知物种的鉴定上具有优势）。



##### 丰度后处理 (Post-processing):

例如在使用 Kraken2 之后，自动调用 Bracken 将分配到高分类层级（如“属”或“科”）的 reads 丰度，利用贝叶斯算法重新分配到“种”水平。



##### 数据标准化 (Standardisation):

不同的 Profiler 输出表格格式千奇百怪。流程使用 Taxpasta 将所有软件的结果统一格式化为标准化的 BIOM 或 TSV 表格，极大地降低了下游进行统计分析的难度。



##### 可视化 (Visualization):

整合 Krona 生成交互式的嵌套饼图，直观展示样本内的物种层级结构；结合 MultiQC 输出整体质控报告。



##### 如何运行


不同于 mag，taxprofiler 除了样本表 samplesheet.csv，还需要一个数据库配置表 databases.csv，用于指定你下载好的 Kraken2、MetaPhlAn 等软件对应的参考数据库路径。


##### 运行示例命令

nextflow run nf-core/taxprofiler -profile docker 
  --input samplesheet.csv 
  --databases databases.csv 
  --outdir ./taxprofiler_results 
  --run_kraken2 
  --run_metaphlan 
  --run_krona

### 四： 总结与应用场景选择


选择 nf-core/mag (依赖组装) 

如果： 你的目标是挖掘未知的“微生物暗物质”、获得新菌株的完整基因组 (MAGs)、探究微生物之间的代谢互养关系，或者你的样本来自未被深入研究的极端环境。

选择 nf-core/taxprofiler (不依赖组装) 

如果： 你的目标是快速得知已知环境（如人类肠道、常见土壤）中的物种群落组成及其相对丰度，进行大规模临床样本的病原体筛查，或者用于疾病组与健康组的标志物差异比较研究。
