/Work_bio/lhn_work/cellranger/Demo2/work_for_ath.sh
同事 Demo 流程结果详解
整个流程分两步，结果也分两部分：QC 结果（01.QC/）+ Cell Ranger VDJ 分析结果（N2-SI_TT_D1_23575YLT4/outs/）。

一、QC 阶段输出（Demo2/01.QC/N2-SI_TT_D1_23575YLT4/）
由 fastp 生成：

文件	含义
*_R1_001.fastq.gz / *_R2_001.fastq.gz	过滤后的 R1/R2（去除低质量、N过多、过短的 reads，去接头后）的测序数据，作为下游 Cell Ranger 的输入
*.json	结构化的 QC 统计（reads数、碱基数、Q20/Q30占比、重复率等），供程序解析
*.html	可视化版 QC 报告，人工查看用
从 test.e 日志可以看到关键质控指标：总 reads ~7700万对，通过过滤 7511万对（98%通过），Q30碱基占比 ~95%+，重复率 22.98%——说明原始测序数据质量很好。

二、Cell Ranger VDJ 分析输出（outs/ 目录）

/Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/web_summary.html

这是核心结果，按"概览报告 → 细胞/克隆型层面 → contig（重建出的抗体序列）层面 → 比对/参考文件"几个层次理解：

1. 概览报告
web_summary.html：最重要的总览报告，网页打开即可看到所有关键指标、图表（细胞数、UMI分布、V/J基因使用情况等）
metrics_summary.csv：上面报告里所有数字指标的 CSV 版（详见下表）
2. 克隆型（Clonotype）层面 —— 免疫组分析的核心结果
clonotypes.csv：把携带相同 CDR3 序列组合的细胞归为同一个"克隆型"（clonotype），每行代表一种克隆型：

frequency：该克隆型对应的细胞数
proportion：占总细胞的比例
cdr3s_aa / cdr3s_nt：该克隆型的重链(IGH)+轻链(IGK/IGL) CDR3 氨基酸/核酸序列（CDR3 是抗体识别抗原特异性的关键区域）
例如 clonotype1 有 52 个细胞共享同一对 IGH+IGK CDR3 序列，说明这是一个被显著扩增的 B 细胞克隆（很可能是针对某抗原产生免疫反应的优势克隆）。

consensus.fasta / consensus_annotations.csv：每个克隆型的"代表序列"（取该克隆型内多条 contig 的共识序列），及其 V/D/J/C 基因注释

consensus.bam(.bai)`：原始 reads 比对回 consensus 序列的结果，用于验证组装质量

3. Contig（单条重建出来的转录本）层面
"Contig" = 从每个细胞的 reads 组装出的一条完整抗体转录本序列（重链或轻链）。

filtered_contig_annotations.csv / all_contig_annotations.csv：每条 contig 的详细注释表，逐细胞逐contig列出：

barcode：细胞的10x条形码（细胞ID）
chain：链类型（IGH重链 / IGK κ轻链 / IGL λ轻链）
v_gene/d_gene/j_gene/c_gene：该contig比对到的V(D)J/恒定区基因片段
full_length/productive：是否全长、是否为有功能的（可翻译、无终止码、符合阅读框）抗体序列
fwr1~4/cdr1~3：抗体序列的框架区(FWR)和互补决定区(CDR)氨基酸及核酸序列——CDR3 列就是上面 clonotypes.csv 用来分组的依据
reads/umis：支持这条contig的read数/UMI数（反映表达量及置信度）
raw_clonotype_id：该contig所属的克隆型编号
filtered_ 版本只保留高置信度、判定为真实细胞的contig；all_ 版本包含所有（含背景/低质量）。

*.fasta / *.fastq（all_contig、filtered_contig）：对应的序列文件本身

*.bam(.bai)（all_contig.bam）：原始测序reads比对到组装出的contig上的结果，用于评估组装准确性

all_contig_annotations.json/.bed：JSON/BED 格式的同一份注释信息，方便程序化读取或在基因组浏览器中查看

4. 细胞判定与参考相关
cell_barcodes.json：被判定为"真实细胞"（而非空液滴/背景噪音）的===============条形码列表
concat_ref.fasta(.fai)/concat_ref.bam(.bai)：把比对用的V(D)J参考基因片段拼接起来，及contig比对到参考上的结果——用来看每条contig具体匹配到哪些已知基因
donor_regions.fa：推断出的"种系（germline）序列"，即去除体细胞高频突变后、该样本个体本身的V(D)J基因型序列
vdj_reference/：本次分析实际用到的参考库副本（reference.json + fasta），保证结果可追溯
vdj_contig_info.pb：所有contig信息的二进制(ProtoBuf)版本，供 Cell Ranger 内部/Loupe使用，不用人工查看
airr_rearrangement.tsv：按照免疫组学界通用的 AIRR 标准格式输出的重排序列表，方便用其他免疫组学分析工具（如 immcantation）继续分析
vloupe.vloupe：10x官方可视化软件 Loupe V(D)J Browser 的专用文件，可交互查看每个克隆型、每条contig、序列比对等（需要下载Loupe软件打开）
三、metrics_summary.csv 关键指标解读
指标	数值	含义
Estimated Number of Cells	2,267	判定为真实细胞的数量
Mean Read Pairs per Cell	16,568	平均每个细胞测了多少reads
Number of Read Pairs	37,559,174	总测序read pair数
Valid Barcodes	83.0%	条形码有效率
Reads Mapped to Any V(D)J Gene	96.6%	reads比对到V(D)J基因的比例（高=数据质量好）
Reads Mapped to IGH/IGK/IGL	73.9% / 22.2% / 51.7%	分别比对到重链/κ轻链/λ轻链的比例
Median IGH/IGK/IGL UMIs per Cell	60.0 / 1302.0 / 0.0	每细胞重链/κ链/λ链的UMI中位数——IGL为0很反常，说明这批细胞几乎都用κ轻链（IGK），没有λ轻链(IGL)，符合很多免疫库（如人源化小鼠/大鼠模型常偏向使用κ链）的特点
Cells With Productive V-J Spanning Pair	76.0%	同时有功能性重链+轻链配对的细胞比例（这是能形成完整抗体的细胞）
Paired Clonotype Diversity	217.06	克隆型多样性指数，越大说明免疫库越多样化（未被单一克隆主导）
简单总结：这个 demo 跑出来的结果显示——从约2267个B细胞中，重建出了每个细胞的抗体重链+轻链序列，按CDR3序列把它们分成了不同克隆型（其中最大的克隆型clonotype1占了2.3%的细胞），最终目的通常是为后续抗体序列分析/抗体药物候选筛选提供数据（找出扩增明显的优势克隆，对应序列可用于表达验证）。


clonotypes.csv  中 clonotype3	43	0.018967798853109837	IGK:CQKYNSAPWTF, 

也就是只有轻链单独序列也可以成为一个克隆型？ ===================  ????

Bash Find contigs belonging to clonotype3 to verify chain composition
IN
awk -F',' '$NF=="clonotype3" || $0 ~ /clonotype3,/' /Work_bio/lhn_work/cellranger/Demo2/N2-SI_TT_D1_23575YLT4/outs/filtered_contig_annotations.csv | grep ",clonotype3," | head -10

OUT
N2-SI_TT_D1_23575YLT4,AACGAGCCAATACGTC-1,true,AACGAGCCAATACGTC-1_contig_1,true,677,IGK,IGKV1-27,,IGKJ1,IGKC-01,true,true,DIQMTQSPSSLSASVGDRVTITC,GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGC,RASQGISNYLA,CGGGCGAGTCAGGGCATTAGCAATTATTTAGCC,WYQQKPGKVPKLLIY,TGGTATCAGCAGAAACCAGGGAAAGTTCCTAAGCTCCTGATCTAT,AASTLQS,GCTGCATCCACTTTGCAATCA,GVPSRFSGSGSGTDFTLTISSLQPEDVATYY,GGGGTCCCATCTCGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGCCTGCAGCCTGAAGATGTTGCAACTTATTAC,CQKYNSAPWTF,TGTCAAAAGTATAACAGTGCCCCGTGGACGTTC,GQGTKVEIK,GGCCAAGGGACCAAGGTGGAAATCAAAC,14446,4424,clonotype3,clonotype3_consensus_1,1

验证了一下：clonotype3 下的每个细胞条形码（如 AACGAGCCAATACGTC-1、AACTAGGCACATGAAG-1 等）在 filtered_contig_annotations.csv 里都只有一条 contig（_contig_1，IGK链），完全没有对应的 IGH contig。所以确认：

是的，Cell Ranger VDJ 的克隆型定义不要求重链+轻链都存在，单独的轻链（或单独的重链）也可以独立构成一个克隆型。

原理：

Cell Ranger 的克隆型分组规则是：把该细胞内所有被判定为"高置信度(high_confidence)+有功能性(productive)"的链的 CDR3 序列组合作为分组依据。
如果一个细胞只成功组装/判定出一条有效链（比如只有 IGK，没有 IGH），那它的"组合"里就只有这一条链，Cell Ranger 依然会把它当作一个独立的克隆型，用这一条链的 CDR3 作为标识。
这也解释了为什么整体指标里 Cells With Productive V-J Spanning Pair 只有 76.0%——意味着约 24% 的细胞没有同时拿到重链+轻链的完整配对，clonotype3 这种"单轻链"细胞就是这 24% 里的一部分。
常见原因：

测序深度/捕获效率问题：该细胞的重链转录本UMI/reads数太少，没能组装出 full-length 或未通过 productive 判定阈值（不是真的没有重链，只是没测到/没拼出来）。
生物学原因：部分B细胞确实存在重链或轻链表达异常（如等位基因排斥不完全、轻链编辑、细胞凋亡/降解中的细胞等），属于正常生物学噪音的一部分。
需要注意的是：clonotype3 这种"只有 IGK"的克隆型，严格来说不能代表一个完整、有功能的抗体（缺重链就无法配对形成完整Ig分子），如果你们后续要挑选候选抗体做表达验证，通常应该优先看 Cells With Productive V-J Spanning Pair 为 true 的（即配对完整的）克隆型，比如 clonotype1、clonotype2、clonotype4 这种同时有 IGH+IGK 的。