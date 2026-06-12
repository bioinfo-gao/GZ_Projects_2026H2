# 接受你的最终工作流建议：（DE Analysis）的最终工作流建议：
# 建立索引：
# 使用 hisat2-build 为 DNA Zoo 基因组建立索引。
# 生成注释：
# 运行一次 Liftoff，把南美灰负鼠的基因映射到 DNA Zoo 基因组上，得到 Dv_ref.gtf。
# 定量：
# 使用 FeatureCounts 或 StringTie 进行表达定量。
# 定量后的基因 ID 就会是南美负鼠的基因名（如 TP53, GAPDH）
# 差异分析：
# 将 count 矩阵导入 DESeq2 或 edgeR。　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　这是我的南美灰负鼠的refence 文件列表： ls /Work_bio/references/Monodelphis_domestica/MonDom5/ncbi_refseq
# didelphis_virginiana_genome.zip  mDidVir1.annotation.gff  mDidVir1.genome.fa.0123  mDidVir1.genome.fa.bwt.2bit.64  ncbi_dataset_all_versions
# md5sum.txt                       mDidVir1.annotation.gtf  mDidVir1.genome.fa.amb   mDidVir1.genome.fa.fai          README.md
# mDidVir1.annotation.agat.log     mDidVir1.genome.fa       mDidVir1.genome.fa.ann   mDidVir1.genome.fa.pac          star_index   　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　这是我的北美负鼠的refence 文件列表：ls /Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/
# dv-2k.fasta       dv-2k.fasta.amb  dv-2k.fasta.bwt.2bit.64  dv-2k.fasta.pac    star_index
# dv-2k.fasta.0123  dv-2k.fasta.ann  dv-2k.fasta.fai          run_fixed_star.sh                              请给我各个步骤的详细命令并使用合理的参数，我的服务器硬盘足够大8T, 128G内存


# 既然你拥有 128G 内存 和 8T 硬盘，且确定了使用 Liftoff 同源迁移注释 + HISAT2 比对 的工作流，这是目前处理“高质量基因组+无注释”最专业且最高效的方案。
# 由于你的 dv-2k.fasta 有 50 万条序列，HISAT2 比 STAR 更适合处理这种极其破碎的基因组。
# 以下是使用 Vim 编写的全套详细命令：
# 第一阶段：构建 HISAT2 索引（针对北美负鼠）
# HISAT2 对碎片化基因组的内存压力远小于 STAR，128G 内存绰绰有余。
# 进入目标目录：
# code
# Bash
# cd /Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/
# 创建索引脚本：
# code
# Bash
# vim run_hisat2_index.sh
# 在 Vim 中按 i 插入：
# code
# Bash
# #!/bin/bash
# # 构建 HISAT2 索引，耗时约 1-2 小时
# hisat2-build -p 16 dv-2k.fasta dv_index
# 按 Esc 输入 :wq 保存，然后运行：
# code
# Bash
# bash run_hisat2_index.sh
# 第二阶段：利用 Liftoff 迁移注释（从南美负鼠到北美负鼠）
# 这是最关键的一步，通过同源性把南美负鼠的基因位置“画”到北美负鼠的基因组上。
# 安装 Liftoff (如果尚未安装)：
# code
# Bash
# pip install liftoff
# 编写迁移脚本：
# code
# Bash


# vim run_liftoff.sh
# 插入以下内容（路径已根据你的列表调整）：
# code
# Bash
# #!/bin/bash
# # 目标基因组 (Didelphis)
# TARGET="/Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/dv-2k.fasta"
# # 参考基因组 (Monodelphis)
# REF_FA="/Work_bio/references/Monodelphis_domestica/MonDom5/ncbi_refseq/mDidVir1.genome.fa"
# # 参考注释 (Monodelphis)
# REF_GTF="/Work_bio/references/Monodelphis_domestica/MonDom5/ncbi_refseq/mDidVir1.annotation.gtf"

# liftoff -g $REF_GTF \
#         -o Didelphis_v.liftoff.gtf \
#         -u unmapped_features.txt \
#         -p 16 \
#         -sc 0.85 \
#         $TARGET $REF_FA
# 参数 -sc 0.85：表示序列识别的最小阈值（由于物种不同，建议设在 0.85 左右）。
# 第三阶段：RNA-seq 数据比对 (HISAT2)
# 假设你的 RNA-seq 数据存放在 ~/data/rnaseq/ 目录下。
# 编写比对脚本：
# code
# Bash
# vim align_samples.sh
# 插入以下内容：
# code
# Bash
# #!/bin/bash
# INDEX="/Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/dv_index"
# # 以单样本为例，建议使用循环处理多样本
# hisat2 -p 16 --dta \
#        -x $INDEX \
#        -1 sample_R1.fastq.gz \
#        -2 sample_R2.fastq.gz | samtools view -bS - | samtools sort -@ 8 -o sample_sorted.bam
# 第四阶段：表达量定量 (featureCounts)
# 使用 Liftoff 生成的 GTF 文件进行定量，这样你的基因 ID 就会自动与南美负鼠（RefSeq）对应，方便做后续的功能分析。
# 运行定量命令：
# code
# Bash
# featureCounts -p -T 16 \
#               -t exon \
#               -g gene_id \
#               -a Didelphis_v.liftoff.gtf \
#               -o counts_matrix.txt \
#               *.bam
# 第五阶段：差异分析 (DESeq2)
# 你需要在服务器上运行 R 脚本。
# 编写 R 脚本：
# code
# Bash
# vim de_analysis.R
# 插入核心代码（示例）：
# code
# R
# library(DESeq2)
# counts <- read.table("counts_matrix.txt", header=T, row.names=1)
# # 过滤掉前 5 列非 count 信息
# countData <- counts[, 6:ncol(counts)]
# # 构建分组信息（根据你的实验设计修改）
# colData <- data.frame(row.names=colnames(countData), condition=c("control","control","treat","treat"))
# dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
# dds <- DESeq(dds)
# res <- results(dds)
# write.csv(as.data.frame(res), "DE_results_Didelphis.csv")

# 特别注意：
# 内存监控： 运行 Liftoff 时，由于它需要同时比对两套基因组，内存占用可能会瞬间冲到 80G+，你的 128G 内存刚好可以完美覆盖。
# 基因组 ID： Liftoff 运行结束后，记得检查一下 Didelphis_v.liftoff.gtf 的内容，确保基因 ID 是你熟悉的南美负鼠名称。
# 一切准备就绪，按照脚本顺序运行即可！