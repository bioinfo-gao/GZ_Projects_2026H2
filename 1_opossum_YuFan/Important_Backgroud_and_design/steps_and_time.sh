# 样品存放在嵌套的子目录中（每个样品一个文件夹），我们需要对之前的脚本进行调整，使其能够自动进入每个文件夹寻找 .fq.gz 文件。
# 以下是适配你文件结构的完整工作流脚本，全部建议使用 Vim 进行编辑。
# 第一阶段：准备工作（Index 与 Liftoff）
# 预计耗时：4 - 6 小时

mamba install -c bioconda hisat2 -y
mamba install -c bioconda liftoff -y

which hisat2-build
which liftoff
which samtools

vim stage1_prep.sh

# # 在 Vim 中输入：
# #!/bin/bash
# # 1. 构建 HISAT2 索引
# echo "Building HISAT2 index..."
# hisat2-build -p 16 dv-2k.fasta dv_index

# # 2. 运行 Liftoff 迁移注释
# echo "Running Liftoff..."
# REF_FA="/Work_bio/references/Monodelphis_domestica/MonDom5/ncbi_refseq/mDidVir1.genome.fa"
# REF_GTF="/Work_bio/references/Monodelphis_domestica/MonDom5/ncbi_refseq/mDidVir1.annotation.gtf"
# TARGET_FA="/Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/dv-2k.fasta"

# liftoff -g $REF_GTF -o Didelphis_v.liftoff.gtf -u unmapped.txt -p 16 -sc 0.85 $TARGET_FA $REF_FA

# tmux attach -t sftp_work

# # 启动命令：
# bash stage1_prep.sh

# 1. 确保在南美灰负鼠的工作目录
cd /Work_bio/references/Monodelphis_domestica/MonDom5/ncbi_refseq

# 2. 自动提取 GFF3 中所有的特征类型（比如 mRNA, exon, CDS 等），生成 features.txt
awk -F'\t' '$1 !~ /^#/ {print $3}' clean_ref_annotation.gff3 | sort | uniq > feature_types.txt

# 3. 彻底删除刚才运行中断产生的数据库缓存
rm -f *.db *.db-shm *.db-wal

# 4. 使用 vim 修改运行脚本
vim run_liftoff.sh





















# 第二阶段：批量比对脚本（适配嵌套目录）
# 预计耗时：10 - 12 小时
# 这个脚本会自动遍历 NC_1 到 pi5_4 的子文件夹。
# code
# Bash
# vim run_mapping.sh
# 在 Vim 中输入：
# code
# Bash
# #!/bin/bash
# DATA_DIR="/home/gao/Dropbox/Quote_2605011001"
# INDEX="/Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/dv_index"
# OUT_DIR="./bam_files"
# mkdir -p $OUT_DIR

# # 样品列表
# SAMPLES=("NC_1" "NC_2" "NC_3" "NC_4" "pi5_1" "pi5_2" "pi5_3" "pi5_4")

# for smp in "${SAMPLES[@]}"
# do
#     echo "Processing sample: $smp"
#     # 使用 find 自动查找该目录下的 R1 和 R2 文件
#     R1=$(find $DATA_DIR/$smp -name "*_1.fq.gz")
#     R2=$(find $DATA_DIR/$smp -name "*_2.fq.gz")
    
#     # HISAT2 比对 -> 排序 -> 转 BAM
#     hisat2 -p 16 --dta -x $INDEX -1 $R1 -2 $R2 | \
#     samtools view -@ 4 -bS - | \
#     samtools sort -@ 4 -o ${OUT_DIR}/${smp}_sorted.bam
    
#     # 建立索引
#     samtools index ${OUT_DIR}/${smp}_sorted.bam
# done
# 第三阶段：定量与差异分析（DESeq2）
# 1. 定量 (featureCounts)
# code
# Bash
# vim run_counts.sh
# 输入内容：
# code
# Bash
# #!/bin/bash
# featureCounts -p -T 16 -t exon -g gene_id \
#               -a Didelphis_v.liftoff.gtf \
#               -o counts_matrix.txt \
#               ./bam_files/*_sorted.bam
# 2. R 差异分析（适配你的样本名）
# code
# Bash
# vim run_deseq2.R
# 在 Vim 中输入核心代码：
# code
# R
# library(DESeq2)
# # 读取数据，check.names=F 防止样本名中的横杠被改掉
# counts <- read.table("counts_matrix.txt", header=T, row.names=1, check.names=F)
# # 提取 count 部分（从第 6 列开始，因为前面的列是 GeneID, Chr, Start 等）
# countData <- counts[, 6:ncol(counts)]

# # 设置分组：前 4 个是 NC (Control)，后 4 个是 pi5 (Treat)
# colData <- data.frame(
#     condition = factor(c(rep("NC", 4), rep("pi5", 4))),
#     row.names = colnames(countData)
# )

# dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
# dds <- DESeq(dds)
# # 提取 pi5 相对于 NC 的差异结果
# res <- results(dds, contrast=c("condition", "pi5", "NC"))
# # 排序并保存
# resOrdered <- res[order(res$padj),]
# write.csv(as.data.frame(resOrdered), "DE_results_NC_vs_pi5.csv")
# 💡 管理建议与结果交付
# 监控进程：
# 你现在的任务是在后台跑的，如果想看完成了多少，可以查 BAM 文件夹：
# code
# Bash
# ls -lh ./bam_files/ | grep ".bam" | wc -l
# (如果显示 8，说明比对全部完成)
# 客户结果交付：
# 分析完成后，直接将最终的 CSV 结果发给客户。利用你已经设好的 Bore.pub 隧道：
# code
# Bash
# sudo cp DE_results_NC_vs_pi5.csv /srv/external_data/client_sftp/files/
# 客户操作：
# 客户连接 bore.pub:54833（或其他当前端口）。
# 在 files 目录下双击下载 DE_results_NC_vs_pi5.csv 即可在 Excel 中查看。
# 权限确认：
# 由于你的原始数据在 /home/gao/Dropbox，请确保该路径对 root 是可读的，因为你运行 hisat2 时通常带有 sudo 或者是以管理员身份运行。
# 现在你可以启动第一阶段了！ 如果今晚开始运行，预计明天下午你就能拿到所有的差异分析报表。