# 1. 启动终端复用会话，确保后台稳定运行
tmux new -s download_opossum

# 2. 激活你的常规生信环境，以调用相应的索引工具
mamba activate regular_bioinfo

# 3. 安装 ncbi-datasets-cli 工具 (如果之前没有安装过)
mamba install -c conda-forge -c bioconda ncbi-datasets-cli -y

# 3. 创建北美负鼠的标准参考目录并进入
mkdir -p /Work_bio/references/Didelphis_virginiana/mDidVir1/ncbi_refseq
cd /Work_bio/references/Didelphis_virginiana/mDidVir1/ncbi_refseq

# 4. 下载最新 Didelphis virginiana 基因组 + 注释（GTF）
datasets download genome taxon 9265 --include genome,gtf,protein,seq-report --filename didelphis_virginiana_genome.zip

# 5. 解压并重命名文件，对齐人类目录中的命名规范
unzip didelphis_virginiana_genome.zip -d didelphis_virginiana

# 进入目录查看
cd didelphis_virginiana/ncbi_dataset/data

# /mnt/ex_8T_SSD/references/Didelphis_virginiana/mDidVir1/ncbi_refseq/didelphis_virginiana/ncbi_dataset/data]:
# gao@us1 $ ll
# total 80
# -rw------- 1 gao gao 36542 Jun  8 21:05 assembly_data_report.jsonl
# -rw------- 1 gao gao  4509 Jun  8 21:15 dataset_catalog.json
# drwxrwxr-x 2 gao gao  4096 Jun  8 21:25 GCA_000002295.1
# drwxrwxr-x 2 gao gao  4096 Jun  8 21:25 GCA_016433145.1
# drwxrwxr-x 2 gao gao  4096 Jun  8 21:25 GCA_027887165.1
# drwxrwxr-x 2 gao gao  4096 Jun  8 21:25 GCA_027887165.2
# drwxrwxr-x 2 gao gao  4096 Jun  8 21:25 GCA_027917375.1
# drwxrwxr-x 2 gao gao  4096 Jun  8 21:25 GCA_904810665.1
# drwxrwxr-x 2 gao gao  4096 Jun  8 21:25 GCF_000002295.2
# drwxrwxr-x 2 gao gao  4096 Jun  8 21:25 GCF_016433145.1
# drwxrwxr-x 2 gao gao  4096 Jun  8 21:25 GCF_027887165.1

# 3. 理解 GCF vs. GCA 的命名玄学
# NCBI 的命名本身就带有强烈的质量导向：
# GCF (RefSeq)：代表 NCBI 参考序列数据库。这是 NCBI 经过审核、清洗、并进行标准化注释的基因组，通常质量最高，被全球科研界视为“标准参考”。
# GCA (GenBank)：代表 提交者的原始提交。可能是某个测序中心直接上传的原始数据，质量参差不齐。
# 首选 GCF_027887165.1 (数字越新通常越好)：
# 如果 GCF_027887165.1 的 BUSCO 分数更高且 Contig N50 更大，请直接使用它。这是因为

# 简单来说：GCF_016433145.1 是经典的“骨架版”，而 GCF_027887165.1 通常是基于更新的测序技术（如长读长测序 PacBio/Oxford Nanopore）得到的“高精度版”。

# 1. 为什么会有多个 GCF 版本？
# NCBI 的 RefSeq 数据库会根据提交的组装质量进行更新。通常：

# GCF_016433145.1：可能是基于较早的 Illumina 短读长数据组装的。

# GCF_027887165.1：极有可能是利用了近年来流行的长读长测序（Long-read sequencing）技术重新进行的组装。

# 2. 如何判断谁更适合你的分析？
# 你可以运行下面这行命令，对比两个版本在报告中的差异：


# /mnt/ex_8T_SSD/references/Didelphis_virginiana/mDidVir1/ncbi_refseq/didelphis_virginiana/ncbi_dataset/data
# # 使用 jq 工具提取关键指标（若未安装，可用 cat 查看 assembly_data_report.jsonl） #jq - commandline JSON processor [version 1.8.1]
# # 查看 GCF_016433145.1 的元数据
# grep "GCF_016433145.1" assembly_data_report.jsonl | jq '.assembly_stats'
# # 查看 GCF_027887165.1 的元数据
# grep "GCF_027887165.1" assembly_data_report.jsonl | jq '.assembly_stats'

# 6.cp and keep the same directory structure as Home genome directory
cd /Work_bio/references/Didelphis_virginiana/mDidVir1/ncbi_refseq/didelphis_virginiana/data

# 6.1. 选取最完整的版本 (GCF_027887165.1) 并提取关键文件
# 注意：NCBI 下载的文件通常在子目录下
cp GCF_027887165.1/*.fna ../../mDidVir1.genome.fa
cp GCF_027887165.1/*.gtf ../../mDidVir1.annotation.gff

# 7. 生成 Samtools 的 faidx 索引 (对应你人类目录里的 .fai 文件)
samtools faidx mDidVir1.genome.fa

# 8. 生成 bwa-mem2 索引 (对应你人类目录里的 .bwt.2bit.64 等文件)
# 注意：该步骤极其消耗内存，请确保服务器处于空闲状态
bwa-mem2 index mDidVir1.genome.fa

# 9. 为下一步的单细胞/转录组比对创建 star_index 空目录
mkdir star_index
# 如果你要接着直接构建 STAR 索引，请根据你的测序读长(如 --sjdbOverhang 149)执行以下命令 
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir ./star_index \
     --genomeFastaFiles mDidVir1.genome.fa \
     --sjdbGTFfile mDidVir1.annotation.gtf \
     --sjdbOverhang 149
     
# 注意：这个参数在 STAR 构建索引时的定义是：读长减去 1 (Read Length - 1)