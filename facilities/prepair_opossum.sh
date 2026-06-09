# 1. 启动终端复用会话，确保后台稳定运行
tmux new -s download_opossum

# 2. 激活你的常规生信环境，以调用相应的索引工具
mamba activate regular_bioinfo


# 3. 安装 ncbi-datasets-cli 工具 (如果之前没有安装过)
mamba install -c conda-forge -c bioconda ncbi-datasets-cli -y

# 3. 创建北美负鼠的标准参考目录并进入
mkdir -p /Work_bio/references/Didelphis_virginiana/mDidVir1/ncbi_refseq
cd /Work_bio/references/Didelphis_virginiana/mDidVir1/ncbi_refseq


# 5. 自动检索并下载北美负鼠（Didelphis virginiana）的基因组及相关注释文件
datasets download genome taxon "Didelphis virginiana" \
    --include gff3,rna,cds,protein,genome,seq-report \
    --filename D_virginiana_genome.zip

# 5. 解压并重命名文件，对齐人类目录中的命名规范
gunzip GCF_011100635.1_mDidVir1.pri_genomic.fna.gz
gunzip GCF_011100635.1_mDidVir1.pri_genomic.gtf.gz

mv GCF_011100635.1_mDidVir1.pri_genomic.fna mDidVir1.genome.fa
mv GCF_011100635.1_mDidVir1.pri_genomic.gtf mDidVir1.annotation.gtf

# 6. 生成 Samtools 的 faidx 索引 (对应你人类目录里的 .fai 文件)
samtools faidx mDidVir1.genome.fa

# 7. 生成 bwa-mem2 索引 (对应你人类目录里的 .bwt.2bit.64 等文件)
# 注意：该步骤极其消耗内存，请确保服务器处于空闲状态
bwa-mem2 index mDidVir1.genome.fa

# 8. 为下一步的单细胞/转录组比对创建 star_index 空目录
mkdir star_index

# (可选) 9. 如果你要接着直接构建 STAR 索引，请根据你的测序读长(如 --sjdbOverhang 149)执行以下命令
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir ./star_index \
     --genomeFastaFiles mDidVir1.genome.fa \
     --sjdbGTFfile mDidVir1.annotation.gtf \
     --sjdbOverhang 149