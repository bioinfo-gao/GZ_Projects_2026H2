# /home/gao/Dropbox/Keqiang

# CHO genome wgs data
# Confirmation of CHO cell lines from Chinese hamster ovaries

# ls -lh /home/gao/Dropbox/Keqiang/*R1* /home/gao/Dropbox/Keqiang/*r1* /home/gao/Dropbox/Keqiang/*.fq.gz /home/gao/Dropbox/Keqiang/*.fastq.gz 2>/dev/null | head -30
ls -lh /home/gao/Dropbox/Keqiang/

which bwa samtools bedtools 2>/dev/null && echo "===" && ls /home/gao/Dropbox/Keqiang/ && echo "===" && df -h /home/gao/projects_2026H2/ | tail -1
which bwa-mem2 2>/dev/null && echo "bwa-mem2 found" || echo "bwa-mem2 NOT found"

# mkdir -p /home/gao/projects_2026H2/3_cho_wgs_species_confirm/{ref/cho,qc,align,results} && cd /home/gao/projects_2026H2/3_cho_wgs_species_confirm/ref/cho && echo "Downloading CHO genome (CriGri-PICR)..." && wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/GCF_003668045.3_CriGri-PICR/GCF_003668045.3_CriGri-PICR_genomic.fna.gz && wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/GCF_003668045.3_CriGri-PICR/GCF_003668045.3_CriGri-PICR_genomic.gff.gz && echo "Download complete!"


# 好，只有一个样本（wt1），数据量适中（33GB）。服务器目前 STAR 还在跑（~80GB RAM），但 BWA 内存需求小得多，可以同时做。

# 让我设计完整的实施脚本：

# Run
# Step 1: 快速 QC
# Step 2: 下载 CHO 参考基因组 + 建索引
# Step 3: 取子集做多物种比对 → 确认是 CHO
# Step 4: 检查 DHFR 基因座 → 确认株系

# 创建目录 + 下载 + 解压 + 建索引

cd /Work_bio/references/Cricetulus_griseus/CriGri-PICR/ncbi_refseq && \
echo "=== Downloading CHO genome ===" && \
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/GCF_003668045.1_CriGri-PICR/GCF_003668045.1_CriGri-PICR_genomic.fna.gz && \
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/668/045/GCF_003668045.1_CriGri-PICR/GCF_003668045.1_CriGri-PICR_genomic.gff.gz && \
echo "=== Decompressing ===" && \
gunzip -k GCF_003668045.1_CriGri-PICR_genomic.fna.gz && \
echo "=== Building samtools faidx ===" && \
samtools faidx GCF_003668045.1_CriGri-PICR_genomic.fna && \
echo "=== Building bwa-mem2 index (~30-60 min) ===" && \
bwa-mem2 index GCF_003668045.1_CriGri-PICR_genomic.fna && \
echo "=== Done! ===" && \
ls -lh