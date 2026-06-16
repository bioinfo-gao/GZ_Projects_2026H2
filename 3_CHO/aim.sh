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

# bwa-mem2 index 被 OOM Killed 了 — 因为 STAR 正在占用 ~80GB RAM，bwa-mem2 index 也需要大量内存（~20-30GB），超出了剩余内存。

# 两个方案：

# 方案	操作	等待时间
# A：用常规 bwa index	内存需求小得多（~5-8GB），现在就能跑	约 30 分钟
# B：等 STAR 跑完再建 bwa-mem2 index	NC_2/NC_3 快完成了	等 1-2 小时
# 建议用方案 A，常规 bwa 功能完全一样，只是速度稍慢，对你这个任务没有影响：

cd /Work_bio/references/Cricetulus_griseus/CriGri-PICR/ncbi_refseq && \
echo "=== Building bwa index ===" && \
bwa index GCF_003668045.1_CriGri-PICR_genomic.fna && \
echo "=== Done! ===" && \
ls -lh
# 建好索引后，比对时用 bwa mem 代替 bwa-mem2 mem 即可，分析结果完全相同。