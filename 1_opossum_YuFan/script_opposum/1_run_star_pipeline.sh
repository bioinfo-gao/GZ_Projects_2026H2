#!/bin/bash

# ==================== 参数配置 ====================
DATA_DIR="/home/gao/Dropbox/Quote_2605011001_OP"
GENOME_DIR="/Work_bio/references/Didelphis_virginiana/mDidVir1/STAR_index"
OUT_DIR="/home/gao/projects_2026H2/1_opossum_YuFan/star_alignment"
# 补上由 Liftoff 刚刚生成的全新 GTF 绝对路径
GTF_FILE="/Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/Didelphis_v.liftoff.gtf"

mkdir -p "$OUT_DIR"
samples=("NC_1" "NC_2" "NC_3" "NC_4" "pi5_1" "pi5_2" "pi5_3" "pi5_4")
MAX_JOBS=2

# ==================== 核心控流循环 ====================
for sample in "${samples[@]}"; do
    
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 5
    done
    
    R1=$(ls ${DATA_DIR}/${sample}/*_1.fq.gz)
    R2=$(ls ${DATA_DIR}/${sample}/*_2.fq.gz)
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 🚀 启动样本比对: ${sample}"
    
    # 加入了 --sjdbGTFfile 参数，确保 STAR 吃进 Liftoff 的注释
    STAR --runThreadN 16 \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand zcat \
         --sjdbGTFfile "$GTF_FILE" \
         --outFileNamePrefix "${OUT_DIR}/${sample}_" \
         --outSAMtype BAM SortedByCoordinate \
         --twopassMode Basic \
         --quantMode GeneCounts > "${OUT_DIR}/${sample}_star_run.log" 2>&1 &

done


