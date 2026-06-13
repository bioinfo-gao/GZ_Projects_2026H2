#!/bin/bash

# ==================== 参数配置 ====================
DATA_DIR="/home/gao/Dropbox/Quote_2605011001_OP"
# 请确保此路径指向你为北美负鼠建立的 STAR Index 文件夹
GENOME_DIR="/Work_bio/references/Didelphis_virginiana/mDidVir1/STAR_index"
OUT_DIR="/home/gao/Code/Bioinfo_Analysis_Projects/star_alignment"

# 创建输出文件夹
mkdir -p "$OUT_DIR"

# 定义待处理的 8 个样本名
samples=("NC_1" "NC_2" "NC_3" "NC_4" "pi5_1" "pi5_2" "pi5_3" "pi5_4")

# 最大并发样本数限制
MAX_JOBS=2

# ==================== 核心控流循环 ====================
for sample in "${samples[@]}"; do
    
    # 动态检测当前正在运行的后台任务数 (jobs -r)
    # 如果运行中的任务数达到了 MAX_JOBS (2个)，则暂停脚本，每 5 秒检查一次
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 5
    done
    
    # 精准获取当前样本的 R1 和 R2 压缩文件绝对路径（自动匹配长文件名）
    R1=$(ls ${DATA_DIR}/${sample}/*_1.fq.gz)
    R2=$(ls ${DATA_DIR}/${sample}/*_2.fq.gz)
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 🚀 启动样本比对: ${sample}"
    echo " -> R1: $R1"
    echo " -> R2: $R2"
    
    # 运行 STAR 并将其投入后台 (&)
    # 每个样本分配 16 线程，日志和报错信息重定向到对应样本的 log 中
    STAR --runThreadN 16 \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${OUT_DIR}/${sample}_" \
         --outSAMtype BAM SortedByCoordinate \
         --twopassMode Basic \
         --quantMode GeneCounts > "${OUT_DIR}/${sample}_star_run.log" 2>&1 &

done

# 等待最后留在后台的样本全部运行完毕
wait
echo "[$(date '+%Y-%m-%d %H:%M:%S')] 🎉 所有样本 STAR 2-pass 比对全部平稳完成！"