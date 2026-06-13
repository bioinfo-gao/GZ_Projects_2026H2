#!/bin/bash
source ~/.bashrc
mamba activate regular_bioinfo

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ⏳ 监控启动：等待 Liftoff 任务完成..."

# 每 30 秒检查一次，只要系统里还有属于 gao 的 liftoff 或 gffread 进程在跑，就继续等
while pgrep -u gao -x liftoff > /dev/null || pgrep -u gao -x gffread > /dev/null; do
    sleep 30
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✅ 进程监测显示 Liftoff 已结束，验证产物..."

# 检查最终的 GTF 文件是否存在且大小不为 0
GTF_FILE="/Work_bio/references/Monodelphis_domestica/MonDom5/ncbi_refseq/Didelphis_v.liftoff.gtf"

if [ -s "$GTF_FILE" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 🎉 找到非空注释文件，自动拉起 STAR 并行管线！"
    cd /home/gao/Code/Bioinfo_Analysis_Projects
    bash run_star_pipeline.sh
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ❌ 错误：未检测到生成的 GTF 文件，STAR 未触发。请检查 Liftoff 是否报错中断。"
fi