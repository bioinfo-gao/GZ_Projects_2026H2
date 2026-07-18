#!/usr/bin/env bash
# 编排：等 Bowtie2 宿主索引（idx17）建好再启动 taxprofiler（复用共享索引）；
# 最多等 90 分钟，超时则不再等，让 taxprofiler 从 fasta 自建索引兜底。
set -uo pipefail
PROJ=/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics
IDXDIR=/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/bowtie2_index
for i in $(seq 1 90); do
    if ls "$IDXDIR"/GRCm39.*.bt2l >/dev/null 2>&1 || ls "$IDXDIR"/GRCm39.*.bt2 >/dev/null 2>&1; then
        echo "[$(date)] host index ready after ${i} min"; break
    fi
    sleep 60
done
cd "$PROJ"
export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity'
export NXF_ANSI_LOG=false
# 直接在本 tmux session 内以 TMUX 已存在的方式跑（3_ 脚本检测到 TMUX 则前台执行）
exec bash "$PROJ/scripts/3_run_taxprofiler.sh"
