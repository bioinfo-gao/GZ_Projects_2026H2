#!/bin/bash
# ============================================================================
# Study A / Step 2b — 优雅停止 + 明确带 -resume 重启（一次性运维操作，已执行于 2026-07-08）
#   背景：Study A 跑了10h13m只完成34/72比对子任务，排查发现 local_resources.config
#   的 queueSize=2（全局并发上限，非"2样本并行"）+ BWAMEM2_MEM 内存申报50GB（实测仅
#   ~18GB）双重限制导致本可并发的任务被挤到近乎串行。精调配置后（queueSize→3，
#   BWAMEM2_MEM内存→24GB，见 0_common/local_resources.config 当前版本）需要重启生效。
#
#   ⚠️ 关键点：必须明确带 -resume，不能走 A2_run_sarek_somatic.sh 默认的
#   "先跑一次不带resume，失败了再resume"逻辑——那样第一次不带resume的跑法会
#   丢弃 work_A 里已缓存的完成任务，白白浪费已跑的10小时。
#
#   本脚本记录当时实际执行的操作，供以后复用/审阅，不是自动化的常规脚本
#   （不会自动判断"现在是否需要重启"——需要时人工确认后手动跑）。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
SHEET="$PROJ/scripts/0_common/A_somatic.csv"
REF="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
NF="/home/gao/.conda/envs/regular_bioinfo/bin/nextflow"
CONFIG="$PROJ/scripts/0_common/local_resources.config"
SESSION="p14_sarek_A"

echo ">> Step 1: 优雅停止当前运行（发 Ctrl-C，非强杀）"
tmux send-keys -t "$SESSION" C-c
sleep 8
echo ">>   确认 nextflow 干净退出（应看到 'Execution complete -- Goodbye'）:"
tail -5 "$PROJ/.nextflow.log" 2>/dev/null
echo ">>   确认无残留 bwa-mem2/nextflow 进程:"
ps aux | grep -E "bwa-mem2|nextflow" | grep -v grep || echo "   （干净，无残留）"

echo
echo ">> Step 2: 带 -resume 显式重启（新 tmux session，避免依赖脚本默认的先跑后resume逻辑）"
mkdir -p "$PROJ/logs"
tmux new-session -d -s "$SESSION" "cd '$PROJ' && export NXF_OPTS='-Xms512m -Xmx2g' && export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity' && '$NF' run nf-core/sarek -r 3.8.1 -profile singularity -c '$CONFIG' --input '$SHEET' --outdir '$PROJ/output_A' -work-dir '$PROJ/work_A' --fasta '$REF' --fasta_fai '$REF.fai' --igenomes_ignore --genome null --aligner bwa-mem2 --skip_tools baserecalibrator --tools mutect2,tiddit --wes false -resume 2>&1 | tee -a '$PROJ/logs/sarek_A.log'"

echo "已在 tmux '$SESSION' 里带 -resume 重启。"
echo "验证缓存命中: grep 'cached' $PROJ/logs/sarek_A.log | head -3"
