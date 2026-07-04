#!/bin/bash
# 监控 Nextflow 运行状态（可选）
# 用法: bash 3_Optional_work_monitor.sh

LOG="/home/gao/projects_2026H2/10_Yue_Liu/scripts/nextflow_run.log"

echo "=== tmux session 'rnaseq_yueliu' last 20 lines ==="
tmux capture-pane -t rnaseq_yueliu -p 2>/dev/null | tail -20 || echo "(session not found)"

echo ""
echo "=== Nextflow log tail ==="
tail -30 "$LOG" 2>/dev/null || echo "Log not found: $LOG"

echo ""
echo "=== Output results structure ==="
ls -lh /home/gao/projects_2026H2/10_Yue_Liu/output_results/ 2>/dev/null || echo "(not yet created)"

echo ""
echo "=== Key output files ==="
find /home/gao/projects_2026H2/10_Yue_Liu/output_results/star_salmon/ \
     -name "*.tsv" 2>/dev/null | head -10
