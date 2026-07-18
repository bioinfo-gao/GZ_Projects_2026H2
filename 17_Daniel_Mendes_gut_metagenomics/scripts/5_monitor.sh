#!/usr/bin/env bash
# 快速状态检查：各 tmux 日志末 15 行 + .nextflow.log mtime + 输出目录。
PROJ=/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics
cd "$PROJ"
for s in tax17 mag17 db17 idx17; do
    if tmux has-session -t "$s" 2>/dev/null; then
        echo "===== tmux $s (last 15) ====="
        tmux capture-pane -t "$s" -p 2>/dev/null | tail -n 15
    fi
done
echo; echo "===== .nextflow.log mtime ====="
for f in .nextflow.log work_taxprofiler/../.nextflow.log; do :; done
ls -la --time-style=+%H:%M:%S .nextflow.log 2>/dev/null; date +%H:%M:%S
echo; echo "===== outputs ====="
ls output_results/ 2>/dev/null; ls output_results_mag/ 2>/dev/null
echo; echo "===== load ====="; uptime; free -h | head -2