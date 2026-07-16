#!/usr/bin/env bash
# Ad-hoc status check for the sarek run (last 15 lines convention).
PROJ=/home/gao/projects_2026H2/16_Wenliang_Pan_wgs
echo "=== tmux sessions ==="; tmux ls 2>/dev/null | grep pan_wgs || echo "(no pan_wgs session)"
echo "=== .nextflow.log mtime ==="; ls -l --time-style=+%H:%M:%S "$PROJ"/.nextflow.log 2>/dev/null
echo "=== last 15 lines of run log ==="; tail -n 15 "$PROJ/logs/sarek_run.log" 2>/dev/null
echo "=== load / mem ==="; uptime; free -h | head -2
