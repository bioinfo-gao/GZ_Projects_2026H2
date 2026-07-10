#!/bin/bash
# ============================================================================
# nextflow_watchdog.sh — 事件驱动看门狗，监控 1~N 个在 tmux 里跑的 nextflow 作业。
#
# 设计见 facilities/Server/concurrent_nextflow_resource_lessons_0710.md 与 skill /corun。
# 只监控 + 告警，绝不自行 kill/重启（那是需要判断的运维决策，交给人/agent）。
#
# 用法：
#   nextflow_watchdog.sh <SESSION>:<PROJECT_DIR>:<LOGFILE> [ <SESSION>:<PROJECT_DIR>:<LOGFILE> ... ]
# 例：
#   nextflow_watchdog.sh \
#     ellen_sarek:/home/gao/projects_2026H2/13_Ellen_knockin_wgs:/home/gao/projects_2026H2/13_Ellen_knockin_wgs/logs/sarek_run.log \
#     p14_sarek_A:/home/gao/projects_2026H2/14_geneedit_lats12_wgs:/home/gao/projects_2026H2/14_geneedit_lats12_wgs/logs/sarek_A.log
#
# 环境变量（可选）：
#   WD_LOG           监控日志路径（默认 /tmp/nextflow_watchdog.log）
#   WD_POLL          巡检间隔秒（默认 180）
#   WD_HEARTBEAT     心跳唤醒轮数（默认 30，即 30*180s=90min）
#   WD_LOW_MEM_MB    低内存阈值 MB（默认 2500，连续 2 次触发）
#
# 退出码：10=某 session 结束  11=低内存  12=健康心跳  （均意在触发上层复检）
# ============================================================================
set -uo pipefail
LOG="${WD_LOG:-/tmp/nextflow_watchdog.log}"
POLL="${WD_POLL:-180}"
HEARTBEAT_POLLS="${WD_HEARTBEAT:-30}"
LOW_MEM_MB="${WD_LOW_MEM_MB:-2500}"
[ "$#" -ge 1 ] || { echo "用法: $0 <session>:<projdir>:<logfile> [...]"; exit 2; }

classify_end() {  # $1=logfile -> COMPLETED|FAILED
  if tail -n 40 "$1" 2>/dev/null | grep -qiE 'Pipeline completed successfully|resume 后完成|完成（无需 resume）|Succeeded'; then
    echo COMPLETED; else echo FAILED; fi
}

low_streak=0; poll=0
echo "[$(date '+%F %T')] watchdog 启动: $# 个作业, poll=${POLL}s heartbeat=$((HEARTBEAT_POLLS*POLL/60))min" >> "$LOG"
while true; do
  poll=$((poll+1)); ts=$(date '+%F %T')
  avail=$(free -m | awk '/^Mem:/{print $7}'); swap=$(free -m | awk '/^Swap:/{print $3}')
  line="[$ts] avail=${avail}MB swap=${swap}MB"; ended=""
  for spec in "$@"; do
    sess=${spec%%:*}; rest=${spec#*:}; proj=${rest%%:*}; lf=${rest#*:}
    st=UP; tmux has-session -t "$sess" 2>/dev/null || st=GONE
    nl="$proj/.nextflow.log"
    age=$(( ($(date +%s) - $(stat -c %Y "$nl" 2>/dev/null || echo 0)) / 60 ))
    line="$line | $sess:$st nflog+${age}m"
    if [ "$st" = GONE ]; then ended="$ended $sess=$(classify_end "$lf")"; fi
  done
  echo "$line" >> "$LOG"

  if [ -n "$ended" ]; then echo "[$ts] EVENT=SESSION_END$ended" >> "$LOG"; exit 10; fi
  if [ "${avail:-999999}" -lt "$LOW_MEM_MB" ]; then low_streak=$((low_streak+1)); else low_streak=0; fi
  if [ "$low_streak" -ge 2 ]; then echo "[$ts] EVENT=LOW_MEM avail=${avail}MB" >> "$LOG"; exit 11; fi
  if [ "$poll" -ge "$HEARTBEAT_POLLS" ]; then echo "[$ts] EVENT=HEARTBEAT_HEALTHY" >> "$LOG"; exit 12; fi
  sleep "$POLL"
done
