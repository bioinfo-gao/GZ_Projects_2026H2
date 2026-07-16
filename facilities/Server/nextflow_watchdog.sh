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
# ⚠ 过夜/无人值守必须这样跑（2026-07-16 proj16 事故教训，见下）：
#   tmux new-session -d -s <proj>_watch \
#     'WD_PERSIST=1 bash .../nextflow_watchdog.sh <sess>:<projdir>:<logfile>'
#
# 环境变量（可选）：
#   WD_LOG           巡检流水日志（默认 /tmp/nextflow_watchdog.log）
#   WD_ALERT_LOG     告警日志（默认 <WD_LOG 同目录>/nextflow_watchdog_ALERTS.log）— 异常与流水分离，一眼可见
#   WD_POLL          巡检间隔秒（默认 180）
#   WD_HEARTBEAT     心跳唤醒轮数（默认 30，即 30*180s=90min）
#   WD_LOW_MEM_MB    低内存阈值 MB（默认 2500，连续 2 次触发）
#   WD_STALE_MIN     .nextflow.log 停滞告警阈值 分钟（默认 45；0=禁用）        ← 2026-07-16 新增
#   WD_LOW_LOAD      load 低于此值视为"没在干活"（默认 10；0=禁用）             ← 2026-07-16 新增
#   WD_LOW_LOAD_POLLS 低 load 连续几轮才告警（默认 10 轮，即 180s*10=30min）    ← 2026-07-16 新增
#   WD_PERSIST       =1 则告警后**继续巡检不退出**（过夜必须；默认 0=沿用旧的退出语义）
#
# 退出码：10=某 session 结束  11=低内存  12=健康心跳  13=流程停滞  14=CPU 空转
#         （均意在触发上层复检；WD_PERSIST=1 时除 SESSION_END 外均不退出）
#
# ---------------------------------------------------------------------------
# 2026-07-16 新增规则 13/14 的由来（proj16 事故，净亏 ~15.4h，其中 6.4h 是空转）：
#   `--use_gatk_spark markduplicates` 触发单线程 GATK4_ESTIMATELIBRARYCOMPLEXITY 空转 6h21m，
#   状态是「tmux 会话活着 + 内存充足 + 单线程磨洋工 + ~53 核空转」——
#   **原有 SESSION_END/LOW_MEM/HEARTBEAT 三条规则一条都不命中**，看门狗只会默默记下 nflog+300m 却不告警。
#   且 HEARTBEAT 时 exit 12 等上层复检，**夜里没有上层 = 等于没监控**（故有 WD_PERSIST）。
#   教训：绝大多数监控只查"死没死"，而最贵的故障是**「活着但空转」**——它不触发任何传统告警。
#   详见 facilities/Server/nextflow_pipeline/resume缓存失效与空转监控_教训_0716.md
# ============================================================================
set -uo pipefail
LOG="${WD_LOG:-/tmp/nextflow_watchdog.log}"
ALERT_LOG="${WD_ALERT_LOG:-$(dirname "$LOG")/nextflow_watchdog_ALERTS.log}"
POLL="${WD_POLL:-180}"
HEARTBEAT_POLLS="${WD_HEARTBEAT:-30}"
LOW_MEM_MB="${WD_LOW_MEM_MB:-2500}"
STALE_MIN="${WD_STALE_MIN:-45}"
LOW_LOAD="${WD_LOW_LOAD:-10}"
LOW_LOAD_POLLS="${WD_LOW_LOAD_POLLS:-10}"
PERSIST="${WD_PERSIST:-0}"
[ "$#" -ge 1 ] || { echo "用法: $0 <session>:<projdir>:<logfile> [...]"; exit 2; }

classify_end() {  # $1=logfile -> COMPLETED|FAILED
  if tail -n 40 "$1" 2>/dev/null | grep -qiE 'Pipeline completed successfully|resume 后完成|完成（无需 resume）|Succeeded'; then
    echo COMPLETED; else echo FAILED; fi
}
alert() {  # $1=事件名 $2=详情 ; 写告警日志 + 流水日志
  local ts; ts=$(date '+%F %T')
  echo "[$ts] 🚨 $1 — $2" | tee -a "$ALERT_LOG" >> "$LOG"
}
# 事件处置：PERSIST=1 时告警后继续；否则沿用旧语义退出触发上层复检
handle() {  # $1=事件名 $2=详情 $3=退出码
  alert "$1" "$2"
  [ "$PERSIST" = "1" ] || exit "$3"
}

low_streak=0; poll=0; stale_alerted=0; load_alerted=0
echo "[$(date '+%F %T')] watchdog 启动: $# 个作业, poll=${POLL}s heartbeat=$((HEARTBEAT_POLLS*POLL/60))min persist=${PERSIST} stale>${STALE_MIN}m low_load<${LOW_LOAD}x${LOW_LOAD_POLLS}" >> "$LOG"
while true; do
  poll=$((poll+1)); ts=$(date '+%F %T')
  avail=$(free -m | awk '/^Mem:/{print $7}'); swap=$(free -m | awk '/^Swap:/{print $3}')
  load=$(awk '{printf "%d", $1}' /proc/loadavg)
  line="[$ts] avail=${avail}MB swap=${swap}MB load=${load}"; ended=""; max_age=0
  for spec in "$@"; do
    sess=${spec%%:*}; rest=${spec#*:}; proj=${rest%%:*}; lf=${rest#*:}
    st=UP; tmux has-session -t "$sess" 2>/dev/null || st=GONE
    nl="$proj/.nextflow.log"
    age=$(( ($(date +%s) - $(stat -c %Y "$nl" 2>/dev/null || echo 0)) / 60 ))
    line="$line | $sess:$st nflog+${age}m"
    [ "$st" = UP ] && [ "$age" -gt "$max_age" ] && max_age=$age
    if [ "$st" = GONE ]; then ended="$ended $sess=$(classify_end "$lf")"; fi
  done
  echo "$line" >> "$LOG"

  # --- 10: 某 session 结束（PERSIST 下也退出——没有可监控对象了） ---
  if [ -n "$ended" ]; then
    case "$ended" in *FAILED*) alert "SESSION_END" "$ended —— 无成功标志，需排查" ;;
                     *) echo "[$ts] EVENT=SESSION_END$ended" >> "$LOG" ;; esac
    exit 10
  fi
  # --- 11: 低内存 ---
  if [ "${avail:-999999}" -lt "$LOW_MEM_MB" ]; then low_mem_streak=$(( ${low_mem_streak:-0} + 1 )); else low_mem_streak=0; fi
  if [ "${low_mem_streak:-0}" -ge 2 ]; then handle "LOW_MEM" "avail=${avail}MB —— OOM 风险" 11; low_mem_streak=0; fi

  # --- 13: 流程停滞（.nextflow.log 无 task 事件）=== 2026-07-16 新增 ---
  if [ "$STALE_MIN" -gt 0 ] && [ "$max_age" -ge "$STALE_MIN" ]; then
    if [ "$stale_alerted" -eq 0 ]; then
      stale_alerted=1
      handle "STALLED" ".nextflow.log 已 ${max_age} 分钟无更新（阈值 ${STALE_MIN}m）—— 疑似单任务卡死/空转，查 tmux capture-pane -t <sess> -p | tail -15" 13
    fi
  else stale_alerted=0; fi

  # --- 14: CPU 空转（会话活着但 load 长期低）=== 2026-07-16 新增 ---
  if [ "$LOW_LOAD" -gt 0 ] && [ "$load" -lt "$LOW_LOAD" ]; then low_streak=$((low_streak+1)); else low_streak=0; load_alerted=0; fi
  if [ "$LOW_LOAD" -gt 0 ] && [ "$low_streak" -ge "$LOW_LOAD_POLLS" ] && [ "$load_alerted" -eq 0 ]; then
    load_alerted=1
    handle "LOW_CPU_UTILIZATION" "load=${load} 已持续 $((low_streak*POLL/60)) 分钟低于 ${LOW_LOAD} —— 会话活着但核在空转（ELC 型单线程阻塞？）。这不是'慢'，是设计出了问题，需人工判断" 14
  fi

  # --- 12: 健康心跳（放最后：真事件优先） ---
  if [ "$poll" -ge "$HEARTBEAT_POLLS" ]; then
    echo "[$ts] EVENT=HEARTBEAT_HEALTHY" >> "$LOG"
    [ "$PERSIST" = "1" ] && poll=0 || exit 12
  fi
  sleep "$POLL"
done
