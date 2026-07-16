#!/usr/bin/env bash
# ============================================================================
# proj16 常驻看门狗 —— 必须在 tmux 里跑（会话名 pan_watch），不受 agent/SSH 断开影响。
#
#   tmux new-session -d -s pan_watch 'bash scripts/10_watchdog.sh'
#
# 为什么不直接用 facilities/Server/nextflow_watchdog.sh：
#   那个通用看门狗只有 SESSION_END / LOW_MEM / HEARTBEAT 三条触发规则，会在 heartbeat 时
#   exit 12 等上层复检 —— 夜里没有上层，就等于没监控。且 2026-07-16 的 ELC 空转事故属于
#   「会话活着、内存充足、单线程磨洋工、~53 核空转」的**健康但病态慢**状态，三条规则一条不命中。
# 本脚本：(a) 永不退出、自我循环；(b) 补两条停滞检测规则（nflog 停滞 / CPU 长期低占用）。
# 只告警、绝不自行 kill —— 杀不杀是需要判断的决策，留给人。
# ============================================================================
set -uo pipefail

PROJ=/home/gao/projects_2026H2/16_Wenliang_Pan_wgs
SESS=pan_wgs                      # 被监控的 sarek 会话
NFLOG=$PROJ/.nextflow.log
WD=$PROJ/logs/watchdog.log
ALERT=$PROJ/logs/watchdog_ALERTS.log

POLL=${WD_POLL:-600}              # 巡检间隔 10 min（CLAUDE.md Phase 2: 10–20 min）
STALE_MIN=${WD_STALE_MIN:-45}     # .nextflow.log 停滞超过 N 分钟 → 告警
LOW_LOAD=${WD_LOW_LOAD:-10}       # load 低于此值视为"没在干活"
LOW_LOAD_STREAK=${WD_LOW_STREAK:-3}   # 连续 3 次（30 min）低 load 才告警，避阶段切换误报

mkdir -p "$PROJ/logs"
alert() {  # $1=事件名 $2=详情
  local ts; ts=$(date '+%F %T')
  echo "[$ts] 🚨 $1 — $2" | tee -a "$ALERT" >> "$WD"
}
note() { echo "[$(date '+%F %T')] $*" >> "$WD"; }

note "watchdog 启动 (poll=${POLL}s stale_thresh=${STALE_MIN}m low_load<${LOW_LOAD}x${LOW_LOAD_STREAK})"
low_streak=0; stale_alerted=0; load_alerted=0

while true; do
  # --- 1) 会话还在吗 ---
  if ! tmux has-session -t "$SESS" 2>/dev/null; then
    if tail -n 40 "$PROJ/logs/sarek_run.log" 2>/dev/null | grep -qiE 'Pipeline completed successfully|Succeeded'; then
      note "✅ SESSION_END=COMPLETED — sarek 正常结束，看门狗退出（下游由 pan_down 接手）"
    else
      alert "SESSION_END=FAILED" "pan_wgs 会话消失但日志无成功标志 —— 需排查 logs/sarek_run.log"
    fi
    exit 0
  fi

  # --- 2) 指标采集 ---
  age=$(( ( $(date +%s) - $(stat -c %Y "$NFLOG" 2>/dev/null || echo 0) ) / 60 ))
  load=$(awk '{printf "%d", $1}' /proc/loadavg)
  avail=$(free -m | awk '/^Mem:/{print $7}')
  run=$(grep -c "status: RUNNING" <(tail -n 200 "$NFLOG" 2>/dev/null) 2>/dev/null || echo 0)
  note "nflog+${age}m load=${load} avail=${avail}MB"

  # --- 3) 规则 A：nextflow 日志停滞（无 task 事件）= 流程没在推进 ---
  if [ "$age" -ge "$STALE_MIN" ]; then
    [ "$stale_alerted" -eq 0 ] && alert "STALLED" ".nextflow.log 已 ${age} 分钟无更新（阈值 ${STALE_MIN}m）—— 疑似单任务卡死/空转，查 tmux capture-pane -t $SESS -p | tail -15"
    stale_alerted=1
  else stale_alerted=0; fi

  # --- 4) 规则 B：会话活着但 CPU 长期低占用 = 健康但病态慢（2026-07-16 ELC 事故的特征）---
  if [ "$load" -lt "$LOW_LOAD" ]; then low_streak=$((low_streak+1)); else low_streak=0; load_alerted=0; fi
  if [ "$low_streak" -ge "$LOW_LOAD_STREAK" ] && [ "$load_alerted" -eq 0 ]; then
    alert "LOW_CPU_UTILIZATION" "load=${load} 已持续 $((low_streak*POLL/60)) 分钟低于 ${LOW_LOAD} —— 核在空转（ELC 型单线程阻塞？）。这不是"慢"，是设计出了问题，需人工判断"
    load_alerted=1
  fi

  # --- 5) 规则 C：低内存 ---
  if [ "${avail:-999999}" -lt 2500 ]; then alert "LOW_MEM" "avail=${avail}MB —— OOM 风险"; fi

  sleep "$POLL"
done
