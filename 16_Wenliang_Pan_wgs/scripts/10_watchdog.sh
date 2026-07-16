#!/usr/bin/env bash
# ============================================================================
# proj16 看门狗 —— 薄封装，直接复用共享实现（不再维护重复副本）。
#   共享实现: facilities/Server/nextflow_watchdog.sh（2026-07-16 已并入 rule 13/14/15 + WD_PERSIST）
#
# 启动（必须 tmux，过夜不受 agent/SSH 断开影响）：
#   tmux new-session -d -s pan_watch 'bash scripts/10_watchdog.sh'
#
# 阈值标定（2026-07-16 15:10 实测修正，第一版两条规则都标错了）：
#   - WD_STALE_MIN=45 ：nflog mtime 是【存活心跳】(ELC 空转时仍每 2.5min 一写)，
#                       故此规则只能抓「nextflow 本体死了」，抓不到「流程没进展」。保留但别指望它。
#   - WD_LOW_LOAD_MIN=240(4h)：常规 GATK4_MARKDUPLICATES 近单线程、合法跑 ~3.5h/样本时 load 仅 ~4-6。
#                       阈值必须长于最长合法单线程阶段，否则每次 markdup 都狼来了 →「告警狼来了 = 没有告警」。
#   - WD_MAX_TASK_H=5  ：唯一能区分「ELC 空转 6.4h(无界)」与「markdup 3.5h(有界合法)」的信号。最有用的一条。
# 详见：analysis_plan_0715.md §7 反面教训 3、facilities/Server/nextflow_pipeline/resume缓存失效与空转监控_教训_0716.md
# ============================================================================
set -uo pipefail
PROJ=/home/gao/projects_2026H2/16_Wenliang_Pan_wgs
SHARED=/home/gao/projects_2026H2/facilities/Server/nextflow_watchdog.sh

exec env \
  WD_PERSIST=1 \
  WD_POLL=600 \
  WD_LOG="$PROJ/logs/watchdog.log" \
  WD_ALERT_LOG="$PROJ/logs/watchdog_ALERTS.log" \
  WD_STALE_MIN=45 \
  WD_LOW_LOAD=10 \
  WD_LOW_LOAD_MIN=240 \
  WD_MAX_TASK_H=5 \
  bash "$SHARED" "pan_wgs:$PROJ:$PROJ/logs/sarek_run.log"
