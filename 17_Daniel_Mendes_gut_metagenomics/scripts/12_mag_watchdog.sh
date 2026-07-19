#!/usr/bin/env bash
# Phase 2 常驻看门狗 — 监控 mag17：每 15min 查存活 + .nextflow.log mtime 是否推进 + load/内存。
# 健康但慢绝不擅自杀(见 corun playbook)；只记录/告警。mag17 结束即退出。
set -uo pipefail
PROJ=/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics
LOG="$PROJ/logs/mag_watchdog.log"
NFLOG="$PROJ/.nextflow.log"
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
log "MAG watchdog 启动，监控 mag17。"
STALL=0
while tmux has-session -t mag17 2>/dev/null; do
  sleep 900
  tmux has-session -t mag17 2>/dev/null || break
  now=$(date +%s); mt=$(stat -c %Y "$NFLOG" 2>/dev/null || echo "$now"); age=$(( (now-mt)/60 ))
  ld=$(uptime | sed 's/.*load average: //'); mem=$(free -g | awk '/Mem:/{print $3"/"$2"G used"}')
  last=$(grep -aE "Submitted process|Pipeline completed|ERROR|error" "$NFLOG" 2>/dev/null | tail -1 | sed 's/\x1b\[[0-9;]*m//g' | cut -c1-110)
  if [ "$age" -ge 40 ]; then
    STALL=$((STALL+1))
    log "WARN STALLED: .nextflow.log 停更 ${age}min (连续${STALL}); load=$ld; mem=$mem; last: $last"
    [ "$STALL" -ge 2 ] && log "ALERT: 疑似卡死(≥40min×2) — 人工核查(futex 死锁/OOM/容器拉取超时)。看门狗继续观察不擅自杀。"
  else
    STALL=0
    log "MAG healthy: nflog 更新于 ${age}min 前; load=$ld; mem=$mem; last: $last"
  fi
done
# 收尾判定
if grep -qa "Pipeline completed successfully" "$NFLOG" 2>/dev/null; then
  log "MAG 完成 ✅ (Pipeline completed successfully)。核查 output_results_mag/ + MultiQC + GTDB 分类。"
else
  log "mag17 结束但未见 'Pipeline completed successfully' — 可能失败/中断，核查 logs/mag_run.log 末尾。"
fi
log "MAG watchdog 退出。"
