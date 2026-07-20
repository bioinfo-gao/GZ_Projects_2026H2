#!/usr/bin/env bash
# Phase 2 编排 — 等 HUMAnN 全部跑完(机器空出)后，自动启动 nf-core/mag(脚本 4)。
# 用户 2026-07-18 已批准 Phase 2 MAG，要求"机器可执行时开始"→ 不与 HUMAnN(54线程满载)叠跑。
# 完成判据：logs/humann_run.log 出现 HUMANN_ALL_DONE 且 10 个 *_pathabundance.tsv 齐。
set -uo pipefail
PROJ=/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics
LOG="$PROJ/logs/mag_orchestrate.log"
HRUN="$PROJ/logs/humann_run.log"
POLL=300   # 5 min

log(){ echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

log "orchestrator started; waiting for HUMAnN to finish before launching MAG."

# ---- 1) 等 HUMAnN 完成 ----
while :; do
  ndone=$(find "$PROJ/humann_out" -maxdepth 2 -name "*_pathabundance.tsv" 2>/dev/null | wc -l)
  if grep -q "HUMANN_ALL_DONE" "$HRUN" 2>/dev/null && [ "$ndone" -ge 10 ]; then
    log "HUMAnN DONE (marker + ${ndone}/10 pathabundance). proceeding to MAG."
    break
  fi
  # 若 humann tmux 已死但没到 ALL_DONE，告警但继续等(可能在 join/renorm 收尾)
  if ! tmux has-session -t humann17 2>/dev/null && ! grep -q "HUMANN_ALL_DONE" "$HRUN" 2>/dev/null; then
    log "WARN: humann17 tmux gone but no HUMANN_ALL_DONE yet (done=${ndone}/10). still waiting ${POLL}s; check if stalled."
  fi
  log "waiting… HUMAnN ${ndone}/10 done; load=$(uptime | sed 's/.*load average: //')"
  sleep "$POLL"
done

# ---- 2) 确认 CheckM2 库(非阻塞：最多再等 10 分钟，到点就不带 checkm2 直接上 MAG) ----
for i in $(seq 1 20); do
  if ls /Work_bio/references/Metagenomics/checkm2/*.dmnd >/dev/null 2>&1; then
    log "CheckM2 DB ready → MAG 将带 --run_checkm2."; break
  fi
  [ "$i" -eq 20 ] && log "CheckM2 DB 仍未就位 → MAG 用 BUSCO-only(脚本4自动跳过 checkm2)."
  sleep 30
done

# ---- 3) 启动 MAG(脚本4 自建 tmux mag17) ----
log "launching MAG (scripts/4_run_mag.sh → tmux mag17)…"
bash "$PROJ/scripts/4_run_mag.sh" 2>&1 | tee -a "$LOG"

# ---- 4) 启动自检(Phase-1)：90s 后确认 mag17 活着且日志在动 ----
sleep 90
if tmux has-session -t mag17 2>/dev/null; then
  log "mag17 alive. recent log:"; tail -n 8 "$PROJ/logs/mag_run.log" 2>/dev/null | tee -a "$LOG"
  if grep -qiE "error|exception|command not found|no such file" "$PROJ/logs/mag_run.log" 2>/dev/null; then
    log "WARN: early error keyword in mag_run.log — 需人工核查。"
  fi
else
  log "ERROR: mag17 tmux 未起来 — MAG 启动失败，需人工核查 logs/mag_run.log。exiting watchdog."
  exit 1
fi

# ---- 5) Phase-2 常驻看门狗：每 15min 查 mag17 存活 + .nextflow.log mtime 是否推进 ----
NFLOG="$PROJ/.nextflow.log"
STALL=0
while tmux has-session -t mag17 2>/dev/null; do
  sleep 900
  tmux has-session -t mag17 2>/dev/null || break
  now=$(date +%s); mt=$(stat -c %Y "$NFLOG" 2>/dev/null || echo "$now"); age=$(( (now-mt)/60 ))
  ld=$(uptime | sed 's/.*load average: //')
  last=$(tail -n 1 "$NFLOG" 2>/dev/null | cut -c1-120)
  if [ "$age" -ge 30 ]; then
    STALL=$((STALL+1))
    log "WARN STALLED: .nextflow.log 停更 ${age}min (连续${STALL}次); load=$ld; last: $last"
    [ "$STALL" -ge 2 ] && log "ALERT: 疑似卡死(≥30min×2)—人工核查 mag17/是否 futex 死锁/OOM。看门狗继续观察不擅自杀。"
  else
    STALL=0
    log "MAG healthy: nflog 更新于 ${age}min 前; load=$ld"
  fi
done
log "mag17 session ended — MAG 收尾/结束。检查 output_results_mag/ 与 MultiQC。orchestrator+watchdog done."
