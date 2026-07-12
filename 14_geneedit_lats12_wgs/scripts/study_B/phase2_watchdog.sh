#!/bin/bash
# ============================================================================
# Phase 2 看门狗 —— 在 tmux 中运行,生命周期与用户 SSH 连接解耦(可断开睡觉)。
#   职责:等 Study A(p14_sarek_A)完成 → 把 Study B(p14_sarek_B)从过渡配置
#         (cpus=40)升级到满机配置(cpus=56),带 -resume 复用 work_B、零进度丢失。
#   保守原则:任何一步不确定就【不动 B】,让 B 保持过渡配置继续安全跑(过渡配置本身
#   就是完整可跑的,升级只是提速,不升不会坏事)。绝不在旧 nextflow 未确认退出时启新的
#   (防两进程同写 work_B 损坏)。
# ============================================================================
set -uo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
LOGA="$PROJ/logs/sarek_A.log"
LOGB="$PROJ/logs/sarek_B.log"
WLOG="$PROJ/logs/phase2_watchdog.log"
FULL="$PROJ/scripts/0_common/local_resources_full_machine.config"
SHEET="$PROJ/scripts/0_common/B_germline.csv"
REF="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
NF="/home/gao/.conda/envs/regular_bioinfo/bin/nextflow"
log(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$WLOG"; }

log "===== phase2 watchdog 启动:等 Study A 完成后升级 Study B 到满机配置 ====="

# --- 1) 等 Study A 完成(只认成功标记;session消失但无成功标记=可能失败→保守退出)---
while true; do
  if grep -q "Pipeline completed successfully" "$LOGA" 2>/dev/null; then
    log "Study A 已成功完成。"; break
  fi
  if ! tmux has-session -t p14_sarek_A 2>/dev/null; then
    sleep 10
    if grep -q "Pipeline completed successfully" "$LOGA" 2>/dev/null; then log "Study A 已完成。"; break; fi
    log "⚠ Study A session 消失但无成功标记 → 可能失败/异常。不升级 B,保持过渡配置安全跑。退出。"; exit 0
  fi
  sleep 120
done

# --- 2) 确认 Study B 仍在跑(要升级的对象)---
if ! tmux has-session -t p14_sarek_B 2>/dev/null; then
  log "⚠ Study B session 不在(已自行结束?)→ 不操作。退出。"; exit 0
fi

# --- 3) 优雅停止 Study B 的 nextflow(Ctrl-C,让其 checkpoint 干净退出)---
log "向 p14_sarek_B 发送 C-c,等待 nextflow 干净退出..."
tmux send-keys -t p14_sarek_B C-c
CLEAN=0
for i in $(seq 1 24); do   # 最多等 ~4 分钟
  sleep 10
  if grep -qE "Goodbye|Execution cancelled|Execution complete|Pipeline completed" "$LOGB" 2>/dev/null; then CLEAN=1; break; fi
done
sleep 5
tmux kill-session -t p14_sarek_B 2>/dev/null || true
sleep 5

# --- 4) 硬确认没有残留的 B nextflow/工作进程(防两进程同写 work_B)---
# 只要还有引用 work_B 的进程,就不敢启新的 → 保守退出(B 已停,用户重连后手动升级)
if pgrep -af "work_B|output_B" >/dev/null 2>&1; then
  log "⚠ 仍检出引用 work_B 的残留进程,不敢启新 nextflow(防 work-dir 损坏)。"
  log "  已停止旧 B。用户重连后请手动:带 -resume + $FULL 重启 B。退出。"
  pgrep -af "work_B|output_B" | tee -a "$WLOG"
  exit 0
fi
[ "$CLEAN" = 1 ] && log "旧 B nextflow 已干净退出。" || log "未见明确退出标记,但无残留进程,谨慎继续。"

# --- 5) 满机配置 + -resume 重启 Study B(命令与 B2 run() 完全一致,仅换 config 并加 -resume)---
log "以满机配置 $FULL 带 -resume 重启 Study B(复用 work_B)..."
tmux new-session -d -s p14_sarek_B "cd '$PROJ' && \
  export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity' NXF_OPTS='-Xms512m -Xmx2g' && \
  '$NF' run nf-core/sarek -r 3.8.1 -profile singularity -c '$FULL' \
    --input '$SHEET' --outdir '$PROJ/output_B' -work-dir '$PROJ/work_B' \
    --fasta '$REF' --fasta_fai '$REF.fai' --igenomes_ignore --genome null \
    --aligner bwa-mem2 --skip_tools baserecalibrator --tools haplotypecaller,tiddit --wes false \
    -resume 2>&1 | tee -a '$LOGB'"
sleep 8
if tmux has-session -t p14_sarek_B 2>/dev/null; then
  log "✅ Study B 已在满机配置(cpus=56)下带 -resume 重启。Phase 2 完成。"
else
  log "❌ 重启后 session 未起来!请用户重连后手动检查 $LOGB。"
fi
log "===== watchdog 退出 ====="
