#!/bin/bash
# ============================================================================
# Step 3 — 运行状态速览 + 早期失败检测（CLAUDE.md：前 3 分钟主动轮询）
#   用法:
#     bash 3_work_monitor.sh            # 一次性快照
#     bash 3_work_monitor.sh watch      # 前 3 分钟每 30s 轮询查 error/成功标志
# ============================================================================
PROJ="/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
LOG="$PROJ/logs/sarek_run.log"
SESSION="ellen_sarek"

snapshot() {
    echo "===== tmux 会话 ====="; tmux ls 2>/dev/null | grep "$SESSION" || echo "  (无 $SESSION 会话)"
    echo "===== 日志末尾 ====="; tail -n 25 "$LOG" 2>/dev/null || echo "  (暂无日志)"
    echo "===== 输出目录 ====="; find "$PROJ/output_results" -maxdepth 2 -type d 2>/dev/null | sed "s#$PROJ/##" | head -30
}

if [ "${1:-}" = "watch" ]; then
    for i in 1 2 3 4 5 6; do
        sleep 30
        if grep -qiE "error|halted|failed|traceback|cannot|no such|exit status" "$LOG" 2>/dev/null; then
            echo "!! 可能失败（${i}x30s）——日志相关行:"; grep -iE "error|halted|failed|traceback|cannot|no such|exit status" "$LOG" | tail -15; exit 1
        fi
        if grep -qiE "Pipeline completed successfully" "$LOG" 2>/dev/null; then echo "✓ 完成"; exit 0; fi
        echo "--- 仍在运行 (${i}x30s) ---"; tail -n 3 "$LOG" 2>/dev/null
    done
    echo "前 3 分钟未见明显报错，转入中程监控（每 10-20 分钟看一次即可）。"
else
    snapshot
fi
