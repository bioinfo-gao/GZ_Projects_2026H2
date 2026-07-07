#!/bin/bash
# ============================================================================
# Step 2 — nf-core/sarek on hybrid reference（比对 + QC + 去重 + SV）
#   self-relaunch 进 tmux + 失败自动 -resume（强制模式）。
#
#   用法:
#     bash 2_run_sarek.sh                         # 默认试跑 RAGH（samplesheet_trial_RAGH.csv）
#     bash 2_run_sarek.sh samplesheet_full.csv    # 明天 CD1A 到位后跑全部 6 样
#
#   引擎决策见 docs/analysis_plan_0706.md：
#     - 指向合并 hybrid 参考（--fasta）
#     - --aligner bwa-mem2
#     - --skip_tools baserecalibrator（无 known-sites + NovaSeq 收益极小）
#     - --tools manta,tiddit（BND 断点 = 整合结合部；模式A 可另加 germline caller）
#     - 跳过 snpEff/VEP 人类注释（断点自行注释到 GENCODE vM35）
# ============================================================================
set -euo pipefail

PROJ="/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
SCRIPT="$PROJ/scripts/2_run_sarek.sh"
SHEET="${1:-$PROJ/scripts/samplesheet_trial_RAGH.csv}"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_constructs.fa"
NF="/home/gao/.conda/envs/regular_bioinfo/bin/nextflow"
SESSION="ellen_sarek"

[ -f "$HYBRID" ]     || { echo "ERROR: hybrid 参考不存在，先跑 0_build_hybrid_ref.sh"; exit 1; }
[ -f "$HYBRID.fai" ] || { echo "ERROR: 缺 .fai，重跑 0_build_hybrid_ref.sh"; exit 1; }
[ -f "$SHEET" ]      || { echo "ERROR: samplesheet 不存在: $SHEET（先跑 1_produce_samplesheet.py）"; exit 1; }

# 不在 tmux → 建 session 把自身投入，主 shell 退出
if [ -z "${TMUX:-}" ]; then
    mkdir -p "$PROJ/logs"   # 必须在 tmux/tee 启动前建好，否则 tee 启动即失败、无日志
    tmux kill-session -t "$SESSION" 2>/dev/null || true
    tmux new-session -d -s "$SESSION" "bash '$SCRIPT' '$SHEET' 2>&1 | tee '$PROJ/logs/sarek_run.log'"
    echo "已在 tmux '$SESSION' 后台启动。samplesheet: $SHEET"
    echo "看日志: tail -f $PROJ/logs/sarek_run.log   或  tmux capture-pane -t $SESSION -p"
    exit 0
fi

# 在 tmux 内 → 执行 sarek
mkdir -p "$PROJ/logs" "$PROJ/output_results" "$PROJ/work"
export NXF_OPTS='-Xms512m -Xmx2g'
export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity'

run_sarek() {
    "$NF" run nf-core/sarek \
        -r 3.8.1 \
        -profile singularity \
        -c "$PROJ/scripts/local_resources.config" \
        --input "$SHEET" \
        --outdir "$PROJ/output_results" \
        -work-dir "$PROJ/work" \
        --fasta "$HYBRID" \
        --fasta_fai "$HYBRID.fai" \
        --igenomes_ignore \
        --genome null \
        --aligner bwa-mem2 \
        --skip_tools baserecalibrator \
        --tools manta,tiddit \
        --wes false \
        "$@"
}

echo "===== sarek 首次运行 $(date) ====="
if run_sarek; then echo "===== 完成（无需 resume）====="; exit 0; fi
echo "===== 失败，自动 -resume $(date) ====="
if run_sarek -resume; then echo "===== resume 后完成 ====="; exit 0; fi
echo "===== 两次均失败，请查日志 ====="; exit 1
