#!/bin/bash
# ============================================================================
# Study A / Step 2 — nf-core/sarek SOMATIC（RO 谱系，配对 vs origin）
#   patient=RO；origin=normal(0)，B1TP/B2TP/3tumor=tumor(1)
#   sarek 对每个 tumor 跑 Mutect2 vs normal + TIDDIT SV。跑在合并混合参考上。
#   self-relaunch tmux + 失败自动 -resume（复用 proj13）。Manta 弃用（proj13 教训）。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
SCRIPT="$PROJ/scripts/study_A/A2_run_sarek_somatic.sh"
SHEET="$PROJ/scripts/0_common/A_somatic.csv"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_Cas9_iHPV.fa"
NF="/home/gao/.conda/envs/regular_bioinfo/bin/nextflow"
SESSION="p14_sarek_A"

[ -f "$SHEET" ]  || { echo "ERROR: 缺 $SHEET（先跑 1_make_samplesheets.py）"; exit 1; }
[ -f "$HYBRID" ] || { echo "ERROR: 缺混合参考（先跑 0c）"; exit 1; }

if [ -z "${TMUX:-}" ]; then
    mkdir -p "$PROJ/logs"
    tmux kill-session -t "$SESSION" 2>/dev/null || true
    tmux new-session -d -s "$SESSION" "bash '$SCRIPT' 2>&1 | tee '$PROJ/logs/sarek_A.log'"
    echo "已在 tmux '$SESSION' 启动。日志: $PROJ/logs/sarek_A.log"; exit 0
fi

mkdir -p "$PROJ/logs" "$PROJ/output_A" "$PROJ/work_A"
export NXF_OPTS='-Xms512m -Xmx2g'
export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity'
run() {
    "$NF" run nf-core/sarek -r 3.8.1 -profile singularity \
        -c "$PROJ/scripts/0_common/local_resources.config" \
        --input "$SHEET" --outdir "$PROJ/output_A" -work-dir "$PROJ/work_A" \
        --fasta "$HYBRID" --fasta_fai "$HYBRID.fai" --igenomes_ignore --genome null \
        --aligner bwa-mem2 --skip_tools baserecalibrator \
        --tools mutect2,tiddit --wes false "$@"
}
echo "== sarek A 首跑 $(date) =="; if run; then exit 0; fi
echo "== 失败 -resume $(date) =="; run -resume
