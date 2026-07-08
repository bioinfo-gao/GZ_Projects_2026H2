#!/bin/bash
# ============================================================================
# Study B / Step 2 — nf-core/sarek GERMLINE（6 样 Lats1/2 品系）
#   每样 germline SNV/indel（HaplotypeCaller）+ TIDDIT SV。跑在合并混合参考上
#   （L1L2H 在 TG_iHPV contig 有覆盖 → B3 整合位点）。self-relaunch tmux + resume。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
SCRIPT="$PROJ/scripts/study_B/B2_run_sarek_germline.sh"
SHEET="$PROJ/scripts/0_common/B_germline.csv"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_Cas9_iHPV.fa"
NF="/home/gao/.conda/envs/regular_bioinfo/bin/nextflow"
SESSION="p14_sarek_B"

[ -f "$SHEET" ]  || { echo "ERROR: 缺 $SHEET（先跑 1_make_samplesheets.py）"; exit 1; }
[ -f "$HYBRID" ] || { echo "ERROR: 缺混合参考（先跑 0c）"; exit 1; }

if [ -z "${TMUX:-}" ]; then
    mkdir -p "$PROJ/logs"
    tmux kill-session -t "$SESSION" 2>/dev/null || true
    tmux new-session -d -s "$SESSION" "bash '$SCRIPT' 2>&1 | tee '$PROJ/logs/sarek_B.log'"
    echo "已在 tmux '$SESSION' 启动。日志: $PROJ/logs/sarek_B.log"; exit 0
fi

mkdir -p "$PROJ/logs" "$PROJ/output_B" "$PROJ/work_B"
export NXF_OPTS='-Xms512m -Xmx2g'
export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity'
run() {
    "$NF" run nf-core/sarek -r 3.8.1 -profile singularity \
        -c "$PROJ/scripts/0_common/local_resources.config" \
        --input "$SHEET" --outdir "$PROJ/output_B" -work-dir "$PROJ/work_B" \
        --fasta "$HYBRID" --fasta_fai "$HYBRID.fai" --igenomes_ignore --genome null \
        --aligner bwa-mem2 --skip_tools baserecalibrator \
        --tools haplotypecaller,tiddit --wes false "$@"
}
echo "== sarek B 首跑 $(date) =="; if run; then exit 0; fi
echo "== 失败 -resume $(date) =="; run -resume
