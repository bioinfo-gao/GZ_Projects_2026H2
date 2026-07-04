#!/bin/bash
# 用法: bash 2_nextflow.sh
# 行为: 自动创建 tmux session → 首次运行 nextflow → 失败则自动 -resume 重试一次

# ── 不在 tmux 时自动创建 session 并把自身投入 ──────────────────
if [ -z "$TMUX" ]; then
    SCRIPT="$(realpath "$0")"
    cd /home/gao/projects_2026H2/12_Zhen_Yan/scripts
    tmux kill-session -t rnaseq_zhenyan 2>/dev/null || true
    tmux new-session -d -s rnaseq_zhenyan "bash '${SCRIPT}' 2>&1 | tee nextflow_run.log"
    echo "Nextflow 已在 tmux session 'rnaseq_zhenyan' 启动"
    echo "查看: tmux attach -t rnaseq_zhenyan"
    echo "日志: tail -f /home/gao/projects_2026H2/12_Zhen_Yan/scripts/nextflow_run.log"
    exit 0
fi

# ── 以下在 tmux 内部执行 ───────────────────────────────────────
export NXF_OPTS="-Xms512m -Xmx2g"
export NXF_SINGULARITY_CACHEDIR="/home/gao/.singularity/nf-core"
cd /home/gao/projects_2026H2/12_Zhen_Yan/scripts

run_nextflow() {
    nextflow run nf-core/rnaseq \
        -r 3.15.1 \
        -profile singularity \
        -c local_optimized.config \
        "$@" \
        --input nf_core_samplesheet.csv \
        --outdir ../output_results \
        --fasta /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/GRCh38.primary_assembly.genome.fa \
        --gtf /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/gencode.v45.annotation.gtf \
        --star_index '/Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/star_index' \
        --gencode \
        --aligner star_salmon \
        --max_cpus 28 \
        --max_memory '108.GB'
}

ts() { date '+%Y-%m-%d %H:%M:%S'; }

echo "[$(ts)] ══ First run ══════════════════════════════════"
if run_nextflow; then
    echo "[$(ts)] ✅ SUCCESS"
    exit 0
fi

echo ""
echo "[$(ts)] ⚠️  First run failed — auto-resuming with -resume ..."
echo "[$(ts)] ══ Resume run ═════════════════════════════════"
if run_nextflow -resume; then
    echo "[$(ts)] ✅ SUCCESS after auto-resume"
    exit 0
fi

echo ""
echo "[$(ts)] ❌ Failed after auto-resume. Manual check required."
echo "       Check: less nextflow_run.log"
echo "       Or:    nextflow log"
exit 1
