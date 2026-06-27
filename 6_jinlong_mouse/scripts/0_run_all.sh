#!/bin/bash
# Full automated pipeline — launch and leave
# Steps: 1) samplesheet → 2) nf-core/rnaseq (auto-resume) → 3) DE/PCA → 4) enrichment

PROJECT_DIR="/home/gao/projects_2026H2/6_jinlong_mouse"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

# ── Self-relaunch into tmux ───────────────────────────────────
if [ -z "$TMUX" ]; then
    SCRIPT="$(realpath "$0")"
    tmux kill-session -t rnaseq 2>/dev/null || true
    tmux new-session -d -s rnaseq \
        "bash '${SCRIPT}' 2>&1 | tee '${SCRIPTS_DIR}/pipeline_all.log'"
    echo "Pipeline started in tmux 'rnaseq'"
    echo "Attach: tmux attach -t rnaseq"
    echo "Log:    tail -f ${SCRIPTS_DIR}/pipeline_all.log"
    exit 0
fi

# ── Inside tmux ──────────────────────────────────────────────
ts()   { date '+%Y-%m-%d %H:%M:%S'; }
fail() { echo "[$(ts)] ❌ $1 — stopping."; exit 1; }
SEP="══════════════════════════════════════════════"

cd "$SCRIPTS_DIR" || fail "Cannot cd to $SCRIPTS_DIR"

echo "[$(ts)] $SEP"
echo "[$(ts)]  JINLONG MOUSE — FULL PIPELINE"
echo "[$(ts)] $SEP"
echo ""

# ── Step 1: Samplesheet ──────────────────────────────────────
echo "[$(ts)] ── Step 1 / 4 : Samplesheet ──────────────────"
conda run -n regular_bioinfo python 1_produce_nf-core_Samplesheet.py \
    || fail "Samplesheet generation failed"
echo "[$(ts)] ✅ Step 1 done"
echo ""

# ── Step 2: nf-core/rnaseq (auto-resume on failure) ──────────
export NXF_OPTS="-Xms512m -Xmx2g"
export NXF_SINGULARITY_CACHEDIR="/home/gao/.singularity/nf-core"

run_nextflow() {
    nextflow run nf-core/rnaseq \
        -r 3.15.1 \
        -profile singularity \
        -c local_optimized.config \
        "$@" \
        --input nf_core_samplesheet.csv \
        --outdir ../output_results \
        --fasta /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa \
        --gtf /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/gencode.vM35.annotation.gtf \
        --star_index '/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/star_index' \
        --gencode \
        --aligner star_salmon \
        --max_cpus 28 \
        --max_memory '108.GB'
}

echo "[$(ts)] ── Step 2 / 4 : nf-core/rnaseq (first run) ───"
if run_nextflow; then
    echo "[$(ts)] ✅ Step 2 done (first run succeeded)"
else
    echo ""
    echo "[$(ts)] ⚠️  First run failed — auto-resuming ..."
    echo "[$(ts)] ── Step 2 / 4 : nf-core/rnaseq (resume) ────"
    run_nextflow -resume || fail "nf-core/rnaseq failed after auto-resume"
    echo "[$(ts)] ✅ Step 2 done (succeeded after resume)"
fi
echo ""

# ── Step 3: DE + PCA ─────────────────────────────────────────
echo "[$(ts)] ── Step 3 / 4 : DESeq2 + PCA ─────────────────"
conda run -n DE_R45 Rscript 4_run_DE_PCA.R \
    || fail "DE analysis failed"
echo "[$(ts)] ✅ Step 3 done"
echo ""

# ── Step 4: Enrichment ───────────────────────────────────────
echo "[$(ts)] ── Step 4 / 4 : GO / KEGG / GSEA / StemCell ──"
conda run -n DE_R45 Rscript 5_run_enrichment.R \
    || fail "Enrichment analysis failed"
echo "[$(ts)] ✅ Step 4 done"
echo ""

echo "[$(ts)] $SEP"
echo "[$(ts)]  ALL DONE 🎉"
echo "[$(ts)]  Results: ${PROJECT_DIR}/Data_Analysis/"
echo "[$(ts)] $SEP"
