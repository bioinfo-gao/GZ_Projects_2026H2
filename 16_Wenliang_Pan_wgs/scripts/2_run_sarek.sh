#!/usr/bin/env bash
# Project 16 — standard human germline WGS via nf-core/sarek (Mode A).
# Self-relaunch into tmux + auto-resume on failure. Run: bash scripts/2_run_sarek.sh
set -uo pipefail

PROJ=/home/gao/projects_2026H2/16_Wenliang_Pan_wgs
SCRIPT="$PROJ/scripts/2_run_sarek.sh"
LOG="$PROJ/logs/sarek_run.log"

# --- relaunch self inside tmux so the job outlives the ssh/agent session ---
if [ -z "${TMUX:-}" ]; then
    tmux new-session -d -s pan_wgs "bash '$SCRIPT' 2>&1 | tee '$LOG'"
    echo "launched tmux session 'pan_wgs' (log: $LOG)"
    exit 0
fi

cd "$PROJ"
export NXF_SINGULARITY_CACHEDIR="$PROJ/singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

run_sarek() {
    conda run -n regular_bioinfo nextflow run nf-core/sarek -r 3.8.1 \
        -profile singularity \
        -c scripts/local_resources.config \
        --input scripts/samplesheet.csv \
        --outdir output_results \
        -w work \
        --genome GATK.GRCh38 \
        --aligner bwa-mem2 \
        --trim_fastq \
        --skip_tools baserecalibrator \
        --tools haplotypecaller,manta,tiddit,cnvkit,vep \
        --download_cache \
        --max_memory 120.GB --max_cpus 56 \
        "$@"
}

echo "=== sarek start $(date) ==="
if run_sarek; then echo "=== sarek OK (fresh) $(date) ==="; exit 0; fi
echo "=== fresh run failed, retrying with -resume $(date) ==="
if run_sarek -resume; then echo "=== sarek OK (resume) $(date) ==="; exit 0; fi
echo "=== sarek FAILED even after resume $(date) ==="
exit 1
