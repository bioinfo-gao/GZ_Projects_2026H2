#!/usr/bin/env bash
# Phase 1 — assembly-free taxonomic profiling (nf-core/taxprofiler 2.0.1)
# Kraken2 + Bracken (属/种级丰度) + MetaPhlAn (marker-gene 种级) 双工具交叉验证。
# self-relaunch into tmux + auto -resume（同 rnaseq/wgs 模式）。
set -uo pipefail
PROJ=/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics
SCRIPT="$PROJ/scripts/3_run_taxprofiler.sh"
HOST_FA=/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa
HOST_IDX=/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/bowtie2_index

if [ -z "${TMUX:-}" ]; then
    cd "$PROJ"
    tmux new-session -d -s tax17 "bash '$SCRIPT' 2>&1 | tee $PROJ/logs/taxprofiler_run.log"
    echo "launched tmux tax17"; exit 0
fi

cd "$PROJ"
export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity'
export NXF_ANSI_LOG=false

# 若预建索引就绪则用它，否则退回让 pipeline 从 fasta 现建
IDX_ARG=(--hostremoval_reference "$HOST_FA")
if ls "$HOST_IDX"/GRCm39.*.bt2l >/dev/null 2>&1 || ls "$HOST_IDX"/GRCm39.*.bt2 >/dev/null 2>&1; then
    IDX_ARG+=(--shortread_hostremoval_index "$HOST_IDX")
fi

run() {
    nextflow run nf-core/taxprofiler -r 2.0.1 -profile singularity \
        -c scripts/local_resources_taxprofiler.config "$@" \
        --input samplesheet_taxprofiler.csv --databases databases.csv \
        --outdir output_results -work-dir work_taxprofiler \
        --perform_shortread_qc \
        --perform_shortread_hostremoval "${IDX_ARG[@]}" \
        --run_kraken2 --run_bracken --run_metaphlan \
        --run_profile_standardisation --run_krona
    # 注：taxpasta_add_name/rank/lineage 需 --taxpasta_taxonomy_dir(NCBI taxdump)，本地未备；
    #     下游 R 用 MetaPhlAn lineage + Bracken/Kraken2 per-sample report 的 name 列解析物种名。
}
# 首跑即带 -resume：无缓存时等价于全新跑，有缓存时复用已完成 process
if run -resume; then echo "taxprofiler OK"; exit 0; fi
sleep 10
if run -resume; then echo "taxprofiler OK (retry)"; exit 0; fi
echo "taxprofiler FAILED after resume"; exit 1
