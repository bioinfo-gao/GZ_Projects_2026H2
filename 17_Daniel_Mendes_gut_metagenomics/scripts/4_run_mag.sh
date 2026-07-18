#!/usr/bin/env bash
# Phase 2 — assembly-based MAG recovery (nf-core/mag 5.4.2)
# 分组共组装（--coassemble_group: AL, IF）→ MEGAHIT（--skip_spades，避 125GB 内存墙）
# → MetaBAT2 + MaxBin2 + SemiBin2 分箱 → DAS Tool 精炼 → BUSCO + CheckM2 质控
# → GTDB-Tk 分类（split-tree）→ Prokka 注释 → MultiQC。
# 依赖：GTDB-Tk r226 + CheckM2 库已由 scripts/1_predownload_dbs.sh 备好。
set -uo pipefail
PROJ=/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics
SCRIPT="$PROJ/scripts/4_run_mag.sh"
NEXTFLOW=/Work_bio/gao/configs/.conda/envs/mag_biobakery/bin/nextflow   # 必须绝对路径 (≥25.04.2)
export JAVA_HOME=/Work_bio/gao/configs/.conda/envs/mag_biobakery
HOST_FA=/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa
GTDB_DB=/Work_bio/references/Metagenomics/gtdbtk/release226
CHECKM2_DB=/Work_bio/references/Metagenomics/checkm2
BUSCO_DB=/Work_bio/references/Metagenomics/busco/bacteria_odb10.2024-01-08.tar.gz

if [ -z "${TMUX:-}" ]; then
    cd "$PROJ"
    tmux new-session -d -s mag17 "bash '$SCRIPT' 2>&1 | tee $PROJ/logs/mag_run.log"
    echo "launched tmux mag17"; exit 0
fi

cd "$PROJ"
export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity'
export NXF_ANSI_LOG=false

# GTDB 库未就位则先跳过分类（真实分类需库；库好后 -resume 补跑那步不现实，届时单独跑）
GTDB_ARG=(--skip_gtdbtk)
if [ -d "$GTDB_DB" ] && [ -n "$(ls -A "$GTDB_DB" 2>/dev/null)" ]; then
    GTDB_ARG=(--gtdb_db "$GTDB_DB")
fi
CHECKM2_ARG=()
if [ -d "$CHECKM2_DB" ] && ls "$CHECKM2_DB"/*.dmnd >/dev/null 2>&1; then
    CHECKM2_ARG=(--run_checkm2 --checkm2_db "$CHECKM2_DB")
fi

run() {
    "$NEXTFLOW" run nf-core/mag -r 5.4.2 -profile singularity \
        -c scripts/local_resources_mag.config "$@" \
        --input samplesheet_mag.csv --outdir output_results_mag -work-dir work_mag \
        --host_fasta "$HOST_FA" \
        --coassemble_group --skip_spades \
        --skip_concoct --skip_comebin --skip_metabinner \
        --refine_bins_dastool \
        --busco_db "$BUSCO_DB" --busco_db_lineage bacteria_odb10 \
        "${CHECKM2_ARG[@]}" "${GTDB_ARG[@]}"
}
if run;          then echo "mag OK"; exit 0; fi
if run -resume;  then echo "mag OK (resume)"; exit 0; fi
echo "mag FAILED after resume"; exit 1
