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

# 所有 pipeline 参数(布尔+路径)走 params_mag.yaml —— 新版 Nextflow 把 CLI bare `--flag`
# 解析成字符串 "true"，被 nf-schema 严格校验拒绝；YAML 里布尔=真布尔、路径=字符串，绕过该坑。
# checkm2_db 必须指 .dmnd 文件(非目录)。GTDB/CheckM2 库已就位(见 params_mag.yaml)。
run() {
    "$NEXTFLOW" run nf-core/mag -r 5.4.2 -profile singularity \
        -c scripts/local_resources_mag.config \
        -params-file scripts/params_mag.yaml \
        -work-dir work_mag "$@"
}
if run;          then echo "mag OK"; exit 0; fi
if run -resume;  then echo "mag OK (resume)"; exit 0; fi
echo "mag FAILED after resume"; exit 1
