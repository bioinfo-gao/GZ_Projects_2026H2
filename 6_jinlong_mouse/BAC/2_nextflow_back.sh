#!/bin/bash
# 运行步骤:
#   1. 先运行 python 1_produce_nf-core_Samplesheet.py 生成 nf_core_samplesheet.csv
#   2. 在普通终端中运行本脚本: bash 2_nextflow.sh
#      (脚本会自动在 tmux 'rnaseq' 会话中启动 Nextflow)

cd /home/gao/projects_2026H2/6_jinlong_mouse/scripts

tmux kill-session -t rnaseq 2>/dev/null || true

export NXF_OPTS="-Xms512m -Xmx2g"
export NXF_SINGULARITY_CACHEDIR="/home/gao/.singularity/nf-core"

tmux new-session -d -s rnaseq "
  export NXF_OPTS='-Xms512m -Xmx2g'
  export NXF_SINGULARITY_CACHEDIR='/home/gao/.singularity/nf-core'
  cd /home/gao/projects_2026H2/6_jinlong_mouse/scripts

  nextflow run nf-core/rnaseq \
      -r 3.15.1 \
      -profile singularity \
      -c local_optimized.config \
      --input nf_core_samplesheet.csv \
      --outdir ../output_results \
      --fasta /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa \
      --gtf /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/gencode.vM35.annotation.gtf \
      --star_index '/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/star_index' \
      --gencode \
      --aligner star_salmon \
      --max_cpus 28 \
      --max_memory '108.GB' \
      2>&1 | tee nextflow_run.log
"

echo "Nextflow started in tmux session 'rnaseq'"
echo "Monitor: tmux attach -t rnaseq"
echo "Or tail:  tail -f /home/gao/projects_2026H2/6_jinlong_mouse/scripts/nextflow_run.log"
