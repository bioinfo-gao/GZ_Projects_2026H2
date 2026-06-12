#!/bin/bash

# 1. 进入您的脚本工作目录
cd /home/gao/projects_2026H1/2026_Item16_ZhenYan/scripts

# 2. 清理可能存在的旧会话并创建新的 rnaseq 会话
tmux kill-session -t rnaseq 2>/dev/null || true

# 1. 创建并进入 rnaseq 会话
tmux new -s rnaseq                    # 创建名为 'rnaseq' 的新会话并立即进入（前台运行，您会直接看到终端切换到该会话）
# tmux new-session -d -s rnaseq         # 创建名为 'rnaseq' 的新会话但在后台运行（-d 参数表示 detached，不会自动进入会话）

# 注意：这两行不能同时存在，因为它们会创建两个同名会话导致冲突
# 如果您想创建并立即进入会话，请使用第一行（取消注释第一行，删除或注释第二行）
# 如果您想在后台创建会话稍后再连接，请使用第二行（取消注释第二行，删除或注释第一行）

# . 如果您在普通终端中
# 您看到的是正常的 shell 提示符（如 gao@us1:~$）
# 这表示您在主终端，而不是在 Tmux 会话中
# 2. 如果您在 Tmux 会话中
# 终端底部通常会显示绿色/蓝色的状态栏
# 或者提示符无颜色

# 3. 启动 nf-core 3.15.1 在 tmux 会话中

# 限制 Nextflow 自身的内存开销，确保它不被 Killed
export NXF_OPTS="-Xms512m -Xmx2g"

# WARN: Singularity cache directory has not been defined -- Remote image will be stored in the path: /Work_bio/gao/projects_2026H1/2026_Item16_ZhenYan/scripts/work/singularity
# -- Use the environment variable NXF_SINGULARITY_CACHEDIR to specify a different location

# 这样做的好处：

# 避免重复下载：nf-core 容器镜像只会下载一次并缓存在指定位置
# 节省磁盘空间：不会在每次运行时都在 work/singularity 目录中存储重复的镜像
# 提高性能：后续运行可以直接使用缓存的镜像，加快启动速度
# 消除警告信息：不会再看到关于 Singularity 缓存目录未定义的警告
# 修改后的脚本现在会将 Singularity 容器镜像缓存在 /home/gao/.singularity/nf-core 目录中，这是一个标准且合理的缓存位置。

# 设置 Singularity 缓存目录以避免重复下载容器镜像
export NXF_SINGULARITY_CACHEDIR="/home/gao/.singularity/nf-core"


# 在新的 tmux 会话中运行 Nextflow
#  删除 --remove_ribo_rna 因为是polyA 测序，所以不需要这行

# WARN: Singularity cache directory has not been defined -- Remote image will be stored in the path: /Work_bio/gao/projects_2026H1/2026_Item16_ZhenYan/scripts/work/singularity -- Use the environment variable NXF_SINGULARITY_CACHEDIR to specify a different location
nextflow run nf-core/rnaseq \
    -r 3.15.1 \
    -profile singularity \
    -c local_optimized.config \
    --input nf_core_samplesheet.csv \
    --outdir ../output_results \
    --fasta /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/GRCh38.primary_assembly.genome.fa \
    --gtf /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/gencode.v45.annotation.gtf \
    --star_index '/Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/star_index' \
    --gencode \
    --aligner star_salmon \
    --max_cpus 28 \
    --max_memory '108.GB' 
    # \
    # -resume
# 首次运行，不resume，避免引入无关历史记录
echo "Nextflow 已在 tmux 会话 'rnaseq' 中启动"
echo "使用 'tmux a' 连接查看输出"

# Ctrl + B , then D 退出 tmux 会话