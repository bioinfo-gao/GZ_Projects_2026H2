#!/bin/bash
# 目的：用 nf-core/sarek 官方内置的 test profile（自带极小测试 fastq + 参考序列）
# 验证 nf-core/sarek 在本机的安装、容器拉取、执行器配置是否全部正常工作。
# 不需要真实人类 WGS 数据或本地参考基因组，全部由 pipeline 自动下载。

set -euo pipefail

PROJECT_DIR="/home/gao/projects_2026H2/4_wgs_human_immu"
cd "$PROJECT_DIR/scripts"

tmux kill-session -t sarek_test 2>/dev/null || true

# 注意：export 必须写在 tmux 命令字符串内部执行，而不是在外层脚本里 export ——
# 如果 tmux server 已经在后台常驻运行，新建的 session 会沿用 server 启动时的环境，
# 不会继承外层脚本临时 export 的变量（这是调试时踩过的一个坑，NXF_SINGULARITY_CACHEDIR
# 没生效，容器被缓存进了 work 目录）。

tmux new -d -s sarek_test "
  export NXF_OPTS='-Xms512m -Xmx2g';
  export NXF_SINGULARITY_CACHEDIR='/Work_bio/gao/configs/.singularity';
  /home/gao/.conda/envs/regular_bioinfo/bin/nextflow run nf-core/sarek \
    -r 3.8.1 \
    -profile test,singularity \
    -c $PROJECT_DIR/configs/local_resources.config \
    --outdir $PROJECT_DIR/test_run/outdir \
    -work-dir $PROJECT_DIR/test_run/work \
    -resume \
    2>&1 | tee $PROJECT_DIR/logs/sarek_test_run.log
"

echo "已在 tmux 会话 'sarek_test' 中后台启动 nf-core/sarek 测试运行。"
echo "查看实时日志：tmux attach -t sarek_test   (或 tail -f $PROJECT_DIR/logs/sarek_test_run.log)"
