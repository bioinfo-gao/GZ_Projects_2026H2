# conda config --add channels conda-forge
# conda config --add channels bioconda # Nextflow is available via the Bioconda channel.
# mamba install nextflow 
nextflow -version

# tmux ls
# tmux kill-server                     # Close All Sessions
# tmux kill-session -t session_name    # Close a specific session by name

tmux new -s RNA 
#tmux a

# Check if rRNA databases exist locally
# RIBO_MANIFEST="/Work_bio/gao/projects/2026_Item12_rRNA/scripts/rRNA_databases/sortmerna_database_manifest.txt"
# RIBO_MANIFEST="/Work_bio/gao/projects/rRNA_databases/sortmerna_database_manifest.txt"
RIBO_MANIFEST="/Work_bio/references/rRNA_databases/sortmerna_database_manifest.txt"

if [ ! -f "$RIBO_MANIFEST" ]; then
    echo "ERROR: rRNA database manifest not found at $RIBO_MANIFEST!"
    echo "Please run download_rRNA_databases.sh on a machine with internet access,"
    echo "then transfer the rRNA_databases folder to this directory."
    exit 1
fi
echo "All rRNA database files verified. Starting pipeline..."


# cd /Work_bio/gao/projects/2026_Item12_rRNA/scripts
cd /home/gao/projects_2026H2/5_mouse-rRNA/scripts

# 限制 Nextflow 自身的内存开销，确保它不被 Killed
export NXF_OPTS="-Xms512m -Xmx2g"
# 设置 Singularity 缓存目录以避免重复下载容器镜像
export NXF_SINGULARITY_CACHEDIR="/home/gao/.singularity/nf-core"


# CC: 已切换为小鼠基因组(GRCm39/gencode M35)正式运行，下面这段曾经的预估讨论已落地到 local_optimized.config：
# CC: 机器实测 32 物理核心(64线程超线程)/125GB 内存 -> executor.cpus=28(按物理核心算，超线程对STAR这种计算密集型任务收益有限)，executor.memory=110GB
# CC: STAR_ALIGN 内存由人类的 32-35GB 下调为 30GB(小鼠基因组比人类小约12%，且按用户要求多留缓冲避免swap)
# CC: STAR_ALIGN 是内存瓶颈而非CPU瓶颈，故按内存反推并行数：cpus=9, maxForks=3 (3*9=27 CPUs, 3*30=90GB)，比原来cpus=14,maxForks=2多跑50%样本
# CC: 如果你的索引是用非默认参数建的（比如调大了 --sjdbOverhang 或加了大量 GTF 注释/SNP 信息），实际内存会更高，建议先跑一次用 /usr/bin/time -v 核实

nextflow run nf-core/rnaseq \
    -r 3.15.1 \
    -profile singularity \
    -c local_optimized.config \
    -c avoid_download.config \
    --input /home/gao/projects_2026H2/5_mouse-rRNA/scripts/nf_core_samplesheet.csv \
    --outdir ../output_results \
    --fasta /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa \
    --gtf /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/gencode.vM35.annotation.gtf \
    --star_index /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/star_index \
    --gencode \
    --aligner star_salmon \
    --remove_ribo_rna \
    --ribo_database_manifest "$RIBO_MANIFEST" \
    --save_non_ribo_reads \
    --max_cpus 28 \
    --max_memory '110.GB'

# resume

#   不要使用 -resume 参数首次运行：我已经在脚本中移除了它，确保首次运行使用全新配置,  否则是有可能引用之前的错误配置和缓存， 或者旧配置或缓存 ，第一次意外中断自恨才需要 -resume 

#　还需要去除上一行的末尾的　\
# 缓存目录：Singularity 容器会被缓存到 /home/gao/.singularity/nf-core/，避免重复下载
# 工作目录：Nextflow 的临时工作文件会在 /home/gao/projects_2026H2/5_mouse-rRNA/scripts/work/ 目录中生成
# 输出目录：最终结果将保存在 /home/gao/projects_2026H2/5_mouse-rRNA/output_results/



