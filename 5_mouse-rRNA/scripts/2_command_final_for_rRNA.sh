# conda config --add channels conda-forge
# conda config --add channels bioconda # Nextflow is available via the Bioconda channel.
# mamba install nextflow 
nextflow -version

#tmux ls
tmux new -s RNA  
#tmux a

# Check if rRNA databases exist locally
#RIBO_MANIFEST="/Work_bio/gao/projects/2026_Item12_rRNA/scripts/rRNA_databases/sortmerna_database_manifest.txt"
RIBO_MANIFEST="/Work_bio/gao/projects/rRNA_databases/sortmerna_database_manifest.txt"
if [ ! -f "$RIBO_MANIFEST" ]; then
    echo "ERROR: rRNA database manifest not found at $RIBO_MANIFEST!"
    echo "Please run download_rRNA_databases.sh on a machine with internet access,"
    echo "then transfer the rRNA_databases folder to this directory."
    exit 1
fi
echo "All rRNA database files verified. Starting pipeline..."


cd /Work_bio/gao/projects/2026_Item12_rRNA/scripts

# 限制 Nextflow 自身的内存开销，确保它不被 Killed
export NXF_OPTS="-Xms512m -Xmx2g"
# 设置 Singularity 缓存目录以避免重复下载容器镜像
export NXF_SINGULARITY_CACHEDIR="/home/gao/.singularity/nf-core"


# 2026-05-05 添加  改回同时跑两个样本，MEM 90G,  在Item 11 里面的是常规RNAseq 没有rRNA ，但没有样本有两对fastq， 并行曾经降低到一个样本
# 迄今所有的rRNA 都是两个并行

nextflow run nf-core/rnaseq \
    -r 3.15.1 \
    -profile singularity \
    -c local_optimized.config \
    -c avoid_download.config \
    --input /home/gao/projects/2026_Item12_rRNA/scripts/nf_core_samplesheet.csv \
    --outdir ../output_results \
    --fasta /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/GRCh38.primary_assembly.genome.fa \
    --gtf /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/gencode.v45.annotation.gtf \
    --star_index /Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/star_index \
    --gencode \
    --aligner star_salmon \
    --remove_ribo_rna \
    --ribo_database_manifest "$RIBO_MANIFEST" \
    --save_non_ribo_reads \
    --max_cpus 28 \
    --max_memory '90.GB'



# resume
#   不要使用 -resume 参数首次运行：我已经在脚本中移除了它，确保首次运行使用全新配置,  否则是有可能引用之前的错误配置和缓存， 或者旧配置或缓存 ，第一次意外中断自恨才需要 -resume 
#　还需要去除上一行的末尾的　\
# 缓存目录：Singularity 容器会被缓存到 /home/gao/.singularity/nf-core/，避免重复下载
# 工作目录：Nextflow 的临时工作文件会在 /home/gao/projects/2026_Item12_rRNA/scripts/work/ 目录中生成
# 输出目录：最终结果将保存在 /home/gao/projects/2026_Item12_rRNA/output_results/



