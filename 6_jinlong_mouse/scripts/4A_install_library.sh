#!/bin/bash
# 安装 DE 分析和富集分析所需 R 包 (DE_R45 环境)
# 运行方法: bash 4A_install_library.sh

# 确保频道优先级正确
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict

# 安装 DE 分析基础包 + 鼠富集分析包
mamba install -c conda-forge -c bioconda \
  "r-base=4.5.*" \
  bioconductor-deseq2 \
  bioconductor-ashr \
  r-ggplot2 \
  r-pheatmap \
  r-dplyr \
  r-readr \
  r-tidyr \
  r-ggrepel \
  bioconductor-clusterprofiler \
  bioconductor-enrichplot \
  bioconductor-org.mm.eg.db \
  bioconductor-fgsea \
  r-msigdbr \
  r-ggupset \
  r-patchwork

# 安装 rtracklayer（用于 GTF 解析，可选）
mamba install -c conda-forge -c bioconda bioconductor-rtracklayer

echo "All packages installed. Activate with: conda activate DE_R45"
echo "Test in R: library(clusterProfiler); library(org.Mm.eg.db); library(msigdbr)"
