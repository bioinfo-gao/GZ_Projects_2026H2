#!/bin/bash
# 安装 DE 分析和富集分析所需 R 包 (DE_R45 环境, Human)
# 运行方法: bash 4A_install_library.sh

conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict

mamba install -n DE_R45 -c conda-forge -c bioconda \
  "r-base=4.5.*" \
  bioconductor-deseq2 \
  bioconductor-ashr \
  r-ggplot2 \
  r-pheatmap \
  r-dplyr \
  r-readr \
  r-tidyr \
  r-readxl \
  r-ggrepel \
  bioconductor-clusterprofiler \
  bioconductor-enrichplot \
  bioconductor-org.hs.eg.db \
  bioconductor-fgsea \
  r-msigdbr \
  r-ggupset \
  r-patchwork

echo "All packages installed. Activate with: conda activate DE_R45"
echo "Test in R: library(clusterProfiler); library(org.Hs.eg.db); library(msigdbr)"
