# 1. 激活你的环境
conda activate DE_R45

# 报错的核心原因是 Conda/Mamba 的包命名规则：
# r- 前缀仅用于 CRAN 官方包（如 ggplot2, dplyr）
# bioconductor- 前缀用于 Bioconductor 生态包（如 clusterProfiler, ggtree）
# 你提供的包全部属于 Bioconductor，因此 mamba 在默认频道找不到 r-clusterprofiler 等名称。

# 2. 使用 mamba 安装（严格锁定 R 4.5.*，自动拉取编译好的二进制包）
# 显式声明版本约束可防止 mamba 求解器为了兼容性自动升级 R。星号 * 允许安装 4.5.x 的任意小版本（如 4.5.1/4.5.2），确保环境稳定。
# 1. 确保频道优先级正确（关键！）
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict

# 2. 使用正确前缀安装（锁定 R 4.5.* 防止版本变动）
mamba install -c conda-forge -c bioconda \
  "r-base=4.5.*" \
  bioconductor-clusterprofiler \
  bioconductor-enrichplot \
  bioconductor-ggtree \
  bioconductor-org.at.tair.db \
  r-ggiraph \
  r-gdtools 


# 1. 确保已激活环境
conda activate DE_R45

# 2. 使用 mamba 安装 (注意前缀是 bioconductor-)
mamba install -c conda-forge -c bioconda bioconductor-rtracklayer

#   r-ggrepel \
#   r-pheatmap


# 避坑指南
# 频道优先级：如果安装时提示冲突，先运行：
# conda config --add channels conda-forge
# conda config --add channels bioconda
# conda config --set channel_priority strict
# 网络问题：如果下载卡在 Solving environment，可配置清华镜像：
# bash
# 12