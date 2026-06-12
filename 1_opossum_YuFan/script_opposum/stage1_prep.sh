#!/bin/bash
# 进入工作目录
cd /Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/

echo "Step 1: Building HISAT2 index..."
# 使用绝对路径运行 hisat2-build
hisat2-build -p 16 ./dv-2k.fasta dv_index

echo "Step 2: Running Liftoff..."
REF_FA="/Work_bio/references/Monodelphis_domestica/MonDom5/ncbi_refseq/mDidVir1.genome.fa"
REF_GTF="/Work_bio/references/Monodelphis_domestica/MonDom5/ncbi_refseq/mDidVir1.annotation.gtf"
TARGET_FA="/Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/dv-2k.fasta"

# 运行 Liftoff
liftoff -g $REF_GTF -o Didelphis_v.liftoff.gtf -u unmapped.txt -p 16 -sc 0.85 $TARGET_FA $REF_FA