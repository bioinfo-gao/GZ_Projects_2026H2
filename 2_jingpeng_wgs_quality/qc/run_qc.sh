#!/bin/bash
set -e

DIR="/home/gao/Dropbox/JinPeng"
OUTDIR="${DIR}/fastqc_results"
THREADS=8

mkdir -p "${OUTDIR}"

echo "=============================="
echo "Step 1: 基本统计 (reads数量 & 长度)"
echo "=============================="
for fq in ${DIR}/*.fq.gz ${DIR}/*.fastq.gz; do
    [ -f "$fq" ] || continue
    fname=$(basename "$fq")
    echo "--- $fname ---"
    # 统计前10000条reads的长度分布
    zcat "$fq" | head -40000 | awk 'NR%4==2{len=length($0); total++; sum+=len; if(len<min||min==0)min=len; if(len>max)max=len} END{printf "  Sampled reads: %d\n  Min len: %d\n  Max len: %d\n  Avg len: %.1f\n", total, min, max, sum/total}'
    # 统计总reads数（通过行数/4）
    echo -n "  Total reads: "
    zcat "$fq" | wc -l | awk '{printf "%d (%.1fM)\n", $1/4, $1/4/1000000}'
done 2>&1 | tee "${OUTDIR}/basic_stats.txt"

echo ""
echo "=============================="
echo "Step 2: FastQC (全部样品)"
echo "=============================="
fastqc ${DIR}/*.fq.gz ${DIR}/*.fastq.gz \
    -o "${OUTDIR}" \
    -t ${THREADS} \
    --noextract

echo ""
echo "=============================="
echo "Step 3: MultiQC 汇总报告"
echo "=============================="
multiqc "${OUTDIR}" -o "${OUTDIR}/multiqc" --force

echo ""
echo "=============================="
echo "完成！请查看:"
echo "  基本统计: ${OUTDIR}/basic_stats.txt"
echo "  汇总报告: ${OUTDIR}/multiqc/multiqc_report.html"
echo "=============================="
