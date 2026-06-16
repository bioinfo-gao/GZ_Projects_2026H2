#!/bin/bash
set -e

DIR="/home/gao/Dropbox/JinPeng"
OUTDIR="${DIR}/fastqc_results"
SAMPDIR="${DIR}/sampled"
THREADS=8
NREADS=1000000  # 抽样100万条reads = 400万行

mkdir -p "${OUTDIR}" "${SAMPDIR}"

echo "=============================="
echo "Step 1: 每个文件抽样 ${NREADS} 条 reads"
echo "=============================="
LINES=$((NREADS * 4))
for fq in ${DIR}/*.fq.gz ${DIR}/*.fastq.gz; do
    [ -f "$fq" ] || continue
    fname=$(basename "$fq")
    out="${SAMPDIR}/${fname%.gz}"
    echo "  Sampling: $fname"
    zcat "$fq" | head -${LINES} > "$out"
done

echo ""
echo "=============================="
echo "Step 2: FastQC on sampled reads"
echo "=============================="
fastqc ${SAMPDIR}/*.fq ${SAMPDIR}/*.fastq \
    -o "${OUTDIR}" \
    -t ${THREADS} \
    --noextract 2>/dev/null

echo ""
echo "=============================="
echo "Step 3: MultiQC 汇总"
echo "=============================="
multiqc "${OUTDIR}" -o "${OUTDIR}/multiqc" --force

echo ""
echo "=============================="
echo "Step 4: 快速长度分布统计"
echo "=============================="
for fq in ${SAMPDIR}/*.fq ${SAMPDIR}/*.fastq; do
    [ -f "$fq" ] || continue
    fname=$(basename "$fq")
    echo "--- $fname ---"
    awk 'NR%4==2{print length($0)}' "$fq" | sort -n | uniq -c | sort -rn -k2 | head -10
done 2>&1 | tee "${OUTDIR}/length_distribution.txt"

echo ""
echo "=============================="
echo "完成！预计总耗时 10-15 分钟"
echo "  汇总报告: ${OUTDIR}/multiqc/multiqc_report.html"
echo "  长度分布: ${OUTDIR}/length_distribution.txt"
echo "=============================="
