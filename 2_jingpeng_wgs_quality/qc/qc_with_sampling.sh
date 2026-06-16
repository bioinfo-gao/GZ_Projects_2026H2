#!/bin/bash

INDIR="/home/gao/Dropbox/JinPeng"
OUTDIR="/home/gao/projects_2026H2/2_jingpeng_wgs_quality/qc_out"
SAMPDIR="${OUTDIR}/sampled"
THREADS=8
NREADS=1000000
LINES=$((NREADS * 4))

mkdir -p "${OUTDIR}" "${SAMPDIR}"

echo "=============================="
echo "Step 1: 抽样 ${NREADS} 条 reads / 样品"
echo "=============================="
for fq in $(find "${INDIR}" -maxdepth 1 -name "*.fq.gz" -o -name "*.fastq.gz" | sort); do
    fname=$(basename "$fq")
    # 统一输出为 .fastq 后缀
    outname="${fname%.fq.gz}"
    outname="${outname%.fastq.gz}.fastq"
    out="${SAMPDIR}/${outname}"
    echo "  Sampling: $fname -> $outname"
    zcat "$fq" | head -${LINES} > "$out"
    echo "    Done: $(wc -l < "$out") lines"
done

echo ""
echo "=============================="
echo "Step 2: FastQC 逐个运行"
echo "=============================="
for fq in $(find "${SAMPDIR}" -name "*.fastq" | sort); do
    fname=$(basename "$fq")
    echo "  Running FastQC: $fname"
    fastqc "$fq" -o "${OUTDIR}" -t ${THREADS} --noextract
    echo "    Done: $fname"
done

echo ""
echo "=============================="
echo "Step 3: MultiQC 汇总"
echo "=============================="
multiqc "${OUTDIR}" -o "${OUTDIR}/multiqc" --force

echo ""
echo "=============================="
echo "Step 4: 长度分布统计"
echo "=============================="
for fq in $(find "${SAMPDIR}" -name "*.fastq" | sort); do
    fname=$(basename "$fq")
    echo "--- $fname ---"
    awk 'NR%4==2{print length($0)}' "$fq" | sort -n | uniq -c
done > "${OUTDIR}/length_distribution.txt"
echo "  已保存: ${OUTDIR}/length_distribution.txt"

echo ""
echo "=============================="
echo "全部完成！"
echo "  FastQC 报告: ${OUTDIR}/"
echo "  MultiQC 汇总: ${OUTDIR}/multiqc/multiqc_report.html"
echo "  长度分布: ${OUTDIR}/length_distribution.txt"
echo "=============================="
