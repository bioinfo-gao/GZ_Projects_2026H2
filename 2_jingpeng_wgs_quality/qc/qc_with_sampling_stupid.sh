#!/bin/bash

INDIR="/home/gao/Dropbox/JinPeng"
OUTDIR="/home/gao/projects_2026H2/2_jingpeng_wgs_quality/qc_out"
SAMPDIR="${OUTDIR}/sampled"
THREADS=8
NREADS=1000000

mkdir -p "${OUTDIR}" "${SAMPDIR}"

echo "=============================="
echo "Step 1: seqtk 抽样 ${NREADS} 条 reads / 样品"
echo "=============================="
for fq in $(find "${INDIR}" -maxdepth 1 \( -name "*.fq.gz" -o -name "*.fastq.gz" \) | sort); do
    fname=$(basename "$fq")
    outname="${fname%.fq.gz}"
    outname="${outname%.fastq.gz}.sampled.fastq"
    out="${SAMPDIR}/${outname}"
    echo "  Sampling: $fname"
    seqtk sample -s42 "$fq" ${NREADS} > "$out"
    reads_out=$(($(wc -l < "$out") / 4))
    echo "    Done: ${reads_out} reads"
done

echo ""
echo "=============================="
echo "Step 2: FastQC"
echo "=============================="
for fq in $(find "${SAMPDIR}" -name "*.sampled.fastq" | sort); do
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
for fq in $(find "${SAMPDIR}" -name "*.sampled.fastq" | sort); do
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

