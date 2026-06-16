#!/bin/bash

INDIR="/home/gao/Dropbox/JinPeng"
OUTDIR="/home/gao/projects_2026H2/2_jingpeng_wgs_quality/qc_out"
mkdir -p "${OUTDIR}"

for r1 in $(find "${INDIR}" -maxdepth 1 -name "*_R1.fq.gz" -o -name "*_R1.fastq.gz" | sort); do
    fname=$(basename "$r1")
    r2_fname=$(echo "$fname" | sed 's/_R1\./_R2\./')
    r2="${INDIR}/${r2_fname}"

    if [ ! -f "$r2" ]; then
        echo "WARNING: R2 not found for $fname, skipping"
        continue
    fi

    sample=$(echo "$fname" | sed 's/_R1\.f.*//')
    echo "=== Processing: $sample ==="

    fastp -i "$r1" -I "$r2" \
      -j "${OUTDIR}/${sample}_fastp.json" \
      -h "${OUTDIR}/${sample}_fastp.html" \
      -w 4 \
      --reads_to_process 1000000

    echo "  Done: $sample"
    echo ""
done

echo "=============================="
echo "MultiQC 汇总"
echo "=============================="
multiqc "${OUTDIR}" -o "${OUTDIR}/multiqc" --force

echo ""
echo "全部完成！"
echo "  单样品报告: ${OUTDIR}/*_fastp.html"
echo "  汇总报告:   ${OUTDIR}/multiqc/multiqc_report.html"
