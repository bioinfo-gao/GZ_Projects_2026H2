#!/bin/bash
set -euo pipefail

# ============================================================
# CHO WGS 物种确认 + 株系鉴定
# ============================================================

# === Configuration ===
PROJDIR="/home/gao/projects_2026H2/3_cho_wgs_species_confirm"
R1="/home/gao/Dropbox/Keqiang/23TF7FLT4_3_0469165296_wt1_S8_L003_R1_001.fastq.gz"
R2="/home/gao/Dropbox/Keqiang/23TF7FLT4_3_0469165296_wt1_S8_L003_R2_001.fastq.gz"
CHO_REF="/Work_bio/references/Cricetulus_griseus/CriGri-PICR/ncbi_refseq/GCF_003668045.1_CriGri-PICR_genomic.fna"
CHO_GFF="/Work_bio/references/Cricetulus_griseus/CriGri-PICR/ncbi_refseq/GCF_003668045.1_CriGri-PICR_genomic.gff.gz"
THREADS=8
SUBSAMPLE=5000000  # 5M read pairs for quick analysis

mkdir -p ${PROJDIR}/{qc,align,results}
cd ${PROJDIR}

echo "============================================================"
echo "  CHO WGS Species Confirmation & Strain Identification"
echo "  Sample: wt1"
echo "  Started: $(date)"
echo "============================================================"
echo ""

# ============================================================
# Step 1: Quick QC with fastp
# ============================================================
echo ">>> Step 1: fastp QC (sampling ${SUBSAMPLE} reads)"
echo "------------------------------------------------------------"
fastp -i ${R1} -I ${R2} \
  -j ${PROJDIR}/qc/wt1_fastp.json \
  -h ${PROJDIR}/qc/wt1_fastp.html \
  -w 4 --reads_to_process ${SUBSAMPLE}

echo ""
echo ">>> QC Summary:"
python3 -c "
import json
d = json.load(open('${PROJDIR}/qc/wt1_fastp.json'))
s = d['summary']
bf = s['before_filtering']
af = s['after_filtering']
print(f\"  Total reads (sampled): {bf['total_reads']:,}\")
print(f\"  Q20 rate: {bf['q20_rate']}\")
print(f\"  Q30 rate: {bf['q30_rate']}\")
print(f\"  GC content: {bf['gc_content']}\")
print(f\"  Duplication rate: {d.get('duplication',{}).get('rate','N/A')}\")
print(f\"  Insert size peak: {d.get('insert_size',{}).get('peak','N/A')}\")
print(f\"  After filter Q30: {af['q30_rate']}\")
"
echo ""

# ============================================================
# Step 2: Extract subset reads
# ============================================================
echo ">>> Step 2: Extracting ${SUBSAMPLE} read pairs for alignment"
echo "------------------------------------------------------------"
fastp -i ${R1} -I ${R2} \
  --reads_to_process ${SUBSAMPLE} \
  -o ${PROJDIR}/align/sub_R1.fq.gz \
  -O ${PROJDIR}/align/sub_R2.fq.gz \
  -A -G -Q -L \
  -j ${PROJDIR}/align/sub_fastp.json \
  -h /dev/null 2>/dev/null
echo "  Done: sub_R1.fq.gz and sub_R2.fq.gz"
echo ""

# ============================================================
# Step 3: Align to CHO genome
# ============================================================
echo ">>> Step 3: Aligning to CHO genome (CriGri-PICR)"
echo "------------------------------------------------------------"
bwa mem -t ${THREADS} -R "@RG\tID:wt1\tSM:wt1\tPL:ILLUMINA" \
  ${CHO_REF} \
  ${PROJDIR}/align/sub_R1.fq.gz \
  ${PROJDIR}/align/sub_R2.fq.gz 2>${PROJDIR}/align/bwa_cho.log | \
  samtools sort -@ 4 -o ${PROJDIR}/align/wt1_cho.bam

samtools index ${PROJDIR}/align/wt1_cho.bam
echo ""

echo ">>> CHO Alignment Statistics:"
echo "------------------------------------------------------------"
samtools flagstat ${PROJDIR}/align/wt1_cho.bam | tee ${PROJDIR}/results/cho_flagstat.txt
echo ""

# Extract key mapping rate
MAPPING_RATE=$(samtools flagstat ${PROJDIR}/align/wt1_cho.bam | grep "mapped (" | head -1 | grep -oP '\(\K[0-9.]+')
echo "============================================================"
echo "  >>> CHO MAPPING RATE: ${MAPPING_RATE}%"
echo "============================================================"
echo ""

if (( $(echo "$MAPPING_RATE > 85" | bc -l) )); then
  echo "  RESULT: Mapping rate > 85% => CONFIRMED CHO genome"
elif (( $(echo "$MAPPING_RATE > 50" | bc -l) )); then
  echo "  RESULT: Mapping rate 50-85% => POSSIBLY CHO but with issues"
else
  echo "  RESULT: Mapping rate < 50% => NOT CHO or severely contaminated"
fi
echo ""

# ============================================================
# Step 4: Find DHFR gene and check coverage for strain ID
# ============================================================
echo ">>> Step 4: DHFR gene analysis (strain identification)"
echo "------------------------------------------------------------"

# Find DHFR coordinates from GFF
echo "  Searching DHFR in annotation..."
zgrep -i "dhfr" ${CHO_GFF} | grep -P "\tgene\t" | head -5 > ${PROJDIR}/results/dhfr_gff_entries.txt
cat ${PROJDIR}/results/dhfr_gff_entries.txt

# Extract DHFR region coordinates
DHFR_CHR=$(awk '{print $1}' ${PROJDIR}/results/dhfr_gff_entries.txt | head -1)
DHFR_START=$(awk '{print $4}' ${PROJDIR}/results/dhfr_gff_entries.txt | head -1)
DHFR_END=$(awk '{print $5}' ${PROJDIR}/results/dhfr_gff_entries.txt | head -1)

if [ -n "${DHFR_CHR}" ] && [ -n "${DHFR_START}" ] && [ -n "${DHFR_END}" ]; then
  DHFR_REGION="${DHFR_CHR}:${DHFR_START}-${DHFR_END}"
  echo ""
  echo "  DHFR region: ${DHFR_REGION}"
  echo ""

  # DHFR coverage
  echo "  DHFR locus coverage:"
  DHFR_COV=$(samtools depth -r ${DHFR_REGION} ${PROJDIR}/align/wt1_cho.bam | \
    awk '{sum+=$3; n++} END {if(n>0) printf "%.2f", sum/n; else print "0"}')
  DHFR_BASES=$(samtools depth -r ${DHFR_REGION} ${PROJDIR}/align/wt1_cho.bam | wc -l)
  DHFR_LEN=$((DHFR_END - DHFR_START))
  echo "    Mean depth at DHFR: ${DHFR_COV}x"
  echo "    Bases with coverage: ${DHFR_BASES} / ${DHFR_LEN}"
  echo ""

  # Flanking region coverage (as control) - 500kb upstream and downstream
  FLANK_START=$((DHFR_START - 500000))
  FLANK_END=$((DHFR_END + 500000))
  if [ ${FLANK_START} -lt 1 ]; then FLANK_START=1; fi
  FLANK_REGION="${DHFR_CHR}:${FLANK_START}-${FLANK_END}"

  FLANK_COV=$(samtools depth -r ${FLANK_REGION} ${PROJDIR}/align/wt1_cho.bam | \
    awk '{sum+=$3; n++} END {if(n>0) printf "%.2f", sum/n; else print "0"}')
  echo "    Flanking region (±500kb) mean depth: ${FLANK_COV}x"
  echo ""

  # Interpretation
  echo "  >>> DHFR Strain Interpretation:"
  echo "  ============================================================"
  python3 -c "
dhfr_cov = float('${DHFR_COV}')
flank_cov = float('${FLANK_COV}')
dhfr_bases = int('${DHFR_BASES}')
dhfr_len = int('${DHFR_LEN}')
coverage_pct = dhfr_bases / dhfr_len * 100 if dhfr_len > 0 else 0

print(f'    DHFR depth: {dhfr_cov:.2f}x')
print(f'    Flanking depth: {flank_cov:.2f}x')
print(f'    DHFR bases covered: {coverage_pct:.1f}%')
print()

if dhfr_cov < 0.01 and coverage_pct < 5:
    print('    VERDICT: DHFR DELETED (homozygous)')
    print('    => Likely CHO-DG44')
elif flank_cov > 0 and dhfr_cov / max(flank_cov, 0.001) < 0.3:
    print('    VERDICT: DHFR PARTIALLY DELETED (hemizygous)')
    print('    => Likely CHO-DXB11')
elif flank_cov > 0 and dhfr_cov / max(flank_cov, 0.001) > 0.7:
    print('    VERDICT: DHFR INTACT')
    print('    => Likely CHO-K1 or CHO-S')
else:
    print('    VERDICT: DHFR status inconclusive')
    print('    => Need full-depth alignment for confirmation')
"
  echo "  ============================================================"
else
  echo "  WARNING: Could not find DHFR in GFF annotation"
  echo "  Trying alternative search..."
  zgrep -i "dihydrofolate" ${CHO_GFF} | grep -P "\tgene\t" | head -5
fi

echo ""

# ============================================================
# Step 5: Genome-wide coverage summary
# ============================================================
echo ">>> Step 5: Genome-wide coverage estimate"
echo "------------------------------------------------------------"
python3 -c "
import json
d = json.load(open('${PROJDIR}/qc/wt1_fastp.json'))
total_bases = d['summary']['before_filtering']['total_bases']
# CHO genome ~2.4Gb
genome_size = 2.4e9
# Extrapolate from sampled reads to full dataset
sampled_reads = d['summary']['before_filtering']['total_reads']
# Estimate total reads from file size ratio
import os
r1_size = os.path.getsize('${R1}')
# rough estimate: 17GB compressed ~ 70GB uncompressed ~ 230M reads
est_total_reads = 230000000  # rough estimate for 17GB R1
est_coverage = est_total_reads * 2 * 150 / genome_size
print(f'  Sampled reads (QC): {sampled_reads:,}')
print(f'  Estimated total reads: ~{est_total_reads:,} (from file size)')
print(f'  Estimated genome coverage: ~{est_coverage:.0f}x')
print(f'  CHO genome size: ~2.4 Gb')
"
echo ""

# ============================================================
# Summary
# ============================================================
echo "============================================================"
echo "  ANALYSIS COMPLETE"
echo "  $(date)"
echo "============================================================"
echo ""
echo "  Output files:"
echo "    QC report:   ${PROJDIR}/qc/wt1_fastp.html"
echo "    BAM file:    ${PROJDIR}/align/wt1_cho.bam"
echo "    Flagstat:    ${PROJDIR}/results/cho_flagstat.txt"
echo "    DHFR info:   ${PROJDIR}/results/dhfr_gff_entries.txt"
echo "============================================================"

multiqc /home/gao/projects_2026H2/3_cho_wgs_species_confirm/qc/ \
  -o /home/gao/projects_2026H2/3_cho_wgs_species_confirm/qc/multiqc \
  --force