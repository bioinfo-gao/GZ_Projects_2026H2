#!/bin/bash
# ============================================================================
# Study A / Step 4 — 体细胞变异过滤 + 计数 + Trp53 LOH 初查
#   sarek 只出了 raw Mutect2 VCF(FILTER=".")→ 需 FilterMutectCalls。
#   5 对 tumor-vs-origin：过滤→PASS SNV/indel 计数；Trp53 位点 tumor/origin 深度比(LOH线索)。
#   前置：A2 完成(output_A/variant_calling/mutect2/*)。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
GRCM39="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
OUT="$PROJ/analysis_A/somatic"; mkdir -p "$OUT"
RUN(){ conda run -n regular_bioinfo "$@"; }
M2="$PROJ/output_A/variant_calling/mutect2"
TRP53="chr11:69471185-69482699"
PAIRS=(RO_B1TP RO_B2TP RO_tumor1 RO_tumor2 RO_tumor3)

echo -e "pair\tPASS_total\tPASS_SNV\tPASS_indel" > "$OUT/somatic_counts.tsv"
for t in "${PAIRS[@]}"; do
  d="$M2/${t}_vs_RO_origin"; raw="$d/${t}_vs_RO_origin.mutect2.vcf.gz"
  [ -f "$raw" ] || { echo "缺 $raw"; continue; }
  filt="$OUT/${t}.mutect2.filtered.vcf.gz"
  ob=""; [ -f "$d/${t}_vs_RO_origin.mutect2.artifactprior.tar.gz" ] && ob="--ob-priors $d/${t}_vs_RO_origin.mutect2.artifactprior.tar.gz"
  echo ">> FilterMutectCalls $t"
  RUN gatk FilterMutectCalls -R "$GRCM39" -V "$raw" --stats "$raw.stats" $ob -O "$filt" 2>"$OUT/${t}.filter.log" \
    || { echo "  FilterMutectCalls 失败(见 ${t}.filter.log)"; continue; }
  tot=$(RUN bcftools view -f PASS "$filt" 2>/dev/null | grep -vc '^#')
  snv=$(RUN bcftools view -f PASS -v snps "$filt" 2>/dev/null | grep -vc '^#')
  ind=$(RUN bcftools view -f PASS -v indels "$filt" 2>/dev/null | grep -vc '^#')
  echo -e "${t}\t${tot}\t${snv}\t${ind}" >> "$OUT/somatic_counts.tsv"
  echo "   $t: PASS=$tot SNV=$snv indel=$ind"
done
echo "== 体细胞计数 =="; column -t "$OUT/somatic_counts.tsv"

echo ">> Trp53 位点深度(tumor/edited vs origin，LOH/缺失线索)"
echo -e "sample\tTrp53_meandepth" > "$OUT/trp53_depth.tsv"
for s in RO_origin RO_B1TP RO_B2TP RO_tumor1 RO_tumor2 RO_tumor3; do
  cram="$PROJ/output_A/preprocessing/markduplicates/$s/$s.md.cram"
  dp=$(RUN samtools depth -r "$TRP53" "$cram" 2>/dev/null | awk '{s+=$3;n++} END{if(n)printf "%.1f",s/n; else print "NA"}')
  echo -e "${s}\t${dp}" >> "$OUT/trp53_depth.tsv"
done
column -t "$OUT/trp53_depth.tsv"
echo "DONE A4 → $OUT （Trp53 LOH 精判需 het-SNP VAF，见报告阶段）"
