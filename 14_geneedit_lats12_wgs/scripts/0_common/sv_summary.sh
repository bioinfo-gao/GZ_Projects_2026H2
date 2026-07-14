#!/bin/bash
# ============================================================================
# SV 汇总 — TIDDIT 结构变异（两组），PASS 过滤 + 按类型计数 + 大小分布
#   Study A: 5 tumor/edited + origin（somatic SV = tumor 有 origin 无，脚本仅计数，
#            精细 somatic 差异在报告阶段比对）；Study B: 6 样 L1L2 vs L1L2H / 随龄。
#   ⚠ 背景过滤(MGP/DGV)未接入（0d 未下）→ 此处为 TIDDIT PASS 原始计数，
#     背景过滤作为后续精化（报告中注明）。
#   前置：A2/B2 完成（output_*/variant_calling/tiddit/*）。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
OUT="$PROJ/analysis_common/sv"; mkdir -p "$OUT"
RUN(){ conda run -n regular_bioinfo "$@"; }

echo -e "study\tsample\ttotal\tPASS\tPASS_DEL\tPASS_DUP\tPASS_INV\tPASS_BND\tPASS_other" > "$OUT/sv_counts.tsv"
summarize(){
  local study="$1" sample="$2" vcf="$3"
  [ -f "$vcf" ] || { echo "缺 $vcf"; return; }
  conda run -n regular_bioinfo bash -c "
    tot=\$(bcftools view '$vcf' 2>/dev/null | grep -vc '^#')
    pass=\$(bcftools view -f PASS '$vcf' 2>/dev/null | grep -vc '^#')
    # 按 SVTYPE 计数(PASS)
    types=\$(bcftools view -f PASS '$vcf' 2>/dev/null | bcftools query -f '%INFO/SVTYPE\n' 2>/dev/null | sort | uniq -c)
    del=\$(echo \"\$types\" | awk '\$2==\"DEL\"{print \$1}'); del=\${del:-0}
    dup=\$(echo \"\$types\" | awk '\$2==\"DUP\"{print \$1}'); dup=\${dup:-0}
    inv=\$(echo \"\$types\" | awk '\$2==\"INV\"{print \$1}'); inv=\${inv:-0}
    bnd=\$(echo \"\$types\" | awk '\$2==\"BND\"{print \$1}'); bnd=\${bnd:-0}
    oth=\$((pass-del-dup-inv-bnd))
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' '$study' '$sample' \"\$tot\" \"\$pass\" \"\$del\" \"\$dup\" \"\$inv\" \"\$bnd\" \"\$oth\"
  " >> "$OUT/sv_counts.tsv"
}

# Study A somatic：normal(origin) 用 germline VCF；tumor 用 *_vs_RO_origin.tiddit.tumor.vcf.gz
summarize A RO_origin "$PROJ/output_A/variant_calling/tiddit/RO_origin/RO_origin.tiddit.vcf.gz"
for s in RO_B1TP RO_B2TP RO_tumor1 RO_tumor2 RO_tumor3; do
  summarize A "$s" "$PROJ/output_A/variant_calling/tiddit/${s}_vs_RO_origin/${s}_vs_RO_origin.tiddit.tumor.vcf.gz"
done
for s in L1L2_3M L1L2H_3M L1L2_12M L1L2H_12M L1L2_18M L1L2H_18M; do
  summarize B "$s" "$PROJ/output_B/variant_calling/tiddit/$s/$s.tiddit.vcf.gz"
done
echo "== SV 计数 =="; column -t "$OUT/sv_counts.tsv"
echo "DONE SV → $OUT/sv_counts.tsv"
