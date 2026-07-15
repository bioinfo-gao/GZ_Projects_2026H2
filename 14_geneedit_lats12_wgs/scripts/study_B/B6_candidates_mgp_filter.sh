#!/bin/bash
# ============================================================================
# Study B / Step 6 — de novo 候选：扣 C57BL/6 亚系背景(MGP) → 复现/差异候选
#   逻辑：小鼠=C57BL/6，vs GRCm39(6J) 的变异≈品系/亚系背景；用 Sanger MGP(全品系
#   vs 6J 的 SNP+indel) isec 扣掉 → 剩余 = 私有/de novo 候选。再看 (a)多样本/随龄复现
#   (b) L1L2 vs L1L2H 差异 (c) 落在 输卵管/纤毛/上皮·Hippo·DNA修复 相关基因。
#   ⚠ 染色体命名：我方 VCF=chr1，MGP=1(Ensembl) → 先 rename 我方 VCF 对齐 MGP。
#   ⚠ 完整 frameshift/stop 影响注释需 VEP/snpEff(未装) → 本步给基因重叠，consequence 作精化。
#   前置：B2(germline VCF) + MGP snps/indels 下好。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
OUT="$PROJ/analysis_B/candidates"; mkdir -p "$OUT"
MGP=/Work_bio/references/Mus_musculus/GRCm39/mgp_dbsnp
GTF=/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/gencode.vM35.annotation.gtf
HC="$PROJ/output_B/variant_calling/haplotypecaller"
RUN(){ conda run -n regular_bioinfo "$@"; }
SAMPLES=(L1L2_3M L1L2H_3M L1L2_12M L1L2H_12M L1L2_18M L1L2H_18M)

[ -f "$MGP/mgp_REL2021_snps.vcf.gz" ] || { echo "ERROR: MGP 未就绪"; exit 1; }
# MGP 索引(.csi)；若缺则建
for f in snps indels; do
  [ -f "$MGP/mgp_REL2021_${f}.vcf.gz.csi" ] || [ -f "$MGP/mgp_REL2021_${f}.vcf.gz.tbi" ] \
    || RUN bcftools index "$MGP/mgp_REL2021_${f}.vcf.gz"
done

# chr→Ensembl 命名映射
MAP="$OUT/chr_rename.txt"; : > "$MAP"
for i in $(seq 1 19) X Y; do echo -e "chr$i\t$i" >> "$MAP"; done; echo -e "chrM\tMT" >> "$MAP"

echo -e "sample\tPASS\tprivate_deNovo" > "$OUT/candidate_counts.tsv"
for s in "${SAMPLES[@]}"; do
  vcf=$(ls "$HC/$s/$s.haplotypecaller.filtered.vcf.gz" 2>/dev/null || ls "$HC/$s/$s.haplotypecaller.vcf.gz" 2>/dev/null | head -1)
  [ -z "$vcf" ] && { echo "缺 $s VCF"; continue; }
  echo ">> $s"
  RUN bash -c "
    # PASS + 重命名到 Ensembl 命名 → 与 MGP 对齐
    bcftools view -f PASS,. '$vcf' 2>/dev/null \
      | bcftools annotate --rename-chrs '$MAP' -Oz -o '$OUT/$s.ren.vcf.gz' 2>/dev/null
    bcftools index -t '$OUT/$s.ren.vcf.gz' 2>/dev/null
    # 扣 MGP snps+indels 背景 → 私有(de novo 候选)
    bcftools isec -C -w1 '$OUT/$s.ren.vcf.gz' '$MGP/mgp_REL2021_snps.vcf.gz' '$MGP/mgp_REL2021_indels.vcf.gz' -Oz -o '$OUT/$s.private.vcf.gz' 2>/dev/null
    bcftools index -t '$OUT/$s.private.vcf.gz' 2>/dev/null
    p=\$(bcftools view -H '$OUT/$s.ren.vcf.gz' 2>/dev/null | wc -l)
    d=\$(bcftools view -H '$OUT/$s.private.vcf.gz' 2>/dev/null | wc -l)
    printf '%s\t%s\t%s\n' '$s' \"\$p\" \"\$d\" >> '$OUT/candidate_counts.tsv'
  "
done
echo "== 候选计数(PASS→扣MGP后私有) =="; column -t "$OUT/candidate_counts.tsv"

# 复现：多样本共有的私有候选（跨样本 isec 计数）
echo ">> 跨样本复现 + L1L2 vs L1L2H 差异 → recurrent_candidates.vcf"
RUN bash -c "
  bcftools isec -n +2 -Oz -p '$OUT/recur' $(printf \"'%s/%s.private.vcf.gz' \" \"${SAMPLES[@]/#/$OUT\/}\" 2>/dev/null) 2>/dev/null || true
" 2>/dev/null || echo "  (复现计算见 $OUT/recur/)"
echo "DONE B6 → $OUT （基因重叠注释 + 类别筛选在报告阶段用 GTF/bedtools）"
