#!/bin/bash
# ============================================================================
# Study A / Step 3 — 编辑验证：Brca1/Brca2/Pten CRISPR KO 是否发生 + 肿瘤谱系反推
#   在 sgRNA 切点比 编辑细胞/肿瘤 vs RO_origin，查 indel/frameshift、等位比例。
#   sgRNA: Pten×3(B1TP+B2TP)、Brca1×3(B1TP)；Brca2 未提供 → 扫 Brca2 全基因 indel。
#   前置：A2 sarek 完成（output_A CRAM）。GRCm39 参考。
#
#   流程：
#   (1) seqkit locate 每条 spacer(+反向互补) → 限定在靶基因内 → 切点=PAM上游3bp
#   (2) 各切点±window 多样本 bcftools 联合 call → 提 indel，比各样本 vs origin 的 GT/AD
#   (3) Brca2 全基因扫 indel（B2TP/肿瘤有、origin 无）
#   (4) 汇总 KO 判定 + 用 Brca1 vs Brca2 KO 反推每个 tumor 源自 B1TP/B2TP
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
GRCM39="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
GUIDES="$PROJ/docs/client_materials/sgRNA_guides.md"
OUT="$PROJ/analysis_A/edit_verification"; mkdir -p "$OUT"
RUN(){ conda run -n regular_bioinfo "$@"; }
WIN=60   # 切点±window(bp) 做局部 call

# 样本(顺序固定，origin 放第一列作对照)
SAMPLES=(RO_origin RO_B1TP RO_B2TP RO_tumor1 RO_tumor2 RO_tumor3)
CRAMS=(); for s in "${SAMPLES[@]}"; do CRAMS+=("$PROJ/output_A/preprocessing/markduplicates/$s/$s.md.cram"); done

# 靶基因坐标(GRCm39, GENCODE vM35) 与 spacer→基因映射
declare -A GENE_LOCUS=( [Pten]="chr19:32734897-32803560" [Brca1]="chr11:101379590-101442781" [Brca2]="chr5:150446095-150493794" )

echo "==== (1) 定位 sgRNA 切点（仅在靶基因区搜索，快） ===="
# 从 guides.md 抽 spacer 并带基因标签（表格行: | 基因 | ... | `SPACER` | ...）
grep -E '^\| *(Pten|Brca1|Brca2)' "$GUIDES" | while IFS='|' read -r _ gene _ spacer _; do
  gene=$(echo "$gene" | tr -d ' '); sp=$(echo "$spacer" | grep -oE '[ACGT]{20}' || true)
  [ -n "$sp" ] && echo -e "${gene}\t${sp}"
done > "$OUT/spacer_gene.tsv"
echo "spacer→基因:"; cat "$OUT/spacer_gene.tsv"
awk -F'\t' '{print ">"$1"_g"NR"\n"$2}' "$OUT/spacer_gene.tsv" > "$OUT/spacers.fa"

# 只提取靶基因区做搜索底库（header 形如 chr19:32734897-32803560）
: > "$OUT/target_genes.fa"
for g in "${!GENE_LOCUS[@]}"; do RUN samtools faidx "$GRCM39" "${GENE_LOCUS[$g]}" >> "$OUT/target_genes.fa"; done
RUN seqkit locate --bed -f "$OUT/spacers.fa" "$OUT/target_genes.fa" > "$OUT/spacer_hits.bed" 2>/dev/null || true
echo "靶基因区命中:"; wc -l < "$OUT/spacer_hits.bed"

# 局部坐标→基因组坐标，算切点(protospacer 近 PAM 端内 3bp)
: > "$OUT/cut_sites.tsv"
echo -e "gene\tchr\tcut_pos\tstrand\tspacer_name" >> "$OUT/cut_sites.tsv"
while read -r region ls le name score strand; do
  base=${name%%_g*}
  rchr=${region%%:*}; rstart=${region#*:}; rstart=${rstart%-*}   # 区间起点(1-based)
  gs=$((rstart + ls)); ge=$((rstart + le - 1))                    # 命中的基因组 1-based 起止
  # 仅保留落在正确靶基因染色体的命中
  locus=${GENE_LOCUS[$base]:-}; [ -z "$locus" ] && continue
  [ "$rchr" != "${locus%%:*}" ] && continue
  if [ "$strand" = "+" ]; then cut=$((ge-3)); else cut=$((gs+3)); fi
  echo -e "${base}\t${rchr}\t${cut}\t${strand}\t${name}" >> "$OUT/cut_sites.tsv"
done < "$OUT/spacer_hits.bed"
echo "靶基因内切点:"; column -t "$OUT/cut_sites.tsv"

SAMP_HDR=$(printf '%s\t' "${SAMPLES[@]}")
CRAMSTR="${CRAMS[*]}"

echo "==== (2) 各切点 ±${WIN}bp 多样本联合 call indel ===="
# ⚠ bcftools -R 需 BED(tab分隔)，不能用 chr:s-e 字符串；重叠窗先 bedtools merge
: > "$OUT/_regions.bed"
tail -n +2 "$OUT/cut_sites.tsv" | while read -r gene chr cut strand name; do
  s=$((cut-WIN)); [ "$s" -lt 0 ] && s=0
  printf '%s\t%s\t%s\n' "$chr" "$s" "$((cut+WIN))" >> "$OUT/_regions.bed"
done
# ⚠ pipe 必须在单个 conda run 内（conda run 不转发管道 stdin）
RUN bash -c "sort -k1,1 -k2,2n '$OUT/_regions.bed' | bedtools merge > '$OUT/regions.bed'"
if [ -s "$OUT/regions.bed" ]; then
  # ⚠ 整条管线放进单个 conda run（双 conda-run 管道会静默失败）
  conda run -n regular_bioinfo bash -c "
    bcftools mpileup -f '$GRCM39' -R '$OUT/regions.bed' -a AD,DP -q 20 -Q 15 --ignore-RG -Ou $CRAMSTR 2>/dev/null \
      | bcftools call -mv -Oz -o '$OUT/cutsite_variants.vcf.gz' 2>/dev/null
    bcftools index -t '$OUT/cutsite_variants.vcf.gz'
    { printf 'chr\tpos\tref\talt\ttype\t$SAMP_HDR\n';
      bcftools view -v indels '$OUT/cutsite_variants.vcf.gz' 2>/dev/null \
        | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\tINDEL[\t%GT:%AD]\n'; } > '$OUT/cutsite_indels.tsv'
  "
  echo "切点 indel:"; column -t -s $'\t' "$OUT/cutsite_indels.tsv" | head -40
else
  echo "⚠ 无靶基因内切点命中(检查 spacer 定位)"
fi

echo "==== (3) Brca2 全基因扫 indel (B2TP/肿瘤 vs origin) ===="
conda run -n regular_bioinfo bash -c "
  bcftools mpileup -f '$GRCM39' -r '${GENE_LOCUS[Brca2]}' -a AD,DP -q 20 -Q 15 --ignore-RG -Ou $CRAMSTR 2>/dev/null \
    | bcftools call -mv -Oz -o '$OUT/brca2_variants.vcf.gz' 2>/dev/null
  bcftools index -t '$OUT/brca2_variants.vcf.gz'
  { printf 'chr\tpos\tref\talt\t$SAMP_HDR\n';
    bcftools view -v indels '$OUT/brca2_variants.vcf.gz' 2>/dev/null \
      | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%AD]\n'; } > '$OUT/brca2_indels.tsv'
"
echo "Brca2 indel 数(含背景):"; tail -n +2 "$OUT/brca2_indels.tsv" | wc -l

echo "DONE A3 → $OUT （判读见 cutsite_indels.tsv / brca2_indels.tsv；(4)谱系反推在 A3b_summary.py）"
