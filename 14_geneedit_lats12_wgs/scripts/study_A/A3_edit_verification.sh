#!/bin/bash
# ============================================================================
# Study A / Step 3 — 编辑验证：Brca1/Brca2/Pten CRISPR KO 是否发生 + 肿瘤谱系反推
#   在 sgRNA 切点比 编辑细胞/肿瘤 vs RO_origin，查 indel/frameshift、等位比例。
#   sgRNA: Pten×3(B1TP+B2TP)、Brca1×3(B1TP)、Brca2×3(B2TP)。
#   ⚠ 2026-07-16 更新：客户补齐 Brca2×3 guide（切点 chr5:150452957/961/989，均已核 NGG PAM）
#     → Brca2 从「无 guide、只能扫全基因」升级为与 Brca1/Pten 同级的**定向切点验证**。
#     (3) 的全基因扫描保留为补充（捕捉切点窗口外的大 indel），不再是主证据。
#   前置：A2 sarek 完成（output_A CRAM）。GRCm39 参考。
#
#   流程：
#   (1) seqkit locate 每条 spacer(+反向互补) → 限定在靶基因内 → 切点=PAM上游3bp
#   (2) 各切点±window 多样本 bcftools 联合 call → 提 indel，比各样本 vs origin 的 GT/AD
#   (3) Brca2 全基因扫 indel（补充：切点窗口外的大 indel）
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
# ⚠ 基因名必须**整格精确匹配**(`^| Brca2 |`)，不能只匹配前缀：guides.md 里另有说明性表格
#   （行首形如 `| Brca2+150529497 |`）会被前缀式 grep 误抓成 guide 行 → 该行无 20-nt spacer
#   → 配合 set -e 令整个脚本静默死在第(1)步。2026-07-16 已被此 bug 咬过一次。
# ⚠ 用 if 而非 `[ -n ]&&echo`：后者在最后一行不匹配时返回 1 → while 退出码 1 → set -e/pipefail 杀脚本。
grep -E '^\| *(Pten|Brca1|Brca2) *\|' "$GUIDES" | while IFS='|' read -r _ gene _ spacer _; do
  gene=$(echo "$gene" | tr -d ' '); sp=$(echo "$spacer" | grep -oE '[ACGT]{20}' || true)
  if [ -n "$sp" ]; then echo -e "${gene}\t${sp}"; fi
done > "$OUT/spacer_gene.tsv"

# 硬断言：必须恰好抓到 9 条 guide (Pten×3 + Brca1×3 + Brca2×3)，否则解析出错，宁可炸不可静默少跑
n_sp=$(wc -l < "$OUT/spacer_gene.tsv")
if [ "$n_sp" -ne 9 ]; then
  echo "❌ spacer 解析异常：期望 9 条，实得 $n_sp 条 —— 检查 $GUIDES 表格式" >&2; exit 1
fi
for g in Pten Brca1 Brca2; do
  n=$(awk -F'\t' -v g="$g" '$1==g' "$OUT/spacer_gene.tsv" | wc -l)
  if [ "$n" -ne 3 ]; then echo "❌ $g 期望 3 条 guide，实得 $n 条" >&2; exit 1; fi
done
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

echo "==== (3) Brca2 全基因扫 indel — 补充，捕捉切点窗口外大 indel (B2TP/肿瘤 vs origin) ===="
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
