#!/bin/bash
# ============================================================================
# Study A / Step 6 — 样本 identity 指纹：tumor3 到底是不是 RO_origin 那只小鼠的细胞？
#   (2026-07-16 新增；A3b 把 Brca2 定死之后，tumor3 剩下的唯一开放项已不是"编辑位点"
#    而是"**样本身份**")
#
# 背景：tumor3 在 Pten/Brca1/Brca2 **三个靶基因全 WT**，但它是三个瘤里非整倍体最重的。
#   两种解释：(a) 亲本的**未编辑逃逸亚克隆**（指纹应与 origin 一致）
#            (b) **样本混淆/调包**（指纹应与 origin 不一致）
#
# ⚠ 度量选择（关键）：**不能用朴素的 genotype concordance**——
#   tumor 自带大量 CNV/LOH，会把 het→hom，即使是同一只鼠也会压低一致性；tumor3 恰是
#   非整倍体最重的那个 → 朴素一致性会**系统性冤枉 tumor3**，制造假的"调包"结论。
#   改用 **private allele 计数**：统计"该样本带 ALT、而 origin 在同位点**一条 ALT read 都没有**"
#   的位点数。LOH 只会丢等位、**不会凭空造出 origin 没有的新等位** → 该度量对 LOH 免疫。
#   判读：同一只鼠 → 各样本 private 数量应在同一量级（仅体细胞突变）；
#         不同个体/品系 → private 数量高出 1-2 个数量级（遗传性多态）。
#   对照组：B1TP/B2TP/tumor1/tumor2 已知源自 origin → 它们的 private 数就是"同鼠"基线。
#
# ⚠ 已知局限（必须如实告知客户）：近交系内不同个体之间基因组几乎相同，
#   故本法**能排除"跨品系/跨来源调包"，但无法排除"同品系另一只鼠"**。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
GRCM39="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
OUT="$PROJ/analysis_A/identity_fingerprint"; mkdir -p "$OUT"
MINDP=10          # origin 侧最低深度，低于此不判 private（防覆盖空洞造假阳性）

SAMPLES=(RO_origin RO_B1TP RO_B2TP RO_tumor1 RO_tumor2 RO_tumor3)
CRAMS=(); for s in "${SAMPLES[@]}"; do CRAMS+=("$PROJ/output_A/preprocessing/markduplicates/$s/$s.md.cram"); done
CRAMSTR="${CRAMS[*]}"

# ---- 抽样区间：19 条常染色体各取 2 个 500kb 窗（避开端粒/着丝粒侧），合计 ~19Mb ----
: > "$OUT/regions.bed"
for c in $(seq 1 19); do
  len=$(awk -v C="chr$c" '$1==C{print $2}' "${GRCM39}.fai")
  for frac in 30 60; do
    s=$(( len * frac / 100 )); e=$(( s + 500000 ))
    if [ "$e" -lt "$len" ]; then printf 'chr%s\t%s\t%s\n' "$c" "$s" "$e" >> "$OUT/regions.bed"; fi
  done
done
echo "抽样窗口数: $(wc -l < "$OUT/regions.bed")  总长: $(awk '{s+=$3-$2}END{printf "%.1f Mb", s/1e6}' "$OUT/regions.bed")"

echo "==== 六样本联合 call（抽样区间） ===="
conda run -n regular_bioinfo bash -c "
  bcftools mpileup -f '$GRCM39' -R '$OUT/regions.bed' -a AD,DP -q 20 -Q 15 --ignore-RG -Ou $CRAMSTR 2>/dev/null \
    | bcftools call -m -Oz -o '$OUT/fingerprint.vcf.gz' 2>/dev/null
  bcftools index -t '$OUT/fingerprint.vcf.gz'
  bcftools view -m2 -M2 -v snps '$OUT/fingerprint.vcf.gz' -Oz -o '$OUT/fingerprint.snps.vcf.gz' 2>/dev/null
  bcftools index -t '$OUT/fingerprint.snps.vcf.gz'
  bcftools query -f '%CHROM\t%POS[\t%GT:%AD]\n' '$OUT/fingerprint.snps.vcf.gz' > '$OUT/gt_ad.tsv'
"
echo "双等位 SNP 位点数: $(wc -l < "$OUT/gt_ad.tsv")"

echo "==== private allele 计数（vs RO_origin；对 LOH 免疫） ===="
# 列: chr pos, 然后每样本 GT:AD(ref,alt)。origin 是第 3 列。
awk -F'\t' -v MINDP="$MINDP" '
function alt_ad(f,  p,a){ p=index(f,":"); a=substr(f,p+1); sub(/^[^,]*,/,"",a); sub(/,.*$/,"",a); return a+0 }
function ref_ad(f,  p,a){ p=index(f,":"); a=substr(f,p+1); sub(/,.*$/,"",a); return a+0 }
function has_alt(f,  g){ g=substr(f,1,index(f,":")-1); return (g ~ /1/) }
{
  o=$3; o_dp=ref_ad(o)+alt_ad(o);
  if (o_dp < MINDP) next;            # origin 覆盖不足 → 不判
  n_eval++;
  for(i=4;i<=NF;i++){
    if (has_alt($i) && alt_ad($i)>=3 && alt_ad(o)==0) priv[i]++;
    tot[i]++;
  }
}
END{
  split("RO_B1TP RO_B2TP RO_tumor1 RO_tumor2 RO_tumor3", nm, " ");
  printf "评估位点数(origin DP>=%d): %d\n\n", MINDP, n_eval;
  printf "%-11s %14s %14s\n", "sample", "private_alleles", "per_10k_sites";
  for(i=4;i<=8;i++) printf "%-11s %14d %14.1f\n", nm[i-3], priv[i]+0, (priv[i]+0)/n_eval*10000;
}' "$OUT/gt_ad.tsv" | tee "$OUT/private_alleles.txt"

echo "DONE A6 → $OUT"
