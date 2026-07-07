#!/bin/bash
# ============================================================================
# Step 5 — 拷贝数估算（mosdepth / samtools coverage 深度比值）
#   拷贝数 ≈ 构建体 contig 唯一比对平均深度 / 样本自身单拷贝小鼠基线深度
#   本批常无纯 WT 对照 → 用样本自身全基因组常染色体中位深度作基线。
#   同源/交叉比对风险：用 MAPQ 过滤唯一比对（-q 20），避免虚高。
#
#   用法: bash 5_copy_number.sh RAGH_153
# ============================================================================
set -euo pipefail
SAMPLE="${1:?用法: bash 5_copy_number.sh <sample>}"

PROJ="/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_constructs.fa"
OUTDIR="$PROJ/analysis/copy_number/$SAMPLE"; mkdir -p "$OUTDIR"
RUN(){ conda run -n regular_bioinfo "$@"; }
MAPQ=20

CRAM=$(ls "$PROJ"/output_results/preprocessing/markduplicates/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && CRAM=$(ls "$PROJ"/output_results/preprocessing/recalibrated/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && { echo "ERROR: 找不到 $SAMPLE 的 CRAM"; exit 1; }

# 1) 基线：全基因组常染色体(chr1..chr19)中位深度（唯一比对）
echo ">> mosdepth 全基因组深度（MAPQ>=$MAPQ）..."
cd "$OUTDIR"
RUN mosdepth -t 6 -n --fast-mode -Q $MAPQ -f "$HYBRID" "$SAMPLE" "$CRAM"
# summary 里每条 contig 一行；取 chr1..chr19 加权平均作基线
BASELINE=$(awk '$1 ~ /^chr([1-9]|1[0-9])$/ {bp+=$2*$4; L+=$2} END{if(L>0) printf "%.4f", bp/L; else print "NA"}' "${SAMPLE}.mosdepth.summary.txt")
echo "  常染色体基线平均深度: $BASELINE x"

# 2) 每个构建体 contig 平均深度 + 覆盖广度，算拷贝数比值
echo -e "construct\tmean_depth\tbreadth_pct\tbaseline\tcopy_number_ratio" > "$OUTDIR/copy_number.tsv"
for TG in $(grep '^>TG_' "$HYBRID" | tr -d '>' | awk '{print $1}'); do
    read -r md br <<< $(RUN samtools coverage -q $MAPQ -r "$TG" "$CRAM" --reference "$HYBRID" \
                          | awk 'NR==2{print $7" "$6}')  # $7=meandepth $6=coverage%
    cn=$(awk -v m="$md" -v b="$BASELINE" 'BEGIN{ if(b>0) printf "%.2f", m/b; else print "NA"}')
    echo -e "${TG}\t${md}\t${br}\t${BASELINE}\t${cn}" >> "$OUTDIR/copy_number.tsv"
done

echo ">> 拷贝数结果："
column -t "$OUTDIR/copy_number.tsv"
echo "解读提示：本样本自己品系对应的构建体应有可观深度（ratio 反映拷贝数，单敲入约 ~1×/等位，"
echo "          纯合双等位 ~1× 全深度）；其他品系构建体应 ~0（阴性对照，验证特异性）。"
echo "DONE → $OUTDIR/copy_number.tsv"
