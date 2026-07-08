#!/bin/bash
# ============================================================================
# Study B / Step 5 — 拷贝数/非整倍体/倍性（tumor-only，无配对正常）
#   C57BL/6 背景 → GRCm39≈正常；Control-FREEC tumor-only + mosdepth 分箱 →
#   全基因组拷贝数谱、染色体臂级增删、倍性。比较 L1L2 vs L1L2H、3M→12M→18M。
#   前置：B2 完成；0a 装好 control-freec。
#   用法: bash B5_cnv_ploidy.sh L1L2H_18M
# ============================================================================
set -euo pipefail
SAMPLE="${1:?用法: bash B5_cnv_ploidy.sh <sample>}"
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_Cas9_iHPV.fa"
OUT="$PROJ/analysis_B/cnv_ploidy/$SAMPLE"; mkdir -p "$OUT"
RUN(){ conda run -n regular_bioinfo "$@"; }
CRAM=$(ls "$PROJ"/output_B/preprocessing/markduplicates/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && { echo "ERROR: 无 $SAMPLE CRAM（先跑 B2）"; exit 1; }

echo ">> mosdepth 全基因组分箱深度（500kb 窗，快速倍性/CNV 预览）"
cd "$OUT"; RUN mosdepth -t 6 -n --fast-mode --by 500000 -f "$HYBRID" "$SAMPLE" "$CRAM"
echo "   → $OUT/${SAMPLE}.regions.bed.gz（按 chr 归一化中位深度画拷贝数谱）"

cat <<NOTE
>> Control-FREEC tumor-only（更正式的 CNV/非整倍体）：
   [general] ploidy=2; window=50000; chrLenFile=<GRCm39.fai>; 无 [control]
   [sample] mateFile=$CRAM
   → *_CNVs, *_ratio.txt；R 画谱。跨样本比 L1L2 vs L1L2H、随龄。
NOTE
echo "DONE B5 → $OUT"
