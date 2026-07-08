#!/bin/bash
# ============================================================================
# Study B / Step 4 — 工程等位状态 + 体细胞/渗漏重组（第二假设）
#   (1) Lats1/2 floxed 外显子覆盖：是否完整（无 Cre 应完整）；若某样本/老龄局部覆盖
#       下降 = 渗漏/体细胞 loxP 重组删了外显子 → 局部 Lats 失活。
#   (2) iHPV 的 EGFP-pA "stop" 盒覆盖：是否有体细胞丢失 = E6/E7 被渗漏激活的证据。
#   随龄(18M)尤其看这类嵌合事件。
#   前置：B2 完成；需 Lats1/2 loxP/floxed 外显子坐标（待 0b/文献补）+ iHPV contig 上 stop 盒坐标。
#   用法: bash B4_engineered_alleles.sh L1L2_18M
# ============================================================================
set -euo pipefail
SAMPLE="${1:?用法: bash B4_engineered_alleles.sh <sample>}"
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_Cas9_iHPV.fa"
OUT="$PROJ/analysis_B/engineered_alleles/$SAMPLE"; mkdir -p "$OUT"
RUN(){ conda run -n regular_bioinfo "$@"; }
CRAM=$(ls "$PROJ"/output_B/preprocessing/markduplicates/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && { echo "ERROR: 无 $SAMPLE CRAM（先跑 B2）"; exit 1; }

# TODO(审阅后填)：GRCm39 上 Lats1/Lats2 floxed 外显子坐标；TG_iHPV 上 EGFP-pA stop 盒坐标
LATS_EXONS_BED="$PROJ/refs/lats12_floxed_exons.bed"   # 待建
STOP_REGION="TG_iHPV:START-END"                        # 待填（EGFP-pA）

echo "== Lats1/2 floxed 外显子深度（应完整；局部骤降=体细胞重组）=="
if [ -f "$LATS_EXONS_BED" ]; then
  RUN mosdepth -t 4 -n -b "$LATS_EXONS_BED" -f "$HYBRID" "$OUT/${SAMPLE}_lats" "$CRAM"
  zcat "$OUT/${SAMPLE}_lats.regions.bed.gz" | column -t
else echo "  （待建 lats12_floxed_exons.bed）"; fi

echo "== iHPV stop 盒(EGFP-pA)深度（体细胞丢失=E6/E7 渗漏激活）=="
echo "  待填 STOP_REGION 后: samtools coverage -r $STOP_REGION"
echo "DONE B4（框架；待坐标补齐运行）→ $OUT"
