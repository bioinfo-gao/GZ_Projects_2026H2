#!/bin/bash
# ============================================================================
# Step 7（CD1A 专项）— Neo 盒状态核查
#   `CD1A KI.dna` 里 NeoR/KanR 完整存在（不像 RAGH/MTTH 标注"Neo Deleted"）——
#   不能假设实际小鼠是否已删除，直接查 CD1A_B125 在 Neo 坐标处的 WGS 覆盖深度。
#
#   Neo 坐标不硬编码，现场重新解析 "CD1A KI.dna" 取 NeoR/KanR feature 坐标——
#   TG_CD1A.fa 是从同一份 .dna 直接导出，坐标系一致，可直接用于 hybrid 参考。
#
#   用法: bash 8_cd1a_neo_status.sh CD1A_B125
# ============================================================================
set -euo pipefail
SAMPLE="${1:?用法: bash 8_cd1a_neo_status.sh <sample>}"

PROJ="/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_constructs.fa"
DNA="$PROJ/refs/constructs/CD1A KI.dna"
OUTDIR="$PROJ/analysis/cd1a_neo_status/$SAMPLE"; mkdir -p "$OUTDIR"
RUN(){ conda run -n regular_bioinfo "$@"; }

CRAM=$(ls "$PROJ"/output_results/preprocessing/markduplicates/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && CRAM=$(ls "$PROJ"/output_results/preprocessing/recalibrated/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && { echo "ERROR: 找不到 $SAMPLE 的 CRAM"; exit 1; }
[ -f "$DNA" ] || { echo "ERROR: 缺 $DNA"; exit 1; }

echo ">> 现场解析 '$DNA' 定位 NeoR/KanR 坐标（不用记忆里的数字）"
read -r NEO_S NEO_E <<< $(RUN python3 -c "
from Bio import SeqIO
rec = SeqIO.read('$DNA', 'snapgene')
for ft in rec.features:
    lab = ft.qualifiers.get('label', [''])[0]
    if 'NeoR' in lab or 'KanR' in lab:
        print(int(ft.location.start), int(ft.location.end)); break
")
[ -z "${NEO_S:-}" ] && { echo "ERROR: 未在构建体里找到 NeoR/KanR feature"; exit 1; }
echo "  NeoR/KanR 坐标（构建体内）: $NEO_S-$NEO_E"

echo ">> 样本在 Neo 坐标处的覆盖深度（不加 MAPQ 过滤——Neo 序列是载体特有，无小鼠同源风险）"
RUN samtools coverage -r "TG_CD1A:$((NEO_S+1))-${NEO_E}" "$CRAM" --reference "$HYBRID" | tee "$OUTDIR/neo_coverage.tsv"

# 对比：同一样本在 CD1A 人源区（应有信号）的深度，作为"有覆盖时长什么样"的参照
REGIONS="$PROJ/refs/constructs/construct_regions.tsv"
read -r HS HE <<< $(awk -v c="TG_CD1A" '$1==c && $1!~/^#/ {print $2" "$3}' "$REGIONS")
REF_DEPTH=$(RUN samtools coverage -q 20 -r "TG_CD1A:$HS-$HE" "$CRAM" --reference "$HYBRID" | awk 'NR==2{print $7}')
NEO_DEPTH=$(awk 'NR==2{print $7}' "$OUTDIR/neo_coverage.tsv")

echo
echo "CD1A 人源区参照深度: ${REF_DEPTH}x   Neo 盒深度: ${NEO_DEPTH}x"
VERDICT=$(awk -v n="$NEO_DEPTH" -v r="$REF_DEPTH" 'BEGIN{
    if(r<=0){print "参照深度为0，无法判断（该样本可能非CD1A系或构建体缺失）"}
    else if(n/r < 0.15){print "Neo 已删除（深度远低于参照，与 RAGH/MTTH 一致）"}
    else {print "Neo 仍然存在（深度接近参照）——与 RAGH/MTTH 不同，需在报告中如实说明"}
}')
echo "判定: $VERDICT"
echo "$VERDICT" > "$OUTDIR/verdict.txt"
echo "DONE 8 → $OUTDIR"
