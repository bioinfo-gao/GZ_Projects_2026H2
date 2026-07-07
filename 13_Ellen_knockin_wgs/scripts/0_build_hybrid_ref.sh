#!/bin/bash
# ============================================================================
# Step 0 — 构建合并混合参考 (hybrid reference)
#   GRCm39 (小鼠) + 每个构建体各作独立 contig  →  faidx + dict
#   bwa-mem2 索引交给 sarek 在流程内构建（一次，随参考变化重建）。
#
#   每个样本只在自己构建体 contig 上有覆盖，在其他构建体上零覆盖
#   = 自带跨品系阴性对照（特异性验证）。
#
#   试跑阶段：refs/constructs/ 下现有 TG_RAGH.fa + TG_MTTH.fa。
#   明天 CD1A 真实序列到位后，把 TG_CD1A.fa 放进 refs/constructs/ 重跑本脚本即可。
# ============================================================================
set -euo pipefail

PROJ="/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
GRCM39="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
CONSTRUCT_DIR="$PROJ/refs/constructs"
OUTDIR="$PROJ/refs/hybrid"
OUT="$OUTDIR/GRCm39_plus_constructs.fa"

mkdir -p "$OUTDIR"
cd "$PROJ"

# 收集构建体 contig（只取 TG_*.fa，排除 homology_arms / exons 等辅助文件）
CONSTRUCTS=$(ls "$CONSTRUCT_DIR"/TG_*.fa 2>/dev/null | sort)
if [ -z "$CONSTRUCTS" ]; then echo "ERROR: no TG_*.fa in $CONSTRUCT_DIR"; exit 1; fi
echo "构建体 contigs:"; echo "$CONSTRUCTS" | sed 's/^/  /'

# 1) 拼接（小鼠在前，构建体在后）
echo ">> cat GRCm39 + constructs -> $OUT"
cat "$GRCM39" $CONSTRUCTS > "$OUT"

# 2) faidx + dict
echo ">> samtools faidx / dict"
conda run -n regular_bioinfo samtools faidx "$OUT"
conda run -n regular_bioinfo samtools dict "$OUT" > "${OUT%.fa}.dict"

# 3) 校验：确认构建体 contig 名进入参考
echo ">> hybrid 参考尾部 contigs（应见 TG_*）："
grep '^>TG_' "$OUT" || { echo "ERROR: 构建体 contig 未进入参考"; exit 1; }
echo ">> 参考序列数 / 总长："
conda run -n regular_bioinfo samtools faidx "$OUT" 2>/dev/null
awk '{n++; L+=$2} END{printf "  contigs=%d  total=%.2f Gb\n", n, L/1e9}' "${OUT}.fai"

echo "DONE. hybrid 参考: $OUT"
echo "（bwa-mem2 索引由 sarek 在流程内构建，无需在此预建）"
