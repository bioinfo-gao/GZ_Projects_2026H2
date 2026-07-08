#!/bin/bash
# ============================================================================
# Step 0c — 建合并混合参考：GRCm39 + TG_Cas9 + TG_iHPV（复用 proj13 做法）
#   一次建、12 样共用。A 样本在 Cas9 contig 有覆盖、B L1L2H 在 iHPV contig 有覆盖，
#   互为阴性对照。bwa-mem2 索引交给 sarek 在流程内建。
#   前置：0b 已把 refs/constructs/TG_*.fa 落盘。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
GRCM39="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
CDIR="$PROJ/refs/constructs"
OUTDIR="$PROJ/refs/hybrid"; OUT="$OUTDIR/GRCm39_plus_Cas9_iHPV.fa"
mkdir -p "$OUTDIR"

CONSTRUCTS=$(ls "$CDIR"/TG_*.fa 2>/dev/null | sort)
[ -z "$CONSTRUCTS" ] && { echo "ERROR: 无 TG_*.fa，先跑 0b 落盘构建体序列"; exit 1; }
echo "构建体 contigs:"; echo "$CONSTRUCTS" | sed 's/^/  /'

echo ">> cat GRCm39 + constructs → $OUT"
cat "$GRCM39" $CONSTRUCTS > "$OUT"
conda run -n regular_bioinfo samtools faidx "$OUT"
conda run -n regular_bioinfo samtools dict "$OUT" > "${OUT%.fa}.dict"

echo ">> 校验（应见 TG_ contig）"; grep '^>TG_' "$OUT" || { echo "ERROR: 构建体未进参考"; exit 1; }
awk '{n++; L+=$2} END{printf "  contigs=%d total=%.2f Gb\n", n, L/1e9}' "${OUT}.fai"
echo "DONE 0c → $OUT （bwa-mem2 索引由 sarek 内部构建）"
