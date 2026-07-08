#!/bin/bash
# ============================================================================
# Step 6 — 合子型（纯合 vs 杂合 KI）
#   方法：KI 等位会把"两条同源臂之间的小鼠原生序列"整段切除替换。若杂合，样本还
#   保留一条未编辑的小鼠等位 → 该"被切除区段"在小鼠基因组坐标上应仍有约一半基线
#   深度的覆盖；若纯合 KI，该区段应接近零覆盖（两条等位都被替换/无处比对）。
#
#   "被切除区段"坐标不查记忆/硬编码，而是【现场】用 minimap2 把该构建体的同源臂
#   （从 regions.tsv 的 arm5/arm3 坐标里现取序列）重新比对到 GRCm39 定位 ——
#   避免依赖可能记错的历史坐标（教训见 memory: reading-images-carefully 的姊妹原则，
#   对"记忆里的坐标"同样要重新验证，不能凭印象硬编码）。
#
#   前置：先跑 5_copy_number.sh（复用 mosdepth 基线）。
#   用法: bash 7_zygosity_analysis.sh RAGH_153
# ============================================================================
set -euo pipefail
SAMPLE="${1:?用法: bash 7_zygosity_analysis.sh <sample>}"

PROJ="/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_constructs.fa"
GRCM39="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
REGIONS="$PROJ/refs/constructs/construct_regions.tsv"
OUTDIR="$PROJ/analysis/zygosity/$SAMPLE"; mkdir -p "$OUTDIR"
RUN(){ conda run -n regular_bioinfo "$@"; }
MAPQ=20

CRAM=$(ls "$PROJ"/output_results/preprocessing/markduplicates/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && CRAM=$(ls "$PROJ"/output_results/preprocessing/recalibrated/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && { echo "ERROR: 找不到 $SAMPLE 的 CRAM"; exit 1; }
SUMM="$PROJ/analysis/copy_number/$SAMPLE/${SAMPLE}.mosdepth.summary.txt"
[ -f "$SUMM" ] || { echo "ERROR: 缺 mosdepth summary，先跑 5_copy_number.sh"; exit 1; }
BASELINE=$(awk '$1 ~ /^chr([1-9]|1[0-9])$/ {bp+=$2*$4; L+=$2} END{if(L>0) printf "%.3f", bp/L}' "$SUMM")

echo -e "construct\tdeleted_mouse_span\tspan_depth\tbaseline\tratio\tzygosity_call" > "$OUTDIR/zygosity_summary.tsv"

while read -r TG HS HE A5S A5E A3S A3E ENDO REST; do
    [[ "$TG" =~ ^# || -z "$TG" ]] && continue
    MD=$(RUN samtools coverage -q $MAPQ -r "$TG:$HS-$HE" "$CRAM" --reference "$HYBRID" | awk 'NR==2{print $7}')
    if awk -v m="$MD" 'BEGIN{exit !(m<0.5)}'; then
        echo "== $TG == 深度过低（阴性对照构建体），跳过"; continue
    fi
    echo "== $TG == 现场比对同源臂到 GRCm39 定位被切除区段..."

    # 现取 arm5/arm3 序列（构建体 contig 内坐标，来自 regions.tsv）
    RUN samtools faidx "$HYBRID" "${TG}:$((A5S+1))-${A5E}" > "$OUTDIR/${TG}_arm5.fa" 2>/dev/null
    RUN samtools faidx "$HYBRID" "${TG}:$((A3S+1))-${A3E}" > "$OUTDIR/${TG}_arm3.fa" 2>/dev/null
    A5_HIT=$(RUN minimap2 -x asm5 --secondary=no "$GRCM39" "$OUTDIR/${TG}_arm5.fa" 2>/dev/null | awk '{print $6"\t"$8"\t"$9; exit}')
    A3_HIT=$(RUN minimap2 -x asm5 --secondary=no "$GRCM39" "$OUTDIR/${TG}_arm3.fa" 2>/dev/null | awk '{print $6"\t"$8"\t"$9; exit}')
    if [ -z "$A5_HIT" ] || [ -z "$A3_HIT" ]; then
        echo "  WARN: 同源臂比对失败，跳过 $TG 合子型判定"; continue
    fi
    A5_CHR=$(echo "$A5_HIT" | cut -f1); A5_END=$(echo "$A5_HIT" | cut -f3)
    A3_CHR=$(echo "$A3_HIT" | cut -f1); A3_START=$(echo "$A3_HIT" | cut -f2)
    if [ "$A5_CHR" != "$A3_CHR" ] || [ "$A5_END" -ge "$A3_START" ]; then
        echo "  WARN: 两臂比对不构成合理区间(${A5_CHR}:${A5_END} / ${A3_CHR}:${A3_START})，跳过"; continue
    fi
    DEL_SPAN="${A5_CHR}:${A5_END}-${A3_START}"
    echo "  被切除区段（现场定位）: $DEL_SPAN"

    SPAN_DEPTH=$(RUN samtools coverage -q $MAPQ -r "$DEL_SPAN" "$CRAM" --reference "$HYBRID" | awk 'NR==2{print $7}')
    RATIO=$(awk -v m="$SPAN_DEPTH" -v b="$BASELINE" 'BEGIN{if(b>0)printf "%.3f",m/b; else print 0}')
    CALL=$(awk -v r="$RATIO" 'BEGIN{ if(r<0.15){print "homozygous_KI(~0_residual_WT_allele)"} else if(r<0.75){print "heterozygous(~1_WT_allele_remains)"} else {print "unexpected_high(check:_incomplete_deletion_or_mosaicism)"} }')
    echo "  被切除区段深度=${SPAN_DEPTH}x  基线=${BASELINE}x  ratio=${RATIO}  → ${CALL}"
    echo -e "${TG}\t${DEL_SPAN}\t${SPAN_DEPTH}\t${BASELINE}\t${RATIO}\t${CALL}" >> "$OUTDIR/zygosity_summary.tsv"
done < "$REGIONS"

echo; echo "===== 汇总 ($SAMPLE) ====="; column -t "$OUTDIR/zygosity_summary.tsv"
echo "解读：ratio~0=纯合KI（两条等位都被替换，原生区无处比对）；ratio~0.5=杂合（一条WT等位仍在）；"
echo "      ratio~1.0 或更高=异常，需核查是否切除不完全/嵌合。"
echo "DONE 7 → $OUTDIR"
