#!/bin/bash
# ============================================================================
# Step 3.5（对应客户目标3）— KI 序列完整性核查
#   在人源特异区内做滑动窗口深度扫描，找局部深度骤降(内部缺失)/骤升(重复/重排)。
#   不预设基因边界（避免硬编码可能过时的坐标）：纯深度均匀性扫描，通用于三系。
#   若发现异常窗口，再用 split-read 局部验证（结合脚本4的 split 逻辑，人工核查该区段）。
#
#   前置：先跑 5_copy_number.sh（复用其 mosdepth 基线 summary）。
#   用法: bash 6_ki_integrity_check.sh RAGH_153
# ============================================================================
set -euo pipefail
SAMPLE="${1:?用法: bash 6_ki_integrity_check.sh <sample>}"

PROJ="/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_constructs.fa"
REGIONS="$PROJ/refs/constructs/construct_regions.tsv"
OUTDIR="$PROJ/analysis/ki_integrity/$SAMPLE"; mkdir -p "$OUTDIR"
RUN(){ conda run -n regular_bioinfo "$@"; }
MAPQ=20; WINDOW=500          # 滑动窗口大小(bp)
LOW_FRAC=0.30; HIGH_FRAC=2.5  # 相对该构建体人源区中位深度的异常阈值

CRAM=$(ls "$PROJ"/output_results/preprocessing/markduplicates/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && CRAM=$(ls "$PROJ"/output_results/preprocessing/recalibrated/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && { echo "ERROR: 找不到 $SAMPLE 的 CRAM"; exit 1; }
[ -f "$REGIONS" ] || { echo "ERROR: 缺 $REGIONS"; exit 1; }

for TG in $(grep '^>TG_' "$HYBRID" | tr -d '>' | awk '{print $1}'); do
    read -r HS HE <<< $(awk -v c="$TG" '$1==c && $1!~/^#/ {print $2" "$3}' "$REGIONS")
    [ -z "${HS:-}" ] && { echo "WARN: $TG 无人源区定义，跳过"; continue; }

    # 门控：先看该构建体是否存在（复用脚本5的判断逻辑，避免对阴性构建体做无意义的完整性扫描）
    MD=$(RUN samtools coverage -q $MAPQ -r "$TG:$HS-$HE" "$CRAM" --reference "$HYBRID" | awk 'NR==2{print $7}')
    if awk -v m="$MD" 'BEGIN{exit !(m<0.5)}'; then
        echo "== $TG == 深度 ${MD}x 过低（可能是阴性对照构建体），跳过完整性扫描"; continue
    fi
    echo "== $TG == 人源区 $HS-$HE，深度 ${MD}x，滑动窗口(${WINDOW}bp)扫描"

    BED="$OUTDIR/${TG}.windows.bed"
    awk -v c="$TG" -v s="$HS" -v e="$HE" -v w="$WINDOW" \
        'BEGIN{for(i=s;i<e;i+=w) print c"\t"i"\t"((i+w<e)?i+w:e)}' > "$BED"
    RUN mosdepth -t 4 -n --by "$BED" -Q $MAPQ -f "$HYBRID" "$OUTDIR/${TG}_win" "$CRAM"

    OUT="$OUTDIR/${TG}.integrity_flags.tsv"
    zcat "$OUTDIR/${TG}_win.regions.bed.gz" | awk -v med="$MD" -v lo="$LOW_FRAC" -v hi="$HIGH_FRAC" \
        'BEGIN{OFS="\t"; print "chrom","start","end","depth","flag"}
         {ratio=(med>0)?$4/med:0;
          flag=(ratio<lo)?"LOW(possible_deletion)":(ratio>hi)?"HIGH(possible_dup/rearr)":"ok";
          if(flag!="ok") print $1,$2,$3,$4,flag}' > "$OUT"

    N_FLAG=$(($(wc -l < "$OUT") - 1))
    echo "  异常窗口数: $N_FLAG（阈值：<${LOW_FRAC}x 或 >${HIGH_FRAC}x 中位深度）"
    [ "$N_FLAG" -gt 0 ] && { head -20 "$OUT" | column -t; }   # head先于column,避免 SIGPIPE 在 pipefail 下中断脚本
done
echo "DONE 6 → $OUTDIR （异常区段用 IGV 或脚本4的 split-read 逻辑局部核查）"
