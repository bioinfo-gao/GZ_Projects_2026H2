#!/bin/bash
# ============================================================================
# Step 5 — 拷贝数估算（深度比值，只在人源特异区）
#   拷贝数 ≈ 构建体【人源特异区】唯一比对深度 / 样本自身常染色体基线深度
#   关键：遮蔽小鼠同源臂——同源臂与内源位点全同→多重比对→MAPQ 过滤后为 0，
#         若按全长 contig 平均会严重低估拷贝数（试跑实测 0.08 vs 真实 0.60）。
#   人源区坐标读自 refs/constructs/construct_regions.tsv；MAPQ>=20 去交叉比对虚高；
#   本批无纯 WT → 样本自身常染色体(chr1-19)中位深度作基线。
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

# 2) 拷贝数：只在【人源特异区】算（遮蔽小鼠同源臂——同源臂与内源位点全同→多重比对→
#    被 MAPQ 过滤掉，若按全长 contig 平均会严重低估；见 docs/试跑经验与教训_0707.md §六）
REGIONS="$PROJ/refs/constructs/construct_regions.tsv"
[ -f "$REGIONS" ] || { echo "ERROR: 缺 $REGIONS"; exit 1; }

echo -e "construct\thuman_region\thuman_meandepth\thuman_breadth_pct\twholecontig_meandepth\twholecontig_breadth_pct\tbaseline\tcopy_number_ratio\tinterpretation" > "$OUTDIR/copy_number.tsv"
for TG in $(grep '^>TG_' "$HYBRID" | tr -d '>' | awk '{print $1}'); do
    # 从 regions 表取人源区坐标（跳过注释/表头行）
    read -r HS HE <<< $(awk -v c="$TG" '$1==c && $1!~/^#/ {print $2" "$3}' "$REGIONS")
    if [ -z "${HS:-}" ]; then
        echo "  WARN: $TG 无人源区定义，回退全长（拷贝数会偏低，请补 construct_regions.tsv）"
        REGION="$TG"
    else
        REGION="$TG:$HS-$HE"
    fi
    # 人源区（真实拷贝数）
    read -r mdh brh <<< $(RUN samtools coverage -q $MAPQ -r "$REGION" "$CRAM" --reference "$HYBRID" | awk 'NR==2{print $7" "$6}')
    # 全长 contig（仅供参照，会被同源臂稀释）
    read -r mdw brw <<< $(RUN samtools coverage -q $MAPQ -r "$TG" "$CRAM" --reference "$HYBRID" | awk 'NR==2{print $7" "$6}')
    cn=$(awk -v m="$mdh" -v b="$BASELINE" 'BEGIN{if(b>0) printf "%.2f", m/b; else print "NA"}')
    interp=$(awk -v c="$cn" 'BEGIN{ if(c=="NA"){print "NA"} else if(c<0.1){print "absent(neg_control)"} else if(c<0.75){print "~het_single_copy"} else if(c<1.5){print "~hom_single_copy"} else {print "multi_copy/concatemer_check"} }')
    echo -e "${TG}\t${REGION}\t${mdh}\t${brh}\t${mdw}\t${brw}\t${BASELINE}\t${cn}\t${interp}" >> "$OUTDIR/copy_number.tsv"
done

echo ">> 拷贝数结果（human_* 为真实值；wholecontig_* 仅参照）："
column -t "$OUTDIR/copy_number.tsv"
echo "解读：ratio ~0.5=杂合单拷贝 / ~1.0=纯合单拷贝或双等位 / >1.5=多拷贝(查串联)；本样自己品系构建体应有可观 ratio，"
echo "      其他品系构建体应 ~0（阴性对照）。全长列若远低于 human 列，即同源臂稀释的证据。"
echo "DONE → $OUTDIR/copy_number.tsv"
