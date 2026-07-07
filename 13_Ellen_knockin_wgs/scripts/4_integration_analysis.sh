#!/bin/bash
# ============================================================================
# Step 4 — 整合位点检测（定向嵌合读段法，精修版 2026-07-07）
#   ⚠ 重要认识（试跑 RAGH_153 揭示）：定点打靶构建体的同源臂=内源序列，臂读段在
#     构建体与内源靶点间多重比对(MAPQ0)。所以"配偶落回内源靶点"是【同源臂共定位】，
#     不是真正的跨结合部读段——on-target 不能靠 MAPQ≥20 的唯一读段来"检出"。
#   因此本脚本分两路，语义不同：
#   (1) ON-TARGET 桥接（不加 MAPQ）：构建体读段配偶落自身内源靶点的计数，佐证定点整合；
#       但真正证明"KI 存在+拷贝数"的是脚本5（人源区深度 0.60=杂合单拷贝）+ 定点设计。
#   (2) OFF-TARGET 筛查（MAPQ≥20 唯一比对）：人源特异读段落到【非】预期位点=潜在脱靶随机
#       整合。RAGH_153 结果为 0 = 无脱靶（好结果）。这才是 MAPQ≥20 该回答的问题。
#
#   试跑教训（docs/试跑经验与教训_0707.md）带来的三处精修：
#   1) MAPQ>=20：滤掉同源臂多重比对读段（同源臂与内源位点全同→MAPQ0），
#      否则内源直系同源位点会冒出假整合信号（MTTH 样本假落 chr5 Htt）。
#   2) 高深度 artifact 黑名单：丢弃深度 > 基线×5 的候选位点（chr1:78.58Mb 达 29×）。
#   3) "构建体是否存在"门控：某构建体人源区几乎无覆盖(ratio<0.1)时，其候选位点是
#      交叉比对噪声 → 跳过（阴性构建体不产整合位点）。
#   并对每个位点注释 on-target（本构建体内源靶点）/ cross-map（他构建体靶点，噪声）/
#   novel（潜在脱靶，重点核查）。
#   Manta 已弃用（试跑证明它对这类结合部零检出且 7.5h/样，见 2_run_sarek.sh --tools tiddit）。
#
#   前置：先跑 5_copy_number.sh（提供 mosdepth 基线 + 人源区深度）。
#   用法: bash 4_integration_analysis.sh RAGH_153
# ============================================================================
set -euo pipefail
SAMPLE="${1:?用法: bash 4_integration_analysis.sh <sample>}"

PROJ="/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_constructs.fa"
REGIONS="$PROJ/refs/constructs/construct_regions.tsv"
OUTDIR="$PROJ/analysis/integration/$SAMPLE"; mkdir -p "$OUTDIR"
RUN(){ conda run -n regular_bioinfo "$@"; }
MAPQ=20; PRESENCE_MIN=0.10; ARTIFACT_FOLD=5; MIN_SUPPORT=3

CRAM=$(ls "$PROJ"/output_results/preprocessing/markduplicates/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && CRAM=$(ls "$PROJ"/output_results/preprocessing/recalibrated/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && { echo "ERROR: 找不到 $SAMPLE 的 CRAM"; exit 1; }
[ -f "$REGIONS" ] || { echo "ERROR: 缺 $REGIONS"; exit 1; }
echo "CRAM: $CRAM"

# 基线：优先复用 5_copy_number.sh 的 mosdepth summary，否则快速估计
SUMM="$PROJ/analysis/copy_number/$SAMPLE/${SAMPLE}.mosdepth.summary.txt"
if [ -f "$SUMM" ]; then
    BASELINE=$(awk '$1 ~ /^chr([1-9]|1[0-9])$/ {bp+=$2*$4; L+=$2} END{if(L>0) printf "%.3f", bp/L; else print "NA"}' "$SUMM")
else
    echo "  (无 mosdepth summary，用 chr1 估基线；建议先跑 5_copy_number.sh)"
    BASELINE=$(RUN samtools coverage -q $MAPQ -r chr1 "$CRAM" --reference "$HYBRID" | awk 'NR==2{print $7}')
fi
ARTIFACT_THRESH=$(awk -v b="$BASELINE" -v f="$ARTIFACT_FOLD" 'BEGIN{printf "%.2f", b*f}')
echo "  基线深度=${BASELINE}x  → artifact 阈值(深度>基线×${ARTIFACT_FOLD})=${ARTIFACT_THRESH}x"

# 各构建体的内源靶点（用于 on-target/cross-map 注释）：contig -> "chrom:start-end"
declare -A ENDO
while read -r c hs he a5s a5e a3s a3e endo rest; do
    [[ "$c" =~ ^# || -z "$c" ]] && continue
    ENDO["$c"]="$endo"
done < "$REGIONS"

CONTIGS=$(grep '^>TG_' "$HYBRID" | tr -d '>' | awk '{print $1}')
SUMMARY="$OUTDIR/integration_summary.tsv"
echo -e "construct\tpresent\tcopy_ratio\ton_target_bridging_reads\tn_offtarget_sites\toff_target_loci\tverdict" > "$SUMMARY"

for TG in $CONTIGS; do
    echo "===== $TG ====="
    # 取人源区坐标 + 内源靶点
    read -r HS HE <<< $(awk -v c="$TG" '$1==c && $1!~/^#/ {print $2" "$3}' "$REGIONS")
    ENDO_LOCUS="${ENDO[$TG]:-NA}"
    [ -z "${HS:-}" ] && { echo "  WARN: $TG 无人源区定义，跳过"; continue; }

    # 门控：本构建体是否存在（人源区唯一比对深度/基线）
    HD=$(RUN samtools coverage -q $MAPQ -r "$TG:$HS-$HE" "$CRAM" --reference "$HYBRID" | awk 'NR==2{print $7}')
    RATIO=$(awk -v m="$HD" -v b="$BASELINE" 'BEGIN{if(b>0)printf "%.3f",m/b; else print 0}')
    PRESENT=$(awk -v r="$RATIO" -v t="$PRESENCE_MIN" 'BEGIN{print (r>=t)?"yes":"no"}')
    echo "  人源区深度=${HD}x  ratio=${RATIO}  present=${PRESENT}  (内源靶点 ${ENDO_LOCUS})"
    if [ "$PRESENT" = "no" ]; then
        echo "  → 构建体在本样本近乎缺失（阴性对照），跳过（其位点均为交叉比对噪声）"
        echo -e "${TG}\tno\t${RATIO}\t0\t0\tNA\tabsent(negative_control)" >> "$SUMMARY"
        continue
    fi

    # (ON-TARGET 桥接读段，同源臂辅助、不加 MAPQ)：落在本构建体、配偶落在自身内源靶点
    #   ±20kb 的读段计数。定点打靶的同源臂=内源序列，臂读段在两处多重比对(MAPQ0)，
    #   配偶落回内源靶点——这正是"整合在设计靶点"的桥接证据（不能用 MAPQ≥20，否则被滤空）。
    #   注意：这不能独立区分"正确整合"vs"臂本就同源"，需与 present/拷贝数(脚本5)合起来判读。
    ec="${ENDO_LOCUS%%:*}"; er="${ENDO_LOCUS#*:}"; es="${er%-*}"; ee="${er#*-}"
    ONTGT=$(RUN samtools view -T "$HYBRID" "$CRAM" "$TG" \
      | awk -v ec="$ec" -v es="$((es-20000))" -v ee="$((ee+20000))" \
            '$7==ec && $8>=es && $8<=ee {n++} END{print n+0}')
    echo "  on-target 桥接读段 → 自身内源靶点 ${ENDO_LOCUS}(±20kb)：${ONTGT} 条（臂辅助，佐证定点整合）"

    # ===== off-target 筛查（MAPQ>=20，只留唯一比对的人源读段）=====
    # (A1) discordant：MAPQ>=20，配偶在别的 contig（滤掉同源臂 MAPQ0 多重比对）
    RUN samtools view -q $MAPQ -T "$HYBRID" "$CRAM" "$TG" \
      | awk -v tg="$TG" '$7!="=" && $7!=tg && $7!="*" {print $7"\t"$8}' \
      | sort -k1,1 -k2,2n > "$OUTDIR/${TG}.discordant_mate_pos.tsv"
    # (A2) split：SA 补充比对落回非构建体
    RUN samtools view -q $MAPQ -T "$HYBRID" "$CRAM" "$TG" \
      | awk '/\tSA:Z:/ {for(i=12;i<=NF;i++) if($i~/^SA:Z:/){split($i,a,":"); print a[3]}}' \
      | tr ';' '\n' | awk -F',' 'NF>=2 && $1!~/^TG_/ {print $1"\t"$2}' \
      | sort -k1,1 -k2,2n > "$OUTDIR/${TG}.split_partner_pos.tsv"
    # (A3) 5kb 窗口聚类
    cat "$OUTDIR/${TG}.discordant_mate_pos.tsv" "$OUTDIR/${TG}.split_partner_pos.tsv" \
      | sort -k1,1 -k2,2n \
      | awk 'BEGIN{OFS="\t"} {win=int($2/5000); key=$1"\t"win; c[key]++; if(!(key in mn)||$2<mn[key])mn[key]=$2; if($2>mx[key])mx[key]=$2}
             END{for(k in c) print k, mn[k], mx[k], c[k]}' \
      | sort -k5,5nr | awk -v m="$MIN_SUPPORT" '$5>=m' > "$OUTDIR/${TG}.raw_sites.tsv"

    # (A4) artifact 黑名单（site 深度>基线×5 丢弃）+ on-target/cross/novel 注释
    OUT="$OUTDIR/${TG}.offtarget_screen.tsv"
    echo -e "chrom\tstart\tend\tsupport_reads\tsite_depth\tclass" > "$OUT"
    N_OFF=0; OFF_LOCI=""
    while read -r chrom win start end support; do
        [ -z "${chrom:-}" ] && continue
        depth=$(RUN samtools coverage -q $MAPQ -r "${chrom}:${start}-${end}" "$CRAM" --reference "$HYBRID" | awk 'NR==2{print $7}')
        isart=$(awk -v d="$depth" -v t="$ARTIFACT_THRESH" 'BEGIN{print (d>t)?"1":"0"}')
        [ "$isart" = "1" ] && { echo "    [artifact 丢弃] ${chrom}:${start}-${end} depth=${depth}x (>${ARTIFACT_THRESH})"; continue; }
        # 分类：own-endo(自身靶点，臂辅助) / cross-map(他构建体靶点，噪声) / novel(潜在脱靶)
        cls="novel_offtarget?"
        for oc in "${!ENDO[@]}"; do
            el="${ENDO[$oc]}"; [ "$el" = "NA" ] && continue
            ec2="${el%%:*}"; er2="${el#*:}"; es2="${er2%-*}"; ee2="${er2#*-}"
            hit=$(awk -v c="$chrom" -v e="$end" -v s="$start" -v ec="$ec2" -v es="$es2" -v ee="$ee2" \
                  'BEGIN{print (c==ec && e>=es-20000 && s<=ee+20000)?"1":"0"}')
            if [ "$hit" = "1" ]; then
                if [ "$oc" = "$TG" ]; then cls="own_endo(arm-assisted)"; else cls="cross-map_${oc}(noise)"; fi
                break
            fi
        done
        echo -e "${chrom}\t${start}\t${end}\t${support}\t${depth}\t${cls}" >> "$OUT"
        if [ "$cls" = "novel_offtarget?" ]; then N_OFF=$((N_OFF+1)); OFF_LOCI="${OFF_LOCI}${chrom}:${start}(${support}r);"; fi
    done < "$OUTDIR/${TG}.raw_sites.tsv"

    [ -z "$OFF_LOCI" ] && OFF_LOCI="none"
    verdict=$(awk -v o="$N_OFF" 'BEGIN{print (o==0)?"present_ontarget_no_offtarget":"OFFTARGET_FOUND_investigate"}')
    echo "  off-target 筛查：${N_OFF} 个 novel 位点（唯一比对 MAPQ≥${MAPQ}，已过 artifact 黑名单）"
    column -t "$OUT"
    echo -e "${TG}\tyes\t${RATIO}\t${ONTGT}\t${N_OFF}\t${OFF_LOCI}\t${verdict}" >> "$SUMMARY"
done

echo; echo "===== 汇总 ($SAMPLE) ====="; column -t -s $'\t' "$SUMMARY"
echo "解读："
echo "  copy_ratio + on_target_bridging：证明构建体存在且整合在设计靶点（present + 桥接读段落自身内源靶点）"
echo "  n_offtarget_sites=0 = 无脱靶随机整合（好结果）；>0 = novel 位点需重点核查"
echo "  阴性构建体(present=no)已跳过 = 特异性验证通过"
echo "DONE → $OUTDIR"
