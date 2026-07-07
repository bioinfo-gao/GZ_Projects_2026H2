#!/bin/bash
# ============================================================================
# Step 4 — 整合位点 / 结合部检测（两路交叉验证）
#   (A) 定向嵌合读段：抓落在构建体 contig 上的读段 → 看其配偶(discordant)
#       和软剪切补充比对(SA/split) 落回小鼠何处 → 聚类成候选整合位点
#   (B) 通用 SV caller：从 sarek 的 Manta VCF 里提取涉及 TG_ contig 的 BND
#
#   用法: bash 4_integration_analysis.sh RAGH_153
# ============================================================================
set -euo pipefail
SAMPLE="${1:?用法: bash 4_integration_analysis.sh <sample>}"

PROJ="/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_constructs.fa"
OUTDIR="$PROJ/analysis/integration/$SAMPLE"; mkdir -p "$OUTDIR"
RUN(){ conda run -n regular_bioinfo "$@"; }

# 定位 sarek 的 CRAM（markduplicates 输出）
CRAM=$(ls "$PROJ"/output_results/preprocessing/markduplicates/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && CRAM=$(ls "$PROJ"/output_results/preprocessing/recalibrated/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && { echo "ERROR: 找不到 $SAMPLE 的 CRAM"; exit 1; }
echo "CRAM: $CRAM"

# 构建体 contig 列表
CONTIGS=$(grep '^>TG_' "$HYBRID" | tr -d '>' | awk '{print $1}')
echo "构建体 contigs: $CONTIGS"

for TG in $CONTIGS; do
    echo "===== $TG ====="
    # (A1) discordant：落在 TG contig、但配偶在别的染色体（RNEXT != TG 且 != '=')
    #      输出配偶所在 mouse 染色体+位置 → 候选整合位点
    RUN samtools view -T "$HYBRID" "$CRAM" "$TG" \
      | awk -v tg="$TG" '$7!="=" && $7!=tg && $7!="*" {print $7"\t"$8}' \
      | sort -k1,1 -k2,2n > "$OUTDIR/${TG}.discordant_mate_pos.tsv"

    # (A2) split：软剪切 + 有 SA 补充比对标签的读段（跨结合部）
    RUN samtools view -T "$HYBRID" "$CRAM" "$TG" \
      | awk '/\tSA:Z:/ {for(i=12;i<=NF;i++) if($i~/^SA:Z:/){split($i,a,":"); print a[3]}}' \
      | tr ';' '\n' | awk -F',' 'NF>=2 && $1!~/^TG_/ {print $1"\t"$2}' \
      | sort -k1,1 -k2,2n > "$OUTDIR/${TG}.split_partner_pos.tsv"

    # (A3) 聚类：把 discordant 配偶位置按 5 kb 窗口聚合计数 → 候选整合位点表
    cat "$OUTDIR/${TG}.discordant_mate_pos.tsv" "$OUTDIR/${TG}.split_partner_pos.tsv" \
      | sort -k1,1 -k2,2n \
      | awk 'BEGIN{OFS="\t"} {win=int($2/5000); key=$1"\t"win; c[key]++; if(!(key in mn)||$2<mn[key])mn[key]=$2; if($2>mx[key])mx[key]=$2}
             END{for(k in c) print k, mn[k], mx[k], c[k]}' \
      | sort -k5,5nr \
      | awk 'BEGIN{OFS="\t"; print "chrom","win5kb","start","end","support_reads"} $5>=3' \
      > "$OUTDIR/${TG}.candidate_integration_sites.tsv"

    echo "  候选整合位点(支持读段≥3)："
    column -t "$OUTDIR/${TG}.candidate_integration_sites.tsv" | head -15
    echo "  构建体 contig 覆盖概览:"
    RUN samtools coverage -r "$TG" "$CRAM" --reference "$HYBRID" | sed 's/^/    /'
done

# (B) Manta BND：提取涉及 TG_ contig 的断点
MANTA=$(ls "$PROJ"/output_results/variant_calling/manta/"$SAMPLE"/*.diploid_sv.vcf.gz 2>/dev/null | head -1)
if [ -n "$MANTA" ]; then
    echo "===== Manta BND 涉及构建体 ====="
    RUN bcftools view "$MANTA" | grep -vE '^##' | awk -F'\t' '$5 ~ /TG_/ || $1 ~ /^TG_/ {print $1"\t"$2"\t"$5"\t"$8}' \
      | tee "$OUTDIR/manta_BND_involving_constructs.tsv" | head -20
else
    echo "（未找到 Manta VCF，跳过 B 路）"
fi

echo "DONE → $OUTDIR"
