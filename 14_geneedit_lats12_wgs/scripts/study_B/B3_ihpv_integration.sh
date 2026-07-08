#!/bin/bash
# ============================================================================
# Study B / Step 3 — iHPV 构建体整合位点（★首要假设：插入突变致输卵管表型）
#   方法同 proj13 嵌合读段法：抓 TG_iHPV contig 上读段的 discordant/split 落点 →
#   聚类整合位点 → 注释落在哪个小鼠基因（是否输卵管相关）。MAPQ≥20 + artifact 黑名单。
#   仅 L1L2H(3M/12M/18M) 有此构建体；L1L2 为阴性对照（应零覆盖）。
#   前置：B2 完成（跑在含 TG_iHPV 的混合参考上）。
#   用法: bash B3_ihpv_integration.sh L1L2H_12M
# ============================================================================
set -euo pipefail
SAMPLE="${1:?用法: bash B3_ihpv_integration.sh <sample>}"
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_Cas9_iHPV.fa"
GTF="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/gencode.vM35.annotation.gtf"
OUT="$PROJ/analysis_B/ihpv_integration/$SAMPLE"; mkdir -p "$OUT"
RUN(){ conda run -n regular_bioinfo "$@"; }; MAPQ=20
CRAM=$(ls "$PROJ"/output_B/preprocessing/markduplicates/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && { echo "ERROR: 无 $SAMPLE CRAM（先跑 B2）"; exit 1; }

TG="TG_iHPV"
echo "== $SAMPLE @ $TG 覆盖（L1L2 应≈0）=="; RUN samtools coverage -q $MAPQ -r "$TG" "$CRAM" --reference "$HYBRID" | sed 's/^/  /'
echo "== discordant 配偶落点 → 候选整合位点 =="
RUN samtools view -q $MAPQ -T "$HYBRID" "$CRAM" "$TG" \
  | awk -v tg="$TG" '$7!="=" && $7!=tg && $7!="*"{print $7"\t"$8}' \
  | sort -k1,1 -k2,2n \
  | awk '{w=int($2/5000); c[$1"\t"w]++; if(!(($1"\t"w)in mn)||$2<mn[$1"\t"w])mn[$1"\t"w]=$2}
         END{for(k in c) if(c[k]>=3){split(k,a,"\t"); print a[1]"\t"mn[k]"\t"c[k]}}' \
  | sort -k3,3nr | awk 'BEGIN{print "chrom\tpos\tsupport"}1' | tee "$OUT/ihpv_integration_sites.tsv"

echo ">> 注释：整合位点落在哪个小鼠基因（下一步用 GTF overlap；重点看输卵管/上皮/纤毛相关）"
echo "DONE B3 → $OUT （如需 artifact黑名单/on-off-target 注释，套 proj13 脚本4 完整版）"
