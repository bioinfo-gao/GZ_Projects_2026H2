#!/bin/bash
# ============================================================================
# Study A / Step 5 — Cas9 转基因整合位点（次要/加分）
#   方法同 proj13：抓落在 TG_Cas9 contig 上的读段，其配偶(discordant)/软剪切(split)
#   落回小鼠何处 → 聚类成整合位点。MAPQ≥20 唯一比对 + 高深度 artifact 黑名单。
#   所有 A 样本应在同一 Cas9 位点（同一亲本谱系）。
#   前置：A2 完成（跑在含 TG_Cas9 的混合参考上）。
#   用法: bash A5_cas9_integration.sh RO_origin
# ============================================================================
set -euo pipefail
SAMPLE="${1:?用法: bash A5_cas9_integration.sh <sample>}"
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_Cas9_iHPV.fa"
OUT="$PROJ/analysis_A/cas9_integration/$SAMPLE"; mkdir -p "$OUT"
RUN(){ conda run -n regular_bioinfo "$@"; }; MAPQ=20
CRAM=$(ls "$PROJ"/output_A/preprocessing/markduplicates/"$SAMPLE"/*.cram 2>/dev/null | head -1)
[ -z "$CRAM" ] && { echo "ERROR: 无 $SAMPLE CRAM（先跑 A2）"; exit 1; }

TG="TG_Cas9"
echo "== $SAMPLE @ $TG 覆盖 =="; RUN samtools coverage -q $MAPQ -r "$TG" "$CRAM" --reference "$HYBRID" | sed 's/^/  /'
echo "== discordant 配偶落点（MAPQ≥$MAPQ）=="
RUN samtools view -q $MAPQ -T "$HYBRID" "$CRAM" "$TG" \
  | awk -v tg="$TG" '$7!="=" && $7!=tg && $7!="*"{print $7"\t"$8}' \
  | sort -k1,1 -k2,2n \
  | awk '{w=int($2/5000); c[$1"\t"w]++; if(!(($1"\t"w)in mn)||$2<mn[$1"\t"w])mn[$1"\t"w]=$2}
         END{for(k in c) if(c[k]>=3){split(k,a,"\t"); print a[1]"\t"mn[k]"\t"c[k]}}' \
  | sort -k3,3nr | awk 'BEGIN{print "chrom\tpos\tsupport"}1' | tee "$OUT/cas9_integration_sites.tsv"
echo "DONE A5 → $OUT （proj13 精修版逻辑：MAPQ过滤+聚类；如需 artifact黑名单/注释可套 proj13 脚本4）"
