#!/bin/bash
# ============================================================================
# Study B / Step 3 — iHPV 构建体 存在筛查 + 候选整合位点（首要假设：插入突变）
#   sarek 对齐的是纯 GRCm39 → 构建体来源读段为 unmapped。做法：
#   (1) 一次 CRAM 扫描抽出 "unmapped 或 mate-unmapped" 读段 → 小 BAM（省 I/O）
#   (2) unmapped 读段 → 比对 iHPV 标记参考(HPV16+EGFP+Luc) → 构建体读段计数
#       → L1L2H(带iHPV) 应阳性、L1L2 应阴性（内部特异性对照）
#   (3) 构建体读段的"基因组锚定"(mapped, mate-unmapped) 落点聚类 = 候选整合位点
#   ⚠ 局限：TG_iHPV 为标记序列(HPV16/EGFP/Luc)，非 PMC4662542 完整载体图 →
#     可定"落在哪个基因组区"，碱基级结合部/被打断基因精判需完整载体（报告注明）。
#   前置：B2 完成；refs/constructs/TG_iHPV.fa 已 bwa-mem2 index。
#   ⚡ 6 样 3 路并行(xargs -P3)：I/O 密集，负载低时充分用盘。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
GRCM39="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
TG="$PROJ/refs/constructs/TG_iHPV.fa"
OUT="$PROJ/analysis_B/ihpv_integration"; mkdir -p "$OUT"
export PROJ GRCM39 TG OUT

run_one(){
  local s="$1"
  local cram="$PROJ/output_B/preprocessing/markduplicates/$s/$s.md.cram"
  local sdir="$OUT/$s"; mkdir -p "$sdir"
  echo "[$(date +%H:%M:%S)] >> $s"
  conda run -n regular_bioinfo bash -c "
    # (1) 一次扫描：unmapped 或 mate-unmapped → 小 BAM
    samtools view -@ 4 --reference '$GRCM39' -b -e 'flag.unmap || flag.munmap' '$cram' 2>/dev/null > '$sdir/sub.bam'
    # (2) unmapped 读段 → 构建体
    samtools fastq -f 4 '$sdir/sub.bam' 2>/dev/null > '$sdir/unmapped.fq'
    bwa-mem2 mem -t 6 -v 1 '$TG' '$sdir/unmapped.fq' 2>/dev/null \
      | samtools view -b -F 4 -q 20 -o '$sdir/construct_hits.bam' - 2>/dev/null
    samtools view '$sdir/construct_hits.bam' | cut -f1 | sort -u > '$sdir/hit_names.txt'
    hpv=\$(samtools view '$sdir/construct_hits.bam' | awk '\$3==\"TG_HPV16\"' | wc -l)
    egfp=\$(samtools view '$sdir/construct_hits.bam' | awk '\$3==\"TG_EGFP\"' | wc -l)
    luc=\$(samtools view '$sdir/construct_hits.bam' | awk '\$3==\"TG_Luc\"' | wc -l)
    tot=\$(wc -l < '$sdir/hit_names.txt')
    printf '%s\t%s\t%s\t%s\t%s\n' '$s' \"\$tot\" \"\$hpv\" \"\$egfp\" \"\$luc\" > '$sdir/presence.tsv'
    # (3) 基因组锚定：mapped 且 mate-unmapped，名字∈构建体读段 → 落点聚类
    samtools view -f 8 -F 4 '$sdir/sub.bam' 2>/dev/null \
      | awk 'NR==FNR{h[\$1];next}(\$1 in h){print \$3\"\t\"\$4}' '$sdir/hit_names.txt' - \
      | sort -k1,1 -k2,2n > '$sdir/anchors.tsv'
    awk '{bin=int(\$2/5000); key=\$1\"\t\"bin*5000; c[key]++} END{for(k in c) print k\"\t\"c[k]}' '$sdir/anchors.tsv' \
      | sort -k3,3nr | head -15 > '$sdir/candidate_loci.tsv'
  "
  echo "[$(date +%H:%M:%S)] DONE $s"
}
export -f run_one

printf '%s\n' L1L2_3M L1L2H_3M L1L2_12M L1L2H_12M L1L2_18M L1L2H_18M \
  | xargs -P 3 -I{} bash -c 'run_one "$@"' _ {}

# 汇总
echo -e "sample\tconstruct_reads\tHPV16\tEGFP\tLuc" > "$OUT/construct_presence.tsv"
for s in L1L2_3M L1L2H_3M L1L2_12M L1L2H_12M L1L2_18M L1L2H_18M; do
  cat "$OUT/$s/presence.tsv" 2>/dev/null >> "$OUT/construct_presence.tsv"
done
echo "== 构建体存在(L1L2H应阳/L1L2应阴) =="; column -t "$OUT/construct_presence.tsv"
echo "ALL_B3_DONE $(date +%H:%M:%S)"
