#!/bin/bash
# ============================================================================
# Study B / Step 6 — de novo 候选：扣 C57BL/6 亚系背景 → 高影响/复现变异
#   germline VCF(样本 vs GRCm39) 用 Sanger MGP + 小鼠 dbSNP 扣掉背景 germline →
#   剩余高影响(frameshift/stop/剪接) + 多样本或随龄复现的变异 = 候选。
#   L1L2 vs L1L2H 差异集单列。关联输卵管发育/纤毛/上皮、DNA修复、基因组稳定性基因。
#   前置：B2 完成（germline VCF）；0d 下好 MGP/dbSNP。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
OUT="$PROJ/analysis_B/candidates"; mkdir -p "$OUT"
RUN(){ conda run -n regular_bioinfo "$@"; }
MGP="/Work_bio/references/Mus_musculus/GRCm39/mgp_dbsnp"   # 0d 下载目标

cat <<NOTE
流程（待 B2 germline VCF + 0d 背景库就绪后按实际路径运行）：
  1) 归一化：bcftools norm -f GRCm39 每样 VCF
  2) 扣背景：bcftools isec -C 样本.vcf MGP.vcf dbSNP.vcf → 剩余非背景变异
  3) 注释影响：SnpEff/VEP（小鼠 GRCm39）→ 取 HIGH/MODERATE
  4) 复现/差异：
       - 多样本或随龄(3M→12M→18M)复现 → 候选
       - L1L2 vs L1L2H 差异集（bcftools isec）→ iHPV 相关
  5) 关联：候选是否落在 输卵管发育/纤毛/上皮 · Hippo/YAP · DNA修复 · 基因组稳定性 基因
  输出：$OUT/candidates.tsv + 按类别分组
NOTE
echo "MGP/dbSNP 目录: $MGP  （$( [ -d "$MGP" ] && echo 已建 || echo 待下载 )）"
echo "DONE B6（框架；待 VCF+背景库就绪运行）"
