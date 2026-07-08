#!/bin/bash
# ============================================================================
# Step 0b — 取构建体序列 → ../refs/constructs/（供 0c 建混合参考、下游整合分析）
#   (1) SpCas9：Study A 亲本组成型表达（Rosa26-Cas9 等），定位其整合位点用
#   (2) iHPV 构建体：Study B L1L2H，CAG-loxP-EGFP-pA-loxP-E6E7-IRES-Luc
#       E6/E7 源 Addgene #13712 (pB-actin E6E7)；完整 LSL 载体见 PMC4662542
#
#   ⚠ 需人工/联网补齐（审阅时确认来源）：
#   - SpCas9 CDS：公开（如 Addgene lentiCas9 #52962 序列，或 pX330 的 SpCas9）
#   - Addgene #13712 序列：Addgene 页面 "View all sequences" 可下 GenBank
#   - PMC4662542 补充材料里的 iHPV 完整构建体图/序列
#   到手后放成：refs/constructs/TG_Cas9.fa 、 refs/constructs/TG_iHPV.fa
#   （FASTA，contig 名 >TG_Cas9 / >TG_iHPV；可先只放能拿到的最小可用序列）
# ============================================================================
set -euo pipefail
DEST="/home/gao/projects_2026H2/14_geneedit_lats12_wgs/refs/constructs"
mkdir -p "$DEST"

cat <<'NOTE'
本步骤为"取序列"占位说明，需确认来源后落盘：
  TG_Cas9.fa  ← SpCas9 CDS（公开，Addgene lentiCas9 #52962 / pX330）
  TG_iHPV.fa  ← iHPV 构建体（Addgene #13712 E6/E7 + PMC4662542 CAG-LSL-EGFP-Luc）
建议：优先拿全长 iHPV 载体（整合位点结合部检测需要骨架/接头处序列）；
      Cas9 若只为定位整合，SpCas9 CDS 即够。
落盘后 seqkit 校验：
  conda run -n regular_bioinfo seqkit stats refs/constructs/TG_*.fa
NOTE

echo ">> 现状："; ls -la "$DEST"/TG_*.fa 2>/dev/null || echo "  （尚无 TG_*.fa，待补）"
