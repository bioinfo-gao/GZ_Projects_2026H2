#!/bin/bash
# ============================================================================
# Study A / Step 4 — 拷贝数/倍性（配对，tumor vs RO_origin）
#   Control-FREEC 配对模式：以 origin 为 control，输出各 tumor/编辑细胞的 CNV 谱、
#   染色体臂级增删、倍性。HRD(Brca1/2 缺失) → 预期拷贝数不稳定。
#   前置：A2 完成；0a 装好 control-freec。CRAM→BAM 若 FREEC 需要则先转。
#   本脚本为模板：per-sample 生成 FREEC config 再跑；参数在审阅时按 output_A 实际路径确认。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
HYBRID="$PROJ/refs/hybrid/GRCm39_plus_Cas9_iHPV.fa"
OUT="$PROJ/analysis_A/cnv"; mkdir -p "$OUT"
NORMAL_ID="RO_origin"
TUMORS="RO_B1TP RO_B2TP RO_tumor1 RO_tumor2 RO_tumor3"

cram(){ ls "$PROJ"/output_A/preprocessing/markduplicates/"$1"/*.cram 2>/dev/null | head -1; }
NORMAL_CRAM=$(cram "$NORMAL_ID")

cat <<NOTE
Control-FREEC 配对模式（每个 tumor 一个 config）：
  [general] ploidy=2; window=50000; chrLenFile=<GRCm39.fai>; outputDir=$OUT/<tumor>
  [sample] mateFile=<tumor CRAM/BAM>; inputFormat=BAM (CRAM 需先转 或 用 --reference)
  [control] mateFile=$NORMAL_CRAM
  → FREEC 输出 *_CNVs, *_ratio.txt；R 脚本画 CNV 谱 + 倍性。
备选：cnvkit.py batch <tumors> --normal $NORMAL_ID --fasta $HYBRID （flat/paired）交叉验证。
待 A2 输出后按实际 CRAM 路径填 config 运行；mosdepth 分箱可先出全基因组深度谱预览。
NOTE
echo "NORMAL CRAM: ${NORMAL_CRAM:-未生成}"; echo "tumors: $TUMORS"
echo "DONE A4（模板；待 sarek 输出后填坐标/路径运行）"
