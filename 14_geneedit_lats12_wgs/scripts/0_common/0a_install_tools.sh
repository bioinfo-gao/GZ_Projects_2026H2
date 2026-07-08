#!/bin/bash
# ============================================================================
# Step 0a — 补装 CNV 工具（其余 bwa-mem2/samtools/delly/mosdepth/bcftools 已在 regular_bioinfo）
#   Control-FREEC：tumor-only 与配对 CNV 均支持（Study A 配对、Study B tumor-only）
#   CNVkit：备选/交叉验证
# ============================================================================
set -euo pipefail

echo ">> 确认已有工具"
for t in bwa-mem2 samtools bcftools delly mosdepth; do
  p=$(conda run -n regular_bioinfo bash -lc "command -v $t" 2>/dev/null)
  printf "  %-12s %s\n" "$t" "${p:-MISSING}"
done

echo ">> 补装 Control-FREEC + CNVkit（若已在则跳过）"
if ! conda run -n regular_bioinfo bash -lc "command -v freec" >/dev/null 2>&1; then
  mamba install -n regular_bioinfo -c bioconda -c conda-forge control-freec cnvkit -y
else
  echo "  control-freec 已存在"
fi

echo ">> 验证"
conda run -n regular_bioinfo bash -lc "freec --version 2>&1 | head -1; cnvkit.py version"
echo "DONE 0a"
