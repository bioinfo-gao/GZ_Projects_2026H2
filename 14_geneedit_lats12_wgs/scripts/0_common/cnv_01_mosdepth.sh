#!/bin/bash
# ============================================================================
# CNV step 1 — mosdepth 500kb 分箱深度（全 12 样，GRCm39 参考）
#   拷贝数/非整倍体/倍性天然免对照：基因组内 coverage 比值即拷贝数谱（方案 §3.3）。
#   ⚠ CRAM 对齐的是纯 GRCm39（61 contigs），故 -f 必须用 GRCm39，不是混合参考。
#   输出：analysis_A/cnv/<s>/  和 analysis_B/cnv_ploidy/<s>/  的 *.regions.bed.gz
#   下游由 cnv_02_profile.py 归一化 → 拷贝数谱 + 非整倍体判定 + 图。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
GRCM39="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
RUN(){ conda run -n regular_bioinfo "$@"; }
BIN=500000            # 500kb 窗
THREADS=4             # 每样 mosdepth 4 线程；下方 3 样并发 = 12 线程，安全
export PROJ GRCM39 BIN THREADS

run_one(){
  local study="$1" sample="$2" cram="$3"
  local outdir
  if [ "$study" = A ]; then outdir="$PROJ/analysis_A/cnv/$sample"; else outdir="$PROJ/analysis_B/cnv_ploidy/$sample"; fi
  mkdir -p "$outdir"; cd "$outdir"
  echo "[$(date +%H:%M:%S)] mosdepth $study/$sample"
  conda run -n regular_bioinfo mosdepth -t "$THREADS" -n --fast-mode --by "$BIN" -f "$GRCM39" "$sample" "$cram" \
    && echo "[$(date +%H:%M:%S)] DONE $sample" || echo "[$(date +%H:%M:%S)] FAIL $sample"
}
export -f run_one

# 样本清单：study sample
LIST=$(cat <<EOF
A RO_origin
A RO_B1TP
A RO_B2TP
A RO_tumor1
A RO_tumor2
A RO_tumor3
B L1L2_3M
B L1L2H_3M
B L1L2_12M
B L1L2H_12M
B L1L2_18M
B L1L2H_18M
EOF
)

# 3 样并发（3×4=12 线程）
echo "$LIST" | while read study sample; do
  [ -z "$study" ] && continue
  if [ "$study" = A ]; then cram="$PROJ/output_A/preprocessing/markduplicates/$sample/$sample.md.cram";
  else cram="$PROJ/output_B/preprocessing/markduplicates/$sample/$sample.md.cram"; fi
  echo "$study|$sample|$cram"
done | xargs -P 3 -I{} bash -c 'IFS="|" read s n c <<< "{}"; run_one "$s" "$n" "$c"'

echo "ALL_MOSDEPTH_DONE $(date +%H:%M:%S)"
