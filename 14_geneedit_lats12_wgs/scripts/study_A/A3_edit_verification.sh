#!/bin/bash
# ============================================================================
# Study A / Step 3 — 编辑验证：Brca1/Brca2/Pten CRISPR KO 是否发生
#   在 sgRNA 切点比 编辑细胞/肿瘤 vs RO_origin，查 indel/frameshift、等位比例。
#   sgRNA 序列见 ../../docs/client_materials/sgRNA_guides.md（Pten×3, Brca1×3；Brca2 待补）
#   前置：A2 sarek 完成（CRAM 在 output_A/preprocessing/markduplicates/）。
#
#   流程：
#   (1) 在 GRCm39 定位每条 sgRNA spacer（+反向互补，配 NGG PAM）→ 预测切点±20bp 窗口
#   (2) 对每个窗口，各样本 samtools + bcftools 看 indel；比 vs origin
#   (3) 用 tumor 里 Brca1 vs Brca2 被 KO 反推谱系
#   本脚本先做 (1) 定位并列出窗口；(2)(3) 待 A2 输出后填坐标运行。
# ============================================================================
set -euo pipefail
PROJ="/home/gao/projects_2026H2/14_geneedit_lats12_wgs"
GRCM39="/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa"
GUIDES="$PROJ/docs/client_materials/sgRNA_guides.md"
OUT="$PROJ/analysis_A/edit_verification"; mkdir -p "$OUT"
RUN(){ conda run -n regular_bioinfo "$@"; }

echo ">> (1) 在 GRCm39 定位 sgRNA 切点（seqkit locate spacer + 反向互补）"
# 从 guides.md 抽出 DNA spacer（反引号内 20-mer 大写），写成 fasta
grep -oE '`[ACGT]{20}`' "$GUIDES" | tr -d '`' | awk '{print ">g"NR"\n"$0}' > "$OUT/spacers.fa"
echo "   spacers:"; cat "$OUT/spacers.fa"
RUN seqkit locate -d -f "$OUT/spacers.fa" "$GRCM39" > "$OUT/spacer_hits.tsv" 2>/dev/null || true
echo "   命中（染色体/位置/链）→ $OUT/spacer_hits.tsv"; column -t "$OUT/spacer_hits.tsv" 2>/dev/null | head

cat <<'NOTE'
>> (2)(3) 待 A2 完成后：对每个切点窗口比 tumor/edited vs RO_origin
   for BAM in output_A/.../{RO_B1TP,RO_B2TP,RO_tumor1..3,RO_origin}.md.cram; do
     samtools mpileup / bcftools 在切点±20bp 看 indel 频率
   done
   判读：编辑细胞/肿瘤在靶点有 indel 而 origin 无 → KO 确认；
        tumor 有 Brca1-indel→源自 B1TP，有 Brca2-indel→源自 B2TP。
   关键位点 IGV 目视。
NOTE
echo "DONE A3 (1)；(2)(3) 待 sarek 输出"
