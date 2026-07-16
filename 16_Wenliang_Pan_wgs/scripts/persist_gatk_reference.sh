#!/usr/bin/env bash
# Persist the GATK.GRCh38 (iGenomes analysis-set) reference that sarek downloaded into proj16's
# transient work/ staging, into a PERMANENT shared location so future human WGS reuses it (no re-download).
# COPY only (never move) — the running sarek symlinks to the staging originals.
# Run in tmux: tmux new-session -d -s pan_ref "bash scripts/persist_gatk_reference.sh"
set -euo pipefail

PROJ=/home/gao/projects_2026H2/16_Wenliang_Pan_wgs
DEST=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38     # maps 1:1 to sarek --genome GATK.GRCh38
S=$(find "$PROJ/work" -maxdepth 1 -type d -name "stage-*" | head -1)
[ -z "$S" ] && { echo "no staging dir found"; exit 1; }

FA=$(find "$S" -name "Homo_sapiens_assembly38.fasta" | head -1)
FAI=$(find "$S" -name "Homo_sapiens_assembly38.fasta.fai" | head -1)
DICT=$(find "$S" -name "Homo_sapiens_assembly38.dict" | head -1)
IDX=$(find "$S" -type d -name "BWAmem2Index" | head -1)
BED=$(find "$S" -name "wgs_calling_regions_noseconds.hg38.bed" | head -1)

mkdir -p "$DEST"
echo "=== copying to $DEST (cp -a, preserve) ==="
cp -a "$FA"   "$DEST/"
cp -a "$FAI"  "$DEST/"
cp -a "$DICT" "$DEST/"
cp -a "$BED"  "$DEST/"
cp -a "$IDX"  "$DEST/"          # -> $DEST/BWAmem2Index/

# integrity: sizes must match source
echo "=== verify byte sizes ==="
ok=1
for f in Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict; do
    src=$(find "$S" -name "$f" | head -1)
    [ "$(stat -c%s "$src")" = "$(stat -c%s "$DEST/$f")" ] && echo "OK  $f" || { echo "MISMATCH $f"; ok=0; }
done
for f in Homo_sapiens_assembly38.fasta.0123 Homo_sapiens_assembly38.fasta.bwt.2bit.64; do
    [ "$(stat -c%s "$IDX/$f")" = "$(stat -c%s "$DEST/BWAmem2Index/$f")" ] && echo "OK  $f" || { echo "MISMATCH $f"; ok=0; }
done
[ "$ok" = 1 ] && echo "=== ALL SIZES MATCH ===" || { echo "=== VERIFY FAILED ==="; exit 1; }

cat > "$DEST/README.md" <<EOF
# GATK.GRCh38 — GRCh38 full analysis set (AWS iGenomes)

This is the reference that nf-core/sarek's \`--genome GATK.GRCh38\` resolves to
(\`Homo_sapiens_assembly38\`, the Broad/GATK GRCh38 **full analysis set** with ALT + decoy + HLA contigs).

- **Provenance**: downloaded from AWS iGenomes S3 by project 16 (Wenliang Pan human WGS) and persisted here
  on 2026-07-16 so future human germline WGS reuses it instead of re-downloading (~21 GB) each run.
- **Contents**:
  - \`Homo_sapiens_assembly38.fasta\` (+ \`.fai\`, \`.dict\`)
  - \`BWAmem2Index/\` — pre-built ALT-aware bwa-mem2 index (\`.0123 .bwt.2bit.64 .amb .ann .pac .alt\`)
  - \`wgs_calling_regions_noseconds.hg38.bed\` — GATK WGS calling intervals
- **Known-sites (dbSNP / Mills)** are NOT duplicated here — they live in the sibling
  \`../gatk_bundle/\` (\`dbsnp_146.hg38.vcf.gz\`, \`Mills_and_1000G_gold_standard.indels.hg38.vcf.gz\`).
- **Reuse**: point sarek at these explicitly, e.g.
  \`--fasta $DEST/Homo_sapiens_assembly38.fasta --bwa $DEST/BWAmem2Index\`
  (or configure a local igenomes_base). Contig naming = analysis set, matches gnomAD/ClinVar/dbSNP.
- ⚠ Do NOT rename contigs or mix with the GENCODE \`primary_assembly\` fasta in \`../human_gencode_v45/\`
  (different contig composition; that one is for RNA-seq/annotation, not variant calling).
EOF

echo "=== DONE. contents: ==="
ls -laR "$DEST" | grep -vE "^total" | head -40
du -sh "$DEST"
