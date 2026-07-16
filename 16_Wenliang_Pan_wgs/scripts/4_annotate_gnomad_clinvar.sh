#!/usr/bin/env bash
# Project 16 — post-sarek annotation of germline SNV/indel calls with gnomAD AF + ClinVar.
# sarek already ran VEP (consequence / gene / SIFT / PolyPhen). Here we bolt on population
# frequency (gnomAD) and clinical significance (ClinVar) so downstream rarity+pathogenicity
# filtering (step 5) has everything in the INFO field.
# Run AFTER sarek completes. bash scripts/4_annotate_gnomad_clinvar.sh
set -euo pipefail

PROJ=/home/gao/projects_2026H2/16_Wenliang_Pan_wgs
DB=/Work_bio/references/Homo_sapiens/GRCh38/annotation
OUT=$PROJ/output_results/annotation_gnomad_clinvar
mkdir -p "$OUT"
CR="conda run -n regular_bioinfo"

GNOMAD=$DB/af-only-gnomad.hg38.vcf.gz
CLINVAR=$DB/clinvar_GRCh38.chr.vcf.gz

# Locate the per-sample VEP-annotated haplotypecaller VCFs sarek produced.
mapfile -t VCFS < <(find "$PROJ/output_results/annotation" -name "*haplotypecaller*vep*.vcf.gz" 2>/dev/null | sort)
if [ ${#VCFS[@]} -eq 0 ]; then
    echo "No VEP-annotated haplotypecaller VCFs found; falling back to raw calls."
    mapfile -t VCFS < <(find "$PROJ/output_results/variant_calling/haplotypecaller" -name "*.vcf.gz" ! -name "*.g.vcf.gz" 2>/dev/null | sort)
fi
echo "Annotating ${#VCFS[@]} VCF(s)."

for vcf in "${VCFS[@]}"; do
    s=$(basename "$vcf" | sed -E 's/\..*//')
    echo "--- $s : $vcf ---"
    # 1) gnomAD AF -> INFO/gnomAD_AF
    $CR bcftools annotate -a "$GNOMAD" -c "INFO/gnomAD_AF:=INFO/AF" \
        -h <(echo '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD allele frequency (af-only-gnomad.hg38)">') \
        "$vcf" -Oz -o "$OUT/${s}.gnomad.vcf.gz"
    $CR bcftools index -t "$OUT/${s}.gnomad.vcf.gz"
    # 2) ClinVar CLNSIG / CLNDN
    $CR bcftools annotate -a "$CLINVAR" -c "INFO/CLNSIG,INFO/CLNDN" \
        "$OUT/${s}.gnomad.vcf.gz" -Oz -o "$OUT/${s}.gnomad.clinvar.vcf.gz"
    $CR bcftools index -t "$OUT/${s}.gnomad.clinvar.vcf.gz"
    rm -f "$OUT/${s}.gnomad.vcf.gz" "$OUT/${s}.gnomad.vcf.gz.tbi"
done
echo "Annotated VCFs in $OUT"
ls -lh "$OUT"
