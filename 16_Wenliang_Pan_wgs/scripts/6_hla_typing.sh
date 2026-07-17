#!/usr/bin/env bash
# Project 16 — HLA typing from WGS ("if feasible" deliverable item).
# Tool: T1K (bioconda) — class I + II HLA genotyping from WGS, works from BAM/CRAM by
# extracting MHC-region + unmapped reads. Runs AFTER sarek produces per-sample CRAMs.
# One-time setup builds the IPD-IMGT/HLA reference. Run inside tmux.
set -euo pipefail

PROJ=/home/gao/projects_2026H2/16_Wenliang_Pan_wgs
HLAIDX=/Work_bio/references/Homo_sapiens/GRCh38/t1k_hlaidx
# MUST be the exact reference sarek built the CRAMs with (GATK.GRCh38 / Homo_sapiens_assembly38.fasta).
# Decoding CRAM with any other GRCh38 build (e.g. gencode primary_assembly) fails on @SQ M5
# checksum mismatch — the contig sequences/MD5s differ between builds. (2026-07-17 fix)
FASTA=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/Homo_sapiens_assembly38.fasta
OUT=$PROJ/output_results/hla_typing
mkdir -p "$OUT"

# --- 0) env (create only if absent; do NOT pollute regular_bioinfo) ---
if ! conda env list | grep -qE '^hla\s'; then
    mamba create -y -n hla -c bioconda -c conda-forge t1k samtools
fi
CR="conda run -n hla"

# --- 1) build HLA reference index once ---
if [ ! -s "$HLAIDX/hlaidx_dna_seq.fa" ]; then
    mkdir -p "$HLAIDX"; cd "$HLAIDX"
    T1KDIR=$($CR bash -lc 'dirname $(dirname $(which run-t1k))')/share
    T1KBUILD=$(find "$T1KDIR" -name t1k-build.pl 2>/dev/null | head -1)
    $CR perl "$T1KBUILD" -o "$HLAIDX" -d hla --download IPD-IMGT/HLA
fi

# --- 2) per sample: extract MHC + unmapped reads -> FASTQ -> T1K ---
# MHC region on GRCh38: chr6:28,510,120-33,480,577 (extended xMHC). Add unmapped for HLA reads
# that failed to map to the reference haplotype.
MHC="chr6:28000000-34000000"
for cram in $(find "$PROJ/output_results/preprocessing" -name "*.recal.cram" -o -name "*.md.cram" 2>/dev/null | sort); do
    s=$(basename "$cram" | sed -E 's/\..*//')
    echo "=== HLA $s : $cram ==="
    tmpb="$OUT/${s}.mhc.bam"
    $CR samtools view -b -T "$FASTA" "$cram" $MHC -o "$OUT/${s}.mhc_region.bam"
    $CR samtools view -b -f 4 -T "$FASTA" "$cram" -o "$OUT/${s}.unmapped.bam"
    $CR samtools merge -f "$tmpb" "$OUT/${s}.mhc_region.bam" "$OUT/${s}.unmapped.bam"
    $CR samtools collate -Ou "$tmpb" | $CR samtools fastq -1 "$OUT/${s}_R1.fq" -2 "$OUT/${s}_R2.fq" -0 /dev/null -s /dev/null -n
    $CR run-t1k -1 "$OUT/${s}_R1.fq" -2 "$OUT/${s}_R2.fq" --preset hla -f "$HLAIDX/hlaidx_dna_seq.fa" -t 8 -o "$OUT/${s}_hla" --od "$OUT"
    rm -f "$OUT/${s}.mhc_region.bam" "$OUT/${s}.unmapped.bam" "$tmpb" "$OUT/${s}_R1.fq" "$OUT/${s}_R2.fq"
done
echo "HLA results in $OUT ( *_hla_genotype.tsv )"
