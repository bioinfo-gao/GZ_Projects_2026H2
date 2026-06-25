#!/bin/bash
# =============================================================================
# 5_mouse-rRNA / SortMeRNA_issue — Step 1: subsample raw reads to 10M pairs
#
# Background (for report):
#   Full-depth nf-core/rnaseq run (32M/16M read pairs per sample x 8 rRNA
#   reference databases) caused SortMeRNA 4.3.6 to either exceed the 16h
#   process time limit (exit 143: mouse_32, mouse_41, mouse_48) or get
#   OOM-killed by the kernel (exit 137: mouse_45 at anon-rss ~74GB,
#   mouse_47 at anon-rss ~118GB; confirmed via `journalctl -k`). Only
#   mouse_28 and mouse_29 completed cleanly.
#
#   research_aim.md states the actual goal is QC triage on 9 samples
#   (insert size/dimer content, rRNA contamination %, unique-reads-to-genes
#   mapping) to decide whether libraries are worth re-prepping — NOT a
#   final quantitative dataset. A 10M read-pair subsample is more than
#   sufficient for stable %rRNA and mapping-rate estimates, and is small
#   enough to avoid the memory blow-up entirely while keeping all 8
#   reference databases (mouse_28's full-depth log showed silva-bac-16s
#   matching 30.31% of reads — the dominant rRNA signal — so databases were
#   NOT dropped, to preserve that finding).
#
# Tool: seqtk sample (deterministic with fixed seed -> R1/R2 stay paired)
# Seed: 100 (fixed, for reproducibility)
# Target: 10,000,000 read pairs per sample
# =============================================================================
set -euo pipefail

SEED=100
N_PAIRS=10000000
RAW_SAMPLESHEET=/home/gao/projects_2026H2/5_mouse-rRNA/scripts/nf_core_samplesheet.csv
OUTDIR=/home/gao/projects_2026H2/5_mouse-rRNA/SortMeRNA_issue/subsampled_10M
NEW_SAMPLESHEET=/home/gao/projects_2026H2/5_mouse-rRNA/SortMeRNA_issue/subsampled_10M_samplesheet.csv

mkdir -p "$OUTDIR"

echo "sample,fastq_1,fastq_2,strandedness" > "$NEW_SAMPLESHEET"

tail -n +2 "$RAW_SAMPLESHEET" | while IFS=',' read -r sample fq1 fq2 strandedness; do
    [ -z "$sample" ] && continue
    out1="$OUTDIR/${sample}_sub10M_1.fq.gz"
    out2="$OUTDIR/${sample}_sub10M_2.fq.gz"

    echo "=== Subsampling $sample to ${N_PAIRS} read pairs (seed=$SEED) ==="
    echo "    R1 src: $fq1"
    echo "    R2 src: $fq2"

    seqtk sample -s${SEED} "$fq1" ${N_PAIRS} | gzip -1 > "$out1" &
    seqtk sample -s${SEED} "$fq2" ${N_PAIRS} | gzip -1 > "$out2" &
    wait

    echo "${sample},${out1},${out2},${strandedness}" >> "$NEW_SAMPLESHEET"
    echo "    done -> $out1 / $out2"
done

echo "=== New samplesheet written: $NEW_SAMPLESHEET ==="
cat "$NEW_SAMPLESHEET"
