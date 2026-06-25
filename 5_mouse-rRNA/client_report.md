# RNA-seq Library Quality Assessment Report

**Project**: Mouse RNA-seq Library QC Evaluation  
**Samples**: 9 mouse samples (mouse_28, 29, 32, 41, 42, 45, 46, 47, 48)  
**Date**: 2026-06-25  
**Prepared by**: Zhen Gao, PhD, Principle Bioinformatics Scientist, Athenomics

---

## Executive Summary

We performed a comprehensive quality assessment of 9 mouse RNA-seq libraries to determine their suitability for downstream gene expression analysis. The evaluation addressed three key questions: (1) insert size distribution and adapter dimer content, (2) ribosomal RNA (rRNA) contamination levels, and (3) uniquely mappable reads and gene detection capacity.

**Overall assessment: All 9 libraries fail to meet standard quality thresholds across all three metrics.** The observed defects are systematic and consistent across the batch, indicating a library preparation issue rather than sample-specific degradation. We recommend investigating the root cause of library preparation before proceeding.

---

## Methods Overview

Raw sequencing data (~32 million read pairs per sample, PE150) were subsampled to 10 million read pairs per sample for QC analysis — a depth sufficient to provide stable estimates of all three QC metrics. Reads were processed through the following pipeline:

1. **Adapter trimming**: Trim Galore (minimum length cutoff: 20 bp)
2. **rRNA identification**: SortMeRNA v4.3.6 against 8 reference databases (SILVA 16S/23S bacterial/archaeal, SILVA 18S/28S eukaryotic, Rfam 5S/5.8S)
3. **Genome alignment**: STAR aligner against the mouse reference genome (GRCm39, GENCODE vM35 annotation)
4. **Gene quantification**: Salmon + featureCounts

---

## Results

### 1. Insert Size and Adapter Dimer Content

Adapter dimer contamination is assessed by measuring the fraction of read pairs discarded after trimming due to insert size falling below the minimum usable threshold (20 bp). In a high-quality library, this fraction is typically **< 5–10%**.

| Sample | Read pairs discarded after trimming (%) | Median trimmed read length (bp) |
|--------|----------------------------------------|--------------------------------|
| mouse_28 | **69.7%** | 82 bp |
| mouse_29 | **72.0%** | 147 bp |
| mouse_32 | 26.1% | 42 bp |
| mouse_41 | 37.4% | 52 bp |
| mouse_42 | 22.4% | 42 bp |
| mouse_45 | 23.4% | 42 bp |
| mouse_46 | 29.1% | 147 bp |
| mouse_47 | 42.3% | 42 bp |
| mouse_48 | 22.9% | 42 bp |

**Findings:**

- **mouse_28 and mouse_29** show extreme adapter dimer contamination (~70% of read pairs discarded). The majority of sequenced fragments in these libraries consist of adapter sequences ligated directly to each other, with little to no cDNA insert. These two libraries are essentially non-functional for RNA-seq.
- **The remaining 7 samples** show 22–42% dimer-affected pairs — significantly elevated above acceptable levels, though less severe. The median trimmed read length for mouse_32/42/45/47/48 is only 42 bp, indicating the predominant insert size is well below the PE150 design length.

---

### 2. Ribosomal RNA (rRNA) Contamination

The fraction of reads matching rRNA reference databases (measured before genome alignment, on trimmed reads) reflects the abundance of rRNA in the original sample material. Standard mRNA-seq libraries typically show **< 5–10% rRNA**.

| Sample | Total rRNA (%) | Main rRNA source | Eukaryotic 18S (%) | Eukaryotic 28S (%) | Bacterial 16S (%) |
|--------|---------------|-----------------|---------------------|---------------------|-------------------|
| mouse_28 | 44.6 | Bacterial 16S | 4.1 | 7.9 | **30.4** |
| mouse_29 | 45.1 | Bacterial 16S | 3.6 | 2.7 | **37.7** |
| mouse_32 | 47.2 | Eukaryotic 28S | 7.8 | **31.6** | 3.4 |
| mouse_41 | 37.9 | Eukaryotic 28S | 6.5 | **15.5** | 11.8 |
| mouse_42 | 43.2 | Eukaryotic 28S | 8.7 | **26.7** | 2.8 |
| mouse_45 | 40.4 | Eukaryotic 28S | 11.4 | **22.1** | 1.6 |
| mouse_46 | 27.5 | Bacterial 16S | 10.5 | 0.7 | **15.6** |
| mouse_47 | 44.8 | Eukaryotic 28S | 10.2 | **24.8** | 2.9 |
| mouse_48 | 45.3 | Eukaryotic 28S | 10.4 | **28.0** | 1.9 |

**Findings:**

- **All 9 samples show rRNA contamination in the range of 27.5–47.2%**, far exceeding acceptable levels for mRNA-seq.
- **Two distinct rRNA profiles** are observed across the batch:
  - *Eukaryotic-dominant pattern* (mouse_32/42/45/47/48): rRNA signal is primarily driven by eukaryotic 28S, consistent with incomplete rRNA depletion in standard mammalian tissue libraries.
  - *Bacterial 16S-dominant pattern* (mouse_28/29/46): rRNA signal is dominated by the bacterial 16S database (30–38%). This pattern warrants further investigation — it may reflect true bacterial contamination of the sample material (e.g., from gut-associated tissues), or may include cross-mapping artefacts due to the high taxonomic diversity of bacterial 16S databases. Kraken2 taxonomic classification of unmapped reads from mouse_29 and mouse_46 was performed to help distinguish these possibilities (results in supplementary).

> **Note on rRNA metrics**: The rRNA percentages reported here (from SortMeRNA, pre-alignment) represent the fraction of all sequenced reads that match rRNA sequences. This is the appropriate metric for evaluating total rRNA contamination. A lower "residual rRNA%" sometimes shown in alignment-based QC reports (e.g., MultiQC featureCounts biotype column) reflects only rRNA reads that escaped SortMeRNA removal and still successfully aligned to the genome — a different and less representative statistic.

---

### 3. Uniquely Mappable Reads and Gene Detection

The fraction of reads that uniquely align to the reference genome and the number of genes detected reflect the functional information content of each library. For this analysis, reads were subsampled to 10 million pairs per sample; at this depth, a well-prepared mouse RNA-seq library is expected to achieve a unique mapping rate of **≥ 60–65%** (full-depth runs of ~30 million pairs typically reach 70–90%).

| Sample | Uniquely mapped reads (%) | Reads unmapped — too short (%) | Genes detected |
|--------|--------------------------|-------------------------------|---------------|
| mouse_28 | 19.6% | 66.0% | 10,312 |
| mouse_29 | 12.2% | 79.6% | 6,157 |
| mouse_32 | 28.7% | 51.5% | 19,000 |
| mouse_41 | 22.6% | 59.5% | 16,502 |
| mouse_42 | 26.6% | 54.6% | 20,817 |
| mouse_45 | 24.9% | 54.8% | 20,270 |
| mouse_46 | **2.2%** | **96.2%** | 1,042 |
| mouse_47 | 26.0% | 53.9% | 18,590 |
| mouse_48 | 27.7% | 52.3% | 19,826 |

**Findings:**

- **Unique mapping rates are severely reduced across all samples** (2–29% vs. the expected ≥ 60–65% at this sequencing depth). The predominant failure mode, accounting for 51–96% of unaligned reads across all samples, is classified as "too short" by STAR — reads that are too short after trimming to achieve a reliable alignment.
- This failure mode is a direct downstream consequence of the insert size issue described in Section 1: short inserts produce short trimmed reads that cannot be reliably placed in the genome.
- **mouse_46 is effectively a failed library** with only 2.2% unique mapping and 1,042 genes detected — values consistent with near-complete absence of usable mRNA signal.
- **mouse_28 and mouse_29**, despite having longer remaining reads post-trimming (82 bp and 147 bp respectively), show the lowest gene detection (10,312 and 6,157 genes) because the majority of their read pairs were discarded at the trimming step (~70%), leaving very few reads to align.
- The remaining 6 samples (mouse_32/41/42/45/47/48) detect 16,000–21,000 genes, but this is achieved from only 22–29% of input reads — meaning roughly three-quarters of the sequencing data is unusable, and effective sequencing depth for expression quantification is far below the nominal 10 million pairs.

**Strand specificity**: RSeQC strand inference returned "undetermined" for all 8 analyzable samples (forward: 33–40%, reverse: 54–61%), failing to reach the ≥ 80% threshold required for confident strand assignment. This is most likely a consequence of the fragmented read length rather than an incorrect library preparation protocol — reads that are too short carry insufficient strand-specific signal to allow reliable inference.

---

## Summary Table

| Sample | rRNA (%) | Dimer-affected pairs (%) | Unique mapping (%) | Genes detected | Overall |
|--------|---------|--------------------------|-------------------|----------------|---------|
| mouse_28 | 44.6 | 69.7 | 19.6 | 10,312 | ❌ Fail |
| mouse_29 | 45.1 | 72.0 | 12.2 | 6,157 | ❌ Fail |
| mouse_32 | 47.2 | 26.1 | 28.7 | 19,000 | ❌ Fail |
| mouse_41 | 37.9 | 37.4 | 22.6 | 16,502 | ❌ Fail |
| mouse_42 | 43.2 | 22.4 | 26.6 | 20,817 | ❌ Fail |
| mouse_45 | 40.4 | 23.4 | 24.9 | 20,270 | ❌ Fail |
| mouse_46 | 27.5 | 29.1 | 2.2 | 1,042 | ❌ Fail (废库) |
| mouse_47 | 44.8 | 42.3 | 26.0 | 18,590 | ❌ Fail |
| mouse_48 | 45.3 | 22.9 | 27.7 | 19,826 | ❌ Fail |

**Reference thresholds (healthy mRNA-seq, 10M read pairs)**: rRNA < 10% · Dimer pairs < 10% · Unique mapping ≥ 60–65%

---

## Conclusions and Recommendations

All 9 libraries fail quality standards across all three assessed metrics. The pattern of failure is systematic (affecting all samples simultaneously) rather than sample-specific, pointing to a library preparation protocol issue rather than sample quality.

The core defects are:
1. **High adapter dimer content** — the dominant problem, causing read loss at trimming and cascading into low mapping rates
2. **High rRNA contamination** (27–47%) — indicating insufficient ribosomal RNA depletion
3. **Low unique mapping rates** (2–29%) — a direct consequence of (1)

**We do not recommend simply re-preparing libraries with the same protocol**, as this is likely to reproduce the same defects. Instead, we suggest a root-cause investigation focusing on:

- RNA integrity (RIN values for each sample, if available) — degraded RNA tends to have very short inserts
- Library preparation protocol review — particularly the RNA fragmentation step, rRNA depletion efficiency, and adapter ligation conditions
- Consideration of an alternative rRNA depletion approach (e.g., ribodepletion kit instead of or in addition to the current method)

We are available to discuss these findings in detail and to advise on an optimized re-library strategy.

---

*Report prepared by:*

**Zhen Gao, PhD**  
Principle Bioinformatics Scientist  
Athenomics
