# CHO WGS Species Confirmation & Strain Identification Report

**Project:** CHO Cell Line Genomic Identity Verification  
**Sample:** wt1 (23TF7FLT4_3_0469165296_wt1_S8_L003)  
**Date:** June 17, 2026  
**Analyst:** Zhen Gao, PhD | Principal Bioinformatics Scientist, Athenomics  

---

## 1. Executive Summary

Whole-genome sequencing (WGS) data from sample **wt1** was analyzed to confirm the species origin and identify the cell line strain. The results conclusively demonstrate that:

- **Species:** The sample is confirmed as **Chinese Hamster Ovary (CHO)** cells (*Cricetulus griseus*), with a 99.85% mapping rate to the CHO reference genome.
- **Strain:** The DHFR gene locus is **intact**, ruling out CHO-DG44 and CHO-DXB11. The sample is consistent with **CHO-K1 or CHO-S**.

---

## 2. Methods

### 2.1 Sequencing Data

| Parameter | Value |
|-----------|-------|
| Platform | Illumina |
| Library type | Paired-end WGS |
| Read length | 150 bp |
| Raw data size | R1: 17 GB, R2: 16 GB (compressed) |

### 2.2 Quality Control

Quality assessment was performed using **fastp v1.3.3** on a random subset of 5,000,000 read pairs (10,000,000 reads total). Quality metrics were aggregated using **MultiQC**. No additional filtering was applied for QC evaluation; default fastp quality filters were used to assess post-filtering statistics.

### 2.3 Species Confirmation

A subset of 5,000,000 read pairs was extracted using fastp and aligned to the **CHO-K1 reference genome** (CriGri-PICR, GCF_003668045.1) using **BWA v0.7.17** (Li & Durbin, 2009) with default parameters. Alignment statistics were generated using **SAMtools v1.21** (Danecek et al., 2021). Species identity was confirmed based on the overall mapping rate:

- \>85%: Confirmed CHO
- 50–85%: Ambiguous, requires further investigation
- <50%: Not CHO

### 2.4 Strain Identification

CHO cell line strains were differentiated based on the **DHFR (dihydrofolate reductase) gene** status, which is the primary genetic marker distinguishing major CHO strains:

| Strain | DHFR Status |
|--------|-------------|
| CHO-K1 | Intact (diploid) |
| CHO-S | Intact (diploid) |
| CHO-DG44 | Homozygous deletion |
| CHO-DXB11 | Hemizygous deletion |

DHFR gene coordinates were extracted from the NCBI RefSeq annotation (GCF_003668045.1_CriGri-PICR_genomic.gff). Sequencing depth at the DHFR locus (NW_020822461.1:37,643,639–37,667,418; ~23.8 kb) was compared to the flanking region (±500 kb) using SAMtools depth. The DHFR-to-flanking depth ratio was used for strain classification:

- Ratio ≈ 1.0: DHFR intact → CHO-K1 or CHO-S
- Ratio ≈ 0.5: Hemizygous deletion → CHO-DXB11
- Ratio ≈ 0: Homozygous deletion → CHO-DG44

### 2.5 Reference Genome

| Parameter | Value |
|-----------|-------|
| Species | *Cricetulus griseus* (Chinese hamster) |
| Assembly | CriGri-PICR |
| Accession | GCF_003668045.1 |
| Source | NCBI RefSeq |
| Genome size | ~2.4 Gb |

### 2.6 Software

| Software | Version | Purpose |
|----------|---------|---------|
| fastp | v1.3.3 | Quality control & read filtering |
| MultiQC | v1.27.1 | QC report aggregation |
| BWA | v0.7.17 | Read alignment |
| SAMtools | v1.21 | BAM processing & statistics |

---

## 3. Results

### 3.1 Sequencing Quality

| Metric | Value | Assessment |
|--------|-------|------------|
| Total reads (sampled) | 10,000,000 | — |
| Q20 rate | 98.69% | Excellent |
| Q30 rate | 95.13% | Excellent |
| GC content | 45.25% | Normal for CHO (~45%) |
| Duplication rate | 3.45% | Very low |
| Insert size peak | 175 bp | Normal |
| Reads passing filter | 99.11% | Excellent |
| Adapter-trimmed reads | 18.96% | Acceptable |

**Conclusion:** Sequencing quality is excellent with high base quality scores and low duplication rate. The GC content of 45.25% is consistent with the expected value for the Chinese hamster genome.

### 3.2 Species Confirmation

| Metric | Value |
|--------|-------|
| Total aligned reads | 6,960,739 |
| Mapped reads | 6,950,296 |
| **Mapping rate** | **99.85%** |
| Properly paired | 91.67% |
| Singletons | 0.07% |

**Conclusion:** A mapping rate of **99.85%** to the CHO reference genome (CriGri-PICR) unambiguously confirms that the sample originates from **Chinese Hamster Ovary (CHO) cells** (*Cricetulus griseus*).

### 3.3 Strain Identification

#### DHFR Gene Locus Analysis

| Region | Coordinates | Mean Depth | Bases Covered |
|--------|-------------|------------|---------------|
| **DHFR gene** | NW_020822461.1:37,643,639–37,667,418 | **1.743x** | 4,890 / 23,779 (20.6%) |
| **Flanking region** (±500 kb) | NW_020822461.1:37,143,639–38,167,418 | **1.725x** | 205,256 / 1,023,779 (20.0%) |

| Metric | Value |
|--------|-------|
| **DHFR / Flanking depth ratio** | **1.01** |
| DHFR coverage fraction | 20.6% |
| Flanking coverage fraction | 20.0% |

> Note: The relatively low absolute coverage (~1.7x at covered bases, ~20% bases covered) is expected because only 5,000,000 read pairs (~0.6x genome-wide average) were used for this analysis. This coverage is sufficient for DHFR deletion detection, as a homozygous deletion would show 0x coverage regardless of sequencing depth.

**Conclusion:** The DHFR-to-flanking depth ratio of **1.01** indicates that the DHFR gene is **fully intact** with normal diploid copy number. This result:

- Is consistent with **CHO-K1** or **CHO-S**
- Rules out **CHO-DG44** (would show ratio ≈ 0)
- Rules out **CHO-DXB11** (would show ratio ≈ 0.5)

---

## 4. Summary & Conclusions

| Question | Answer |
|----------|--------|
| Is the sample CHO? | **Yes** (99.85% mapping rate) |
| Which strain? | **CHO-K1 or CHO-S** (DHFR intact) |
| Data quality? | **Excellent** (Q30 > 95%, dup rate 3.45%) |
| Estimated full genome coverage | **~29x** (based on total data volume) |

### Limitations

1. **CHO-K1 vs CHO-S distinction:** These two strains cannot be reliably distinguished by WGS alone, as their genomic differences are minimal. CHO-S was derived from the original CHO line through suspension adaptation, which primarily involves epigenetic and expression-level changes rather than large-scale genomic alterations. Culture conditions (adherent vs. suspension) can confirm the identity.

2. **Subsampling:** This analysis used 5M read pairs (~0.6x coverage) for rapid assessment. Full-depth analysis (~29x) could be performed for comprehensive variant calling, copy number analysis, and karyotyping if needed.

---

## 5. Deliverables

| File | Description |
|------|-------------|
| `Research_Summary.md` | This report: analysis methods and results |
| `qc/multiqc_report.html` | Interactive quality control report |
| `qc/wt1_fastp.html` | Detailed fastp QC report |
| `results/cho_flagstat.txt` | Alignment statistics |
| `results/dhfr_gff.txt` | DHFR gene annotation entry |

---

## 6. References

1. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*. 2018;34(17):i884-i890.
2. Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics*. 2009;25(14):1754-1760.
3. Danecek P, et al. Twelve years of SAMtools and BCFtools. *GigaScience*. 2021;10(2):giab008.
4. Xu X, et al. The genomic sequence of the Chinese hamster ovary (CHO)-K1 cell line. *Nature Biotechnology*. 2011;29(8):735-741.
5. Lewis NE, et al. Genomic landscapes of Chinese hamster ovary cell lines as revealed by the Cricetulus griseus draft genome. *Nature Biotechnology*. 2013;31(8):759-765.
6. Rupp O, et al. A reference genome of the Chinese hamster based on a hybrid assembly strategy. *Biotechnology and Bioengineering*. 2018;115(8):2087-2100.
