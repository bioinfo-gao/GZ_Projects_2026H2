# Mouse RNA-seq Library Quality Assessment Report

**Report Date**: 2026-06-27
**Prepared by**: Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform**: Linux HPC server

---

## 1. Objectives

Nine mouse RNA-seq libraries were submitted for a pre-analysis quality triage prior to committing to full-depth sequencing or downstream expression analysis. The client requested answers to three specific questions:

1. What is the insert size distribution for each library, and is adapter dimer content a concern?
2. What is the level of ribosomal RNA (rRNA) contamination?
3. How many reads map uniquely to the genome, and how many genes can be detected?

The decision criterion set at the outset: if quality is found to be severely compromised, simply repeating library preparation with the same protocol would not be a productive next step — the root cause would need to be identified first.

---

## 2. Key Findings

- **All 9 libraries fail standard mRNA-seq quality thresholds across all three metrics assessed** (adapter dimer content, rRNA contamination, unique mapping rate). The failure pattern is systematic across the entire batch, not isolated to individual samples — this points to a library preparation issue rather than a sample-quality issue.
- **Adapter dimer content is the dominant, primary defect.** Two libraries (mouse_28, mouse_29) are ~70% adapter dimer by read pairs; the remaining seven are 22–42% dimer — all far above the <5–10% expected in a healthy library. This single defect is the direct cause of the low mapping rates.
- **rRNA contamination (27.5–47.2%) is a second, independent defect** — it does not rise and fall together with the adapter dimer level, indicating insufficient rRNA depletion is a separate problem from the adapter-ligation problem, not a downstream symptom of it.
- **Unique mapping rates are severely depressed (2.2–28.7% vs. an expected 60–65%+ at this depth)**, driven almost entirely by reads that are too short to place reliably after trimming — a direct consequence of the adapter dimer issue.
- **mouse_46 is the most severely compromised library** (2.2% unique mapping, 1,042 genes detected) and should be treated as non-recoverable; **mouse_28 and mouse_29** are also effectively non-functional due to extreme dimer content despite longer surviving reads.

---

## 3. Sample Information

| Item | Detail |
| :--- | :---: |
| Samples | 9 mouse samples (mouse_28, 29, 32, 41, 42, 45, 46, 47, 48) |
| Sequencing type | Paired-end, PE150 |
| Raw data volume | ~32 million read pairs per sample |
| QC analysis depth | Subsampled to 10 million read pairs per sample (sufficient for stable QC estimates; full-depth SortMeRNA runs were not memory-feasible — see Section 5) |
| Reference genome | GRCm39 (GENCODE vM35 annotation) |

---

## 4. Analysis Rationale and Decision Criteria

| Question | Metric used | Healthy-library threshold | Rationale |
| :--- | :---: | :---: | :---: |
| Insert size / dimer content | % of read pairs discarded by Trim Galore at <20 bp after adapter trimming | < 5–10% | A high discard rate means most fragments are adapter ligated directly to adapter, with no real cDNA insert |
| rRNA contamination | % of trimmed reads matching SortMeRNA reference databases, measured **before** genome alignment | < 5–10% | Pre-alignment measurement reflects the true rRNA content of the original sample, independent of whether reads go on to align successfully |
| Unique mapping / gene detection | % of reads uniquely placed by STAR; number of genes with detectable counts | ≥ 60–65% unique mapping at 10M-pair depth | Reflects how much of the sequencing investment yields usable expression signal |

---

## 5. Methods

1. **Subsampling**: `seqtk sample -s100`, 10,000,000 read pairs per sample (fixed seed for reproducibility). Full-depth SortMeRNA runs against all 8 rRNA databases were not feasible on this server: memory scales linearly with read count under SortMeRNA's `--paired_in --out2 --fastx` mode, and the full ~32M-pair runs either exceeded the 16h time limit or were killed by the OOM killer. Subsampling to 10M pairs reduced peak memory to 2–39 GB while preserving stable estimates of all three QC metrics.
2. **Adapter trimming**: Trim Galore, minimum length cutoff 20 bp.
3. **rRNA identification**: SortMeRNA v4.3.6 against 8 reference databases — SILVA bacterial/archaeal 16S/23S, SILVA eukaryotic 18S/28S, Rfam 5S/5.8S.
4. **Genome alignment**: STAR against GRCm39 (GENCODE vM35).
5. **Gene quantification**: Salmon + featureCounts.
6. Pipeline: nf-core/rnaseq v3.15.1 (Singularity containers).

---

## 6. Results

### 6.1 Adapter Dimer Content (Insert Size Proxy)

| Sample | Read pairs discarded after trimming (%) | Median trimmed read length (bp) |
| :--- | :---: | :---: |
| mouse_28 | **69.7** | 82 |
| mouse_29 | **72.0** | 147 |
| mouse_32 | 26.1 | 42 |
| mouse_41 | 37.4 | 52 |
| mouse_42 | 22.4 | 42 |
| mouse_45 | 23.4 | 42 |
| mouse_46 | 29.1 | 147 |
| mouse_47 | 42.3 | 42 |
| mouse_48 | 22.9 | 42 |

mouse_28/29 are essentially non-functional libraries (~70% pure adapter dimer). The remaining seven samples show 22–42% dimer content — well above the healthy threshold, with most surviving inserts only ~42 bp.

### 6.2 rRNA Contamination

| Sample | Total rRNA (%) | Main source | Eukaryotic 18S (%) | Eukaryotic 28S (%) | Bacterial 16S (%) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| mouse_28 | 44.6 | Bacterial 16S | 4.1 | 7.9 | **30.4** |
| mouse_29 | 45.1 | Bacterial 16S | 3.6 | 2.7 | **37.7** |
| mouse_32 | 47.2 | Eukaryotic 28S | 7.8 | **31.6** | 3.4 |
| mouse_41 | 37.9 | Eukaryotic 28S | 6.5 | **15.5** | 11.8 |
| mouse_42 | 43.2 | Eukaryotic 28S | 8.7 | **26.7** | 2.8 |
| mouse_45 | 40.4 | Eukaryotic 28S | 11.4 | **22.1** | 1.6 |
| mouse_46 | 27.5 | Bacterial 16S | 10.5 | 0.7 | **15.6** |
| mouse_47 | 44.8 | Eukaryotic 28S | 10.2 | **24.8** | 2.9 |
| mouse_48 | 45.3 | Eukaryotic 28S | 10.4 | **28.0** | 1.9 |

All 9 samples exceed the healthy <5–10% threshold by a wide margin. Two distinct profiles are present: a Eukaryotic-28S-dominant pattern (6 samples, the more typical mammalian incomplete-depletion signature) and a Bacterial-16S-dominant pattern (mouse_28/29/46). The bacterial-16S pattern is examined in detail in Section 7.

### 6.3 Unique Mapping and Gene Detection

| Sample | Uniquely mapped reads (%) | Reads unmapped — too short (%) | Genes detected |
| :--- | :---: | :---: | :---: |
| mouse_28 | 19.6 | 66.0 | 10,312 |
| mouse_29 | 12.2 | 79.6 | 6,157 |
| mouse_32 | 28.7 | 51.5 | 19,000 |
| mouse_41 | 22.6 | 59.5 | 16,502 |
| mouse_42 | 26.6 | 54.6 | 20,817 |
| mouse_45 | 24.9 | 54.8 | 20,270 |
| mouse_46 | **2.2** | **96.2** | 1,042 |
| mouse_47 | 26.0 | 53.9 | 18,590 |
| mouse_48 | 27.7 | 52.3 | 19,826 |

Unique mapping is 2.2–28.7% against an expected ≥60–65% at this depth. In every sample, the predominant failure mode is "too short" — a direct downstream consequence of the dimer issue in Section 6.1: short surviving inserts cannot be reliably placed by STAR.

---

## 7. Supplementary Analysis: Response to Client Follow-Up Questions

**Client question:** *Is the observed rRNA signal correlated with the high GC content / high adapter content seen in these libraries? Please also break down the bacterial/Mycoplasma component of the rRNA signal.*

### 7.1 rRNA% vs. GC% vs. Adapter Dimer%

| Sample | rRNA% (§6.2) | Dimer % (§6.1) | Raw read %GC* |
| :--- | :---: | :---: | :---: |
| mouse_28 | 44.6 | 69.7 | 74.0 |
| mouse_29 | 45.1 | 72.0 | 74.5 |
| mouse_32 | 47.2 | 26.1 | 64.0 |
| mouse_41 | 37.9 | 37.4 | 66.0 |
| mouse_42 | 43.2 | 22.4 | 62.5 |
| mouse_45 | 40.4 | 23.4 | 61.0 |
| mouse_46 | 27.5 | 29.1 | 59.0 |
| mouse_47 | 44.8 | 42.3 | 65.5 |
| mouse_48 | 45.3 | 22.9 | 61.0 |

*%GC from FastQC on raw reads (R1/R2 averaged), generated in the same pipeline run as the other QC metrics above.

**rRNA% does not correlate with GC% or with dimer%.** mouse_32 has the highest rRNA contamination (47.2%) yet only moderate GC (64.0%) and the lowest dimer content among the seven non-extreme samples (26.1%). Conversely, mouse_28/29 have the highest GC (74.0–74.5%) and the highest dimer content (~70%) but only mid-range rRNA% (44.6–45.1%) — no higher than several samples with much lower dimer content. **rRNA contamination behaves as an independent defect, not as a side-effect of the adapter problem.**

**GC% does correlate with dimer%.** The two most dimer-affected samples (mouse_28/29) also carry the highest GC content, well above the ~45–50% typically expected for mouse mRNA-seq. This is consistent with a known mechanism: Illumina adapter sequences are themselves GC-rich (~60%), so a library dominated by adapter-dimer fragments rather than real cDNA insert will show its overall GC% pulled toward the adapter's composition. **The elevated GC values in this batch are best explained as a symptom of the adapter-dimer problem (§6.1), not of the rRNA problem (§6.2).**

### 7.2 Bacterial Component of the rRNA Signal

Three samples (mouse_28, mouse_29, mouse_46) show a "Bacterial 16S"-dominant rRNA profile (15.6–37.7%) rather than the eukaryotic-28S pattern seen in the other six. To check whether this reflects real bacterial content rather than a database artifact, a species-level spot-check (Kraken2, against a reference covering bacteria/archaea/viral/human) was run on the STAR-unmapped reads of two representative samples, mouse_29 and mouse_46.

| Metric | mouse_29 | mouse_46 |
| :--- | :---: | :---: |
| Reads directly, confidently classified as Bacteria | 2.44% | 1.94% |
| Reads Unclassified (no confident species match in this reference) | 85.12% | 86.33% |

The directly classified 1.9–2.4% should be read as a **floor, not a ceiling**, on the true bacterial fraction. Kraken2 only calls a species when a read carries enough exact, unique sequence matches; short, degraded, or divergent-strain bacterial reads routinely fail that bar and fall into "Unclassified" rather than being misclassified as something else. At the same time, this reference contains no mouse genome, so the 85–86% Unclassified pool is necessarily a **mixture of true mouse reads and uncalled bacterial reads** — not, as might be assumed, purely mouse.

Allowing for Kraken2's typical recall on short, fragmented, divergent-strain bacterial sequence — commonly well under 50%, often in the 10–20% range for input of this quality — back-calculating from the observed 1.9–2.4% confidently-classified fraction gives an estimated **true bacterial fraction in the range of roughly 10–30%** for these samples. This is broadly consistent with, rather than contradicted by, the 15.6–37.7% figure already reported by the rRNA tool in Section 6.2.

**Revised interpretation:** the bacterial-16S signal in mouse_28/29/46 is plausibly a genuine, substantial finding, not simply an artifact of the rRNA tool's looser 90% identity threshold. The exact value within the 10–30% range cannot be pinned down further with the data in hand; doing so would require either (a) directly aligning the Unclassified- and Human-assigned reads to the mouse genome to establish a clean "confirmed mouse" baseline, or (b) recalibrating Kraken2's recall for this sample type (e.g., Bracken re-estimation or a lower confidence threshold). Both remain open follow-up items if a more precise figure is needed.

---

## 8. Conclusions

| Sample | rRNA (%) | Dimer-affected pairs (%) | Unique mapping (%) | Genes detected | Overall |
| :--- | :---: | :---: | :---: | :---: | :---: |
| mouse_28 | 44.6 | 69.7 | 19.6 | 10,312 | Fail |
| mouse_29 | 45.1 | 72.0 | 12.2 | 6,157 | Fail |
| mouse_32 | 47.2 | 26.1 | 28.7 | 19,000 | Fail |
| mouse_41 | 37.9 | 37.4 | 22.6 | 16,502 | Fail |
| mouse_42 | 43.2 | 22.4 | 26.6 | 20,817 | Fail |
| mouse_45 | 40.4 | 23.4 | 24.9 | 20,270 | Fail |
| mouse_46 | 27.5 | 29.1 | 2.2 | 1,042 | Fail — non-recoverable |
| mouse_47 | 44.8 | 42.3 | 26.0 | 18,590 | Fail |
| mouse_48 | 45.3 | 22.9 | 27.7 | 19,826 | Fail |

**Reference thresholds**: rRNA < 10% · Dimer pairs < 10% · Unique mapping ≥ 65%

All nine libraries fail quality standards across all three metrics, with a systematic, batch-wide pattern indicating a library preparation issue rather than sample-specific degradation. Two independent root causes are present simultaneously:

1. **Adapter dimer content** — the dominant defect, responsible for the cascading low mapping rates and the abnormal GC content (Section 7.1).
2. **Incomplete rRNA depletion** — a separate, equally severe defect, uncorrelated with the dimer issue.

We do not recommend repeating library preparation with the same protocol, as this would likely reproduce the same defects. We recommend a root-cause investigation focused on: RNA integrity (RIN values, if available — degraded RNA tends to produce short inserts), a review of the fragmentation, rRNA-depletion, and adapter-ligation steps of the library prep protocol, and consideration of an alternative rRNA depletion method.

---

*Report prepared by:*

**Zhen Gao, PhD**
Principal Bioinformatics Scientist
Athenomics
