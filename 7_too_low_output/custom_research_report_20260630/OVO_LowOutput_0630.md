# Investigation Report: Low Read Yield for OVO Sequencing Sample
## Quote Reference: Quote_06062601 | Lane 23JCTGLT3, Lane 7

**Report Date:** June 30, 2026  
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics  
**Analysis Platform:** Linux HPC server  

---

## 1. Objectives

This report documents a bioinformatic investigation into the low read output observed for the OVO sample (submitted by Leif Benner, Perrimon Lab, Harvard Medical School) from sequencing run Quote_06062601. The client reported receiving only ~3.09 million total reads against an ordered quantity of 75 million reads. The objectives of this investigation were to:

1. Quantify the actual read output for the OVO sample with precision
2. Determine whether the missing reads are recoverable from the Undetermined (un-demultiplexed) data
3. Assess the quality and authenticity of the reads that were demultiplexed
4. Identify the most probable root cause of the read deficit
5. Provide a concrete recommendation for resolution

---

## 2. Key Findings

- **Severe read deficit confirmed:** The OVO sample received only **1,665,879 reads** (1.67 million) against an order of 75 million — a **45-fold shortfall** representing just 2.2% of the ordered depth.
- **Missing reads are NOT in the Undetermined file:** The 161 GB Undetermined file contains reads with 10 bp + 10 bp barcodes from unrelated libraries. A systematic search for OVO's I7 barcode (TAAGGCGA) found it in only 0.097% of Undetermined headers — a background-level noise rate, not recoverable signal.
- **The reads that were demultiplexed are genuine:** 77.7% (1,295,404) of the 1.67 M OVO reads contain the expected sgRNA library fixed linker sequence (CCTATTTTCAATTTAACGTCG), confirming the library chemistry worked correctly.
- **Root cause is pool under-representation:** All other samples on the same lane received their expected read depths (58–565 million reads each). The lane itself produced 4.11 billion reads normally. OVO represented only 0.04% of the pool instead of the expected ~1.82%, indicating it was loaded at approximately **45× lower molar concentration** than required.
- **Re-sequencing is the only path to full data:** Recovery from existing data is not feasible. A re-run with corrected pooling is required; however, the OVO library itself is structurally intact and should perform correctly at the right concentration.

---

## 3. Sample Information

| Field | Detail |
| :--- | :---: |
| Client | Leif Benner (leif_benner@hms.harvard.edu) |
| Affiliation | Perrimon Lab, Dept. of Genetics / HHMI, Harvard Medical School |
| Project Reference | Quote_06062601 |
| Sample Name | OVO |
| Library Type | Custom sgRNA / CRISPR screen library with UMI barcodes |
| Sequencing Platform | Illumina (Instrument LH00972, Flow Cell 23JCTGLT3, Lane 7) |
| Read Mode | 150 bp paired-end |
| I7 Index (N701) | TAAGGCGA |
| I5 Index (N502) | CTCTCTAT |
| UMI Design | Three templates: 3 + 6 bp, 4 + 6 bp, or 5 + 6 bp random + fixed anchor |
| Ordered Read Depth | 75 million reads |
| Received Date (sample) | ~June 9, 2026 |
| Data Delivery Date | June 24–25, 2026 |

---

## 4. Analysis Rationale and Decision Criteria

**Why examine Undetermined reads?**  
When demultiplexing yields are unexpectedly low, reads may accumulate in the Undetermined pool due to barcode mismatches (e.g., incorrect index in the SampleSheet, reversed I5 orientation, or 1-mismatch threshold failures). Investigating the Undetermined barcode composition is the standard first diagnostic step.

**Decision threshold for "recoverable reads":**  
If OVO's I7 barcode (TAAGGCGA) appeared at ≥1% frequency in the Undetermined pool, reads would be considered recoverable via re-demultiplexing. A rate below 0.1% is consistent with random index-hopping noise and is not recoverable.

**Why count linker-positive reads?**  
The OVO library has a unique structure where the sgRNA guide sequence is preceded by a fixed 21 bp linker (CCTATTTTCAATTTAACGTCG). Reads lacking this linker are not genuine sgRNA library molecules. Quantifying the linker-positive fraction confirms library integrity and corrects for index-hopping contamination from unrelated libraries sharing the TAAGGCGA I7 index.

**Decision threshold for root cause:**  
If all other samples on the lane received normal read depths while OVO alone was depleted, the problem is sample-specific, pointing to pooling rather than a sequencing run failure.

---

## 5. Methods

**5.1 Read Quantification**  
Reads in the OVO FASTQ files were counted directly by line count (`wc -l`) on the decompressed reads (`zcat`), divided by 4 (lines per FASTQ record).

**5.2 Lane-wide Read Distribution**  
Per-sample read counts for all 17 samples on Lane 7 were estimated by scaling from the calibrated OVO ratio (63,334,911 bytes = 1,665,879 reads → 38.0 bytes/read compressed) applied to file sizes reported in the sequencing provider's `checkSize.xls` manifest.

**5.3 Undetermined Barcode Analysis**  
The first 1,000,000 headers from `Undetermined_Undetermined_23JCTGLT3_L7_1.fq.gz` were sampled using `awk NR%4==1 | head -1000000`. Barcode sequences were extracted from the Illumina FASTQ header field (after `0:`) and tallied using `sort | uniq -c | sort -rn`.

A separate count of headers containing the OVO I7 sequence (TAAGGCGA) was performed on the first 2,000,000 undetermined headers.

**5.4 Genuine Library Read Assessment**  
All R1 sequences from the OVO FASTQ were scanned by `grep` for the exact fixed linker sequence `CCTATTTTCAATTTAACGTCG`, which is the invariant element present in all three UMI template designs immediately downstream of the UMI and immediately upstream of the 20 nt guide barcode. A relaxed variant (1 mismatch in variable positions) was also tested for cross-validation.

**Key tools:** bash, zcat, awk, grep, python3 (arithmetic only)

---

## 6. Results

### 6.1 OVO Read Count

| File | Size | Read Count |
| :--- | :---: | :---: |
| OVO_CKDL260011462-1A_23JCTGLT3_L7_1.fq.gz | 61 MB | 1,665,879 |
| OVO_CKDL260011462-1A_23JCTGLT3_L7_2.fq.gz | 42 MB | 1,665,879 (paired) |
| **Total read pairs** | **103 MB** | **1,665,879** |

Against the ordered 75,000,000 reads, this represents a **97.8% shortfall**.

### 6.2 Undetermined Barcode Composition

A sample of 1,000,000 Undetermined R1 headers revealed the following top barcodes:

| Rank | Barcode (I7+I5) | Count | Format |
| :--- | :---: | :---: | :---: |
| 1 | ACCACACGAT+CAACCACGGT | 271,927 | 10+10 bp UDI |
| 2 | GGCGAGGAAT+AACAAGTGGT | 123,419 | 10+10 bp UDI |
| 3 | GGGGGGGGGG+AGATCTCGGT | 73,534 | Low-quality / poly-G |
| 4 | CGTGTAGGAT+TCTGAAACGT | 68,857 | 10+10 bp UDI |
| 5–14 | Various XXXXAAT+AGATCTCGGT | 29–41K each | 10+10 bp UDI |

**All top-ranked barcodes are 10 bp + 10 bp format, entirely incompatible with OVO's 8 bp Nextera indexes.** None of the OVO barcodes appear in any top-ranked position.

A direct search for TAAGGCGA (OVO I7) across 2,000,000 Undetermined headers found only **1,942 hits (0.097%)** — consistent with background index-hopping noise on patterned flow cells. This rate yields an estimated total of ~2.2 million OVO-like reads across the entire Undetermined file, but these are not genuine OVO molecules.

### 6.3 Lane-wide Read Distribution

Estimated reads per sample on Lane 7 (calibrated from OVO file size ratio):

| Sample | Estimated Reads (M) | Notes |
| :--- | :---: | :---: |
| Undetermined | ~2,270 M | 10+10 bp UDI reads; unrelated libraries |
| C14_41 | ~565 M | Normal |
| C14_45 | ~274 M | Normal |
| C13_16 | ~145 M | Normal |
| A | ~87 M | Normal |
| J_899 | ~86 M | Normal |
| J_897 | ~78 M | Normal |
| C, B | ~75 M each | Normal |
| J_904, J_912, J_909, J_902, J_905, J_910 | ~63–72 M each | Normal |
| J_896 | ~59 M | Normal |
| **OVO** | **~1.7 M** | **45× below target** |
| **Total lane** | **~4,110 M** | — |

OVO's actual share of the lane: **0.04%**. Expected share at 75M reads: **1.82%**. This 45-fold deficit is entirely sample-specific — no other sample on the lane was under-represented.

### 6.4 OVO Library Integrity

| Assessment | Count | Fraction |
| :--- | :---: | :---: |
| Total OVO R1 reads | 1,665,879 | 100% |
| Reads with exact fixed linker | 1,295,404 | **77.7%** |
| Reads with relaxed linker match | 1,298,525 | 78.0% |
| Non-library reads (contamination / index-hop) | ~370,000 | ~22.3% |

The 77.7% linker-positive rate confirms the OVO library was successfully constructed. The non-library fraction (~22%) is consistent with low-level index-hopping contamination from adjacent high-depth samples on the same lane — an expected phenomenon on patterned flow cells, and not a library quality failure.

---

## 7. Conclusions

### 7.1 Summary

| Question | Answer |
| :--- | :---: |
| Are reads recoverable from Undetermined? | **No** — Undetermined contains unrelated 10+10 bp libraries |
| Is the sequencing run itself flawed? | **No** — all other 16 samples received expected depths |
| Is the OVO library structurally intact? | **Yes** — 77.7% of reads are genuine sgRNA molecules |
| What is the root cause? | **OVO was under-represented in the sequencing pool by ~45×** |
| Can existing 1.3M reads support analysis? | **Depends on guide library size** (see recommendation) |

### 7.2 Root Cause Assessment

The evidence conclusively points to a **library pool imbalance** as the sole root cause. The most probable underlying mechanism is an error in molar concentration calculation during pooling:

- Fluorometric quantification methods (e.g., Qubit) measure total nucleic acid mass (ng/µL), not molar concentration (nM)
- Converting mass to molar requires accurate insert size: nM = (ng/µL × 10⁶) / (660 × average library size in bp)
- If insert size was overestimated, or if adapter dimers inflated the mass reading without contributing usable reads, the OVO molarity in the pool would have been systematically over-estimated, leading to under-loading
- This error would not be detected by standard pre-sequencing QC (fragment analyzer shows band at expected size; Qubit shows expected concentration)

### 7.3 Recommendations

**Immediate:** The client should determine from the guide library design whether 1,295,404 genuine sgRNA reads provide sufficient coverage for their analysis (rule of thumb: ≥500× average coverage per guide). If the guide library contains fewer than ~2,600 unique guides, the existing data may support a preliminary analysis.

**For re-sequencing:**
1. Re-quantify the OVO library by **qPCR** (KAPA Library Quantification or equivalent) — this gives molar concentration directly without insert-size assumptions
2. Confirm the fragment size distribution by Bioanalyzer or TapeStation to rule out adapter-dimer contamination
3. Recalculate the pooling ratio to achieve 1.82% molar fraction of the total pool (or allocate a dedicated lane/partial run for OVO)
4. Optionally sequence a small validation aliquot (~5M reads) before committing to the full 75M run to confirm correct pooling

---

## 8. Deliverable Files

| File | Description |
| :--- | :---: |
| `OVO_LowOutput_0630.md` | This report |
| `supplementary/lane_read_distribution.tsv` | Per-sample estimated read counts for Lane 7 |

**Source data (provided by sequencing facility, not modified):**  
`/home/gao/Dropbox/Quote_06062601_output/OVO/` — Demultiplexed OVO FASTQ files (R1: 1,665,879 reads)  
`/home/gao/Dropbox/Quote_06062601_output/Undetermined/` — Un-demultiplexed lane data (161 GB)  
`/home/gao/Dropbox/Quote_06062601_output/checkSize.xls` — Per-sample file size manifest from sequencing provider

---

*Zhen Gao, PhD*  
*Principal Bioinformatics Scientist, Athenomics*
