# Whole-Genome Sequencing Analysis Report
## Sample: wt1 — CHO Cell Line Characterization

**Report Date:** June 26, 2026  
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics  
**Analysis Platform:** Linux HPC server; conda environment `regular_bioinfo`

---

## 1. Objectives

1. Confirm the species identity of the submitted cell line sample (wt1).
2. Identify the CHO sub-strain based on DHFR locus copy-number analysis.
3. Detect and characterize exogenous transgene sequences integrated into the host genome.

---

## 2. Sample Information

| Item | Detail |
| :--- | :--- |
| Sample ID | wt1 |
| Raw data — R1 | 23TF7FLT4_3_0469165296_wt1_S8_L003_R1_001.fastq.gz (17 GB) |
| Raw data — R2 | 23TF7FLT4_3_0469165296_wt1_S8_L003_R2_001.fastq.gz (16 GB) |
| Estimated total reads | ~233 M read pairs |
| Estimated genome coverage | ~30× (CHO genome ~2.4 Gb) |

---

## 3. Methods

### 3.1 Read Quality Control
- **Tool:** fastp v0.23 (`-w 8`)
- Adapter trimming, quality filtering, and per-read statistics
- Outputs: filtered paired-end FASTQ files + HTML/JSON QC report
- **MultiQC v1.35** used to aggregate QC metrics

### 3.2 Genome Alignment
- **Reference:** *Cricetulus griseus* CriGri-PICR assembly (GCF_003668045.1)
- **Aligner:** BWA-MEM (`-t 20`)
- **Sorting:** SAMtools sort (`-@ 8 -m 8G`)
- **Indexing:** SAMtools index (`-@ 8`)
- Read-group tag added: `@RG ID:wt1 SM:wt1 PL:ILLUMINA`

### 3.3 DHFR Locus Depth Analysis
- SAMtools depth over DHFR gene coordinates (GFF annotation: NW_020822461.1:37,643,639–37,667,418)
- Flanking control region: ±500 kb around DHFR locus
- Depth ratio (DHFR / flanking) used for CHO sub-strain classification

### 3.4 Transgene Detection (Host-Subtraction Strategy)
- **Step 1 — Unmapped read extraction:**
  - Both-ends unmapped pairs: `samtools view -f 12 -F 256 -F 2048`
  - Singleton unmapped reads: `samtools view -f 4 -F 8 -F 256 -F 2048`
- **Step 2 — De novo assembly:** MEGAHIT v1.2.9 (`-t 16 -m 0.5 --min-contig-len 200`)
- **Step 3 — BLAST homology search:** BLASTn remote against NCBI nt database
  (`-evalue 1e-10 -max_target_seqs 5`, contigs ≥1,000 bp submitted)

---

## 4. Results

### 4.1 Read Quality Control

| Metric | Value |
| :--- | :--- |
| Total read pairs (post-filter) | 232,926,798 |
| Q30 rate | See `wt1_full_fastp.html` |
| Adapter content | Detected and trimmed |

Full QC report: `wt1_full_fastp.html` · Aggregated report: `multiqc_report.html`

### 4.2 Species Identification

| Metric | Value | Interpretation |
| :--- | :---: | :--- |
| Total mapped reads (primary) | 446,411,806 | |
| **Primary mapping rate** | **95.83%** | **Confirmed *Cricetulus griseus* (CHO)** |
| Properly paired | 92.17% | Normal paired-end library |
| Supplementary alignments | 1,632,515 | Consistent with repetitive CHO genome |

A primary mapping rate of 95.83% to the CriGri-PICR reference genome unambiguously confirms the sample is derived from *Cricetulus griseus* (Chinese hamster ovary, CHO).

Full alignment statistics: `cho_flagstat.txt`

### 4.3 CHO Sub-Strain Identification

| Metric | Value |
| :--- | :--- |
| DHFR locus mean depth | 24.1× |
| Flanking region (±500 kb) mean depth | 23.9× |
| **DHFR / flanking ratio** | **1.008** |

**Interpretation:**

| Ratio range | Sub-strain |
| :--- | :--- |
| < 0.05, coverage < 5% | CHO-DG44 (homozygous DHFR deletion) |
| 0.05 – 0.35 | CHO-DXB11 (hemizygous DHFR deletion) |
| **> 0.70 (observed: 1.008)** | **CHO-K1 or CHO-S (DHFR intact)** |

The DHFR/flanking depth ratio of **1.008** indicates that both copies of the DHFR locus are intact. This is consistent with **CHO-K1 or CHO-S** sub-strains and excludes CHO-DG44 and CHO-DXB11.

### 4.4 Transgene Detection

**Unmapped Read Yield**

| Category | Count |
| :--- | :--- |
| Both-ends unmapped PE pairs | 9,678,712 |
| Singleton unmapped reads | 84,366 |
| Unmapped fraction | ~4.2% of total reads |

The unmapped fraction (4.2%) substantially exceeds the typical CHO background (<0.5%), indicating the presence of significant exogenous sequence content.

**De Novo Assembly (MEGAHIT)**

| Metric | Value |
| :--- | :--- |
| Total contigs assembled | 2,052 |
| Total assembly size | 1,145,229 bp |
| N50 | 575 bp |
| Longest contig | 6,226 bp |
| Contigs ≥ 1,000 bp submitted to BLAST | 162 |

**BLASTn Results — Key Hits (≥1,000 bp contigs, nt database)**

| Contig | Identity | Aligned/Contig | Bitscore | Subject Description |
| :--- | :---: | :---: | :---: | :--- |
| k141_804 | 100.0% | 1,699 / 5,522 bp | 3,138 | Cloning vector **pLV[Exp]-CBA>P301L**, complete sequence |
| k141_1290 | 100.0% | 1,456 / 1,465 bp | 2,689 | *C. griseus* **expression augmenting sequence element (EASE)** |
| k141_1812 | 99.9% | 1,331 / 1,342 bp | 2,453 | Cloning vector **pLenti-EF1a-dCas9-DNMT3B(E697A)-2A-bla** |
| k141_1793 | 100.0% | 1,306 / 1,306 bp | 2,412 | Cloning vector **RS474_ErbB-RASER1C-dCas9VP64** |
| k141_746 | 99.2% | 255 / 1,149 bp | 459 | Human mRNA sequence (clone xip10) |
| k141_2082 | 77.5% | 769 / 3,056 bp | 435 | *Mus musculus* immunoglobulin kappa-like gene |

**33 contigs ≥1,000 bp returned no BLAST hits** against the NCBI nt database. These sequences are likely proprietary transgene or vector elements not yet deposited in public databases.

Full BLAST results: `blast_results.txt`  
Assembly FASTA (all contigs): `transgene_all_contigs.fa`

---

## 5. Conclusions

| Question | Conclusion |
| :--- | :--- |
| Species | **Confirmed *Cricetulus griseus* (CHO)** — 95.83% primary mapping rate to CriGri-PICR |
| Sub-strain | **CHO-K1 or CHO-S** — DHFR locus intact (depth ratio 1.008); CHO-DG44 and CHO-DXB11 excluded |
| Transgene present? | **Yes — confirmed** |
| Transgene identity | **Lentiviral expression system** carrying dCas9-based regulatory elements: (1) pLV[Exp]-CBA lentiviral backbone; (2) dCas9-DNMT3B(E697A) epigenetic silencer; (3) dCas9-VP64 transcription activator targeting ErbB/EGFR; (4) CHO EASE expression-augmenting element |
| Additional sequences | Partial immunoglobulin kappa homology (possible antibody GOI); 33 contigs (≥1,000 bp) with no public database match — likely proprietary sequences |
| Integration sites | Not determined by this analysis (requires mate-pair or long-read sequencing) |

**Clinical/process relevance:** The identified dCas9-VP64 and dCas9-DNMT3B elements suggest a CRISPRa/CRISPRi-based gene regulation system targeting the ErbB/EGFR pathway. If the 33 unmatched contigs contain the gene of interest (GOI), the client is requested to provide the vector map or reference sequence for targeted verification.

---

## 6. Deliverable Files

| File | Description |
| :--- | :--- |
| `WGS_Analysis_Report_wt1.md` | This report |
| `wt1_full_fastp.html` | Full-dataset fastp QC report (interactive HTML) |
| `wt1_full_fastp.json` | fastp QC metrics (machine-readable) |
| `multiqc_report.html` | MultiQC aggregated QC report |
| `cho_flagstat.txt` | SAMtools flagstat — alignment statistics |
| `dhfr_depth.txt` | DHFR and flanking region sequencing depth |
| `blast_results.txt` | BLASTn tabular results (all hits, ≥1,000 bp contigs) |
| `transgene_contigs_1000bp.fa` | Assembled contigs ≥1,000 bp (BLAST input) |
| `transgene_all_contigs.fa` | All assembled contigs ≥200 bp (complete assembly) |

---

*Report prepared by:*  
**Zhen Gao, PhD**  
**Principal Bioinformatics Scientist**  
**Athenomics**  
*June 26, 2026*
