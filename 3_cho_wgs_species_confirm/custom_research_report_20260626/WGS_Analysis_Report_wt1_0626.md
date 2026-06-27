# Whole-Genome Sequencing Analysis Report
## Sample: wt1 — CHO Cell Line Characterization

**Report Date:** June 26, 2026  
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics  
**Analysis Platform:** Linux HPC server

---

## 1. Objectives

1. Confirm the species identity of the submitted cell line sample (wt1).
2. Identify the CHO sub-strain based on DHFR locus copy-number analysis.
3. Detect and characterize exogenous transgene sequences integrated into the host genome.

---

## 2. Key Findings

- **The cell line is confirmed to be CHO** (*Cricetulus griseus*). 95.83% of sequencing reads map to the Chinese hamster reference genome, well above the 85% threshold for confident species confirmation.
- **The sub-strain is CHO-K1 or CHO-S.** The DHFR gene is fully intact (depth ratio 1.008 ≈ 1.0), which rules out CHO-DG44 and CHO-DXB11. CHO-K1 and CHO-S cannot be distinguished by WGS alone.
- **A lentiviral transgene system is present.** 4.2% of reads do not map to the CHO genome — far above normal background. Assembly and database search identified a lentiviral expression system carrying dCas9-based gene regulation elements targeting the ErbB/EGFR pathway, plus a CHO-specific expression enhancer (EASE).
- **33 assembled sequences have no public database match.** These likely represent proprietary vector or GOI sequences. The client is invited to provide the vector map for targeted confirmation.

---

## 3. Sample Information

| Item | Detail |
| :--- | :--- |
| Sample ID | wt1 |
| Raw data — R1 | 23TF7FLT4_3_0469165296_wt1_S8_L003_R1_001.fastq.gz (17 GB) |
| Raw data — R2 | 23TF7FLT4_3_0469165296_wt1_S8_L003_R2_001.fastq.gz (16 GB) |
| Estimated total reads | ~233 M read pairs |
| Estimated genome coverage | ~30× (CHO genome ~2.4 Gb) |

---

## 4. Analysis Rationale and Decision Criteria

This section explains the scientific basis for each analytical step and the criteria used to reach conclusions, so that results can be independently interpreted.

### 4.1 Why Align to a Reference Genome for Species Identification?

Short-read whole-genome sequencing (WGS) produces millions of DNA fragments whose sequences are characteristic of the source organism. By aligning these reads against a high-quality reference genome — in this case the *Cricetulus griseus* CriGri-PICR assembly — we measure what fraction of the reads are consistent with that species. A **primary mapping rate > 85%** is taken as confident confirmation of species identity: by chance, reads from an unrelated organism would align at rates near zero. Rates between 50–85% are ambiguous and require further investigation; rates below 50% indicate the sample is not CHO. The 4.2% unmapped fraction in this sample is not indicative of a different species; rather, it reflects exogenous transgene sequences (see Section 4.4).

### 4.2 Why Use DHFR Locus Copy Number for Sub-Strain Identification?

The dihydrofolate reductase (*DHFR*) gene locus is the defining genomic marker for distinguishing the three major CHO sub-strains used in biopharmaceutical manufacturing:

- **CHO-DG44**: Both copies of *DHFR* are deleted. The locus is absent, yielding near-zero sequencing depth relative to flanking regions.
- **CHO-DXB11**: One copy of *DHFR* is deleted and one remains. Sequencing depth at the locus is approximately half the flanking depth (ratio ~0.5).
- **CHO-K1 / CHO-S**: Both copies of *DHFR* are intact. Sequencing depth at the locus equals the flanking depth (ratio ~1.0).

WGS provides quantitative, locus-specific sequencing depth at single-nucleotide resolution. By computing the mean depth over the annotated *DHFR* exon region and comparing it to a 1 Mb flanking control region (which has no known copy-number variation in any CHO sub-strain), we obtain a robust, normalization-free depth ratio. This approach is insensitive to overall sequencing depth and is unaffected by GC bias because both the locus and its flanking control are on the same chromosome and share similar sequence composition.

### 4.3 Why Use a Host-Subtraction Strategy for Transgene Detection?

Direct de novo assembly of the entire ~30× WGS dataset (~460 M reads, 33 GB) is computationally prohibitive for transgene discovery: it would require >500 GB RAM and would produce a genome-scale assembly dominated by CHO host sequence. The **host-subtraction** strategy bypasses this problem:

1. Reads that fail to align to the CHO reference genome are, by definition, derived from sequences not present in the reference — including integrated transgenes, viral vector elements, and expression cassettes.
2. These non-host reads are isolated and assembled independently. Because the unmapped fraction is small (here ~19 M reads), de novo assembly is fast, memory-efficient, and produces high-quality contigs enriched for transgene content.
3. The strategy detects **any** foreign sequence regardless of prior knowledge of the transgene — it does not require the client to provide a vector map in advance.

**Limitation:** If a transgene sequence is highly similar to an endogenous CHO gene (e.g., a humanized gene whose rodent ortholog is in the reference), the corresponding reads may align to the CHO genome and be excluded from the unmapped pool, potentially reducing sensitivity. In such cases, targeted re-alignment using only the transgene sequence as reference is recommended.

### 4.4 Why De Novo Assembly Rather Than Direct BLAST of Reads?

Individual short reads (150 bp) are too short for reliable BLAST-based identification: many reads match repetitive elements or low-complexity sequences shared across many organisms, generating noisy, ambiguous results. De novo assembly with MEGAHIT reconstructs longer contiguous sequences (contigs) by finding overlaps among reads. Longer contigs (≥200 bp, preferably ≥1,000 bp) provide sufficient sequence context for BLAST to distinguish closely related sequences and assign confident taxonomic or functional identity. Assembly also consolidates coverage: dozens or hundreds of reads supporting the same transgene region collapse into a single representative contig, dramatically reducing the number of BLAST queries while increasing specificity.

### 4.5 Decision Criteria Summary

| Analysis step | Metric | Threshold | Conclusion |
| :--- | :---: | :---: | :--- |
| Species ID | Primary mapping rate | > 85% = confirmed; 50–85% = ambiguous; < 50% = not CHO | Confirmed if > 85% |
| Sub-strain | DHFR / flanking depth ratio | < 0.05 → DG44; 0.05–0.35 → DXB11; > 0.70 → K1/S | Match to ratio range |
| Transgene presence | Unmapped read fraction | > 0.5% above background is significant | Flag for assembly |
| Transgene identity | BLASTn e-value / identity | e-value ≤ 1×10⁻¹⁰, identity ≥ 95% for high-confidence hit | Report top hits per contig |
| Unknown sequences | No BLAST hit | Contigs ≥ 1,000 bp with zero hits | Likely proprietary; request vector map |

---

## 5. Methods

> Key parameters are listed here for reproducibility. For the scientific rationale behind each step, see Section 4.

### 5.1 Read Quality Control
- **Tool:** fastp v0.23 (`-w 8`)
- Adapter trimming, quality filtering, and per-read statistics
- Outputs: filtered paired-end FASTQ files + HTML/JSON QC report
- **MultiQC v1.35** used to aggregate QC metrics

### 5.2 Genome Alignment
- **Reference:** *Cricetulus griseus* CriGri-PICR assembly (GCF_003668045.1)
- **Aligner:** BWA-MEM (`-t 20`)
- **Sorting:** SAMtools sort (`-@ 8 -m 8G`)
- **Indexing:** SAMtools index (`-@ 8`)
- Read-group tag added: `@RG ID:wt1 SM:wt1 PL:ILLUMINA`

### 5.3 DHFR Locus Depth Analysis
- SAMtools depth over DHFR gene coordinates (GFF annotation: NW_020822461.1:37,643,639–37,667,418)
- Flanking control region: ±500 kb around DHFR locus
- Depth ratio (DHFR / flanking) used for CHO sub-strain classification

### 5.4 Transgene Detection (Host-Subtraction Strategy)
- **Step 1 — Unmapped read extraction:**
  - Both-ends unmapped pairs: `samtools view -f 12 -F 256 -F 2048`
  - Singleton unmapped reads: `samtools view -f 4 -F 8 -F 256 -F 2048`
- **Step 2 — De novo assembly:** MEGAHIT v1.2.9 (`-t 16 -m 0.5 --min-contig-len 200`)
- **Step 3 — BLAST homology search:** BLASTn remote against NCBI nt database
  (`-evalue 1e-10 -max_target_seqs 5`, contigs ≥1,000 bp submitted)

---

## 6. Results

### 6.1 Read Quality Control

| Metric | Value |
| :--- | :--- |
| Total read pairs (post-filter) | 232,926,798 |
| Q30 rate | See `wt1_full_fastp.html` |
| Adapter content | Detected and trimmed |

Full QC report: `qc/wt1_full_fastp.html` · Aggregated report: `qc/multiqc_report.html`

### 6.2 Species Identification

| Metric | Value | Interpretation |
| :--- | :---: | :--- |
| Total mapped reads (primary) | 446,411,806 | |
| **Primary mapping rate** | **95.83%** | **Confirmed *Cricetulus griseus* (CHO)** |
| Properly paired | 92.17% | Normal paired-end library |
| Supplementary alignments | 1,632,515 | Consistent with repetitive CHO genome |

A primary mapping rate of 95.83% to the CriGri-PICR reference genome unambiguously confirms the sample is derived from *Cricetulus griseus* (Chinese hamster ovary, CHO).

Full alignment statistics: `species_strain/cho_flagstat.txt`

### 6.3 CHO Sub-Strain Identification

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

### 6.4 Transgene Detection

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

Full BLAST results: `transgene/blast_results.txt`  
Assembly FASTA (all contigs): `transgene/transgene_all_contigs.fa`

---

## 7. Conclusions

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

## 8. Deliverable Files

```
custom_research_report_20260626/
├── WGS_Analysis_Report_wt1_0626.md     ← This report
├── qc/
│   ├── wt1_full_fastp.html             ← Full-dataset fastp QC (interactive HTML)
│   ├── wt1_full_fastp.json             ← fastp QC metrics (machine-readable)
│   └── multiqc_report.html             ← MultiQC aggregated QC report
├── species_strain/
│   ├── cho_flagstat.txt                ← SAMtools flagstat — alignment statistics
│   └── dhfr_depth.txt                  ← DHFR and flanking region sequencing depth
└── transgene/
    ├── blast_results.txt               ← BLASTn tabular results (all hits, ≥1,000 bp contigs)
    ├── transgene_contigs_1000bp.fa     ← Assembled contigs ≥1,000 bp (BLAST input)
    └── transgene_all_contigs.fa        ← All assembled contigs ≥200 bp (complete assembly)
```

---

*Report prepared by:*  
**Zhen Gao, PhD**  
**Principal Bioinformatics Scientist**  
**Athenomics**  
*June 26, 2026*
