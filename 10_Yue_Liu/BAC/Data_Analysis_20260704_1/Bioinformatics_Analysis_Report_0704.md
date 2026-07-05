# Human Bulk RNA-seq Analysis Report — Yue Liu Project (MOLM13 X-ray Radiation Time Course)

**Report Date:** July 04, 2026
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Data folder:** `Data_Analysis_20260704`

## 1. Objectives

Characterise the transcriptomic response of MOLM13 cells to X-ray radiation over a
**Control → 4h → 8h → 16h** time course, integrating a client-supplied legacy dataset
(Control, 4h) with the current sequencing batch (8h, 16h). Specific aims:

- Identify differentially expressed genes (DEGs) at each post-radiation time point vs. untreated Control.
- Perform GO and KEGG pathway enrichment to determine biological processes affected at each time point.
- Run Gene Set Enrichment Analysis (GSEA) for pathway-level signal detection.
- Transparently flag the cross-batch statistical limitation described in Section 4.

## 2. Key Findings

This study integrates **2 sequencing batches** into **4 groups** (Control, 4h, 8h, 16h) across **3 contrasts**, all vs. Control.

**Differentially expressed genes** (criteria: padj ≤ 0.05 AND |log2FC| ≥ 0.263; log2(1.2) = 0.263, equivalent to ≥1.2-fold change):

| Contrast | Total DEGs | Upregulated (log2FC ≥ 0.263) | Downregulated (log2FC ≤ −0.263) |
| :--- | :---: | :---: | :---: |
| 4h_vs_Control | **249** | 139 | 110 |
| 8h_vs_Control | **11980** | 6241 | 5739 |
| 16h_vs_Control | **12054** | 6235 | 5819 |

**Note the sharp jump in DEG count from 4h (same-batch comparison) to 8h/16h (cross-batch comparison)** —
this is a direct symptom of the batch confound described in Section 4, not necessarily biology alone.

**Top pathway findings per comparison:**

- **4h_vs_Control**:
  - GSEA KEGG: Non-alcoholic fatty liver disease (padj=0.0282)
  - GSEA Hallmark: HALLMARK E2F TARGETS (padj=0.000113)

- **8h_vs_Control**:
  - Top GO (BP): proton transmembrane transport (padj=0.0227)
  - Top KEGG: Diabetic cardiomyopathy (padj=0.000162)
  - GSEA KEGG: Ribosome (padj=1.54e-08)
  - GSEA Hallmark: HALLMARK MYC TARGETS V1 (padj=2.5e-09)

- **16h_vs_Control**:
  - Top KEGG: Huntington disease (padj=5.3e-05)
  - GSEA KEGG: Ribosome (padj=1.03e-08)
  - GSEA Hallmark: HALLMARK MYC TARGETS V1 (padj=2.5e-09)

## 3. Sample Information

| Group | Samples | Sequencing Batch | Role |
| :--- | :---: | :---: | :---: |
| Control | Control_1, Control_2, Control_3 | Legacy batch (client-supplied count matrix) | Untreated baseline |
| 4h | 4h_1, 4h_2, 4h_3 (orig. X-Ray-1/2/3) | Legacy batch (client-supplied count matrix) | 4h post X-ray |
| 8h | 8h_1, 8h_2, 8h_3 (orig. "A" group, corrected) | Current NovaSeq X Plus batch | 8h post X-ray |
| 16h | 16h_1, 16h_2, 16h_3 | Current NovaSeq X Plus batch | 16h post X-ray |

Total samples: **12**  |  Comparisons: 4h vs Control, 8h vs Control, 16h vs Control

**Sample relabeling correction (2026-07-04):** the sheet originally labeled the 8h group as "A" —
the client confirmed this was a documentation error; the samples are MOLM13 cells 8h post X-ray
radiation. Raw FASTQ files on disk retain their original names (A_1/A_2/A_3); only the group
label was corrected, no realignment was needed.

## 4. Analysis Rationale and Decision Criteria

### 4.1 ⚠️ Cross-batch statistical limitation (read before interpreting 8h/16h results)

Control and 4h come from a **client-supplied legacy sequencing batch** (pre-quantified count
matrix, quantification method not independently verified by this analysis). 8h and 16h come
from the **current NovaSeq X Plus batch**, quantified with Salmon via nf-core/rnaseq in this project.

**No experimental condition was sequenced in both batches** — there is no way to separate a true
"batch effect" (different sequencing run, possibly different library prep/quantification pipeline)
from the true biological radiation-response effect for the 8h and 16h vs Control contrasts. A
`~Batch + Group` DESeq2 design is **not statistically estimable** here because Batch and Group are
fully aliased (Control/4h only ever appear in the old batch; 8h/16h only ever appear in the new batch).

Per client's explicit decision, Control is used as the common reference for all three contrasts to
produce a unified time-course view. The **4h vs Control** contrast is a clean, same-batch comparison
and its result (see Section 6.1) can be interpreted with normal confidence. The **8h vs Control** and
**16h vs Control** contrasts should be read as *exploratory / batch-confounded* — the very large jump
in DEG count relative to 4h (
249 at 4h vs 11980 at 8h vs 12054 at 16h
) is consistent with batch-driven variance rather than a biologically plausible ~48× increase in
radiation response between 4h and 8h. We recommend treating 8h/16h DEG lists as candidate genes
requiring independent validation (e.g. qPCR), not as confirmed radiation-response findings.

### 4.2 Gene filtering — methodology consistency check vs. client's legacy file

The client's legacy count matrix (`counts-filtered protein coding original file.xlsx`) was pre-filtered to **20,065 protein-coding genes**.
Cross-checking by base Ensembl ID against our GENCODE v45 reference: **20,056 / 20,065 (99.96%)** of the client's genes matched;
**9 genes** were absent from our v45 annotation entirely (IDs in the ENSG00000293xxx range — likely added in a newer Ensembl/GENCODE release than v45, i.e. the client's original pipeline used a slightly newer annotation).
Of the 20,056 genes in common, our own annotation independently classified **20,047 (99.96%) as protein_coding** — confirming our biotype filter is highly consistent with the client's pre-filtering methodology (the remaining genes are classified as lncRNA/pseudogene under GENCODE v45, likely due to biotype reclassification between annotation versions — a known, expected phenomenon, not an error).

### 4.3 Decision criteria

| Step | Decision | Rationale |
| :--- | :---: | :---: |
| Quantification | Salmon (new batch, via nf-core/rnaseq); client-supplied matrix (legacy batch) | Bias-corrected transcript-level quantification for the new batch |
| Gene ID matching | Base Ensembl ID (version suffix stripped) | 99.96% of client's genes matched our annotation this way |
| Biotype gene filter | Keep gene_type == protein_coding (GENCODE v45 annotation) | Matches client's own pre-filtering convention (see 4.2) |
| Low-expression filter | **Not applied** | Matches client's legacy file, which retains all-zero genes (4,746/20,065 such genes present); DESeq2's built-in independent filtering still applies at the padj-correction step |
| DE threshold | padj ≤ 0.05, \|log2FC\| ≥ 0.263 (= 1.2×) | Conservative fold-change avoids over-reporting small but significant changes |
| LFC shrinkage | ashr | Appropriate for unbalanced designs and small sample sizes |
| Cross-batch design | Common Control reference for all 3 contrasts (client decision) | Enables a unified time-course view; batch/time confound disclosed above for 8h/16h |

## 5. Methods

| Tool | Version | Parameters |
| :--- | :---: | :---: |
| nf-core/rnaseq (new batch: 8h/16h) | 3.15.1 | --aligner star_salmon |
| Legacy batch (Control/4h) | client-supplied | pre-quantified count matrix, method not independently verified |
| DESeq2 | R package | design = ~Group, lfcShrink type = ashr |
| clusterProfiler | R package | GO ORA (BP/MF/CC), KEGG ORA, GSEA |
| org.Hs.eg.db | R package | Human gene ID mapping |
| msigdbr | R package | Hallmark gene sets (Homo sapiens) |
| Reference genome (new batch) | GRCh38 / GENCODE v45 | — |

## 6. Results

### 6.1 Differential Expression

| Contrast | Total DEGs | Upregulated | Downregulated | Batch status |
| :--- | :---: | :---: | :---: | :---: |
| 4h_vs_Control | 249 | 139 | 110 | Same-batch (reliable) |
| 8h_vs_Control | 11980 | 6241 | 5739 | Cross-batch (exploratory — see 4.1) |
| 16h_vs_Control | 12054 | 6235 | 5819 | Cross-batch (exploratory — see 4.1) |

### 6.2 Pathway Enrichment (GO BP + KEGG, ALL direction)

| Contrast | GO BP terms | KEGG pathways |
| :--- | :---: | :---: |
| 4h_vs_Control | 0 | 0 |
| 8h_vs_Control | 8 | 8 |
| 16h_vs_Control | 0 | 5 |

## 7. Conclusions

- **4h_vs_Control**: 249 DEGs identified. Predominantly upregulated (139 up vs 110 down).
- **8h_vs_Control** (interpret cautiously — cross-batch comparison, see Section 4.1): 11980 DEGs identified. Predominantly upregulated (6241 up vs 5739 down).
- **16h_vs_Control** (interpret cautiously — cross-batch comparison, see Section 4.1): 12054 DEGs identified. Predominantly upregulated (6235 up vs 5819 down).
- The 4h vs Control result is the only contrast free of batch confound and can anchor confident conclusions
  about the early (4h) transcriptomic response to X-ray radiation in MOLM13 cells.
- The 8h/16h vs Control DEG lists are useful as a candidate/exploratory resource but should be validated
  independently (e.g. qPCR on a shortlist of genes) before being treated as confirmed findings, given the
  batch/time confound.
- Pathway enrichment and GSEA results are available in the `Enrichment/` directory for detailed biological interpretation.

## 8. Deliverable Files

| File / Folder | Contents |
| :--- | :---: |
| `DE_PCA_Results/DEG_*.csv` | Full DEG tables (all genes, with log2FC, padj, raw counts) |
| `DE_PCA_Results/PCA.pdf` | PCA plot, colored by Group and shaped by Batch (visualizes the batch effect) |
| `DE_PCA_Results/Volcano_*.png` | Volcano plots per contrast |
| `DE_PCA_Results/Heatmap_top50_*.pdf` | Heatmaps of top 50 DEGs per contrast |
| `Reads/New_batch_8h_16h_gene_counts.tsv` | Raw count matrix, current NovaSeq batch (8h/16h) |
| `Reads/New_batch_8h_16h_gene_tpm.tsv` | TPM matrix, current NovaSeq batch (8h/16h) |
| `Reads/Old_batch_Control_4h_gene_counts.xlsx` | Client-supplied legacy count matrix (Control/4h), copied as-is |
| `human_Gene_annotation_*.xlsx` | Full human gene annotation with GO/KEGG/UniProt (GENCODE v45) |
| `Enrichment/*/GO/` | GO ORA results (BP/MF/CC, Up/Down/ALL) with dot plots |
| `Enrichment/*/KEGG/` | KEGG ORA results with dot plots |
| `Enrichment/*/GSEA/` | GSEA results (GO BP + KEGG + Hallmark) with ridge/dot plots |
| `QC/multiqc/` | MultiQC report — covers current NovaSeq batch (8h/16h) only; no QC available for the legacy batch |

---
*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*
