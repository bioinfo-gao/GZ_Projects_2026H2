# Human Bulk RNA-seq Analysis Report — Zhen Yan Project

**Report Date:** July 05, 2026
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Data folder:** `Data_Analysis_20260705`

## 1. Objectives

Characterise transcriptomic differences between two experimental groups
(**CHD3**, **ASHL**) and the control group (**CTRL**) using bulk RNA-seq. Specific aims:

- Identify differentially expressed genes (DEGs) for each treatment vs. control contrast.
- Perform GO and KEGG pathway enrichment to determine biological processes affected.
- Run Gene Set Enrichment Analysis (GSEA) for pathway-level signal detection.

## 2. Key Findings

This study comprises **3 groups** (CHD3, ASHL = treatment; CTRL = control) across **2 contrasts**.

**Differentially expressed genes** (criteria: padj ≤ 0.05 AND |log2FC| ≥ 0.263; log2(1.2) = 0.263, equivalent to ≥1.2-fold change):

| Contrast | Total DEGs | Upregulated (log2FC ≥ 0.263) | Downregulated (log2FC ≤ −0.263) |
| :--- | :---: | :---: | :---: |
| CHD3_vs_CTRL | **3571** | 1436 | 2135 |
| ASHL_vs_CTRL | **2406** | 1289 | 1117 |

**Top pathway findings per comparison:**

- **CHD3_vs_CTRL**:
  - Top GO (BP): regulation of transmembrane transport (padj=0.00435)
  - Top KEGG: p53 signaling pathway (padj=0.000836)
  - GSEA KEGG: ECM-receptor interaction (padj=0.0358)

- **ASHL_vs_CTRL**:
  - Top GO (BP): blood vessel morphogenesis (padj=6.36e-09)
  - Top KEGG: Cytoskeleton in muscle cells (padj=0.000101)
  - GSEA KEGG: Cytoskeleton in muscle cells (padj=0.00135)
  - GSEA Hallmark: HALLMARK COAGULATION (padj=0.000109)

## 3. Sample Information

| Group | Samples | Role |
| :--- | :---: | :---: |
| CTRL | CTRL_1, CTRL_2, CTRL_3, CTRL_4 | Control |
| CHD3 | CHD3_1, CHD3_2, CHD3_3, CHD3_4 | Treatment 1 |
| ASHL | ASHL_1, ASHL_2, ASHL_3 | Treatment 2 |

Total samples: **11**  |  Comparisons: CHD3 vs CTRL, ASHL vs CTRL

## 4. Analysis Rationale and Decision Criteria

| Step | Decision | Rationale |
| :--- | :---: | :---: |
| Quantification | Salmon (via nf-core/rnaseq) | Bias-corrected transcript-level quantification |
| Biotype gene filter | Keep gene_type == protein_coding (GENCODE annotation) | 63,187 → 20,049 genes retained (removed 43,138 non-coding/pseudogenes) |
| Low-expression filter | ≥10 counts in ≥9 of 11 samples | 20,049 → 14,081 robustly expressed genes input to DESeq2 |
| DE threshold | padj ≤ 0.05, \|log2FC\| ≥ 0.263 (= 1.2×) | Conservative fold-change avoids over-reporting small but significant changes |
| LFC shrinkage | ashr | Appropriate for unbalanced designs and small sample sizes |
| Control designation | CTRL | Client note: "control vs other groups"; consistent with prior analysis for this client |

## 5. Methods

| Tool | Version | Parameters |
| :--- | :---: | :---: |
| nf-core/rnaseq | 3.15.1 | --aligner star_salmon |
| STAR | 2.7.x | default (2-pass) |
| Salmon | — | default |
| DESeq2 | R package | design = ~Group, lfcShrink type = ashr |
| clusterProfiler | R package | GO ORA (BP/MF/CC), KEGG ORA, GSEA |
| org.Hs.eg.db | R package | Human gene ID mapping |
| msigdbr | R package | Hallmark gene sets (Homo sapiens) |
| Reference genome | GRCh38 / GENCODE v45 | — |

## 6. Results

### 6.1 Differential Expression

| Contrast | Total DEGs | Upregulated | Downregulated |
| :--- | :---: | :---: | :---: |
| CHD3_vs_CTRL | 3571 | 1436 | 2135 |
| ASHL_vs_CTRL | 2406 | 1289 | 1117 |

### 6.2 Pathway Enrichment (GO BP + KEGG, ALL direction)

| Contrast | GO BP terms | KEGG pathways |
| :--- | :---: | :---: |
| CHD3_vs_CTRL | 27 | 1 |
| ASHL_vs_CTRL | 427 | 11 |

## 7. Conclusions

- **CHD3_vs_CTRL**: 3571 DEGs identified. Predominantly downregulated (2135 down vs 1436 up), suggesting suppression of gene expression relative to control.
- **ASHL_vs_CTRL**: 2406 DEGs identified. Predominantly upregulated (1289 up vs 1117 down), suggesting activation of transcriptional programs relative to control.
- Pathway enrichment and GSEA results are available in the `Enrichment/` directory for detailed biological interpretation.

## 8. Deliverable Files

| File / Folder | Contents |
| :--- | :---: |
| `DE_PCA_Results/DEG_*.csv` | Full DEG tables (all genes, with log2FC, padj, raw counts) |
| `DE_PCA_Results/PCA.pdf` | PCA plot |
| `DE_PCA_Results/Volcano_*.png` | Volcano plots per contrast |
| `DE_PCA_Results/Heatmap_top50_*.pdf` | Heatmaps of top 50 DEGs per contrast |
| `Reads/All_sample_gene_counts.tsv` | Raw count matrix |
| `Reads/All_sample_gene_tpm.tsv` | TPM matrix |
| `human_Gene_annotation_*.xlsx` | Full human gene annotation with GO/KEGG/UniProt (GENCODE v45) |
| `Enrichment/*/GO/` | GO ORA results (BP/MF/CC, Up/Down/ALL) with dot plots |
| `Enrichment/*/KEGG/` | KEGG ORA results with dot plots |
| `Enrichment/*/GSEA/` | GSEA results (GO BP + KEGG + Hallmark) with ridge/dot plots |
| `QC/multiqc/` | MultiQC report (sequencing QC, alignment stats) |

---
*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*
