# Human Bulk RNA-seq Analysis Report — Yue Liu Project

**Report Date:** July 04, 2026
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Data folder:** `Data_Analysis_20260704`

## 1. Objectives

Characterise transcriptomic differences between experimental group **16** and
control group **A** using bulk RNA-seq. Specific aims:

- Identify differentially expressed genes (DEGs) for the 16 vs A contrast.
- Perform GO and KEGG pathway enrichment to determine biological processes affected.
- Run Gene Set Enrichment Analysis (GSEA) for pathway-level signal detection.

## 2. Key Findings

This study comprises **2 groups** (16 = treatment; A = control) across **1 contrast**.

**Differentially expressed genes** (criteria: padj ≤ 0.05 AND |log2FC| ≥ 0.263; log2(1.2) = 0.263, equivalent to ≥1.2-fold change):

| Contrast | Total DEGs | Upregulated (log2FC ≥ 0.263) | Downregulated (log2FC ≤ −0.263) |
| :--- | :---: | :---: | :---: |
| 16_vs_A | **3994** | 1906 | 2088 |

**Top pathway findings:**

- **16_vs_A**:
  - Top GO (BP): ribosome biogenesis (padj=6.53e-22)
  - Top KEGG: Ribosome biogenesis in eukaryotes (padj=4.91e-08)
  - GSEA KEGG: Ribosome (padj=3.25e-08)
  - GSEA Hallmark: HALLMARK MITOTIC SPINDLE (padj=1.67e-09)

## 3. Sample Information

| Group | Samples | Role |
| :--- | :---: | :---: |
| A | A_1, A_2, A_3 | Control |
| 16 | 16_1, 16_2, 16_3 | Treatment |

Total samples: **6**  |  Comparison: 16 vs A

## 4. Analysis Rationale and Decision Criteria

| Step | Decision | Rationale |
| :--- | :---: | :---: |
| Quantification | Salmon (via nf-core/rnaseq) | Bias-corrected transcript-level quantification |
| Biotype gene filter | Keep gene_type == protein_coding (GENCODE annotation) | 63,187 → 20,049 genes retained (removed 43,138 non-coding/pseudogenes) |
| Low-expression filter | ≥10 counts in ≥4 of 6 samples | 20,049 → 11,934 robustly expressed genes input to DESeq2 |
| DE threshold | padj ≤ 0.05, \|log2FC\| ≥ 0.263 (= 1.2×) | Conservative fold-change avoids over-reporting small but significant changes |
| LFC shrinkage | ashr | Appropriate for unbalanced designs and small sample sizes |

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
| 16_vs_A | 3994 | 1906 | 2088 |

### 6.2 Pathway Enrichment (GO BP + KEGG, ALL direction)

| Contrast | GO BP terms | KEGG pathways |
| :--- | :---: | :---: |
| 16_vs_A | 32 | 1 |

## 7. Conclusions

- **16_vs_A**: 3994 DEGs identified. Predominantly downregulated (2088 down vs 1906 up), suggesting suppression of gene expression in group 16 relative to control A.
- Pathway enrichment and GSEA results are available in the `Enrichment/` directory for detailed biological interpretation.

## 8. Deliverable Files

| File / Folder | Contents |
| :--- | :---: |
| `DE_PCA_Results/DEG_*.csv` | Full DEG table (all genes, with log2FC, padj, raw counts) |
| `DE_PCA_Results/PCA.pdf` | PCA plot |
| `DE_PCA_Results/Volcano_*.png` | Volcano plot |
| `DE_PCA_Results/Heatmap_top50_*.pdf` | Heatmap of top 50 DEGs |
| `Reads/All_sample_gene_counts.tsv` | Raw count matrix |
| `Reads/All_sample_gene_tpm.tsv` | TPM matrix |
| `human_Gene_annotation_*.xlsx` | Full human gene annotation with GO/KEGG/UniProt (GENCODE v45) |
| `Enrichment/*/GO/` | GO ORA results (BP/MF/CC, Up/Down/ALL) with dot plots |
| `Enrichment/*/KEGG/` | KEGG ORA results with dot plots |
| `Enrichment/*/GSEA/` | GSEA results (GO BP + KEGG + Hallmark) with ridge/dot plots |
| `QC/multiqc/` | MultiQC report (sequencing QC, alignment stats) |

---
*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*
