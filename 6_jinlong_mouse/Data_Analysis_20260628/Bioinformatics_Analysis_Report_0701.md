# Mouse Bulk RNA-seq Analysis Report — Jinlong Project

**Report Date:** July 01, 2026
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Data folder:** `Data_Analysis_20260628`

## 1. Objectives

Characterise transcriptomic differences between three experimental groups (G1, G2, G3)
and the control group (G4) using bulk RNA-seq. Specific aims:

- Identify differentially expressed genes (DEGs) for each treatment vs. control contrast.
- Perform GO and KEGG pathway enrichment to determine biological processes affected.
- Run Gene Set Enrichment Analysis (GSEA) for pathway-level signal detection.
- Assess whether stem cell marker genes are significantly altered.
- Characterise cell differentiation and growth pathway activity across treatment groups.
- Evaluate the Notch signalling pathway (receptors, ligands, effectors, and target genes).

## 2. Key Findings

This study comprises **4 groups** (G1, G2, G3 = treatment; G4 = control) across **3 contrasts**.

**Differentially expressed genes** (criteria: padj ≤ 0.05 AND |log2FC| ≥ 0.263; log2(1.2) = 0.263, equivalent to ≥1.2-fold change):

| Contrast | Total DEGs | Upregulated (log2FC ≥ 0.263) | Downregulated (log2FC ≤ −0.263) |
| :--- | :---: | :---: | :---: |
| G1_vs_G4 | **1340** | 857 | 483 |
| G2_vs_G4 | **279** | 112 | 167 |
| G3_vs_G4 | **1146** | 602 | 544 |

**Top pathway findings per comparison:**

- **G1_vs_G4**:
  - Top GO (BP): oxidative phosphorylation (padj=8.79e-45)
  - Top KEGG: Oxidative phosphorylation (padj=1.97e-46)
  - GSEA KEGG: Oxidative phosphorylation (padj=2.73e-09)
  - GSEA Hallmark: HALLMARK OXIDATIVE PHOSPHORYLATION (padj=7.14e-10)
  - Stem cell markers: 6 significant
  - Cell differentiation & growth markers: 5 significant; top MSigDB term: GOBP_BROWN_FAT_CELL_DIFFERENTIATION (padj=0.0426)
  - Notch pathway genes: 5 significant

- **G2_vs_G4**:
  - Top GO (BP): circulatory system process (padj=0.00278)
  - Top KEGG: Cytokine-cytokine receptor interaction (padj=0.00566)
  - GSEA KEGG: Cytoskeleton in muscle cells (padj=3.28e-08)
  - GSEA Hallmark: HALLMARK TNFA SIGNALING VIA NFKB (padj=6.04e-08)
  - Stem cell markers: 2 significant
  - Cell differentiation & growth markers: 1 significant
  - Notch pathway genes: 0 significant

- **G3_vs_G4**:
  - Top GO (BP): oxidative phosphorylation (padj=1.47e-37)
  - Top KEGG: Oxidative phosphorylation (padj=2.28e-41)
  - GSEA KEGG: Oxidative phosphorylation (padj=2.75e-09)
  - GSEA Hallmark: HALLMARK OXIDATIVE PHOSPHORYLATION (padj=1e-09)
  - Stem cell markers: 3 significant
  - Cell differentiation & growth markers: 8 significant; top MSigDB term: GOBP_MITOTIC_CELL_CYCLE_PROCESS (padj=1.52e-07)
  - Notch pathway genes: 3 significant

## 3. Sample Information

| Group | Samples | Role |
| :--- | :---: | :---: |
| G1 | J_902, J_912, J_896 | Treatment 1 |
| G2 | J_910, J_909, J_905 | Treatment 2 |
| G3 | J_904, J_897, J_899 | Treatment 3 |
| G4 | A, B, C | Control |

Total samples: **12**  |  Comparisons: G1 vs G4, G2 vs G4, G3 vs G4

## 4. Analysis Rationale and Decision Criteria

| Step | Decision | Rationale |
| :--- | :---: | :---: |
| Quantification | Salmon (via nf-core/rnaseq) | Bias-corrected transcript-level quantification |
| Regex gene filter | Remove ribo/noncoding/Gm[0-9] genes | 57,132 → 26,572 genes retained (removed 30,560 noise genes) |
| Low-expression filter | ≥10 counts in ≥10 of 12 samples | 26,572 → 14,090 robustly expressed genes input to DESeq2 |
| DE threshold | padj ≤ 0.05, |log2FC| ≥ 0.263 (= 1.2×) | Conservative fold-change avoids over-reporting small but significant changes |
| LFC shrinkage | ashr | Appropriate for unbalanced designs and small sample sizes |

## 5. Methods

| Tool | Version | Parameters |
| :--- | :---: | :---: |
| nf-core/rnaseq | 3.15.1 | --aligner star_salmon |
| STAR | 2.7.x | --twopassMode None --outFilterMultimapNmax 3 |
| Salmon | — | default |
| DESeq2 | R package | design = ~Group, lfcShrink type = ashr |
| clusterProfiler | R package | GO ORA (BP/MF/CC), KEGG ORA, GSEA |
| org.Mm.eg.db | R package | Mouse gene ID mapping |
| msigdbr | R package | Hallmark gene sets (MM) |
| Reference genome | GRCm39 / GENCODE M35 | — |

## 6. Results

### 6.1 Differential Expression

| Contrast | Total DEGs | Upregulated | Downregulated |
| :--- | :---: | :---: | :---: |
| G1_vs_G4 | 1340 | 857 | 483 |
| G2_vs_G4 | 279 | 112 | 167 |
| G3_vs_G4 | 1146 | 602 | 544 |

### 6.2 Pathway Enrichment (GO BP + KEGG, ALL direction)

| Contrast | GO BP terms | KEGG pathways |
| :--- | :---: | :---: |
| G1_vs_G4 | 121 | 27 |
| G2_vs_G4 | 3 | 1 |
| G3_vs_G4 | 147 | 16 |

### 6.3 Stem Cell Markers

Significant stem cell markers detected: **11**

| Gene | Category | log2FC | padj | Comparison |
| :--- | :---: | :---: | :---: | :---: |
| Kit | Hematopoietic_SC | 0.787 | 0.000111 | G3_vs_G4 |
| Id1 | General_Stemness | 1.055 | 0.000272 | G1_vs_G4 |
| Klf4 | Pluripotency | 0.66 | 0.00332 | G1_vs_G4 |
| Klf4 | Pluripotency | 0.584 | 0.00891 | G3_vs_G4 |
| Prom1 | Neural_SC | -0.526 | 0.0164 | G2_vs_G4 |
| Aldh1a1 | Mesenchymal_SC | 0.608 | 0.0191 | G1_vs_G4 |
| Aldh1a1 | General_Stemness | 0.608 | 0.0191 | G1_vs_G4 |
| Id1 | General_Stemness | 0.389 | 0.0366 | G2_vs_G4 |
| Itgb1 | Mesenchymal_SC | -0.374 | 0.0383 | G3_vs_G4 |
| Procr | Hematopoietic_SC | 0.33 | 0.0394 | G1_vs_G4 |
| Cd24a | General_Stemness | 0.465 | 0.0451 | G1_vs_G4 |

### 6.4 Cell Differentiation & Growth Markers

Significant cell differentiation & growth markers detected: **14**

| Gene | Category | log2FC | padj | Comparison |
| :--- | :---: | :---: | :---: | :---: |
| Vegfa | Growth_Factors | 0.858 | 0.00017 | G1_vs_G4 |
| Vegfa | Growth_Factors | 0.764 | 0.00073 | G3_vs_G4 |
| Klf4 | Differentiation_TF | 0.66 | 0.00332 | G1_vs_G4 |
| Ccnd1 | Proliferation | -0.757 | 0.00356 | G1_vs_G4 |
| Runx1 | Differentiation_TF | -0.563 | 0.00619 | G1_vs_G4 |
| Igf1r | Growth_Factors | -0.418 | 0.00731 | G1_vs_G4 |
| Klf4 | Differentiation_TF | 0.584 | 0.00891 | G3_vs_G4 |
| Cdk6 | Proliferation | -0.81 | 0.00893 | G3_vs_G4 |
| Igf1r | Growth_Factors | -0.382 | 0.0144 | G3_vs_G4 |
| Ccne1 | Proliferation | -0.67 | 0.0166 | G3_vs_G4 |
| Top2a | Proliferation | -0.504 | 0.0306 | G3_vs_G4 |
| Mki67 | Proliferation | -0.513 | 0.0317 | G3_vs_G4 |
| Zeb1 | EMT_Markers | -0.314 | 0.0459 | G2_vs_G4 |
| Fn1 | EMT_Markers | -0.438 | 0.0491 | G3_vs_G4 |

### 6.5 Notch Pathway Genes

Significant Notch pathway genes detected: **8**

| Gene | Category | log2FC | padj | Comparison |
| :--- | :---: | :---: | :---: | :---: |
| Cdkn1a | Target_Genes | -0.551 | 4.98e-07 | G1_vs_G4 |
| Cdkn1a | Target_Genes | -0.467 | 3.28e-05 | G3_vs_G4 |
| Ccnd1 | Target_Genes | -0.757 | 0.00356 | G1_vs_G4 |
| Hes1 | Target_Genes | 0.773 | 0.00408 | G1_vs_G4 |
| Hes1 | Target_Genes | 0.73 | 0.00917 | G3_vs_G4 |
| Hey1 | Target_Genes | 0.605 | 0.0158 | G1_vs_G4 |
| Mfng | Neg_Regulators | 0.475 | 0.0338 | G1_vs_G4 |
| Hey2 | Target_Genes | 0.497 | 0.0412 | G3_vs_G4 |

## 7. Conclusions

- **G1_vs_G4**: 1340 DEGs identified. Predominantly upregulated (857 up vs 483 down), suggesting activation of transcriptional programs in this group.
- **G2_vs_G4**: 279 DEGs identified. Predominantly downregulated (167 down vs 112 up), suggesting suppression of gene expression relative to control.
- **G3_vs_G4**: 1146 DEGs identified. Predominantly upregulated (602 up vs 544 down), suggesting activation of transcriptional programs in this group.
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
| `mouse_Gene_annotation_*.xlsx` | Full mouse gene annotation with GO/KEGG/UniProt (GENCODE M35) |
| `Enrichment_Standard/*/GO/` | GO ORA results (BP/MF/CC, Up/Down/ALL) with dot plots |
| `Enrichment_Standard/*/KEGG/` | KEGG ORA results with dot plots |
| `Enrichment_Standard/*/GSEA/` | GSEA results (KEGG + Hallmark) with ridge/dot plots |
| `Enrichment_Custom_Designed_Pathways/*/StemCell/` | Stem cell marker DE results and bar plots |
| `Enrichment_Custom_Designed_Pathways/*/CellDiff/` | Cell differentiation & growth marker DE results, bar plots, and MSigDB ORA |
| `Enrichment_Custom_Designed_Pathways/*/Notch/` | Notch pathway gene DE results, bar plots, KEGG, and MSigDB ORA |
| `Enrichment_Custom_Designed_Pathways/CellDiff_AllComparisons_Summary.csv` | Cross-comparison cell diff/growth marker summary |
| `Enrichment_Custom_Designed_Pathways/Notch_AllComparisons_Summary.csv` | Cross-comparison Notch pathway gene summary |
| `QC/multiqc/` | MultiQC report (sequencing QC, alignment stats) |

---
*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*
