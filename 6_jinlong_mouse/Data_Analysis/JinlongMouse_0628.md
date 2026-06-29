# Mouse Bulk RNA-seq Analysis Report — Jinlong Project

**Report Date:** June 28, 2026
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

## 2. Key Findings

- **G1_vs_G4**: 1340 DEGs (857 up, 483 down; padj ≤ 0.05, |log2FC| ≥ 0.263)
- **G2_vs_G4**: 279 DEGs (112 up, 167 down; padj ≤ 0.05, |log2FC| ≥ 0.263)
- **G3_vs_G4**: 1146 DEGs (602 up, 544 down; padj ≤ 0.05, |log2FC| ≥ 0.263)
- Top enriched GO biological process in G1_vs_G4: **oxidative phosphorylation** (padj = 8.79e-45)
- Stem cell marker analysis identified **11 significantly altered stem-cell-associated genes** across all comparisons.

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
| Aligner | STAR (1-pass, --outFilterMultimapNmax 3) | Mouse genome has ~14% multi-mapped reads; 2-pass and N=20 would add ~7× CPU cost with minimal accuracy gain |
| Quantification | Salmon (via nf-core/rnaseq) | Bias-corrected transcript-level quantification |
| Gene filtering | Regex (ribo/noncoding/Gm[0-9]) + low-count (≥10 in n−2 samples) | Removes noise genes; retains biologically informative signal |
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
| `Enrichment/*/GO/` | GO ORA results (BP/MF/CC, Up/Down/ALL) with dot plots |
| `Enrichment/*/KEGG/` | KEGG ORA results with dot plots |
| `Enrichment/*/GSEA/` | GSEA results (KEGG + Hallmark) with ridge/dot plots |
| `Enrichment/*/StemCell/` | Stem cell marker DE results and bar plots |
| `QC/multiqc/` | MultiQC report (sequencing QC, alignment stats) |

---
*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*
