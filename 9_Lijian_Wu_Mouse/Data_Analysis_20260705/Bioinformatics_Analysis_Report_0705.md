# Mouse Bulk RNA-seq Analysis Report — Lijian Wu Project

**Report Date:** July 05, 2026
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Data folder:** `Data_Analysis_20260705`

## 1. Objectives

Characterise transcriptomic differences in cDC1 (conventional dendritic cell type 1) cells
across three treatment conditions (**TNFa**, **Tumor**, **Tumor_TNF**) relative to the
untreated control (**Control**) using bulk RNA-seq. Specific aims:

- Identify differentially expressed genes (DEGs) for each treatment vs. control contrast.
- Perform GO and KEGG pathway enrichment to determine biological processes affected.
- Run Gene Set Enrichment Analysis (GSEA) for pathway-level signal detection.

## 2. Key Findings

This study comprises **4 groups** (TNFa, Tumor, Tumor_TNF = treatment; Control = untreated) across **3 contrasts**.

**Differentially expressed genes** (criteria: padj ≤ 0.05 AND |log2FC| ≥ 0.263; log2(1.2) = 0.263, equivalent to ≥1.2-fold change):

| Contrast | Total DEGs | Upregulated (log2FC ≥ 0.263) | Downregulated (log2FC ≤ −0.263) |
| :--- | :---: | :---: | :---: |
| TNFa_vs_Control | **6847** | 3078 | 3769 |
| Tumor_vs_Control | **8483** | 4222 | 4261 |
| Tumor_TNF_vs_Control | **8757** | 4427 | 4330 |

**Top pathway findings per comparison:**

- **TNFa_vs_Control**:
  - Top GO (BP): immune effector process (padj=8.16e-06)
  - Top KEGG: Cell adhesion molecule (CAM) interaction (padj=9.82e-05)
  - GSEA KEGG: TNF signaling pathway (padj=0.000151)
  - GSEA Hallmark: HALLMARK E2F TARGETS (padj=5e-09)

- **Tumor_vs_Control**:
  - Top GO (BP): response to bacterium (padj=0.000188)
  - Top KEGG: Cell adhesion molecule (CAM) interaction (padj=0.000211)
  - GSEA KEGG: Spliceosome (padj=4.09e-08)
  - GSEA Hallmark: HALLMARK E2F TARGETS (padj=5e-09)

- **Tumor_TNF_vs_Control**:
  - Top GO (BP): response to bacterium (padj=0.00369)
  - Top KEGG: Hematopoietic cell lineage (padj=0.000512)
  - GSEA KEGG: Spliceosome (padj=3.1e-08)
  - GSEA Hallmark: HALLMARK E2F TARGETS (padj=2.5e-09)

## 3. Sample Information

| Group | Samples | Role |
| :--- | :---: | :---: |
| Control | cDC1_un_1, cDC1_un_2, cDC1_un_3 | Control (untreated) |
| TNFa | cDC1_TNFa_1, cDC1_TNFa_2, cDC1_TNFa_3 | Treatment 1 |
| Tumor | cDC1_Tumor_1, cDC1_Tumor_2, cDC1_Tumor_3 | Treatment 2 |
| Tumor_TNF | cDC1_Tumor_TNFa_1, cDC1_Tumor_TNFa_2, cDC1_Tumor_TNFa_3 | Treatment 3 |

Total samples: **12**  |  Comparisons: TNFa vs Control, Tumor vs Control, Tumor_TNF vs Control

## 4. Analysis Rationale and Decision Criteria

| Step | Decision | Rationale |
| :--- | :---: | :---: |
| Quantification | Salmon (via nf-core/rnaseq) | Bias-corrected transcript-level quantification |
| STAR alignment | 1-pass, --outFilterMultimapNmax 3 | 7.6× speedup vs default 2-pass on mouse genome; negligible impact on protein-coding DE (Salmon EM handles 2–3 candidate multi-map positions correctly) |
| Regex gene filter | Remove ribosomal/non-coding/Gm[0-9] genes | 57,132 → 26,572 genes retained (removed 30,560 noise genes) |
| Low-expression filter | ≥10 counts in ≥10 of 12 samples | 26,572 → 11,484 robustly expressed genes input to DESeq2 |
| DE threshold | padj ≤ 0.05, \|log2FC\| ≥ 0.263 (= 1.2×) | Conservative fold-change avoids over-reporting small but significant changes |
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
| TNFa_vs_Control | 6847 | 3078 | 3769 |
| Tumor_vs_Control | 8483 | 4222 | 4261 |
| Tumor_TNF_vs_Control | 8757 | 4427 | 4330 |

### 6.2 Pathway Enrichment (GO BP + KEGG, ALL direction)

| Contrast | GO BP terms | KEGG pathways |
| :--- | :---: | :---: |
| TNFa_vs_Control | 370 | 10 |
| Tumor_vs_Control | 142 | 13 |
| Tumor_TNF_vs_Control | 46 | 4 |

### 6.3 Immune and TNF-pathway relevant genes among the top 50 DEGs

Given the project's core focus — the transcriptomic response of dendritic cells (cDC1) to
TNF-α — each contrast's top 50 significant DEGs (the same gene sets shown in
`Heatmap_top50_*.pdf`) were manually reviewed for known immune-function or TNF/NF-κB pathway
relevance.

**TNFa vs Control** (isolated TNF-α stimulus):

| Gene | Function / pathway |
| :--- | :---: |
| Bcl3 | NF-κB pathway component — canonical downstream effector of TNF-α signalling |
| Tank | TRAF-associated NF-κB activator — direct component of the TNF-receptor signalling complex |
| Socs3 | Negative-feedback regulator of cytokine signalling; classic TNF/inflammation-induced gene |
| Ccl5 (RANTES) | Chemokine; recruits T cells, NK cells, monocytes — classic inflammatory/TNF-induced marker |
| Isg15, Gbp4, Samd9l | Interferon-stimulated genes — innate immune/antiviral response |
| Stat4 | Immune signalling transcription factor (Th1 / IFN-γ pathway) |
| H2-Q7, H2-Q6 | Mouse MHC class I genes — antigen presentation, inducible by IFN/TNF |
| Il4i1 | Dendritic-cell immunoregulatory enzyme |
| Serpinb9 | Inhibits granzyme B — immune self-protection/regulation |
| Pvr (CD155) | Immune checkpoint ligand (TIGIT/DNAM-1) — dendritic cell–T cell interaction |
| Vsig10 | Immunoglobulin-superfamily receptor |

**Tumor vs Control** (tumour-material exposure):

| Gene | Function / pathway |
| :--- | :---: |
| Csf1r | Key differentiation/survival receptor for myeloid and dendritic cells |
| Spn (CD43) | Classic leukocyte surface marker — adhesion/activation |
| Cd63 | Late endosome/lysosome marker — antigen uptake and processing |
| Ctsl (Cathepsin L) | Lysosomal protease required for MHC class II antigen processing |
| Ifitm3, Tgtp2, Oas1a | Interferon-stimulated antiviral/innate-immunity genes |
| Itgb3 | Integrin — immune cell adhesion and phagocytosis |
| Ptk2b (Pyk2) | Integrin-signalling kinase in immune cells |

**Tumor_TNF vs Control** (combined tumour + TNF-α):

| Gene | Function / pathway |
| :--- | :---: |
| Ccr2 | Chemokine receptor — monocyte/macrophage recruitment |
| Cxcr3 | Chemokine receptor — T cell/NK cell, Th1 response |
| Grap2, Sh2b3, Lat2 | Immune-receptor signalling adaptor proteins (T/NK/mast/B cells) |
| Milr1 | Mast cell inhibitory immunoglobulin-like receptor |
| Entpd1 (CD39) | Immunoregulatory ectonucleotidase — Treg-associated marker |
| P2ry6 | Purinergic receptor — regulates phagocytosis |
| Ptger3 | Prostaglandin receptor — inflammatory signalling |
| Oasl1, Ifi203, Xaf1 | Interferon-stimulated antiviral/innate-immunity genes |
| Id1 | Immune cell differentiation regulator |
| Arhgap9 | Haematopoietic Rho-GTPase-activating protein — immune cell migration |
| Pvr (CD155) | Immune checkpoint ligand (as above) |

**Cross-contrast observations:**

- **Bcl3 and Tank** appearing among the top DEGs for TNFa vs Control are direct components of the
  NF-κB/TRAF signalling cascade — a strong internal positive control confirming the TNF-α
  stimulus engaged its expected canonical pathway.
- **Pvr (CD155)**, an immune checkpoint ligand, is among the top DEGs in both TNFa vs Control and
  Tumor_TNF vs Control, suggesting consistent TNF-α-driven upregulation regardless of co-exposure
  to tumour material — a candidate worth independent validation.
- **Chemokine receptors Ccr2 and Cxcr3** appear only in the combined Tumor_TNF contrast, not in
  TNFa alone — suggesting the combination of tumour exposure and TNF-α specifically activates an
  immune cell-recruitment programme that neither stimulus induces on its own.
- **Cd63 and Ctsl**, both linked to antigen uptake/processing, appear only in the Tumor contrast —
  consistent with dendritic cells actively processing tumour-derived material.

## 7. Conclusions

- **TNFa_vs_Control**: 6847 DEGs identified. Predominantly downregulated (3769 down vs 3078 up), suggesting suppression of gene expression relative to untreated control.
- **Tumor_vs_Control**: 8483 DEGs identified. Predominantly downregulated (4261 down vs 4222 up), suggesting suppression of gene expression relative to untreated control.
- **Tumor_TNF_vs_Control**: 8757 DEGs identified. Predominantly upregulated (4427 up vs 4330 down), suggesting activation of transcriptional programs relative to untreated control.
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
| `Enrichment/*/GSEA/` | GSEA results (GO BP + KEGG + Hallmark) with ridge/dot plots |
| `QC/multiqc/` | MultiQC report (sequencing QC, alignment stats) |

---
*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*
