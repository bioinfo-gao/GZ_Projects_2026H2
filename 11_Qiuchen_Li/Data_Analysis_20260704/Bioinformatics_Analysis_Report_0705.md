# Human Bulk RNA-seq Analysis Report — Qiuchen Li Project

**Report Date:** July 05, 2026
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Data folder:** `Data_Analysis_20260704`

## 1. Objectives

Characterise transcriptomic differences among three experimental groups
(**Mix**, **NT**, **A5BKO**) using bulk RNA-seq, via all three pairwise contrasts
(client-requested pairwise design — no single designated control). Specific aims:

- Identify differentially expressed genes (DEGs) for each of the 3 pairwise contrasts.
- Perform GO and KEGG pathway enrichment to determine biological processes affected.
- Run Gene Set Enrichment Analysis (GSEA) for pathway-level signal detection.

## 2. Key Findings

This study comprises **3 groups** (Mix, NT, A5BKO) across **3 pairwise contrasts**.

**Differentially expressed genes** (criteria: padj ≤ 0.05 AND |log2FC| ≥ 0.263; log2(1.2) = 0.263, equivalent to ≥1.2-fold change):

| Contrast | Total DEGs | Upregulated (log2FC ≥ 0.263) | Downregulated (log2FC ≤ −0.263) |
| :--- | :---: | :---: | :---: |
| Mix_vs_NT | **5478** | 2670 | 2808 |
| A5BKO_vs_NT | **18** | 12 | 6 |
| A5BKO_vs_Mix | **6212** | 3160 | 3052 |

**Top pathway findings per comparison:**

- **Mix_vs_NT**:
  - Top GO (BP): ribonucleoprotein complex biogenesis (padj=4.81e-09)
  - Top KEGG: Ribosome biogenesis in eukaryotes (padj=0.0072)
  - GSEA KEGG: Spliceosome (padj=1.57e-05)
  - GSEA Hallmark: HALLMARK MYC TARGETS V1 (padj=5e-09)

- **A5BKO_vs_NT**:
  - GSEA KEGG: Ribosome biogenesis in eukaryotes (padj=2.9e-08)
  - GSEA Hallmark: HALLMARK MYC TARGETS V2 (padj=1.67e-09)

- **A5BKO_vs_Mix**:
  - Top GO (BP): ribonucleoprotein complex biogenesis (padj=1.82e-09)
  - Top KEGG: Ribosome biogenesis in eukaryotes (padj=0.00354)
  - GSEA KEGG: Spliceosome (padj=5.32e-06)
  - GSEA Hallmark: HALLMARK MYC TARGETS V1 (padj=2.5e-09)

### 2.1 A recurring signal across all three comparisons: MYC target gene programs

All three pairwise contrasts independently identify a **MYC target gene signature** (MSigDB
Hallmark collection) as their single most significant GSEA hit — Mix vs NT and A5BKO vs Mix
both top out on **MYC_TARGETS_V1**, while A5BKO vs NT tops out on the related **MYC_TARGETS_V2**
set. This is a genuine, independently-computed result in each contrast (confirmed by inspecting
the underlying leading-edge gene lists, which differ between comparisons), not a reporting artefact.

**V1 vs V2 — two related but distinct gene sets:** `HALLMARK_MYC_TARGETS_V1` (200 genes) is the
broader, classic MYC target signature, drawing heavily on ribosomal-protein and translation/
ribosome-biogenesis genes that MYC is known to transcriptionally activate. `HALLMARK_MYC_TARGETS_V2`
(58 genes) is a smaller, more stringent signature enriched for core cell-cycle and DNA-replication
genes more directly bound by MYC. The two sets overlap partially but capture different facets of
MYC-driven transcriptional activity (broad translational output vs. core proliferation machinery).

**Why this is noteworthy:** the direction of enrichment is internally consistent and transitive
across all three contrasts — NT scores higher than Mix, A5BKO scores higher than NT, and A5BKO
scores higher than Mix for their respective MYC signatures. That ordering (**A5BKO > NT > Mix**)
is exactly what would be expected if the same underlying biological axis (MYC-driven proliferation/
translation activity) is the dominant driver of transcriptomic variation across all three groups,
rather than three unrelated, coincidental hits. This is consistent with the very large DEG counts
seen for Mix_vs_NT and A5BKO_vs_Mix (Section 2) alongside the far smaller A5BKO_vs_NT DEG count —
groups closer together on this MYC axis (A5BKO and NT) show fewer DEGs, while groups further apart
(Mix vs either) show many more. We'd recommend treating MYC pathway activity as a candidate primary
axis of biological variation in this experiment, worth confirming with a marker-gene qPCR panel or
by correlating a MYC target module score against the client's expected experimental design.

## 3. Sample Information

| Group | Samples | Role |
| :--- | :---: | :---: |
| Mix | Mix_1, Mix_2, Mix_3 | Experimental group |
| NT | NT_1, NT_2, NT_3 | Experimental group |
| A5BKO | A5BKO_1, A5BKO_2, A5BKO_3 | Experimental group |

Total samples: **9**  |  Comparisons: Mix vs NT, A5BKO vs NT, A5BKO vs Mix

## 4. Analysis Rationale and Decision Criteria

| Step | Decision | Rationale |
| :--- | :---: | :---: |
| Quantification | Salmon (via nf-core/rnaseq) | Bias-corrected transcript-level quantification |
| Biotype gene filter | Keep gene_type == protein_coding (GENCODE annotation) | 63,187 → 20,049 genes retained (removed 43,138 non-coding/pseudogenes) |
| Low-expression filter | ≥10 counts in ≥7 of 9 samples | 20,049 → 11,112 robustly expressed genes input to DESeq2 |
| DE threshold | padj ≤ 0.05, \|log2FC\| ≥ 0.263 (= 1.2×) | Conservative fold-change avoids over-reporting small but significant changes |
| LFC shrinkage | ashr | Appropriate for unbalanced designs and small sample sizes |
| Comparison design | All 3 pairwise contrasts | Client requested full pairwise comparison (两两对比); no single reference group designated |

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
| Mix_vs_NT | 5478 | 2670 | 2808 |
| A5BKO_vs_NT | 18 | 12 | 6 |
| A5BKO_vs_Mix | 6212 | 3160 | 3052 |

### 6.2 Pathway Enrichment (GO BP + KEGG, ALL direction)

| Contrast | GO BP terms | KEGG pathways |
| :--- | :---: | :---: |
| Mix_vs_NT | 31 | 6 |
| A5BKO_vs_NT | 0 | 0 |
| A5BKO_vs_Mix | 18 | 1 |

## 7. Conclusions

- **Mix_vs_NT**: 5478 DEGs identified. Predominantly downregulated (2808 down vs 2670 up).
- **A5BKO_vs_NT**: 18 DEGs identified. Predominantly upregulated (12 up vs 6 down).
- **A5BKO_vs_Mix**: 6212 DEGs identified. Predominantly upregulated (3160 up vs 3052 down).
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
