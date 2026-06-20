# Bioinformatics Analysis Report

Date: 2026-06-19
Project: 1_opossum_YuFan (Didelphis virginiana, NC vs pi5)
Author: Zhen Gao, PhD, Principal Bioinformatics Scientist

## 1. Overview
This report summarizes the differential expression analysis and quality control metrics for the RNA-seq dataset.
- **Analysis Tool**: DESeq2
- **Normalization**: VST (Variance Stabilizing Transformation) for PCA/Heatmap, Median-of-ratios for DE
- **LFC shrinkage**: ashr (Stephens 2016) — adaptive, shrinkage strength depends on each gene's own standard error
- **Significance Thresholds**: padj < 0.05, |log2FoldChange| >= 0.585 (approximately a 1.5-fold change in expression, up or down)

## 2. Upstream Pipeline & Reference Caveats
This experiment uses a non-model organism with a liftoff-transferred annotation, which has
real implications for how the downstream DE results below should be interpreted.

**Genome / annotation**
- Species: *Didelphis virginiana* (opossum). Genome assembly: `dv-2k.fasta` (Hi-C scaffolded, DNA Zoo).
- Gene annotation: `Didelphis_v.liftoff.gtf` — produced by **liftoff** (homology-based annotation
  transfer from a related reference species), **not** a native, experimentally-curated annotation
  for this species. Two characteristics of this annotation should be noted:
  - This liftoff GTF has only `transcript`/`exon`/`CDS` feature rows (no `gene` rows); gene-level
    coordinates were derived by aggregating transcript records.
  - **415 of 27,668 genes (~1.5%)** have their `gene_id` mapped to two different scaffolds
    simultaneously (a known liftoff artifact, e.g. paralog/repeat region mis-mapping). The locus
    with the most supporting transcripts was kept as primary; see the `n_loci` column in
    `Didelphis_virginiana_Gene_Annotation_Client.csv` to identify which genes are affected.

**Alignment (nf-core/rnaseq, aligner = star_salmon)**
- STAR was run with 2-pass mode (`--twopassMode Basic`). Splice junctions discovered across all
  samples in the first pass were used when finalizing alignments in the second pass, which
  partially compensates for the liftoff annotation being an imperfect/incomplete transfer by
  recovering species-specific or novel junctions that the liftoff GTF alone would have missed.
- STAR filtering parameters were tightened for this non-model-organism / imperfect-reference
  scenario: `--outFilterMultimapNmax 8`, `--alignSJoverhangMin 8`,
  `--alignSJDBoverhangMin 1`, `--outFilterMismatchNmax 2`.
- Quantification: Salmon, using a tx2gene mapping built from the same liftoff GTF.

## 3. Quality Control (QC)
- QC reports were generated using MultiQC.
- Mean library size differs by **23.8%** between NC and pi5 groups (NC mean = 43,729,998 reads; pi5 mean = 34,437,340 reads). DESeq2 size-factor normalization corrects for this at the model-fitting level, but a difference of this size is worth confirming against the sequencing/library-prep batch records (see Section 6).

## 4. Differential Expression Analysis Results

### Contrast: pi5_vs_NC
- Genes with padj < 0.05 (regardless of fold-change size): 16
- Of those, passing the |log2FC| >= 0.585 effect-size filter (i.e. at least a 1.5-fold change in expression; final "sig" call): 0 (Up: 0, Down: 0)
- Output File: `DEG_pi5_vs_NC.csv`

## 5. Key Caveats Found During This Analysis

- **PCA shows NC and pi5 samples do not separate** (see `PCA.pdf`) — samples from the two groups
  are interspersed rather than forming distinct clusters, indicating the overall transcriptome is
  highly similar between groups at this sample size.
- **ashr shrinkage compresses nearly all effect sizes toward zero** (most genes' shrunk log2FC fall
  within roughly ±0.1, i.e. well under a 1.1-fold change), consistent with the weak overall signal
  seen in the PCA. This is why the |log2FC| >= 0.585 (~1.5-fold) filtered "sig" column above is much smaller than the raw padj-significant count.
- For each contrast, the padj-significant genes (regardless of fold-change size) were checked for
  whether they are dominated by a single direction and/or show complete (non-overlapping) separation
  between groups despite modest effect sizes — a pattern that, when it affects *all* significant
  genes simultaneously, is more consistent with a shared systematic/technical confound (e.g. the
  library-size difference noted in Section 3, or an unmodeled batch effect correlated with Group)
  than with that many independent gene-specific biological regulation events:
  - pi5_vs_NC: 16 padj-significant genes — 16 up / 0 down (**100% same direction**); 16/16 show complete non-overlapping separation between groups.

- **Recommendation**: before treating any gene from this dataset's DEG list as a confirmed
  biological finding, (1) confirm with the wet-lab team whether NC and pi5 samples were
  prepared/sequenced in the same batch, (2) for genes of interest, inspect per-sample normalized
  counts individually (see `Check_padj_sig_genes_per_sample_dotplot.png` and
  `Sig_padj_genes_manual_check.csv`) rather than relying on padj/log2FC alone.

## 6. Visualizations

### Principal Component Analysis (PCA)
- **File**: `PCA.pdf`
- **Description**: Shows sample clustering based on the top 500 most variable genes. Samples should cluster by biological group if the treatment effect is strong.

### Volcano Plots
- **Files**: `Volcano_*.png`
- **Description**: Displays the relationship between statistical significance (-log10 padj) and magnitude of change (log2FC). Red points indicate upregulated genes, blue points indicate downregulated genes.
  Note: with strong ashr shrinkage and a narrow effect-size range, the two "wings" can look like a
  smooth curve rather than a scattered cloud — this is expected when shrunk LFC and -log10(padj)
  both become near-monotonic functions of the same underlying z-statistic.

### Heatmaps
- No contrast had any gene passing both the padj and |log2FC| thresholds, so no `Heatmap_top50_*.pdf` files were produced.
- **File**: `Heatmap_padj_sig_genes_pi5_vs_NC.pdf`
- **Description**: Heatmap of all genes passing the padj threshold alone (regardless of fold-change size), generated separately for manual review.
- **File**: `Check_padj_sig_genes_per_sample_dotplot.png`
- **Description**: Per-sample normalized counts for the same genes, plotted individually, to verify that significance calls are not driven by a single outlier sample.

## 7. Generated Data Files

| File Name | Description |
| :--- | :--- |
| `Bioinformatics_Analysis_Report.md` | This report. |
| `All_sample_gene_counts.tsv` (in `Reads/`) | Raw count matrix for all samples. |
| `All_sample_gene_tpm.tsv` (in `Reads/`) | TPM (Transcripts Per Million) matrix. |
| `DEG_*.csv` | Differential expression results per contrast, including log2FC, p-values, and base means. |
| `PCA.pdf` | PCA plot showing sample relationships. |
| `Volcano_*.png` | Volcano plot for each contrast. |
| `Heatmap_padj_sig_genes_pi5_vs_NC.pdf` | Heatmap of padj-significant genes (see Section 6). |
| `Check_padj_sig_genes_per_sample_dotplot.png` | Per-sample verification plot for padj-significant genes (see Section 6). |
| `Sig_padj_genes_manual_check.csv` | Per-sample raw and normalized counts plus pre/post-shrinkage log2FC for the padj-significant genes, for manual review. |
| `Didelphis_virginiana_Gene_Annotation_Client.csv` (in `Data_Analysis/`) | Gene-level annotation derived from the liftoff GTF; see the `n_loci` column for genes with multiple candidate loci. |
| `QC/` (in `Data_Analysis/`) | MultiQC and other QC reports. |
