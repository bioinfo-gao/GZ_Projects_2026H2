# Bioinformatics Analysis Report

Project: 1_opossum_YuFan (Didelphis virginiana, NC vs pi5)  
Date: 2026-06-20  
by Zhen Gao, PhD  
Principal Bioinformatics Scientist, Athenomics

## 1. Overview
This report summarizes the differential expression analysis and quality control metrics for the RNA-seq dataset.
- **Analysis Tool**: DESeq2
- **Normalization**: VST (Variance Stabilizing Transformation) for PCA/Heatmap, Median-of-ratios for DE
- **LFC shrinkage**: ashr (Stephens 2016) — adaptive, shrinkage strength depends on each gene's own standard error
- **Significance Thresholds**: padj < 0.05, |log2FoldChange| >= 0.585 (approximately a 1.5-fold change in expression, up or down)

**Key Finding**: Principal component analysis (PCA) shows that the NC and pi5 sample groups do
not form separate, distinguishable clusters (PC1 explains 62% of variance, PC2 explains 15% of variance; samples from both groups overlap along both axes — see `PCA.pdf` and Section 6).
This indicates that, overall, the transcriptomes of the two groups are highly similar at this
sample size, and the differential expression results in Section 4 below should be read with that
context in mind.

## 2. Upstream Pipeline & Reference Caveats
This experiment uses a non-model organism with a liftoff-transferred annotation, which has
real implications for how the downstream DE results below should be interpreted.

**Genome / annotation**
- Species: *Didelphis virginiana* (Virginia opossum). Genome assembly: `mDidVir1`
  (DNA Zoo Consortium, Hi-C-scaffolded; FASTA file `dv-2k.fasta`; total length ~3.42 Gb across
  499,601 sequences, largest scaffold `HiC_scaffold_1` = 428.6 Mb; downloaded 2026-06-11). This
  assembly has no GenBank/RefSeq accession; for citation purposes it should be cited via the DNA
  Zoo resource directly (`https://www.dnazoo.org/assemblies/Didelphis_virginiana`).
- **Why an annotation transfer was needed**: the DNA Zoo `mDidVir1` assembly provides genome
  *sequence* only — there is no native, experimentally-curated gene annotation for this assembly.
  Generating one from scratch (ab initio gene prediction, or full RNA-seq-based annotation) was out
  of scope for this project, so a homology-based annotation transfer from a well-annotated relative
  species was used instead.
- **Why *Monodelphis domestica* was chosen as the source**: *M. domestica* (gray short-tailed
  opossum) is, like *D. virginiana*, a New World marsupial (family Didelphidae) with a mature,
  NCBI RefSeq-curated genome annotation. Within the marsupials with high-quality reference
  annotations, it is the closest well-annotated relative to *D. virginiana*, which gives the
  homology mapping step the best chance of finding conserved synteny and accurately transferring
  exon/intron structure.
- **Source genome/annotation version**: *M. domestica* assembly `MonDom5` (GenBank
  `GCA_000002295.1` / RefSeq `GCF_000002295.2`, released 2007-01-25), gene models from NCBI
  Annotation Release 103 (released 2016-04-14).
- **Method**: gene models were transferred from *M. domestica* (RefSeq annotation, assembly
  `MonDom5`) onto the *D. virginiana* `mDidVir1` assembly using **liftoff v1.6.3** (bioconda;
  Shumate & Salzberg, *Bioinformatics* 2021). Liftoff aligns each annotated gene region (plus a
  flanking margin) from the source genome directly onto the target genome by sequence homology,
  then maps the exact exon/intron structure onto the matching target sequence rather than
  predicting genes from scratch
  — this preserves the source annotation's structure but means liftoff can only transfer genes that
  exist in the source genome; it cannot detect genes that are unique to *D. virginiana*. Parameters
  used: minimum sequence identity `-sc 0.85` (mappings below 85% identity are discarded), search
  window extended `-flank 0.1` (10% of gene length) on each side to accommodate small structural
  rearrangements between the two genomes, restricted to `transcript`/`exon`/`CDS` feature types.
  The `gene` feature type was excluded because including it caused liftoff to fail for this
  particular genome pair (a tool/configuration incompatibility, not a data quality issue); this
  has no material effect on the results, since gene-level coordinates were reconstructed
  afterward by aggregating the transferred transcript records for each gene (see below), which is
  equivalent to the standard gene-level span.
- **Transfer outcome**: of 75,270 transcripts annotated in the *M. domestica* reference, 71,687
  (95.2%) were successfully mapped onto the *D. virginiana* genome; 3,583 (4.8%) could not be
  confidently mapped (below the 85% identity threshold or no syntenic region found) and were
  excluded.
- Two characteristics of the resulting annotation should be noted when interpreting gene-level
  results below:
  - This liftoff GTF has only `transcript`/`exon`/`CDS` feature rows (no `gene` rows, as noted
    above); gene-level coordinates were derived by aggregating transcript records.
  - **415 of 27,668 genes (~1.5%)** have their `gene_id` mapped to two different scaffolds
    simultaneously (a known liftoff artifact: when a gene has a paralog, repeat-region match, or
    ambiguous syntenic candidate elsewhere in the target genome, liftoff can place a copy at more
    than one locus). The locus with the most supporting transcripts was kept as primary; see the
    `n_loci` column in `Didelphis_virginiana_Gene_Annotation_Client.csv` to identify which genes
    are affected.

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
- **Overall sequencing and alignment quality is good and consistent across all 8 samples**: STAR alignment rate 89.57-91.02% (uniquely mapped 86.57-87.8%), sequencing error rate 0.69-0.79%, GC content 41-43% (no outlier samples). This indicates the library size difference noted below reflects differences in sequencing depth/read count between groups, not a difference in data quality.
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
- **Description**: Shows sample clustering based on the top 500 most variable genes.
- **Result for this dataset**: NC and pi5 samples overlap and do not form separate clusters along PC1 (62% variance) and PC2 (15% variance) — see the Key Finding in Section 1 for what this means for interpreting the DE results below.

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
| `Bioinformatics_Analysis_Report.md` | This report (located directly in `Data_Analysis/`). |
| `DE_PCA_Results/` | Output folder for this analysis; contains all files below up to `Reads/`. |
| &nbsp;&nbsp;&nbsp;&nbsp;`DEG_*.csv` | Differential expression results per contrast, including log2FC, p-values, and base means. |
| &nbsp;&nbsp;&nbsp;&nbsp;`PCA.pdf` | PCA plot showing sample relationships. |
| &nbsp;&nbsp;&nbsp;&nbsp;`Volcano_*.png` | Volcano plot for each contrast. |
| &nbsp;&nbsp;&nbsp;&nbsp;`Heatmap_padj_sig_genes_pi5_vs_NC.pdf` | Heatmap of padj-significant genes (see Section 6). |
| &nbsp;&nbsp;&nbsp;&nbsp;`Check_padj_sig_genes_per_sample_dotplot.png` | Per-sample verification plot for padj-significant genes (see Section 6). |
| &nbsp;&nbsp;&nbsp;&nbsp;`Sig_padj_genes_manual_check.csv` | Per-sample raw and normalized counts plus pre/post-shrinkage log2FC for the padj-significant genes, for manual review. |
| `Reads/` | Folder containing copies of the raw count and TPM matrices used as input for this analysis. |
| &nbsp;&nbsp;&nbsp;&nbsp;`All_sample_gene_counts.tsv` | Raw count matrix for all samples. |
| &nbsp;&nbsp;&nbsp;&nbsp;`All_sample_gene_tpm.tsv` | TPM (Transcripts Per Million) matrix. |
| `QC/` | MultiQC and other QC reports. |
| `Didelphis_virginiana_Gene_Annotation_Client.csv` | Gene-level annotation derived from the liftoff GTF (located directly in `Data_Analysis/`, not inside the folders above); see the `n_loci` column for genes with multiple candidate loci. |
