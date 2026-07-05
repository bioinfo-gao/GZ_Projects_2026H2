# Human Bulk RNA-seq Analysis Report — Yue Liu Project (MOLM13 X-ray Radiation Time Course)

**Report Date:** July 04, 2026
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Data folder:** `Data_Analysis_20260704`

> ## ⚠️ READ THIS FIRST: cross-batch sequencing effect
>
> This project integrates two sequencing batches — a **client-supplied legacy batch**
> (Control, 4h post X-ray) and the **current NovaSeq X Plus batch** (8h, 16h post X-ray).
> No condition was sequenced in both batches, so it is statistically impossible to
> separate a true radiation-time effect from a sequencing-batch effect for any
> comparison that spans the two batches (i.e. 8h vs Control, 16h vs Control).
>
> **All results in this report are split into two clearly labeled sets:**
> - **Reliable** (`DE_PCA_Results_Reliable/`): Control vs 4h, and 8h vs 16h — each fit
>   as its own independent, same-batch DESeq2 model. Full GO/KEGG/GSEA enrichment provided.
> - **Unreliable — cross-batch reference only** (`DE_PCA_Results_Unreliable_CrossBatch/`):
>   8h vs Control and 16h vs Control — provided only as a visual/reference resource
>   showing what a naive full-time-course comparison would look like. **Do not use these
>   for conclusions.** No enrichment analysis was run on them.
>
> This split was made after a concrete diagnostic finding (Section 4.1) showed that fitting
> all 4 groups in one shared DESeq2 model contaminates variance estimation badly enough to
> produce false positives even in the same-batch Control-vs-4h comparison.

## 1. Objectives

Characterise the transcriptomic response of MOLM13 cells to X-ray radiation, using two
statistically independent same-batch comparisons: **Control vs 4h** (legacy batch) and
**8h vs 16h** (current NovaSeq batch). Specific aims:

- Identify differentially expressed genes (DEGs) for each reliable, same-batch contrast.
- Perform GO and KEGG pathway enrichment and GSEA for the reliable contrasts.
- Provide a clearly-flagged cross-batch reference view (unreliable) for visual context only.

## 2. Key Findings

**Reliable, same-batch DEGs** (criteria: padj ≤ 0.05 AND |log2FC| ≥ 0.263; log2(1.2) = 0.263, ≥1.2-fold):

| Contrast | Batch | Total DEGs | Upregulated | Downregulated |
| :--- | :---: | :---: | :---: | :---: |
| Control vs 4h | Legacy (n=6) | **0** | 0 | 0 |
| 8h vs 16h | Current NovaSeq (n=6) | **3969** | 1913 | 2056 |

**Top pathway findings (reliable contrasts only):**

- **4h_vs_Control**:
  - GSEA KEGG: Cell cycle (padj=0.0361)
  - GSEA Hallmark: HALLMARK E2F TARGETS (padj=5.1e-09)

- **16h_vs_8h**:
  - Top GO (BP): ribosome biogenesis (padj=2.67e-28)
  - Top KEGG: Ribosome biogenesis in eukaryotes (padj=4.31e-10)
  - GSEA KEGG: Ribosome (padj=3.3e-08)
  - GSEA Hallmark: HALLMARK MITOTIC SPINDLE (padj=1.67e-09)

**For reference only — NOT to be used for conclusions** (cross-batch, see warning above):

| Contrast | Total DEGs | Upregulated | Downregulated |
| :--- | :---: | :---: | :---: |
| 8h vs Control (cross-batch) | 11980 | 6241 | 5739 |
| 16h vs Control (cross-batch) | 12054 | 6235 | 5819 |

## 3. Sample Information

| Group | Samples | Sequencing Batch | Role |
| :--- | :---: | :---: | :---: |
| Control | Control_1, Control_2, Control_3 | Legacy batch (client-supplied count matrix) | Untreated baseline |
| 4h | 4h_1, 4h_2, 4h_3 (orig. X-Ray-1/2/3) | Legacy batch (client-supplied count matrix) | 4h post X-ray |
| 8h | 8h_1, 8h_2, 8h_3 (orig. "A" group, corrected) | Current NovaSeq X Plus batch | 8h post X-ray |
| 16h | 16h_1, 16h_2, 16h_3 | Current NovaSeq X Plus batch | 16h post X-ray |

Total samples: **12**  |  Reliable comparisons: Control vs 4h, 8h vs 16h

**Sample relabeling correction (2026-07-04):** the sheet originally labeled the 8h group as "A" —
the client confirmed this was a documentation error; the samples are MOLM13 cells 8h post X-ray
radiation. Raw FASTQ files on disk retain their original names (A_1/A_2/A_3); only the group
label was corrected, no realignment was needed.

## 4. Analysis Rationale and Decision Criteria

### 4.1 ⚠️ Diagnostic finding that drove this design: shared-model contamination

An initial attempt fit ALL 4 groups (Control/4h/8h/16h, 12 samples) in ONE DESeq2 model
(`~Group`). Client review flagged gene **FSCN1 (ENSG00000075618)** as suspicious: its raw
counts are essentially flat between Control (528, 618, 415) and 4h (593, 512, 434), yet the
merged model reported it as a significant hit (padj = 0.0017, log2FC = 0.49, "Up").

Diagnosis: refitting Control+4h **in isolation** (excluding 8h/16h entirely) gives
padj = 0.9996, log2FC ≈ 0.0005 — i.e. genuinely no effect, consistent with the raw counts.
The merged model's `baseMean` for this gene jumps from 521 (Control+4h only) to 3020 once
8h/16h samples are included (FSCN1 is ~10× higher in 8h/16h), which drags the mean-dispersion
trend curve DESeq2 fits across the whole dataset and corrupts variance shrinkage for genes
throughout the dataset — not just the ones truly affected by 8h/16h. In short: **combining
the two batches into one DESeq2 model produces false positives even for the same-batch
Control vs 4h contrast**, which is more severe than an ordinary batch confound.

**Fix applied:** Control vs 4h and 8h vs 16h are each now fit as their own independent
2-group DESeq2 model, with no shared dispersion estimation between batches. The merged
12-sample model is retained only to generate the cross-batch reference materials in
`DE_PCA_Results_Unreliable_CrossBatch/` (PCA + 8h/16h vs Control DEG lists), which remain
additionally confounded by batch on top of this contamination issue and are not used for
any conclusion or enrichment analysis.

### 4.2 Gene filtering — methodology consistency check vs. client's legacy file

The client's legacy count matrix (`counts-filtered protein coding original file.xlsx`) was pre-filtered to **20,065 protein-coding genes**.
Cross-checking by base Ensembl ID against our GENCODE v45 reference: **20,056 / 20,065 (99.96%)** of the client's genes matched;
**9 genes** were absent from our v45 annotation entirely (IDs in the ENSG00000293xxx range — likely added in a newer Ensembl/GENCODE release than v45).
Of the 20,056 genes in common, our own annotation independently classified **20,047 (99.96%) as protein_coding** — confirming our biotype filter is highly consistent with the client's pre-filtering methodology.

### 4.3 Decision criteria

| Step | Decision | Rationale |
| :--- | :---: | :---: |
| Quantification | Salmon (new batch, via nf-core/rnaseq); client-supplied matrix (legacy batch) | Bias-corrected transcript-level quantification for the new batch |
| Gene ID matching | Base Ensembl ID (version suffix stripped) | 99.96% of client's genes matched our annotation this way |
| Biotype gene filter | Keep gene_type == protein_coding (GENCODE v45 annotation) | Matches client's own pre-filtering convention (see 4.2) |
| Low-expression filter | **Not applied** | Matches client's legacy file, which retains all-zero genes (4,746/20,065 such genes present); DESeq2's built-in independent filtering still applies at the padj-correction step |
| DE threshold | padj ≤ 0.05, \|log2FC\| ≥ 0.263 (= 1.2×) | Conservative fold-change avoids over-reporting small but significant changes |
| LFC shrinkage | ashr | Appropriate for unbalanced designs and small sample sizes |
| Model structure | Two independent same-batch DESeq2 models (NOT one shared 4-group model) | Prevents cross-batch variance contamination (see 4.1) |

## 5. Methods

| Tool | Version | Parameters |
| :--- | :---: | :---: |
| nf-core/rnaseq (new batch: 8h/16h) | 3.15.1 | --aligner star_salmon |
| Legacy batch (Control/4h) | client-supplied | pre-quantified count matrix, method not independently verified |
| DESeq2 | R package | design = ~Group (independent 2-group models), lfcShrink type = ashr |
| clusterProfiler | R package | GO ORA (BP/MF/CC), KEGG ORA, GSEA (reliable contrasts only) |
| org.Hs.eg.db | R package | Human gene ID mapping |
| msigdbr | R package | Hallmark gene sets (Homo sapiens) |
| Reference genome (new batch) | GRCh38 / GENCODE v45 | — |

## 6. Results

### 6.1 Reliable, same-batch Differential Expression

| Contrast | Total DEGs | Upregulated | Downregulated |
| :--- | :---: | :---: | :---: |
| Control vs 4h | 0 | 0 | 0 |
| 8h vs 16h | 3969 | 1913 | 2056 |

### 6.2 Pathway Enrichment (GO BP + KEGG, ALL direction; reliable contrasts only)

| Contrast | GO BP terms | KEGG pathways |
| :--- | :---: | :---: |
| 4h_vs_Control | 0 | 0 |
| 16h_vs_8h | 38 | 2 |

### 6.3 Cross-batch reference (UNRELIABLE — not used for conclusions)

| Contrast | Total DEGs | Upregulated | Downregulated |
| :--- | :---: | :---: | :---: |
| 8h vs Control | 11980 | 6241 | 5739 |
| 16h vs Control | 12054 | 6235 | 5819 |

## 7. Conclusions

- **Control vs 4h** (legacy batch, reliable): 0 DEGs identified at the padj≤0.05 & |log2FC|≥0.263 threshold. With FSCN1's false-positive signal removed, this same-batch comparison shows essentially no detectable bulk transcriptomic change at 4h post X-ray radiation under this threshold (GSEA still finds coordinated but sub-threshold shifts in cell-cycle/E2F-target genes — see Section 2).
- **8h vs 16h** (current batch, reliable): 3969 DEGs identified (1913 up / 2056 down in 16h relative to 8h), indicating substantial continued transcriptomic change between these two later time points.
- **No statistically valid comparison is available linking the legacy batch (Control/4h) to the
  current batch (8h/16h)** — any apparent 8h/16h vs Control signal reflects sequencing batch
  differences at least as much as biology, and should not be interpreted as a radiation-time trend.
- To build a genuine, statistically sound Control→4h→8h→16h time course, we recommend re-sequencing
  a Control condition alongside the 8h/16h samples in the same batch (or re-sequencing all four
  time points together) so that a proper `~Batch + Group` design becomes estimable.
- Pathway enrichment and GSEA results for the two reliable contrasts are available in `Enrichment/`.

## 8. Deliverable Files

| File / Folder | Contents |
| :--- | :---: |
| `DE_PCA_Results_Reliable/PCA_*.pdf` | PCA plots, one per same-batch model |
| `DE_PCA_Results_Reliable/DEG_*.csv` | Full DEG tables for the 2 reliable contrasts |
| `DE_PCA_Results_Reliable/Volcano_*.png`, `Heatmap_top50_*.pdf` | Volcano / heatmap, reliable contrasts |
| `DE_PCA_Results_Unreliable_CrossBatch/` | Cross-batch reference PCA + DEG/volcano/heatmap for 8h/16h vs Control — reference only, not for conclusions |
| `DE_PCA_Results_Unreliable_CrossBatch/All_12samples_gene_counts.tsv` | Merged 12-sample count matrix (all 4 timepoints) |
| `Reads/New_batch_8h_16h_gene_counts.tsv` | Raw count matrix, current NovaSeq batch (8h/16h) |
| `Reads/Old_batch_Control_4h_gene_counts.xlsx` | Client-supplied legacy count matrix (Control/4h), copied as-is |
| `human_Gene_annotation_*.xlsx` | Full human gene annotation with GO/KEGG/UniProt (GENCODE v45) |
| `Enrichment/{4h_vs_Control,16h_vs_8h}/{GO,KEGG,GSEA}/` | Enrichment results, reliable contrasts only |
| `QC/multiqc/` | MultiQC report — covers current NovaSeq batch (8h/16h) only; no QC available for the legacy batch |

---
*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*
