#!/usr/bin/env Rscript
# 运行环境: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript
# 运行方法: cd /home/gao/projects_2026H2/10_Yue_Liu/scripts && Rscript 6_generate_research_report.R
#
# 前置条件: 4_run_DE_PCA.R + 5_run_enrichment.R 均已完成
# (REDO #2: 独立同批次模型 vs 合并跨批次模型，详见脚本4顶部注释)

library(dplyr)
library(readr)

# ================= 1. 路径设置 =================
setwd("/home/gao/projects_2026H2/10_Yue_Liu/scripts/")

data_dirs <- sort(list.dirs("..", full.names = TRUE, recursive = FALSE), decreasing = TRUE)
data_dirs <- data_dirs[grepl("/Data_Analysis_[0-9]{8}$", data_dirs)]
if (length(data_dirs) == 0) stop("No Data_Analysis_YYYYMMDD folder found. Run 4_run_DE_PCA.R first.")
DATA_DIR <- data_dirs[1]
cat("Using:", DATA_DIR, "\n")

REL_DIR   <- file.path(DATA_DIR, "DE_PCA_Results_Reliable")
UNREL_DIR <- file.path(DATA_DIR, "DE_PCA_Results_Unreliable_CrossBatch")
ENR_DIR   <- file.path(DATA_DIR, "Enrichment")

REPORT_DATE <- format(Sys.Date(), "%m%d")
REPORT_FILE <- file.path(DATA_DIR, paste0("Bioinformatics_Analysis_Report_", REPORT_DATE, ".md"))

# ================= 2. 加载 DE 结果 =================
res_reliable   <- readRDS(file.path(REL_DIR, "res_list.rds"))
res_unreliable <- readRDS(file.path(UNREL_DIR, "res_list.rds"))
sig_col  <- "sig (padj<=0.05 & |log2FC|>=0.263)"
fs_file  <- file.path(REL_DIR, "filter_stats.rds")
fstats   <- if (file.exists(fs_file)) readRDS(fs_file) else NULL

deg_summary <- function(res_list) lapply(res_list, function(df) {
  up   <- sum(df[[sig_col]] == "Up",   na.rm = TRUE)
  down <- sum(df[[sig_col]] == "Down", na.rm = TRUE)
  list(up = up, down = down, total = up + down)
})
deg_rel   <- deg_summary(res_reliable)
deg_unrel <- deg_summary(res_unreliable)

# ================= 3. 读取 enrichment 汇总 (仅 reliable 对比) =================
enr_summary <- list()
for (comp_name in names(res_reliable)) {
  comp_dir  <- file.path(ENR_DIR, comp_name)
  go_file   <- file.path(comp_dir, "GO",   "GO_BP_ALL.csv")
  kegg_file <- file.path(comp_dir, "KEGG", "KEGG_ALL.csv")
  enr_summary[[comp_name]] <- list(
    go_terms   = if (file.exists(go_file))   nrow(read_csv(go_file,   show_col_types = FALSE)) else 0,
    kegg_paths = if (file.exists(kegg_file)) nrow(read_csv(kegg_file, show_col_types = FALSE)) else 0
  )
}

read_top <- function(path, desc_col = "Description", padj_col = "p.adjust") {
  tryCatch({
    df <- read_csv(path, show_col_types = FALSE)
    if (nrow(df) == 0) return(NA_character_)
    paste0(df[[desc_col]][1], " (padj=", signif(df[[padj_col]][1], 3), ")")
  }, error = function(e) NA_character_)
}

kf_lines <- c()
for (comp_name in names(deg_rel)) {
  std_dir <- file.path(ENR_DIR, comp_name)
  top_go        <- read_top(file.path(std_dir, "GO",   "GO_BP_ALL.csv"))
  top_kegg      <- read_top(file.path(std_dir, "KEGG", "KEGG_ALL.csv"), padj_col = "p.adjust")
  top_gsea_kegg <- read_top(file.path(std_dir, "GSEA", "GSEA_KEGG.csv"),     padj_col = "p.adjust")
  top_hallmark  <- read_top(file.path(std_dir, "GSEA", "GSEA_Hallmark.csv"), padj_col = "p.adjust")
  top_hallmark  <- if (!is.na(top_hallmark))
    sub("^HALLMARK_", "", gsub("_", " ", sub(",.*", "", top_hallmark))) else NA_character_
  kf_lines <- c(kf_lines,
    paste0("- **", comp_name, "**:"),
    if (!is.na(top_go))        paste0("  - Top GO (BP): ", top_go),
    if (!is.na(top_kegg))      paste0("  - Top KEGG: ",   top_kegg),
    if (!is.na(top_gsea_kegg)) paste0("  - GSEA KEGG: ",  top_gsea_kegg),
    if (!is.na(top_hallmark))  paste0("  - GSEA Hallmark: ", top_hallmark),
    ""
  )
}

report <- c(
  paste0("# Human Bulk RNA-seq Analysis Report — Yue Liu Project (MOLM13 X-ray Radiation Time Course)"),
  "",
  paste0("**Report Date:** ", format(Sys.Date(), "%B %d, %Y")),
  paste0("**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics"),
  paste0("**Analysis Platform:** Linux HPC server"),
  paste0("**Data folder:** `", basename(DATA_DIR), "`"),
  "",

  "> ## ⚠️ READ THIS FIRST: cross-batch sequencing effect",
  ">",
  "> This project integrates two sequencing batches — a **client-supplied legacy batch**",
  "> (Control, 4h post X-ray) and the **current NovaSeq X Plus batch** (8h, 16h post X-ray).",
  "> No condition was sequenced in both batches, so it is statistically impossible to",
  "> separate a true radiation-time effect from a sequencing-batch effect for any",
  "> comparison that spans the two batches (i.e. 8h vs Control, 16h vs Control).",
  ">",
  "> **All results in this report are split into two clearly labeled sets:**",
  "> - **Reliable** (`DE_PCA_Results_Reliable/`): Control vs 4h, and 8h vs 16h — each fit",
  ">   as its own independent, same-batch DESeq2 model. Full GO/KEGG/GSEA enrichment provided.",
  "> - **Unreliable — cross-batch reference only** (`DE_PCA_Results_Unreliable_CrossBatch/`):",
  ">   8h vs Control and 16h vs Control — provided only as a visual/reference resource",
  ">   showing what a naive full-time-course comparison would look like. **Do not use these",
  ">   for conclusions.** No enrichment analysis was run on them.",
  ">",
  "> This split was made after a concrete diagnostic finding (Section 4.1) showed that fitting",
  "> all 4 groups in one shared DESeq2 model contaminates variance estimation badly enough to",
  "> produce false positives even in the same-batch Control-vs-4h comparison.",
  "",

  "## 1. Objectives",
  "",
  "Characterise the transcriptomic response of MOLM13 cells to X-ray radiation, using two",
  "statistically independent same-batch comparisons: **Control vs 4h** (legacy batch) and",
  "**8h vs 16h** (current NovaSeq batch). Specific aims:",
  "",
  "- Identify differentially expressed genes (DEGs) for each reliable, same-batch contrast.",
  "- Perform GO and KEGG pathway enrichment and GSEA for the reliable contrasts.",
  "- Provide a clearly-flagged cross-batch reference view (unreliable) for visual context only.",
  "",

  "## 2. Key Findings",
  "",
  "**Reliable, same-batch DEGs** (criteria: padj ≤ 0.05 AND |log2FC| ≥ 0.263; log2(1.2) = 0.263, ≥1.2-fold):",
  "",
  "| Contrast | Batch | Total DEGs | Upregulated | Downregulated |",
  "| :--- | :---: | :---: | :---: | :---: |",
  paste0("| Control vs 4h | Legacy (n=6) | **", deg_rel[["4h_vs_Control"]]$total, "** | ",
         deg_rel[["4h_vs_Control"]]$up, " | ", deg_rel[["4h_vs_Control"]]$down, " |"),
  paste0("| 8h vs 16h | Current NovaSeq (n=6) | **", deg_rel[["16h_vs_8h"]]$total, "** | ",
         deg_rel[["16h_vs_8h"]]$up, " | ", deg_rel[["16h_vs_8h"]]$down, " |"),
  "",
  "**Top pathway findings (reliable contrasts only):**",
  "",
  kf_lines,

  "**For reference only — NOT to be used for conclusions** (cross-batch, see warning above):",
  "",
  "| Contrast | Total DEGs | Upregulated | Downregulated |",
  "| :--- | :---: | :---: | :---: |",
  paste0("| 8h vs Control (cross-batch) | ", deg_unrel[["8h_vs_Control"]]$total, " | ",
         deg_unrel[["8h_vs_Control"]]$up, " | ", deg_unrel[["8h_vs_Control"]]$down, " |"),
  paste0("| 16h vs Control (cross-batch) | ", deg_unrel[["16h_vs_Control"]]$total, " | ",
         deg_unrel[["16h_vs_Control"]]$up, " | ", deg_unrel[["16h_vs_Control"]]$down, " |"),
  "",

  "## 3. Sample Information",
  "",
  "| Group | Samples | Sequencing Batch | Role |",
  "| :--- | :---: | :---: | :---: |",
  paste0("| Control | Control_1, Control_2, Control_3 | Legacy batch (client-supplied count matrix) | Untreated baseline |"),
  paste0("| 4h | 4h_1, 4h_2, 4h_3 (orig. X-Ray-1/2/3) | Legacy batch (client-supplied count matrix) | 4h post X-ray |"),
  paste0("| 8h | 8h_1, 8h_2, 8h_3 (orig. \"A\" group, corrected) | Current NovaSeq X Plus batch | 8h post X-ray |"),
  paste0("| 16h | 16h_1, 16h_2, 16h_3 | Current NovaSeq X Plus batch | 16h post X-ray |"),
  "",
  "Total samples: **12**  |  Reliable comparisons: Control vs 4h, 8h vs 16h",
  "",
  "**Sample relabeling correction (2026-07-04):** the sheet originally labeled the 8h group as \"A\" —",
  "the client confirmed this was a documentation error; the samples are MOLM13 cells 8h post X-ray",
  "radiation. Raw FASTQ files on disk retain their original names (A_1/A_2/A_3); only the group",
  "label was corrected, no realignment was needed.",
  "",

  "## 4. Analysis Rationale and Decision Criteria",
  "",
  "### 4.1 ⚠️ Diagnostic finding that drove this design: shared-model contamination",
  "",
  "An initial attempt fit ALL 4 groups (Control/4h/8h/16h, 12 samples) in ONE DESeq2 model",
  "(`~Group`). Client review flagged gene **FSCN1 (ENSG00000075618)** as suspicious: its raw",
  "counts are essentially flat between Control (528, 618, 415) and 4h (593, 512, 434), yet the",
  "merged model reported it as a significant hit (padj = 0.0017, log2FC = 0.49, \"Up\").",
  "",
  "Diagnosis: refitting Control+4h **in isolation** (excluding 8h/16h entirely) gives",
  "padj = 0.9996, log2FC ≈ 0.0005 — i.e. genuinely no effect, consistent with the raw counts.",
  "The merged model's `baseMean` for this gene jumps from 521 (Control+4h only) to 3020 once",
  "8h/16h samples are included (FSCN1 is ~10× higher in 8h/16h), which drags the mean-dispersion",
  "trend curve DESeq2 fits across the whole dataset and corrupts variance shrinkage for genes",
  "throughout the dataset — not just the ones truly affected by 8h/16h. In short: **combining",
  "the two batches into one DESeq2 model produces false positives even for the same-batch",
  "Control vs 4h contrast**, which is more severe than an ordinary batch confound.",
  "",
  "**Fix applied:** Control vs 4h and 8h vs 16h are each now fit as their own independent",
  "2-group DESeq2 model, with no shared dispersion estimation between batches. The merged",
  "12-sample model is retained only to generate the cross-batch reference materials in",
  "`DE_PCA_Results_Unreliable_CrossBatch/` (PCA + 8h/16h vs Control DEG lists), which remain",
  "additionally confounded by batch on top of this contamination issue and are not used for",
  "any conclusion or enrichment analysis.",
  "",

  "### 4.2 Gene filtering — methodology consistency check vs. client's legacy file",
  "",
  {
    goc <- if (!is.null(fstats)) fstats$gene_overlap_check else NULL
    if (!is.null(goc)) c(
      sprintf("The client's legacy count matrix (`counts-filtered protein coding original file.xlsx`) was pre-filtered to **%s protein-coding genes**.", format(goc$old_total_genes, big.mark=",")),
      sprintf("Cross-checking by base Ensembl ID against our GENCODE v45 reference: **%s / %s (%.2f%%)** of the client's genes matched;",
              format(goc$common_genes, big.mark=","), format(goc$old_total_genes, big.mark=","), 100*goc$common_genes/goc$old_total_genes),
      sprintf("**%d genes** were absent from our v45 annotation entirely (IDs in the ENSG00000293xxx range — likely added in a newer Ensembl/GENCODE release than v45).", goc$missing_from_new),
      sprintf("Of the %s genes in common, our own annotation independently classified **%s (%.2f%%) as protein_coding** — confirming our biotype filter is highly consistent with the client's pre-filtering methodology.",
              format(goc$common_genes, big.mark=","), format(goc$old_common_protein_coding, big.mark=","), 100*goc$old_common_protein_coding/goc$common_genes)
    ) else "Gene overlap validation data not available."
  },
  "",

  "### 4.3 Decision criteria",
  "",
  "| Step | Decision | Rationale |",
  "| :--- | :---: | :---: |",
  "| Quantification | Salmon (new batch, via nf-core/rnaseq); client-supplied matrix (legacy batch) | Bias-corrected transcript-level quantification for the new batch |",
  "| Gene ID matching | Base Ensembl ID (version suffix stripped) | 99.96% of client's genes matched our annotation this way |",
  "| Biotype gene filter | Keep gene_type == protein_coding (GENCODE v45 annotation) | Matches client's own pre-filtering convention (see 4.2) |",
  "| Low-expression filter | **Not applied** | Matches client's legacy file, which retains all-zero genes (4,746/20,065 such genes present); DESeq2's built-in independent filtering still applies at the padj-correction step |",
  "| DE threshold | padj ≤ 0.05, \\|log2FC\\| ≥ 0.263 (= 1.2×) | Conservative fold-change avoids over-reporting small but significant changes |",
  "| LFC shrinkage | ashr | Appropriate for unbalanced designs and small sample sizes |",
  "| Model structure | Two independent same-batch DESeq2 models (NOT one shared 4-group model) | Prevents cross-batch variance contamination (see 4.1) |",
  "",

  "## 5. Methods",
  "",
  "| Tool | Version | Parameters |",
  "| :--- | :---: | :---: |",
  "| nf-core/rnaseq (new batch: 8h/16h) | 3.15.1 | --aligner star_salmon |",
  "| Legacy batch (Control/4h) | client-supplied | pre-quantified count matrix, method not independently verified |",
  "| DESeq2 | R package | design = ~Group (independent 2-group models), lfcShrink type = ashr |",
  "| clusterProfiler | R package | GO ORA (BP/MF/CC), KEGG ORA, GSEA (reliable contrasts only) |",
  "| org.Hs.eg.db | R package | Human gene ID mapping |",
  "| msigdbr | R package | Hallmark gene sets (Homo sapiens) |",
  "| Reference genome (new batch) | GRCh38 / GENCODE v45 | — |",
  "",

  "## 6. Results",
  "",
  "### 6.1 Reliable, same-batch Differential Expression",
  "",
  "| Contrast | Total DEGs | Upregulated | Downregulated |",
  "| :--- | :---: | :---: | :---: |",
  paste0("| Control vs 4h | ", deg_rel[["4h_vs_Control"]]$total, " | ", deg_rel[["4h_vs_Control"]]$up, " | ", deg_rel[["4h_vs_Control"]]$down, " |"),
  paste0("| 8h vs 16h | ", deg_rel[["16h_vs_8h"]]$total, " | ", deg_rel[["16h_vs_8h"]]$up, " | ", deg_rel[["16h_vs_8h"]]$down, " |"),
  "",
  "### 6.2 Pathway Enrichment (GO BP + KEGG, ALL direction; reliable contrasts only)",
  "",
  "| Contrast | GO BP terms | KEGG pathways |",
  "| :--- | :---: | :---: |",
  {
    sapply(names(enr_summary), function(cn) {
      e <- enr_summary[[cn]]
      paste0("| ", cn, " | ", e$go_terms, " | ", e$kegg_paths, " |")
    })
  },
  "",
  "### 6.3 Cross-batch reference (UNRELIABLE — not used for conclusions)",
  "",
  "| Contrast | Total DEGs | Upregulated | Downregulated |",
  "| :--- | :---: | :---: | :---: |",
  paste0("| 8h vs Control | ", deg_unrel[["8h_vs_Control"]]$total, " | ", deg_unrel[["8h_vs_Control"]]$up, " | ", deg_unrel[["8h_vs_Control"]]$down, " |"),
  paste0("| 16h vs Control | ", deg_unrel[["16h_vs_Control"]]$total, " | ", deg_unrel[["16h_vs_Control"]]$up, " | ", deg_unrel[["16h_vs_Control"]]$down, " |"),
  "",

  "## 7. Conclusions",
  "",
  paste0("- **Control vs 4h** (legacy batch, reliable): ", deg_rel[["4h_vs_Control"]]$total,
         " DEGs identified at the padj≤0.05 & |log2FC|≥0.263 threshold. With FSCN1's false-positive",
         " signal removed, this same-batch comparison shows essentially no detectable bulk",
         " transcriptomic change at 4h post X-ray radiation under this threshold (GSEA still finds",
         " coordinated but sub-threshold shifts in cell-cycle/E2F-target genes — see Section 2)."),
  paste0("- **8h vs 16h** (current batch, reliable): ", deg_rel[["16h_vs_8h"]]$total,
         " DEGs identified (", deg_rel[["16h_vs_8h"]]$up, " up / ", deg_rel[["16h_vs_8h"]]$down,
         " down in 16h relative to 8h), indicating substantial continued transcriptomic change between",
         " these two later time points."),
  "- **No statistically valid comparison is available linking the legacy batch (Control/4h) to the",
  "  current batch (8h/16h)** — any apparent 8h/16h vs Control signal reflects sequencing batch",
  "  differences at least as much as biology, and should not be interpreted as a radiation-time trend.",
  "- To build a genuine, statistically sound Control→4h→8h→16h time course, we recommend re-sequencing",
  "  a Control condition alongside the 8h/16h samples in the same batch (or re-sequencing all four",
  "  time points together) so that a proper `~Batch + Group` design becomes estimable.",
  "- Pathway enrichment and GSEA results for the two reliable contrasts are available in `Enrichment/`.",
  "",

  "## 8. Deliverable Files",
  "",
  "| File / Folder | Contents |",
  "| :--- | :---: |",
  paste0("| `DE_PCA_Results_Reliable/PCA_*.pdf` | PCA plots, one per same-batch model |"),
  paste0("| `DE_PCA_Results_Reliable/DEG_*.csv` | Full DEG tables for the 2 reliable contrasts |"),
  paste0("| `DE_PCA_Results_Reliable/Volcano_*.png`, `Heatmap_top50_*.pdf` | Volcano / heatmap, reliable contrasts |"),
  paste0("| `DE_PCA_Results_Unreliable_CrossBatch/` | Cross-batch reference PCA + DEG/volcano/heatmap for 8h/16h vs Control — reference only, not for conclusions |"),
  paste0("| `DE_PCA_Results_Unreliable_CrossBatch/All_12samples_gene_counts.tsv` | Merged 12-sample count matrix (all 4 timepoints) |"),
  paste0("| `Reads/New_batch_8h_16h_gene_counts.tsv` | Raw count matrix, current NovaSeq batch (8h/16h) |"),
  paste0("| `Reads/Old_batch_Control_4h_gene_counts.xlsx` | Client-supplied legacy count matrix (Control/4h), copied as-is |"),
  paste0("| `human_Gene_annotation_*.xlsx` | Full human gene annotation with GO/KEGG/UniProt (GENCODE v45) |"),
  paste0("| `Enrichment/{4h_vs_Control,16h_vs_8h}/{GO,KEGG,GSEA}/` | Enrichment results, reliable contrasts only |"),
  paste0("| `QC/multiqc/` | MultiQC report — covers current NovaSeq batch (8h/16h) only; no QC available for the legacy batch |"),
  "",
  "---",
  "*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*"
)

writeLines(report, con = REPORT_FILE)
cat("Report saved:", REPORT_FILE, "\n")
