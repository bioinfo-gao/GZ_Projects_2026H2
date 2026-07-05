#!/usr/bin/env Rscript
# 运行环境: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript
# 运行方法: cd /home/gao/projects_2026H2/10_Yue_Liu/scripts && Rscript 6_generate_research_report.R
#
# 前置条件: 4_run_DE_PCA.R + 5_run_enrichment.R 均已完成 (整合版: Control/4h/8h/16h)

library(dplyr)
library(readr)

# ================= 1. 路径设置 =================
setwd("/home/gao/projects_2026H2/10_Yue_Liu/scripts/")

data_dirs <- sort(list.dirs("..", full.names = TRUE, recursive = FALSE), decreasing = TRUE)
data_dirs <- data_dirs[grepl("/Data_Analysis_[0-9]{8}$", data_dirs)]
if (length(data_dirs) == 0) stop("No Data_Analysis_YYYYMMDD folder found. Run 4_run_DE_PCA.R first.")
DATA_DIR <- data_dirs[1]
cat("Using:", DATA_DIR, "\n")

DE_DIR  <- file.path(DATA_DIR, "DE_PCA_Results")
ENR_DIR <- file.path(DATA_DIR, "Enrichment")

REPORT_DATE <- format(Sys.Date(), "%m%d")
REPORT_FILE <- file.path(DATA_DIR, paste0("Bioinformatics_Analysis_Report_", REPORT_DATE, ".md"))

# ================= 2. 加载 DE 结果 =================
res_list <- readRDS(file.path(DE_DIR, "res_list.rds"))
meta     <- readRDS(file.path(DE_DIR, "meta.rds"))
sig_col  <- "sig (padj<=0.05 & |log2FC|>=0.263)"
fs_file  <- file.path(DE_DIR, "filter_stats.rds")
fstats   <- if (file.exists(fs_file)) readRDS(fs_file) else NULL

# 保持时间顺序: 4h, 8h, 16h (而非字母序)
time_order <- c("4h_vs_Control", "8h_vs_Control", "16h_vs_Control")
res_list <- res_list[intersect(time_order, names(res_list))]

deg_summary <- lapply(res_list, function(df) {
  up   <- sum(df[[sig_col]] == "Up",   na.rm = TRUE)
  down <- sum(df[[sig_col]] == "Down", na.rm = TRUE)
  list(up = up, down = down, total = up + down)
})

# ================= 3. 读取 enrichment 汇总 =================
enr_summary <- list()
for (comp_name in names(res_list)) {
  comp_dir  <- file.path(ENR_DIR, comp_name)
  go_file   <- file.path(comp_dir, "GO",   "GO_BP_ALL.csv")
  kegg_file <- file.path(comp_dir, "KEGG", "KEGG_ALL.csv")
  enr_summary[[comp_name]] <- list(
    go_terms   = if (file.exists(go_file))   nrow(read_csv(go_file,   show_col_types = FALSE)) else 0,
    kegg_paths = if (file.exists(kegg_file)) nrow(read_csv(kegg_file, show_col_types = FALSE)) else 0
  )
}

# ================= 4. 生成报告 =================
cat("Writing report:", REPORT_FILE, "\n")

read_top <- function(path, desc_col = "Description", padj_col = "p.adjust") {
  tryCatch({
    df <- read_csv(path, show_col_types = FALSE)
    if (nrow(df) == 0) return(NA_character_)
    paste0(df[[desc_col]][1], " (padj=", signif(df[[padj_col]][1], 3), ")")
  }, error = function(e) NA_character_)
}

kf_lines <- c()
for (comp_name in names(deg_summary)) {
  s   <- deg_summary[[comp_name]]
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

  "## 1. Objectives",
  "",
  "Characterise the transcriptomic response of MOLM13 cells to X-ray radiation over a",
  "**Control → 4h → 8h → 16h** time course, integrating a client-supplied legacy dataset",
  "(Control, 4h) with the current sequencing batch (8h, 16h). Specific aims:",
  "",
  "- Identify differentially expressed genes (DEGs) at each post-radiation time point vs. untreated Control.",
  "- Perform GO and KEGG pathway enrichment to determine biological processes affected at each time point.",
  "- Run Gene Set Enrichment Analysis (GSEA) for pathway-level signal detection.",
  "- Transparently flag the cross-batch statistical limitation described in Section 4.",
  "",

  "## 2. Key Findings",
  "",
  "This study integrates **2 sequencing batches** into **4 groups** (Control, 4h, 8h, 16h) across **3 contrasts**, all vs. Control.",
  "",
  "**Differentially expressed genes** (criteria: padj ≤ 0.05 AND |log2FC| ≥ 0.263; log2(1.2) = 0.263, equivalent to ≥1.2-fold change):",
  "",
  "| Contrast | Total DEGs | Upregulated (log2FC ≥ 0.263) | Downregulated (log2FC ≤ −0.263) |",
  "| :--- | :---: | :---: | :---: |",
  {
    sapply(names(deg_summary), function(cn) {
      s <- deg_summary[[cn]]
      paste0("| ", cn, " | **", s$total, "** | ", s$up, " | ", s$down, " |")
    })
  },
  "",
  "**Note the sharp jump in DEG count from 4h (same-batch comparison) to 8h/16h (cross-batch comparison)** —",
  "this is a direct symptom of the batch confound described in Section 4, not necessarily biology alone.",
  "",
  "**Top pathway findings per comparison:**",
  "",
  kf_lines,

  "## 3. Sample Information",
  "",
  "| Group | Samples | Sequencing Batch | Role |",
  "| :--- | :---: | :---: | :---: |",
  paste0("| Control | Control_1, Control_2, Control_3 | Legacy batch (client-supplied count matrix) | Untreated baseline |"),
  paste0("| 4h | 4h_1, 4h_2, 4h_3 (orig. X-Ray-1/2/3) | Legacy batch (client-supplied count matrix) | 4h post X-ray |"),
  paste0("| 8h | 8h_1, 8h_2, 8h_3 (orig. \"A\" group, corrected) | Current NovaSeq X Plus batch | 8h post X-ray |"),
  paste0("| 16h | 16h_1, 16h_2, 16h_3 | Current NovaSeq X Plus batch | 16h post X-ray |"),
  "",
  paste0("Total samples: **", nrow(meta), "**  |  ",
         "Comparisons: 4h vs Control, 8h vs Control, 16h vs Control"),
  "",
  "**Sample relabeling correction (2026-07-04):** the sheet originally labeled the 8h group as \"A\" —",
  "the client confirmed this was a documentation error; the samples are MOLM13 cells 8h post X-ray",
  "radiation. Raw FASTQ files on disk retain their original names (A_1/A_2/A_3); only the group",
  "label was corrected, no realignment was needed.",
  "",

  "## 4. Analysis Rationale and Decision Criteria",
  "",
  "### 4.1 ⚠️ Cross-batch statistical limitation (read before interpreting 8h/16h results)",
  "",
  "Control and 4h come from a **client-supplied legacy sequencing batch** (pre-quantified count",
  "matrix, quantification method not independently verified by this analysis). 8h and 16h come",
  "from the **current NovaSeq X Plus batch**, quantified with Salmon via nf-core/rnaseq in this project.",
  "",
  "**No experimental condition was sequenced in both batches** — there is no way to separate a true",
  "\"batch effect\" (different sequencing run, possibly different library prep/quantification pipeline)",
  "from the true biological radiation-response effect for the 8h and 16h vs Control contrasts. A",
  "`~Batch + Group` DESeq2 design is **not statistically estimable** here because Batch and Group are",
  "fully aliased (Control/4h only ever appear in the old batch; 8h/16h only ever appear in the new batch).",
  "",
  "Per client's explicit decision, Control is used as the common reference for all three contrasts to",
  "produce a unified time-course view. The **4h vs Control** contrast is a clean, same-batch comparison",
  "and its result (see Section 6.1) can be interpreted with normal confidence. The **8h vs Control** and",
  "**16h vs Control** contrasts should be read as *exploratory / batch-confounded* — the very large jump",
  "in DEG count relative to 4h (",
  {
    if (length(deg_summary) >= 3)
      paste0(deg_summary[[1]]$total, " at 4h vs ", deg_summary[[2]]$total, " at 8h vs ", deg_summary[[3]]$total, " at 16h")
    else ""
  },
  ") is consistent with batch-driven variance rather than a biologically plausible ~48× increase in",
  "radiation response between 4h and 8h. We recommend treating 8h/16h DEG lists as candidate genes",
  "requiring independent validation (e.g. qPCR), not as confirmed radiation-response findings.",
  "",

  "### 4.2 Gene filtering — methodology consistency check vs. client's legacy file",
  "",
  {
    goc <- if (!is.null(fstats)) fstats$gene_overlap_check else NULL
    if (!is.null(goc)) c(
      sprintf("The client's legacy count matrix (`counts-filtered protein coding original file.xlsx`) was pre-filtered to **%s protein-coding genes**.", format(goc$old_total_genes, big.mark=",")),
      sprintf("Cross-checking by base Ensembl ID against our GENCODE v45 reference: **%s / %s (%.2f%%)** of the client's genes matched;",
              format(goc$common_genes, big.mark=","), format(goc$old_total_genes, big.mark=","), 100*goc$common_genes/goc$old_total_genes),
      sprintf("**%d genes** were absent from our v45 annotation entirely (IDs in the ENSG00000293xxx range — likely added in a newer Ensembl/GENCODE release than v45, i.e. the client's original pipeline used a slightly newer annotation).", goc$missing_from_new),
      sprintf("Of the %s genes in common, our own annotation independently classified **%s (%.2f%%) as protein_coding** — confirming our biotype filter is highly consistent with the client's pre-filtering methodology (the remaining genes are classified as lncRNA/pseudogene under GENCODE v45, likely due to biotype reclassification between annotation versions — a known, expected phenomenon, not an error).",
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
  {
    if (!is.null(fstats))
      sprintf("| Low-expression filter | **Not applied** | Matches client's legacy file, which retains all-zero genes (%s such genes present); DESeq2's built-in independent filtering still applies at the padj-correction step |", "4,746/20,065")
    else
      "| Low-expression filter | Not applied | Matches client's legacy file convention |"
  },
  "| DE threshold | padj ≤ 0.05, \\|log2FC\\| ≥ 0.263 (= 1.2×) | Conservative fold-change avoids over-reporting small but significant changes |",
  "| LFC shrinkage | ashr | Appropriate for unbalanced designs and small sample sizes |",
  "| Cross-batch design | Common Control reference for all 3 contrasts (client decision) | Enables a unified time-course view; batch/time confound disclosed above for 8h/16h |",
  "",

  "## 5. Methods",
  "",
  "| Tool | Version | Parameters |",
  "| :--- | :---: | :---: |",
  "| nf-core/rnaseq (new batch: 8h/16h) | 3.15.1 | --aligner star_salmon |",
  "| Legacy batch (Control/4h) | client-supplied | pre-quantified count matrix, method not independently verified |",
  "| DESeq2 | R package | design = ~Group, lfcShrink type = ashr |",
  "| clusterProfiler | R package | GO ORA (BP/MF/CC), KEGG ORA, GSEA |",
  "| org.Hs.eg.db | R package | Human gene ID mapping |",
  "| msigdbr | R package | Hallmark gene sets (Homo sapiens) |",
  "| Reference genome (new batch) | GRCh38 / GENCODE v45 | — |",
  "",

  "## 6. Results",
  "",
  "### 6.1 Differential Expression",
  "",
  "| Contrast | Total DEGs | Upregulated | Downregulated | Batch status |",
  "| :--- | :---: | :---: | :---: | :---: |",
  {
    sapply(names(deg_summary), function(cn) {
      s <- deg_summary[[cn]]
      batch_note <- if (grepl("^4h", cn)) "Same-batch (reliable)" else "Cross-batch (exploratory — see 4.1)"
      paste0("| ", cn, " | ", s$total, " | ", s$up, " | ", s$down, " | ", batch_note, " |")
    })
  },
  "",
  "### 6.2 Pathway Enrichment (GO BP + KEGG, ALL direction)",
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

  "## 7. Conclusions",
  "",
  {
    lines <- c()
    for (comp_name in names(deg_summary)) {
      s <- deg_summary[[comp_name]]
      batch_caveat <- if (grepl("^4h", comp_name)) "" else " (interpret cautiously — cross-batch comparison, see Section 4.1)"
      lines <- c(lines,
        paste0("- **", comp_name, "**", batch_caveat, ": ", s$total, " DEGs identified. ",
               if (s$up > s$down) paste0("Predominantly upregulated (", s$up, " up vs ", s$down, " down).")
               else if (s$down > s$up) paste0("Predominantly downregulated (", s$down, " down vs ", s$up, " up).")
               else paste0("Balanced bidirectional change (", s$up, " up, ", s$down, " down).")))
    }
    lines
  },
  "- The 4h vs Control result is the only contrast free of batch confound and can anchor confident conclusions",
  "  about the early (4h) transcriptomic response to X-ray radiation in MOLM13 cells.",
  "- The 8h/16h vs Control DEG lists are useful as a candidate/exploratory resource but should be validated",
  "  independently (e.g. qPCR on a shortlist of genes) before being treated as confirmed findings, given the",
  "  batch/time confound.",
  "- Pathway enrichment and GSEA results are available in the `Enrichment/` directory for detailed biological interpretation.",
  "",

  "## 8. Deliverable Files",
  "",
  "| File / Folder | Contents |",
  "| :--- | :---: |",
  paste0("| `DE_PCA_Results/DEG_*.csv` | Full DEG tables (all genes, with log2FC, padj, raw counts) |"),
  paste0("| `DE_PCA_Results/PCA.pdf` | PCA plot, colored by Group and shaped by Batch (visualizes the batch effect) |"),
  paste0("| `DE_PCA_Results/Volcano_*.png` | Volcano plots per contrast |"),
  paste0("| `DE_PCA_Results/Heatmap_top50_*.pdf` | Heatmaps of top 50 DEGs per contrast |"),
  paste0("| `Reads/New_batch_8h_16h_gene_counts.tsv` | Raw count matrix, current NovaSeq batch (8h/16h) |"),
  paste0("| `Reads/New_batch_8h_16h_gene_tpm.tsv` | TPM matrix, current NovaSeq batch (8h/16h) |"),
  paste0("| `Reads/Old_batch_Control_4h_gene_counts.xlsx` | Client-supplied legacy count matrix (Control/4h), copied as-is |"),
  paste0("| `human_Gene_annotation_*.xlsx` | Full human gene annotation with GO/KEGG/UniProt (GENCODE v45) |"),
  paste0("| `Enrichment/*/GO/` | GO ORA results (BP/MF/CC, Up/Down/ALL) with dot plots |"),
  paste0("| `Enrichment/*/KEGG/` | KEGG ORA results with dot plots |"),
  paste0("| `Enrichment/*/GSEA/` | GSEA results (GO BP + KEGG + Hallmark) with ridge/dot plots |"),
  paste0("| `QC/multiqc/` | MultiQC report — covers current NovaSeq batch (8h/16h) only; no QC available for the legacy batch |"),
  "",
  "---",
  "*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*"
)

writeLines(report, con = REPORT_FILE)
cat("Report saved:", REPORT_FILE, "\n")
