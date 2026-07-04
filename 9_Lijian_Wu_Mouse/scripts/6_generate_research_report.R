#!/usr/bin/env Rscript
# 运行环境: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript
# 运行方法: cd /home/gao/projects_2026H2/9_Lijian_Wu_Mouse/scripts && Rscript 6_generate_research_report.R
#
# 前置条件: 4_run_DE_PCA.R + 5_run_enrichment.R 均已完成

library(dplyr)
library(readr)

# ================= 1. 路径设置 =================
setwd("/home/gao/projects_2026H2/9_Lijian_Wu_Mouse/scripts/")

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
  paste0("# Mouse Bulk RNA-seq Analysis Report — Lijian Wu Project"),
  "",
  paste0("**Report Date:** ", format(Sys.Date(), "%B %d, %Y")),
  paste0("**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics"),
  paste0("**Analysis Platform:** Linux HPC server"),
  paste0("**Data folder:** `", basename(DATA_DIR), "`"),
  "",

  "## 1. Objectives",
  "",
  "Characterise transcriptomic differences in cDC1 (conventional dendritic cell type 1) cells",
  "across three treatment conditions (**TNFa**, **Tumor**, **Tumor_TNF**) relative to the",
  "untreated control (**Control**) using bulk RNA-seq. Specific aims:",
  "",
  "- Identify differentially expressed genes (DEGs) for each treatment vs. control contrast.",
  "- Perform GO and KEGG pathway enrichment to determine biological processes affected.",
  "- Run Gene Set Enrichment Analysis (GSEA) for pathway-level signal detection.",
  "",

  "## 2. Key Findings",
  "",
  "This study comprises **4 groups** (TNFa, Tumor, Tumor_TNF = treatment; Control = untreated) across **3 contrasts**.",
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
  "**Top pathway findings per comparison:**",
  "",
  kf_lines,

  "## 3. Sample Information",
  "",
  "| Group | Samples | Role |",
  "| :--- | :---: | :---: |",
  paste0("| Control | cDC1_un_1, cDC1_un_2, cDC1_un_3 | Control (untreated) |"),
  paste0("| TNFa | cDC1_TNFa_1, cDC1_TNFa_2, cDC1_TNFa_3 | Treatment 1 |"),
  paste0("| Tumor | cDC1_Tumor_1, cDC1_Tumor_2, cDC1_Tumor_3 | Treatment 2 |"),
  paste0("| Tumor_TNF | cDC1_Tumor_TNFa_1, cDC1_Tumor_TNFa_2, cDC1_Tumor_TNFa_3 | Treatment 3 |"),
  "",
  paste0("Total samples: **", nrow(meta), "**  |  ",
         "Comparisons: TNFa vs Control, Tumor vs Control, Tumor_TNF vs Control"),
  "",

  "## 4. Analysis Rationale and Decision Criteria",
  "",
  "| Step | Decision | Rationale |",
  "| :--- | :---: | :---: |",
  "| Quantification | Salmon (via nf-core/rnaseq) | Bias-corrected transcript-level quantification |",
  "| STAR alignment | 1-pass, --outFilterMultimapNmax 3 | 7.6× speedup vs default 2-pass on mouse genome; negligible impact on protein-coding DE (Salmon EM handles 2–3 candidate multi-map positions correctly) |",
  {
    if (!is.null(fstats))
      sprintf("| Regex gene filter | Remove ribosomal/non-coding/Gm[0-9] genes | %s → %s genes retained (removed %s noise genes) |",
              format(fstats$n_original, big.mark=","),
              format(fstats$n_after_regex, big.mark=","),
              format(fstats$n_original - fstats$n_after_regex, big.mark=","))
    else
      "| Regex gene filter | Remove ribosomal/non-coding/Gm[0-9] genes | Removes ribosomal, non-coding, and predicted genes |"
  },
  {
    if (!is.null(fstats))
      sprintf("| Low-expression filter | ≥%d counts in ≥%d of %d samples | %s → %s robustly expressed genes input to DESeq2 |",
              fstats$low_count_min, fstats$low_count_min_samples, fstats$n_samples,
              format(fstats$n_after_regex, big.mark=","),
              format(fstats$n_final, big.mark=","))
    else
      "| Low-expression filter | ≥10 counts in n−2 samples | Removes unreliably detected genes; improves statistical power |"
  },
  "| DE threshold | padj ≤ 0.05, \\|log2FC\\| ≥ 0.263 (= 1.2×) | Conservative fold-change avoids over-reporting small but significant changes |",
  "| LFC shrinkage | ashr | Appropriate for unbalanced designs and small sample sizes |",
  "",

  "## 5. Methods",
  "",
  "| Tool | Version | Parameters |",
  "| :--- | :---: | :---: |",
  "| nf-core/rnaseq | 3.15.1 | --aligner star_salmon |",
  "| STAR | 2.7.x | --twopassMode None --outFilterMultimapNmax 3 |",
  "| Salmon | — | default |",
  "| DESeq2 | R package | design = ~Group, lfcShrink type = ashr |",
  "| clusterProfiler | R package | GO ORA (BP/MF/CC), KEGG ORA, GSEA |",
  "| org.Mm.eg.db | R package | Mouse gene ID mapping |",
  "| msigdbr | R package | Hallmark gene sets (MM) |",
  "| Reference genome | GRCm39 / GENCODE M35 | — |",
  "",

  "## 6. Results",
  "",
  "### 6.1 Differential Expression",
  "",
  "| Contrast | Total DEGs | Upregulated | Downregulated |",
  "| :--- | :---: | :---: | :---: |",
  {
    sapply(names(deg_summary), function(cn) {
      s <- deg_summary[[cn]]
      paste0("| ", cn, " | ", s$total, " | ", s$up, " | ", s$down, " |")
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
      lines <- c(lines,
        paste0("- **", comp_name, "**: ", s$total, " DEGs identified. ",
               if (s$up > s$down) paste0("Predominantly upregulated (", s$up, " up vs ", s$down, " down), ",
                 "suggesting activation of transcriptional programs relative to untreated control.")
               else if (s$down > s$up) paste0("Predominantly downregulated (", s$down, " down vs ", s$up, " up), ",
                 "suggesting suppression of gene expression relative to untreated control.")
               else paste0("Balanced bidirectional change (", s$up, " up, ", s$down, " down).")))
    }
    lines
  },
  "- Pathway enrichment and GSEA results are available in the `Enrichment/` directory for detailed biological interpretation.",
  "",

  "## 8. Deliverable Files",
  "",
  "| File / Folder | Contents |",
  "| :--- | :---: |",
  paste0("| `DE_PCA_Results/DEG_*.csv` | Full DEG tables (all genes, with log2FC, padj, raw counts) |"),
  paste0("| `DE_PCA_Results/PCA.pdf` | PCA plot |"),
  paste0("| `DE_PCA_Results/Volcano_*.png` | Volcano plots per contrast |"),
  paste0("| `DE_PCA_Results/Heatmap_top50_*.pdf` | Heatmaps of top 50 DEGs per contrast |"),
  paste0("| `Reads/All_sample_gene_counts.tsv` | Raw count matrix |"),
  paste0("| `Reads/All_sample_gene_tpm.tsv` | TPM matrix |"),
  paste0("| `mouse_Gene_annotation_*.xlsx` | Full mouse gene annotation with GO/KEGG/UniProt (GENCODE M35) |"),
  paste0("| `Enrichment/*/GO/` | GO ORA results (BP/MF/CC, Up/Down/ALL) with dot plots |"),
  paste0("| `Enrichment/*/KEGG/` | KEGG ORA results with dot plots |"),
  paste0("| `Enrichment/*/GSEA/` | GSEA results (GO BP + KEGG + Hallmark) with ridge/dot plots |"),
  paste0("| `QC/multiqc/` | MultiQC report (sequencing QC, alignment stats) |"),
  "",
  "---",
  "*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*"
)

writeLines(report, con = REPORT_FILE)
cat("Report saved:", REPORT_FILE, "\n")
