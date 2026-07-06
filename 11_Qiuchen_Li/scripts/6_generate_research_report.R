#!/usr/bin/env Rscript
# 运行环境: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript
# 运行方法: cd /home/gao/projects_2026H2/11_Qiuchen_Li/scripts && Rscript 6_generate_research_report.R
#
# 前置条件: 4_run_DE_PCA.R + 5_run_enrichment.R 均已完成

library(dplyr)
library(readr)

# ================= 1. 路径设置 =================
setwd("/home/gao/projects_2026H2/11_Qiuchen_Li/scripts/")

data_dirs <- sort(list.dirs("..", full.names = TRUE, recursive = FALSE), decreasing = TRUE)
data_dirs <- data_dirs[grepl("/Data_Analysis_[0-9]{8}$", data_dirs)]
if (length(data_dirs) == 0) stop("No Data_Analysis_YYYYMMDD folder found. Run 4_run_DE_PCA.R first.")
DATA_DIR <- data_dirs[1]
cat("Using:", DATA_DIR, "\n")

DE_DIR  <- "rds_cache"  # internal handoff objects (not the client-facing DE_PCA_Results/ folder)
ENR_DIR <- file.path(DATA_DIR, "Enrichment")

REPORT_DATE <- format(Sys.Date(), "%m%d")
REPORT_FILE <- file.path(DATA_DIR, paste0("Bioinformatics_Analysis_Report_", REPORT_DATE, ".md"))

# ================= 2. 加载 DE 结果 =================
res_list <- readRDS(file.path(DE_DIR, "res_list.rds"))
meta     <- readRDS(file.path(DE_DIR, "meta.rds"))
sig_col  <- "sig (padj<=0.05 & |log2FC|>=0.263)"
fs_file  <- file.path(DE_DIR, "filter_stats.rds")
fstats   <- if (file.exists(fs_file)) readRDS(fs_file) else NULL
pca_file <- file.path(DE_DIR, "pca_summary.rds")
pca      <- if (file.exists(pca_file)) readRDS(pca_file) else NULL

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
  paste0("# Human Bulk RNA-seq Analysis Report — Qiuchen Li Project"),
  "",
  paste0("**Report Date:** ", format(Sys.Date(), "%B %d, %Y")),
  paste0("**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics"),
  paste0("**Analysis Platform:** Linux HPC server"),
  paste0("**Data folder:** `", basename(DATA_DIR), "`"),
  "",

  "## 1. Objectives",
  "",
  "Characterise transcriptomic differences among three experimental groups",
  "(**Mix**, **NT**, **A5BKO**) using bulk RNA-seq, via all three pairwise contrasts",
  "(client-requested pairwise design — no single designated control). Specific aims:",
  "",
  "- Identify differentially expressed genes (DEGs) for each of the 3 pairwise contrasts.",
  "- Perform GO and KEGG pathway enrichment to determine biological processes affected.",
  "- Run Gene Set Enrichment Analysis (GSEA) for pathway-level signal detection.",
  "",

  "## 2. Key Findings",
  "",
  "This study comprises **3 groups** (Mix, NT, A5BKO) across **3 pairwise contrasts**.",
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

  "### 2.1 A recurring signal across all three comparisons: MYC target gene programs",
  "",
  "All three pairwise contrasts independently identify a **MYC target gene signature** (MSigDB",
  "Hallmark collection) as their single most significant GSEA hit — Mix vs NT and A5BKO vs Mix",
  "both top out on **MYC_TARGETS_V1**, while A5BKO vs NT tops out on the related **MYC_TARGETS_V2**",
  "set. This is a genuine, independently-computed result in each contrast (confirmed by inspecting",
  "the underlying leading-edge gene lists, which differ between comparisons), not a reporting artefact.",
  "",
  "**V1 vs V2 — two related but distinct gene sets:** `HALLMARK_MYC_TARGETS_V1` (200 genes) is the",
  "broader, classic MYC target signature, drawing heavily on ribosomal-protein and translation/",
  "ribosome-biogenesis genes that MYC is known to transcriptionally activate. `HALLMARK_MYC_TARGETS_V2`",
  "(58 genes) is a smaller, more stringent signature enriched for core cell-cycle and DNA-replication",
  "genes more directly bound by MYC. The two sets overlap partially but capture different facets of",
  "MYC-driven transcriptional activity (broad translational output vs. core proliferation machinery).",
  "",
  "**Why this is noteworthy:** the direction of enrichment is internally consistent and transitive",
  "across all three contrasts — NT scores higher than Mix, A5BKO scores higher than NT, and A5BKO",
  "scores higher than Mix for their respective MYC signatures. That ordering (**A5BKO > NT > Mix**)",
  "is exactly what would be expected if the same underlying biological axis (MYC-driven proliferation/",
  "translation activity) is the dominant driver of transcriptomic variation across all three groups,",
  "rather than three unrelated, coincidental hits. We'd recommend treating MYC pathway activity as a",
  "candidate primary axis of biological variation in this experiment, worth confirming with a",
  "marker-gene qPCR panel or by correlating a MYC target module score against the client's expected",
  "experimental design.",
  "",

  "## 3. Sample Information",
  "",
  "| Group | Samples | Role |",
  "| :--- | :---: | :---: |",
  paste0("| Mix | Mix_1, Mix_2, Mix_3 | Experimental group |"),
  paste0("| NT | NT_1, NT_2, NT_3 | Experimental group |"),
  paste0("| A5BKO | A5BKO_1, A5BKO_3 | Experimental group (n=2, see note below) |"),
  "",
  paste0("Total samples analysed: **", nrow(meta), "**  |  ",
         "Comparisons: Mix vs NT, A5BKO vs NT, A5BKO vs Mix"),
  "",
  "**Sample exclusion:** sample **A5BKO_2** was excluded from this analysis. PCA placed it",
  "essentially on top of the three NT replicates (far from its own A5BKO_1/A5BKO_3 replicates,",
  "which cluster tightly together) — a pattern consistent with sample mislabeling rather than",
  "ordinary biological/technical replicate variation. The A5BKO group therefore has n=2 in this",
  "analysis; all A5BKO-related contrasts below reflect A5BKO_1 and A5BKO_3 only.",
  "",

  "## 4. Analysis Rationale and Decision Criteria",
  "",
  "| Step | Decision | Rationale |",
  "| :--- | :---: | :---: |",
  "| Quantification | Salmon (via nf-core/rnaseq) | Bias-corrected transcript-level quantification |",
  {
    if (!is.null(fstats))
      sprintf("| Biotype gene filter | Keep gene_type == protein_coding (GENCODE annotation) | %s → %s genes retained (removed %s non-coding/pseudogenes) |",
              format(fstats$n_original, big.mark=","),
              format(fstats$n_after_regex, big.mark=","),
              format(fstats$n_original - fstats$n_after_regex, big.mark=","))
    else
      "| Biotype gene filter | Keep protein_coding genes | Removes non-coding RNA, pseudogenes |"
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
  "| Comparison design | All 3 pairwise contrasts | Client requested full pairwise comparison (两两对比); no single reference group designated |",
  "",

  "## 5. Methods",
  "",
  "| Tool | Version | Parameters |",
  "| :--- | :---: | :---: |",
  "| nf-core/rnaseq | 3.15.1 | --aligner star_salmon |",
  "| STAR | 2.7.x | default (2-pass) |",
  "| Salmon | — | default |",
  "| DESeq2 | R package | design = ~Group, lfcShrink type = ashr |",
  "| clusterProfiler | R package | GO ORA (BP/MF/CC), KEGG ORA, GSEA |",
  "| org.Hs.eg.db | R package | Human gene ID mapping |",
  "| msigdbr | R package | Hallmark gene sets (Homo sapiens) |",
  "| Reference genome | GRCh38 / GENCODE v45 | — |",
  "",

  "## 6. Results",
  "",
  "### 6.1 PCA Overview",
  "",
  {
    if (!is.null(pca)) c(
      sprintf("PC1 explains **%d%%** of total variance, PC2 explains **%d%%**.",
              pca$percentVar[1], pca$percentVar[2]),
      "",
      sprintf("Average within-group distance (PC1–PC2 space): **%.2f**; average between-group distance: **%.2f** (separation ratio **%.2fx**).",
              pca$within_dist, pca$between_dist, pca$separation_ratio),
      "",
      if (pca$separation_ratio > 2)
        "Replicates cluster tightly within each group, and the three groups are clearly separated from one another — indicating good replicate reproducibility and a strong, consistent transcriptomic effect between groups (this is with the mislabeled A5BKO_2 sample already excluded — see Section 3)."
      else if (pca$separation_ratio > 1)
        "Groups show a moderate degree of separation, with some overlap between conditions — group differences are detectable but not the dominant source of variance in the dataset."
      else
        "Groups show substantial overlap in PCA space, with within-group variability comparable to or exceeding between-group differences."
    ) else "PCA summary data not available."
  },
  "",
  "### 6.2 Differential Expression",
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
  "### 6.3 Pathway Enrichment (GO BP + KEGG, ALL direction)",
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
               if (s$up > s$down) paste0("Predominantly upregulated (", s$up, " up vs ", s$down, " down).")
               else if (s$down > s$up) paste0("Predominantly downregulated (", s$down, " down vs ", s$up, " up).")
               else paste0("Balanced bidirectional change (", s$up, " up, ", s$down, " down).")))
    }
    lines
  },
  "- Pathway enrichment and GSEA results are available in the `Enrichment/` directory for detailed biological interpretation.",
  "- **Sample A5BKO_2 is highly likely mislabeled or affected by an experimental handling error**",
  "  upstream of sequencing (e.g. sample swap, tube mix-up, or contamination during processing):",
  "  its transcriptome is essentially indistinguishable from the NT group and unrelated to its own",
  "  A5BKO replicates (Section 6.1). We recommend the client cross-check this sample's collection,",
  "  labeling, and processing records before relying on any other data associated with it (e.g. if",
  "  it was part of a larger batch, other samples processed alongside it may warrant a similar check).",
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
  paste0("| `human_Gene_annotation_*.xlsx` | Full human gene annotation with GO/KEGG/UniProt (GENCODE v45) |"),
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
