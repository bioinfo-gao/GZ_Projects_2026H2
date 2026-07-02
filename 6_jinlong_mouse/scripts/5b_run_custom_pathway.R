#!/usr/bin/env Rscript
# 运行环境: DE_R45
# 前置条件: 先运行 5_run_enrichment.R 生成 Enrichment_Standard/ 结果
#
# 分析内容 (客户指定定制通路):
#   A. Stem cell marker gene analysis
#   B. Cell differentiation & growth pathway analysis
#   C. Notch signalling pathway analysis

library(clusterProfiler) 
library(enrichplot)
library(org.Mm.eg.db)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(readr)

# ================= 1. Path setup =================
setwd("/home/gao/projects_2026H2/6_jinlong_mouse/scripts/")

# Automatically pick the most recent Data_Analysis_YYYYMMDD folder
# (sorted descending so [1] is always the latest date)
data_dirs <- sort(list.dirs("..", full.names = TRUE, recursive = FALSE), decreasing = TRUE)
data_dirs <- data_dirs[grepl("/Data_Analysis_[0-9]{8}$", data_dirs)]
if (length(data_dirs) == 0) stop("No Data_Analysis_YYYYMMDD folder found.")
DATA_DIR <- data_dirs[1]
cat("Using:", DATA_DIR, "\n")

DE_DIR  <- file.path(DATA_DIR, "DE_PCA_Results")
STD_DIR <- file.path(DATA_DIR, "Enrichment_Standard")                  # read-only: extract pre-computed KEGG results
ENR_DIR <- file.path(DATA_DIR, "Enrichment_Custom_Designed_Pathways")  # write: all custom pathway outputs
dir.create(ENR_DIR, showWarnings = FALSE, recursive = TRUE)

# ================= 2. Load DE results =================
res_list   <- readRDS(file.path(DE_DIR, "res_list.rds"))
meta       <- readRDS(file.path(DE_DIR, "meta.rds"))
counts_raw <- readRDS(file.path(DE_DIR, "counts_raw.rds"))

# Must match the column name written by 4_run_DE_PCA.R exactly
sig_col <- "sig (padj<=0.05 & |log2FC|>=0.263)"
cat("✅ DE results loaded. Comparisons:", paste(names(res_list), collapse = ", "), "\n")

# ================= 3. Main loop: one comparison at a time =================
for (comp_name in names(res_list)) {
  cat("\n", strrep("=", 60), "\n")
  cat("Processing:", comp_name, "\n")
  cat(strrep("=", 60), "\n")

  comp_dir <- file.path(ENR_DIR, comp_name)
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

  df <- res_list[[comp_name]]
  # Drop rows with NA padj / log2FC (low-count genes filtered by DESeq2)
  df_clean <- df %>% filter(!is.na(padj), !is.na(log2FoldChange))
  # Custom pathway ORA uses all DEGs combined (Up + Down); no directional split needed
  sig_all  <- df_clean %>% filter(.data[[sig_col]] != "NS")

  # Build Ensembl→Entrez map from all expressed genes (used as ORA universe)
  # suppressMessages: bitr prints a warning for unmapped IDs, which is expected and harmless
  id_map <- suppressMessages(
    bitr(sub("\\..*$", "", df_clean$gene_id),
         fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  )
  # Helper closure: look up Entrez IDs for any sub-dataframe using the parent-scope id_map
  get_entrez <- function(sub_df) {
    clean <- sub("\\..*$", "", sub_df$gene_id)
    id_map$ENTREZID[id_map$ENSEMBL %in% clean]
  }
  entrez_all <- get_entrez(sig_all)
  entrez_bg  <- id_map$ENTREZID   # ORA universe: all detected genes mappable to Entrez

  cat("  Total sig DEGs for ORA:", length(entrez_all), "\n")

  # ---- A. Stem cell marker gene analysis ----
  # Strategy: pull DE statistics for a curated list of known stem cell markers,
  # then run MSigDB C8 ORA as an unbiased complement.
  cat("\n--- Stem Cell Marker Analysis ---\n")
  sc_dir <- file.path(comp_dir, "StemCell")
  dir.create(sc_dir, showWarnings = FALSE, recursive = TRUE)

  # Six functional categories covering major stem cell types in mouse
  stemcell_markers <- list(
    Pluripotency     = c("Pou5f1","Sox2","Nanog","Klf4","Myc","Zfp42","Esrrb","Utf1","Sall4","Lin28a"),
    Neural_SC        = c("Nes","Sox1","Sox9","Pax6","Vim","Gfap","Msi1","Prom1","Blbp"),
    Hematopoietic_SC = c("Cd34","Ly6a","Slamf1","Kit","Cd48","Procr","Gata2","Tal1"),
    Mesenchymal_SC   = c("Cd44","Nt5e","Thy1","Eng","Itgb1","Pdgfra","Aldh1a1"),
    Intestinal_SC    = c("Lgr5","Axin2","Olfm4","Ascl2","Lrig1","Troy"),
    General_Stemness = c("Aldh1a1","Aldh1a3","Epcam","Cd24a","Sox4","Id1","Id3","Notch1","Notch2")
  )
  # Flatten and deduplicate: some genes (e.g. Aldh1a1) appear in multiple categories
  all_markers <- unique(unlist(stemcell_markers))

  # Pull DE results for marker genes; join category labels for downstream faceting
  marker_de <- df %>%
    filter(gene_name %in% all_markers) %>%
    select(gene_id, gene_name, log2FoldChange, padj, .data[[sig_col]]) %>%
    arrange(padj)
  marker_category <- data.frame(
    gene_name = unlist(stemcell_markers),
    Category  = rep(names(stemcell_markers), sapply(stemcell_markers, length)),
    stringsAsFactors = FALSE
  )
  marker_de <- left_join(marker_de, marker_category, by = "gene_name")
  write_csv(marker_de, file.path(sc_dir, "StemCell_Markers_DE_Results.csv"))
  cat("  ✅ Stem cell markers found:", nrow(marker_de), "/", length(all_markers), "\n")

  sig_markers <- marker_de %>% filter(.data[[sig_col]] != "NS")
  if (nrow(sig_markers) > 0) {
    cat("  ✅ Significant:", nrow(sig_markers), "\n")
    p_bar <- marker_de %>%
      filter(!is.na(log2FoldChange)) %>%
      arrange(Category, log2FoldChange) %>%
      mutate(gene_name = factor(gene_name, levels = unique(gene_name)),
             Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
      ggplot(aes(x = gene_name, y = log2FoldChange, fill = Direction)) +
      geom_col() +
      geom_hline(yintercept = c(-0.263, 0.263), linetype = "dashed", color = "grey50") +  # significance threshold lines
      facet_wrap(~Category, scales = "free_x", nrow = 2) +   # free_x: each panel shows only its own genes
      scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            strip.text  = element_text(size = 8)) +
      labs(title = paste("Stem Cell Marker Genes:", comp_name), x = NULL, y = "log2 Fold Change")
    tryCatch(
      ggsave(file.path(sc_dir, "StemCell_Markers_log2FC.pdf"), p_bar, width = 14, height = 7, dpi = 300),
      error = function(e) cat("  ⚠️  StemCell bar plot failed:", e$message, "\n")
    )
  }

  # MSigDB C8 ORA: C8 = cell type signature gene sets; filter for stem/progenitor/pluripotency sets
  tryCatch({
    msig_c8 <- msigdbr(species = "Mus musculus", category = "C8") %>%
      filter(grepl("STEM|PROGENITOR|PLURIP", gs_name, ignore.case = TRUE)) %>%
      select(gs_name, entrez_gene) %>% rename(term = gs_name, gene = entrez_gene)
    if (nrow(msig_c8) > 0 && length(entrez_all) >= 5) {
      enr_sc <- enricher(gene = as.character(entrez_all), universe = as.character(entrez_bg),
                         TERM2GENE = msig_c8 %>% mutate(gene = as.character(gene)),
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      if (!is.null(enr_sc) && nrow(as.data.frame(enr_sc)) > 0) {
        df_sc <- as.data.frame(enr_sc)
        write_csv(df_sc, file.path(sc_dir, "MSigDB_StemCell_Signatures.csv"))
        p_sc <- dotplot(enr_sc, showCategory = min(15, nrow(df_sc)),
                        title = paste("MSigDB Stem Cell Signatures:", comp_name)) +
          theme(axis.text.y = element_text(size = 7))
        ggsave(file.path(sc_dir, "MSigDB_StemCell_dotplot.pdf"), p_sc, width = 10, height = 6, dpi = 300)
        cat("  ✅ MSigDB stem cell ORA:", nrow(df_sc), "gene sets\n")
      }
    }
  }, error = function(e) cat("  ⚠️  MSigDB StemCell ORA error:", e$message, "\n"))

  # ---- B. Cell Differentiation & Growth ----
  # Strategy: curated marker bar chart (same approach as Section A) +
  # MSigDB C5 GO:BP ORA restricted to differentiation/proliferation/growth/cycle/development terms.
  cat("\n--- Cell Differentiation & Growth Analysis ---\n")
  diff_dir <- file.path(comp_dir, "CellDiff")
  dir.create(diff_dir, showWarnings = FALSE, recursive = TRUE)

  # Five functional categories covering the major axes of cell growth and fate decisions
  diff_markers <- list(
    Proliferation      = c("Mki67","Pcna","Top2a","Ccnd1","Ccnd2","Ccne1","Cdk4","Cdk6","Cdk2","Myc","E2f1","Mcm2"),
    Growth_Factors     = c("Igf1","Igf1r","Egfr","Fgfr1","Fgfr2","Met","Pdgfra","Pdgfrb","Vegfa","Tgfb1","Tgfb2","Tgfb3"),
    Differentiation_TF = c("Myod1","Myog","Runx1","Runx2","Sox9","Pparg","Cebpa","Gata1","Gata2","Klf4","Klf5","Atoh1"),
    Wnt_Signaling      = c("Ctnnb1","Wnt5a","Wnt3a","Wnt7a","Fzd1","Fzd4","Axin2","Apc","Lef1","Tcf7","Wnt2b"),
    EMT_Markers        = c("Vim","Cdh1","Cdh2","Fn1","Twist1","Snai1","Snai2","Zeb1","Zeb2","Col1a1","Col3a1")
  )
  all_diff <- unique(unlist(diff_markers))

  diff_de <- df %>%
    filter(gene_name %in% all_diff) %>%
    select(gene_id, gene_name, log2FoldChange, padj, .data[[sig_col]]) %>%
    arrange(padj)
  diff_category <- data.frame(
    gene_name = unlist(diff_markers),
    Category  = rep(names(diff_markers), sapply(diff_markers, length)),
    stringsAsFactors = FALSE
  )
  diff_de <- left_join(diff_de, diff_category, by = "gene_name")
  write_csv(diff_de, file.path(diff_dir, "CellDiff_Markers_DE_Results.csv"))
  cat("  ✅ Cell diff/growth markers found:", nrow(diff_de), "/", length(all_diff), "\n")

  sig_diff <- diff_de %>% filter(.data[[sig_col]] != "NS")
  if (nrow(sig_diff) > 0) {
    cat("  ✅ Significant:", nrow(sig_diff), "\n")
    p_diff <- diff_de %>%
      filter(!is.na(log2FoldChange)) %>%
      arrange(Category, log2FoldChange) %>%
      mutate(gene_name = factor(gene_name, levels = unique(gene_name)),
             Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
      ggplot(aes(x = gene_name, y = log2FoldChange, fill = Direction)) +
      geom_col() +
      geom_hline(yintercept = c(-0.263, 0.263), linetype = "dashed", color = "grey50") +  # significance threshold lines
      facet_wrap(~Category, scales = "free_x", nrow = 2) +
      scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            strip.text  = element_text(size = 8)) +
      labs(title = paste("Cell Differentiation & Growth Markers:", comp_name), x = NULL, y = "log2 Fold Change")
    tryCatch(
      ggsave(file.path(diff_dir, "CellDiff_Markers_log2FC.pdf"), p_diff, width = 16, height = 7, dpi = 300),
      error = function(e) cat("  ⚠️  CellDiff bar plot failed:", e$message, "\n")
    )
  } else {
    cat("  ℹ️  No significant cell diff/growth markers\n")
  }

  # MSigDB C5 GO:BP ORA: C5 = ontology gene sets; restrict to differentiation/proliferation-related terms
  # gs_name pattern match avoids running the full 7,000+ GO:BP terms (noise reduction)
  tryCatch({
    msig_diff <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>%
      filter(grepl("DIFFERENTIATION|PROLIFERATION|CELL_GROWTH|CELL_CYCLE|DEVELOPMENT", gs_name)) %>%
      select(gs_name, entrez_gene) %>% rename(term = gs_name, gene = entrez_gene)
    if (nrow(msig_diff) > 0 && length(entrez_all) >= 5) {
      enr_diff <- enricher(gene = as.character(entrez_all), universe = as.character(entrez_bg),
                           TERM2GENE = msig_diff %>% mutate(gene = as.character(gene)),
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      if (!is.null(enr_diff) && nrow(as.data.frame(enr_diff)) > 0) {
        df_d <- as.data.frame(enr_diff)
        write_csv(df_d, file.path(diff_dir, "MSigDB_CellDiff_ORA.csv"))
        p_dd <- dotplot(enr_diff, showCategory = min(20, nrow(df_d)),
                        title = paste("Cell Diff & Growth (MSigDB GO:BP):", comp_name)) +
          theme(axis.text.y = element_text(size = 7))
        ggsave(file.path(diff_dir, "MSigDB_CellDiff_dotplot.pdf"), p_dd,
               width = 12, height = max(5, min(20, nrow(df_d)) * 0.35 + 3), dpi = 300)
        cat("  ✅ MSigDB CellDiff ORA:", nrow(df_d), "gene sets\n")
      } else {
        cat("  ℹ️  No significant CellDiff gene sets (MSigDB)\n")
      }
    }
  }, error = function(e) cat("  ⚠️  MSigDB CellDiff ORA error:", e$message, "\n"))

  # ---- C. Notch Pathway ----
  # Three-layer analysis:
  #   1. Curated marker bar chart (same approach as A/B)
  #   2. Extract mmu04330 entry from pre-computed KEGG ORA in Enrichment_Standard/
  #      (avoids re-running enrichKEGG which requires an internet KEGG API call)
  #   3. MSigDB ORA across all collections filtered for "NOTCH" gene sets
  cat("\n--- Notch Pathway Analysis ---\n")
  notch_dir <- file.path(comp_dir, "Notch")
  dir.create(notch_dir, showWarnings = FALSE, recursive = TRUE)

  # Notch pathway components organised by functional role in the signaling cascade
  notch_markers <- list(
    Receptors      = c("Notch1","Notch2","Notch3","Notch4"),
    Ligands        = c("Dll1","Dll3","Dll4","Jag1","Jag2"),
    Core_Complex   = c("Rbpj","Maml1","Maml2","Maml3","Ep300","Hdac1","Hdac2"),  # transcriptional activation complex
    Target_Genes   = c("Hes1","Hes5","Hes6","Hey1","Hey2","Heyl","Nrarp","Myc","Ccnd1","Cdkn1a"),
    Neg_Regulators = c("Numb","Numbl","Fbxw7","Mfng","Lfng","Rfng","Deltex1","Itch")
  )
  all_notch <- unique(unlist(notch_markers))

  notch_de <- df %>%
    filter(gene_name %in% all_notch) %>%
    select(gene_id, gene_name, log2FoldChange, padj, .data[[sig_col]]) %>%
    arrange(padj)
  notch_category <- data.frame(
    gene_name = unlist(notch_markers),
    Category  = rep(names(notch_markers), sapply(notch_markers, length)),
    stringsAsFactors = FALSE
  )
  notch_de <- left_join(notch_de, notch_category, by = "gene_name")
  write_csv(notch_de, file.path(notch_dir, "Notch_Pathway_DE_Results.csv"))
  cat("  ✅ Notch genes found:", nrow(notch_de), "/", length(all_notch), "\n")

  sig_notch <- notch_de %>% filter(.data[[sig_col]] != "NS")
  if (nrow(sig_notch) > 0) {
    cat("  ✅ Significant:", nrow(sig_notch), "\n")
    print(sig_notch %>% select(gene_name, Category, log2FoldChange, padj, .data[[sig_col]]))  # console inspection
    p_notch <- notch_de %>%
      filter(!is.na(log2FoldChange)) %>%
      arrange(Category, log2FoldChange) %>%
      mutate(gene_name = factor(gene_name, levels = unique(gene_name)),
             Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
      ggplot(aes(x = gene_name, y = log2FoldChange, fill = Direction)) +
      geom_col() +
      geom_hline(yintercept = c(-0.263, 0.263), linetype = "dashed", color = "grey50") +  # significance threshold lines
      facet_wrap(~Category, scales = "free_x", nrow = 2) +
      scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            strip.text  = element_text(size = 8)) +
      labs(title = paste("Notch Pathway Genes:", comp_name), x = NULL, y = "log2 Fold Change")
    tryCatch(
      ggsave(file.path(notch_dir, "Notch_Pathway_log2FC.pdf"), p_notch, width = 14, height = 7, dpi = 300),
      error = function(e) cat("  ⚠️  Notch bar plot failed:", e$message, "\n")
    )
  } else {
    cat("  ℹ️  No significant Notch pathway genes\n")
  }

  # Extract Notch entry from pre-computed KEGG ORA (avoids a redundant enrichKEGG internet call)
  # mmu04330 = KEGG Notch signaling pathway ID for Mus musculus
  kegg_all_csv <- file.path(STD_DIR, comp_name, "KEGG", "KEGG_ALL.csv")
  if (file.exists(kegg_all_csv)) {
    kegg_all_df <- read_csv(kegg_all_csv, show_col_types = FALSE)
    notch_kegg  <- kegg_all_df %>%
      filter(ID == "mmu04330" | grepl("Notch", Description, ignore.case = TRUE))
    if (nrow(notch_kegg) > 0) {
      write_csv(notch_kegg, file.path(notch_dir, "KEGG_Notch_from_Standard.csv"))
      cat("  ✅ Notch in KEGG ORA:", notch_kegg$Description[1],
          "padj =", signif(notch_kegg$p.adjust[1], 3), "\n")
    } else {
      cat("  ℹ️  Notch (mmu04330) not significant in standard KEGG ORA\n")
    }
  } else {
    cat("  ℹ️  Enrichment_Standard KEGG results not found — run 5_run_enrichment.R first\n")
  }

  # MSigDB ORA: search all MSigDB collections (no category filter) for gene sets named with "NOTCH"
  tryCatch({
    msig_notch <- msigdbr(species = "Mus musculus") %>%
      filter(grepl("NOTCH", gs_name)) %>%
      select(gs_name, entrez_gene) %>% rename(term = gs_name, gene = entrez_gene)
    if (nrow(msig_notch) > 0 && length(entrez_all) >= 5) {
      enr_notch <- enricher(gene = as.character(entrez_all), universe = as.character(entrez_bg),
                            TERM2GENE = msig_notch %>% mutate(gene = as.character(gene)),
                            pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      if (!is.null(enr_notch) && nrow(as.data.frame(enr_notch)) > 0) {
        df_n <- as.data.frame(enr_notch)
        write_csv(df_n, file.path(notch_dir, "MSigDB_Notch_ORA.csv"))
        p_nm <- dotplot(enr_notch, showCategory = min(15, nrow(df_n)),
                        title = paste("Notch Signaling (MSigDB):", comp_name)) +
          theme(axis.text.y = element_text(size = 7))
        ggsave(file.path(notch_dir, "MSigDB_Notch_dotplot.pdf"), p_nm,
               width = 10, height = max(5, min(15, nrow(df_n)) * 0.35 + 3), dpi = 300)
        cat("  ✅ MSigDB Notch ORA:", nrow(df_n), "gene sets\n")
      } else {
        cat("  ℹ️  No significant Notch gene sets (MSigDB)\n")
      }
    }
  }, error = function(e) cat("  ⚠️  MSigDB Notch ORA error:", e$message, "\n"))

  cat("✅ Completed custom pathway analysis for:", comp_name, "\n")
}

# ================= 4. Cross-comparison summary =================
# Merge per-comparison marker CSVs into a single file per pathway type,
# adding a "Comparison" column so the reader can filter by G1_vs_G4, G2_vs_G4, etc.
# Output goes to ENR_DIR root (not inside a comparison subfolder) for easy delivery.
cat("\n", strrep("=", 60), "\n")
cat("Cross-comparison summary\n")
cat(strrep("=", 60), "\n")

for (summary_info in list(
  list(name = "StemCell",  csv = "StemCell_Markers_DE_Results.csv",  out = "StemCell_AllComparisons_Summary.csv"),
  list(name = "CellDiff",  csv = "CellDiff_Markers_DE_Results.csv",  out = "CellDiff_AllComparisons_Summary.csv"),
  list(name = "Notch",     csv = "Notch_Pathway_DE_Results.csv",     out = "Notch_AllComparisons_Summary.csv")
)) {
  rows <- list()
  for (comp_name in names(res_list)) {
    f <- file.path(ENR_DIR, comp_name, summary_info$name, summary_info$csv)
    if (file.exists(f)) rows[[comp_name]] <- read_csv(f, show_col_types = FALSE) %>% mutate(Comparison = comp_name)
  }
  if (length(rows) > 0) {
    write_csv(bind_rows(rows), file.path(ENR_DIR, summary_info$out))
    cat("✅ Cross-comparison", summary_info$name, "summary saved\n")
  }
}

cat("\n🎉 Custom pathway analyses complete. Results:", ENR_DIR, "\n")
