#!/usr/bin/env Rscript
# 运行环境: conda activate DE_R45
# 运行方法: cd /home/gao/projects_2026H2/6_jinlong_mouse/scripts && Rscript 5_run_enrichment.R
#
# 分析内容:
#   1. GO enrichment analysis  (BP / MF / CC)
#   2. KEGG pathway analysis
#   3. GSEA (Gene Set Enrichment Analysis)  — GO BP + KEGG + MSigDB Hallmark
#   4. Stem cell marker gene analysis
#
# 前置条件: 先运行 4_run_DE_PCA.R 生成 res_list.rds / meta.rds / counts_raw.rds

library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# ================= 1. 路径设置 =================
setwd("/home/gao/projects_2026H2/6_jinlong_mouse/scripts/")

# 自动找最新的 Data_Analysis_YYYYMMDD 文件夹
data_dirs <- sort(list.dirs("..", full.names = TRUE, recursive = FALSE), decreasing = TRUE)
data_dirs <- data_dirs[grepl("/Data_Analysis_[0-9]{8}$", data_dirs)]
if (length(data_dirs) == 0) stop("No Data_Analysis_YYYYMMDD folder found. Run 4_run_DE_PCA.R first.")
DATA_DIR <- data_dirs[1]
cat("Using:", DATA_DIR, "\n")

DE_DIR   <- file.path(DATA_DIR, "DE_PCA_Results")
ENR_DIR  <- file.path(DATA_DIR, "Enrichment")
dir.create(ENR_DIR, showWarnings = FALSE, recursive = TRUE)

# ================= 2. 加载 DE 结果 =================
res_list   <- readRDS(file.path(DE_DIR, "res_list.rds"))
meta       <- readRDS(file.path(DE_DIR, "meta.rds"))
counts_raw <- readRDS(file.path(DE_DIR, "counts_raw.rds"))

sig_col <- "sig (padj<=0.05 & |log2FC|>=0.263)"

cat("✅ DE results loaded. Comparisons:", paste(names(res_list), collapse=", "), "\n")

# ================= 3. 辅助函数 =================

# Ensembl → Entrez ID 转换（鼠: org.Mm.eg.db）
ensembl_to_entrez <- function(ensembl_ids) {
  # 去掉版本号 (ENSMUSG00000051951.6 → ENSMUSG00000051951)
  clean_ids <- sub("\\..*$", "", ensembl_ids)
  map <- bitr(clean_ids, fromType = "ENSEMBL", toType = "ENTREZID",
              OrgDb = org.Mm.eg.db)
  return(map)
}

# 保存富集结果 CSV + dotplot
save_enrichment <- function(enr_res, prefix, out_dir, n_show = 20) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(enr_res) || nrow(as.data.frame(enr_res)) == 0) {
    cat("  ⚠️  No significant results for:", prefix, "\n")
    return(invisible(NULL))
  }

  df <- as.data.frame(enr_res)
  write_csv(df, file.path(out_dir, paste0(prefix, ".csv")))

  p <- dotplot(enr_res, showCategory = min(n_show, nrow(df)), title = prefix) +
    theme(axis.text.y = element_text(size = 7))
  ggsave(file.path(out_dir, paste0(prefix, "_dotplot.pdf")), p,
         width = 10, height = max(5, min(nrow(df), n_show) * 0.35 + 3), dpi = 300)

  cat("  ✅ Saved:", prefix, "(", nrow(df), "terms )\n")
  return(invisible(df))
}

# ================= 4. 主循环：每个对比 =================
for (comp_name in names(res_list)) {
  cat("\n", strrep("=", 60), "\n")
  cat("Processing:", comp_name, "\n")
  cat(strrep("=", 60), "\n")

  comp_dir <- file.path(ENR_DIR, comp_name)
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

  df <- res_list[[comp_name]]

  # 去除 NA
  df_clean <- df %>% filter(!is.na(padj), !is.na(log2FoldChange))

  # 基因集分类
  sig_up   <- df_clean %>% filter(.data[[sig_col]] == "Up")
  sig_down <- df_clean %>% filter(.data[[sig_col]] == "Down")
  sig_all  <- df_clean %>% filter(.data[[sig_col]] != "NS")

  # Ensembl → Entrez 映射 (基于全部基因)
  id_map <- ensembl_to_entrez(sub("\\..*$", "", df_clean$gene_id))

  # 为差异基因提取 Entrez
  get_entrez <- function(sub_df) {
    clean <- sub("\\..*$", "", sub_df$gene_id)
    id_map$ENTREZID[id_map$ENSEMBL %in% clean]
  }
  entrez_up   <- get_entrez(sig_up)
  entrez_down <- get_entrez(sig_down)
  entrez_all  <- get_entrez(sig_all)
  entrez_bg   <- id_map$ENTREZID  # 背景基因集

  cat("  DEGs — Up:", length(entrez_up), " Down:", length(entrez_down),
      " Total sig:", length(entrez_all), "\n")

  # ---- 4a. GO 富集 ----------------------------------------
  cat("\n--- GO Enrichment ---\n")
  go_dir <- file.path(comp_dir, "GO")

  for (ont in c("BP", "MF", "CC")) {
    for (direction in list(list(ids=entrez_up,   name="UP"),
                           list(ids=entrez_down, name="DOWN"),
                           list(ids=entrez_all,  name="ALL"))) {
      if (length(direction$ids) < 5) next

      tryCatch({
        enr <- enrichGO(gene         = direction$ids,
                        universe      = entrez_bg,
                        OrgDb         = org.Mm.eg.db,
                        ont           = ont,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.2,
                        readable      = TRUE)
        save_enrichment(enr, paste0("GO_", ont, "_", direction$name), go_dir)
      }, error = function(e) cat("  ⚠️  GO", ont, direction$name, "error:", e$message, "\n"))
    }
  }

  # ---- 4b. KEGG 富集 --------------------------------------
  cat("\n--- KEGG Enrichment ---\n")
  kegg_dir <- file.path(comp_dir, "KEGG")

  for (direction in list(list(ids=entrez_up,   name="UP"),
                         list(ids=entrez_down, name="DOWN"),
                         list(ids=entrez_all,  name="ALL"))) {
    if (length(direction$ids) < 5) next

    tryCatch({
      enr <- enrichKEGG(gene         = direction$ids,
                        organism      = "mmu",
                        universe      = entrez_bg,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.2)

      # 将 Entrez ID 转换为基因名称，方便阅读
      if (!is.null(enr) && nrow(as.data.frame(enr)) > 0) {
        enr <- setReadable(enr, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
      }
      save_enrichment(enr, paste0("KEGG_", direction$name), kegg_dir)
    }, error = function(e) cat("  ⚠️  KEGG", direction$name, "error:", e$message, "\n"))
  }

  # ---- 4c. GSEA -------------------------------------------
  cat("\n--- GSEA ---\n")
  gsea_dir <- file.path(comp_dir, "GSEA")
  dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

  # 构建 ranked gene list: sign(log2FC) * (-log10(pvalue))
  df_ranked <- df_clean %>%
    mutate(clean_ens = sub("\\..*$", "", gene_id)) %>%
    left_join(id_map, by = c("clean_ens" = "ENSEMBL")) %>%
    filter(!is.na(ENTREZID), !is.na(pvalue), pvalue > 0) %>%
    mutate(rank_stat = sign(log2FoldChange) * (-log10(pvalue))) %>%
    arrange(desc(rank_stat)) %>%
    distinct(ENTREZID, .keep_all = TRUE)

  gene_ranks <- setNames(df_ranked$rank_stat, df_ranked$ENTREZID)
  cat("  Ranked gene list:", length(gene_ranks), "genes\n")

  if (length(gene_ranks) >= 100) {
    # GSEA — GO BP
    tryCatch({
      gsea_go <- gseGO(geneList     = gene_ranks,
                       OrgDb        = org.Mm.eg.db,
                       ont          = "BP",
                       keyType      = "ENTREZID",
                       minGSSize    = 15,
                       maxGSSize    = 500,
                       pvalueCutoff = 0.05,
                       verbose      = FALSE)
      if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
        gsea_go <- setReadable(gsea_go, OrgDb = org.Mm.eg.db)
        df_g <- as.data.frame(gsea_go)
        write_csv(df_g, file.path(gsea_dir, "GSEA_GO_BP.csv"))

        # Ridge plot (top 20)
        p_ridge <- ridgeplot(gsea_go, showCategory = min(20, nrow(df_g))) +
          labs(title = paste("GSEA GO BP:", comp_name)) +
          theme(axis.text.y = element_text(size = 6))
        ggsave(file.path(gsea_dir, "GSEA_GO_BP_ridgeplot.pdf"), p_ridge,
               width = 12, height = 8, dpi = 300)

        # GSEA plot — top 3 enriched/depleted
        top_ids <- head(df_g$ID[df_g$NES > 0], 3)
        bot_ids <- head(df_g$ID[df_g$NES < 0], 3)
        for (id in c(top_ids, bot_ids)) {
          tryCatch({
            p_gsea <- gseaplot2(gsea_go, geneSetID = id,
                                title = df_g$Description[df_g$ID == id][1])
            ggsave(file.path(gsea_dir, paste0("GSEA_GO_", id, ".pdf")), p_gsea,
                   width = 8, height = 5, dpi = 300)
          }, error = function(e) NULL)
        }
        cat("  ✅ GSEA GO BP:", nrow(df_g), "terms\n")
      }
    }, error = function(e) cat("  ⚠️  GSEA GO BP error:", e$message, "\n"))

    # GSEA — KEGG
    tryCatch({
      gsea_kegg <- gseKEGG(geneList     = gene_ranks,
                           organism     = "mmu",
                           minGSSize    = 15,
                           maxGSSize    = 500,
                           pvalueCutoff = 0.05,
                           verbose      = FALSE)
      if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
        gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
        df_k <- as.data.frame(gsea_kegg)
        write_csv(df_k, file.path(gsea_dir, "GSEA_KEGG.csv"))
        p_dot_k <- dotplot(gsea_kegg, showCategory = min(20, nrow(df_k)), split = ".sign") +
          facet_grid(. ~ .sign) +
          theme(axis.text.y = element_text(size = 7)) +
          labs(title = paste("GSEA KEGG:", comp_name))
        ggsave(file.path(gsea_dir, "GSEA_KEGG_dotplot.pdf"), p_dot_k,
               width = 14, height = max(5, min(20, nrow(df_k)) * 0.4 + 3), dpi = 300)
        cat("  ✅ GSEA KEGG:", nrow(df_k), "pathways\n")
      }
    }, error = function(e) cat("  ⚠️  GSEA KEGG error:", e$message, "\n"))

    # GSEA — MSigDB Hallmark (H collection, mouse)
    tryCatch({
      msig_h <- msigdbr(species = "Mus musculus", category = "H") %>%
        select(gs_name, entrez_gene) %>%
        rename(term = gs_name, gene = entrez_gene) %>%
        mutate(gene = as.character(gene))
      gene_ranks_char <- setNames(gene_ranks, names(gene_ranks))

      gsea_h <- GSEA(geneList  = gene_ranks,
                     TERM2GENE = msig_h,
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose   = FALSE)
      if (!is.null(gsea_h) && nrow(as.data.frame(gsea_h)) > 0) {
        df_h <- as.data.frame(gsea_h)
        write_csv(df_h, file.path(gsea_dir, "GSEA_Hallmark.csv"))
        p_dot_h <- dotplot(gsea_h, showCategory = min(20, nrow(df_h)), split = ".sign") +
          facet_grid(. ~ .sign) +
          theme(axis.text.y = element_text(size = 7)) +
          labs(title = paste("GSEA Hallmark:", comp_name))
        ggsave(file.path(gsea_dir, "GSEA_Hallmark_dotplot.pdf"), p_dot_h,
               width = 14, height = max(5, min(20, nrow(df_h)) * 0.4 + 3), dpi = 300)
        cat("  ✅ GSEA Hallmark:", nrow(df_h), "gene sets\n")
      }
    }, error = function(e) cat("  ⚠️  GSEA Hallmark error:", e$message, "\n"))
  } else {
    cat("  ⚠️  Not enough ranked genes for GSEA, skip\n")
  }

  # ---- 4d. 干细胞 marker 基因分析 -------------------------
  cat("\n--- Stem Cell Marker Analysis ---\n")
  sc_dir <- file.path(comp_dir, "StemCell")
  dir.create(sc_dir, showWarnings = FALSE, recursive = TRUE)

  stemcell_markers <- list(
    Pluripotency    = c("Pou5f1","Sox2","Nanog","Klf4","Myc","Zfp42","Esrrb","Utf1","Sall4","Lin28a"),
    Neural_SC       = c("Nes","Sox1","Sox9","Pax6","Vim","Gfap","Msi1","Prom1","Blbp"),
    Hematopoietic_SC= c("Cd34","Ly6a","Slamf1","Kit","Cd48","Procr","Gata2","Tal1"),
    Mesenchymal_SC  = c("Cd44","Nt5e","Thy1","Eng","Itgb1","Pdgfra","Aldh1a1"),
    Intestinal_SC   = c("Lgr5","Axin2","Olfm4","Ascl2","Lrig1","Troy"),
    General_Stemness= c("Aldh1a1","Aldh1a3","Epcam","Cd24a","Sox4","Id1","Id3","Notch1","Notch2")
  )

  # 在 counts_raw 的 gene_name 列中查找 marker
  all_markers <- unique(unlist(stemcell_markers))
  marker_hits <- counts_raw %>%
    filter(gene_name %in% all_markers) %>%
    select(gene_id, gene_name)

  # 合并 DE 结果
  marker_de <- df %>%
    filter(gene_name %in% all_markers) %>%
    select(gene_id, gene_name, log2FoldChange, padj, .data[[sig_col]]) %>%
    arrange(padj)

  # 添加 marker 类别
  marker_category <- data.frame(
    gene_name = unlist(stemcell_markers),
    Category  = rep(names(stemcell_markers), sapply(stemcell_markers, length)),
    stringsAsFactors = FALSE
  )
  marker_de <- left_join(marker_de, marker_category, by = "gene_name")

  write_csv(marker_de, file.path(sc_dir, "StemCell_Markers_DE_Results.csv"))
  cat("  ✅ Stem cell markers found in data:", nrow(marker_de), "/", length(all_markers), "\n")

  # 显著变化的干细胞 marker 汇总
  sig_markers <- marker_de %>% filter(.data[[sig_col]] != "NS")
  if (nrow(sig_markers) > 0) {
    cat("  ✅ Significant stem cell markers:", nrow(sig_markers), "\n")
    print(sig_markers %>% select(gene_name, Category, log2FoldChange, padj, .data[[sig_col]]))

    # 条形图：所有检测到的 marker 的 log2FC
    p_bar <- marker_de %>%
      filter(!is.na(log2FoldChange)) %>%
      arrange(Category, log2FoldChange) %>%
      mutate(gene_name = factor(gene_name, levels = unique(gene_name)),
             Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
      ggplot(aes(x = gene_name, y = log2FoldChange, fill = Direction)) +
      geom_col() +
      geom_hline(yintercept = c(-0.263, 0.263), linetype = "dashed", color = "grey50") +
      facet_wrap(~Category, scales = "free_x", nrow = 2) +
      scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            strip.text   = element_text(size = 8)) +
      labs(title = paste("Stem Cell Marker Genes:", comp_name),
           x = NULL, y = "log2 Fold Change")

    tryCatch(
      ggsave(file.path(sc_dir, "StemCell_Markers_log2FC.pdf"), p_bar,
             width = 14, height = 7, dpi = 300),
      error = function(e) cat("  ⚠️  StemCell bar plot save failed:", e$message, "\n")
    )
  } else {
    cat("  ℹ️  No significant stem cell markers detected\n")
  }

  # MSigDB 干细胞相关 gene set ORA (C8 cell type signatures)
  tryCatch({
    msig_c8 <- msigdbr(species = "Mus musculus", category = "C8") %>%
      filter(grepl("STEM|PROGENITOR|PLURIP", gs_name, ignore.case = TRUE)) %>%
      select(gs_name, entrez_gene) %>%
      rename(term = gs_name, gene = entrez_gene)

    if (nrow(msig_c8) > 0 && length(entrez_all) >= 5) {
      enr_sc <- enricher(gene     = as.character(entrez_all),
                         universe = as.character(entrez_bg),
                         TERM2GENE = msig_c8 %>% mutate(gene = as.character(gene)),
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.2)
      if (!is.null(enr_sc) && nrow(as.data.frame(enr_sc)) > 0) {
        df_sc <- as.data.frame(enr_sc)
        write_csv(df_sc, file.path(sc_dir, "MSigDB_StemCell_Signatures.csv"))
        p_sc <- dotplot(enr_sc, showCategory = min(15, nrow(df_sc)),
                        title = paste("MSigDB Stem Cell Signatures:", comp_name)) +
          theme(axis.text.y = element_text(size = 7))
        ggsave(file.path(sc_dir, "MSigDB_StemCell_dotplot.pdf"), p_sc,
               width = 10, height = 6, dpi = 300)
        cat("  ✅ MSigDB stem cell ORA:", nrow(df_sc), "gene sets\n")
      }
    }
  }, error = function(e) cat("  ⚠️  MSigDB stem cell ORA error:", e$message, "\n"))

  cat("✅ Completed enrichment for:", comp_name, "\n")
}

# ================= 5. 跨对比汇总 =================
cat("\n", strrep("=", 60), "\n")
cat("Cross-comparison summary\n")
cat(strrep("=", 60), "\n")

# 汇总所有对比的显著干细胞 marker
summary_list <- list()
for (comp_name in names(res_list)) {
  sc_csv <- file.path(ENR_DIR, comp_name, "StemCell", "StemCell_Markers_DE_Results.csv")
  if (file.exists(sc_csv)) {
    df_tmp <- read_csv(sc_csv, show_col_types = FALSE) %>%
      mutate(Comparison = comp_name)
    summary_list[[comp_name]] <- df_tmp
  }
}
if (length(summary_list) > 0) {
  summary_all <- bind_rows(summary_list)
  write_csv(summary_all, file.path(ENR_DIR, "StemCell_AllComparisons_Summary.csv"))
  cat("✅ Cross-comparison stem cell summary saved\n")
}

cat("\n🎉 All enrichment analyses complete. Results:", ENR_DIR, "\n")
