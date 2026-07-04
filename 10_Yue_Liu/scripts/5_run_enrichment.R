#!/usr/bin/env Rscript
# 运行环境: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript
# 运行方法: cd /home/gao/projects_2026H2/10_Yue_Liu/scripts && Rscript 5_run_enrichment.R
#
# 分析内容:
#   1. GO enrichment analysis  (BP / MF / CC)
#   2. KEGG pathway analysis
#   3. GSEA (GO BP + KEGG + MSigDB Hallmark)
#
# 前置条件: 先运行 4_run_DE_PCA.R 生成 res_list.rds / meta.rds / counts_raw.rds
# 无需 StemCell 分析（客户确认：普通细胞类型对比研究）

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# ================= 1. 路径设置 =================
setwd("/home/gao/projects_2026H2/10_Yue_Liu/scripts/")

data_dirs <- sort(list.dirs("..", full.names = TRUE, recursive = FALSE), decreasing = TRUE)
data_dirs <- data_dirs[grepl("/Data_Analysis_[0-9]{8}$", data_dirs)]
if (length(data_dirs) == 0) stop("No Data_Analysis_YYYYMMDD folder found. Run 4_run_DE_PCA.R first.")
DATA_DIR <- data_dirs[1]
cat("Using:", DATA_DIR, "\n")

DE_DIR  <- file.path(DATA_DIR, "DE_PCA_Results")
ENR_DIR <- file.path(DATA_DIR, "Enrichment")
dir.create(ENR_DIR, showWarnings = FALSE, recursive = TRUE)

# ================= 2. 加载 DE 结果 =================
res_list   <- readRDS(file.path(DE_DIR, "res_list.rds"))
meta       <- readRDS(file.path(DE_DIR, "meta.rds"))
counts_raw <- readRDS(file.path(DE_DIR, "counts_raw.rds"))

sig_col <- "sig (padj<=0.05 & |log2FC|>=0.263)"
cat("DE results loaded. Comparisons:", paste(names(res_list), collapse = ", "), "\n")

# ================= 3. 辅助函数 =================
ensembl_to_entrez <- function(ensembl_ids) {
  clean_ids <- sub("\\..*$", "", ensembl_ids)
  bitr(clean_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}

save_enrichment <- function(enr_res, prefix, out_dir, n_show = 20) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(enr_res) || nrow(as.data.frame(enr_res)) == 0) {
    cat("  No significant results for:", prefix, "\n")
    return(invisible(NULL))
  }
  df <- as.data.frame(enr_res)
  write_csv(df, file.path(out_dir, paste0(prefix, ".csv")))
  p <- dotplot(enr_res, showCategory = min(n_show, nrow(df)), title = prefix) +
    theme(axis.text.y = element_text(size = 7))
  ggsave(file.path(out_dir, paste0(prefix, "_dotplot.pdf")), p,
         width = 10, height = max(5, min(nrow(df), n_show) * 0.35 + 3), dpi = 300)
  cat("  Saved:", prefix, "(", nrow(df), "terms )\n")
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
  df_clean <- df %>% filter(!is.na(padj), !is.na(log2FoldChange))

  sig_up   <- df_clean %>% filter(.data[[sig_col]] == "Up")
  sig_down <- df_clean %>% filter(.data[[sig_col]] == "Down")
  sig_all  <- df_clean %>% filter(.data[[sig_col]] != "NS")

  id_map <- ensembl_to_entrez(sub("\\..*$", "", df_clean$gene_id))

  get_entrez <- function(sub_df) {
    clean <- sub("\\..*$", "", sub_df$gene_id)
    id_map$ENTREZID[id_map$ENSEMBL %in% clean]
  }
  entrez_up   <- get_entrez(sig_up)
  entrez_down <- get_entrez(sig_down)
  entrez_all  <- get_entrez(sig_all)
  entrez_bg   <- id_map$ENTREZID

  cat("  DEGs — Up:", length(entrez_up), " Down:", length(entrez_down),
      " Total sig:", length(entrez_all), "\n")

  # ---- 4a. GO 富集 ----------------------------------------
  cat("\n--- GO Enrichment ---\n")
  go_dir <- file.path(comp_dir, "GO")

  for (ont in c("BP", "MF", "CC")) {
    for (direction in list(
      list(ids = entrez_up,   name = "UP"),
      list(ids = entrez_down, name = "DOWN"),
      list(ids = entrez_all,  name = "ALL")
    )) {
      if (length(direction$ids) < 5) next
      tryCatch({
        enr <- enrichGO(gene = direction$ids, universe = entrez_bg,
                         OrgDb = org.Hs.eg.db, ont = ont,
                         pAdjustMethod = "BH", pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2, readable = TRUE)
        save_enrichment(enr, paste0("GO_", ont, "_", direction$name), go_dir)
      }, error = function(e) cat("  GO", ont, direction$name, "error:", e$message, "\n"))
    }
  }

  # ---- 4b. KEGG 富集 --------------------------------------
  cat("\n--- KEGG Enrichment ---\n")
  kegg_dir <- file.path(comp_dir, "KEGG")

  for (direction in list(
    list(ids = entrez_up,   name = "UP"),
    list(ids = entrez_down, name = "DOWN"),
    list(ids = entrez_all,  name = "ALL")
  )) {
    if (length(direction$ids) < 5) next
    tryCatch({
      enr <- enrichKEGG(gene = direction$ids, organism = "hsa",
                         universe = entrez_bg, pAdjustMethod = "BH",
                         pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      if (!is.null(enr) && nrow(as.data.frame(enr)) > 0) {
        enr <- setReadable(enr, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      }
      save_enrichment(enr, paste0("KEGG_", direction$name), kegg_dir)
    }, error = function(e) cat("  KEGG", direction$name, "error:", e$message, "\n"))
  }

  # ---- 4c. GSEA -------------------------------------------
  cat("\n--- GSEA ---\n")
  gsea_dir <- file.path(comp_dir, "GSEA")
  dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

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
      gsea_go <- gseGO(geneList = gene_ranks, OrgDb = org.Hs.eg.db, ont = "BP",
                        keyType = "ENTREZID", minGSSize = 15, maxGSSize = 500,
                        pvalueCutoff = 0.05, verbose = FALSE)
      if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
        gsea_go <- setReadable(gsea_go, OrgDb = org.Hs.eg.db)
        df_g <- as.data.frame(gsea_go)
        write_csv(df_g, file.path(gsea_dir, "GSEA_GO_BP.csv"))

        p_ridge <- ridgeplot(gsea_go, showCategory = min(20, nrow(df_g))) +
          labs(title = paste("GSEA GO BP:", comp_name)) +
          theme(axis.text.y = element_text(size = 6))
        ggsave(file.path(gsea_dir, "GSEA_GO_BP_ridgeplot.pdf"), p_ridge,
               width = 12, height = 8, dpi = 300)

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
        cat("  GSEA GO BP:", nrow(df_g), "terms\n")
      }
    }, error = function(e) cat("  GSEA GO BP error:", e$message, "\n"))

    # GSEA — KEGG
    tryCatch({
      gsea_kegg <- gseKEGG(geneList = gene_ranks, organism = "hsa",
                            minGSSize = 15, maxGSSize = 500,
                            pvalueCutoff = 0.05, verbose = FALSE)
      if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
        gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        df_k <- as.data.frame(gsea_kegg)
        write_csv(df_k, file.path(gsea_dir, "GSEA_KEGG.csv"))
        p_dot_k <- dotplot(gsea_kegg, showCategory = min(20, nrow(df_k)), split = ".sign") +
          facet_grid(. ~ .sign) + theme(axis.text.y = element_text(size = 7)) +
          labs(title = paste("GSEA KEGG:", comp_name))
        ggsave(file.path(gsea_dir, "GSEA_KEGG_dotplot.pdf"), p_dot_k,
               width = 14, height = max(5, min(20, nrow(df_k)) * 0.4 + 3), dpi = 300)
        cat("  GSEA KEGG:", nrow(df_k), "pathways\n")
      }
    }, error = function(e) cat("  GSEA KEGG error:", e$message, "\n"))

    # GSEA — MSigDB Hallmark (H collection, human)
    tryCatch({
      msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
        select(gs_name, entrez_gene) %>%
        rename(term = gs_name, gene = entrez_gene) %>%
        mutate(gene = as.character(gene))

      gsea_h <- GSEA(geneList = gene_ranks, TERM2GENE = msig_h,
                     minGSSize = 10, maxGSSize = 500,
                     pvalueCutoff = 0.05, verbose = FALSE)
      if (!is.null(gsea_h) && nrow(as.data.frame(gsea_h)) > 0) {
        df_h <- as.data.frame(gsea_h)
        write_csv(df_h, file.path(gsea_dir, "GSEA_Hallmark.csv"))
        p_dot_h <- dotplot(gsea_h, showCategory = min(20, nrow(df_h)), split = ".sign") +
          facet_grid(. ~ .sign) + theme(axis.text.y = element_text(size = 7)) +
          labs(title = paste("GSEA Hallmark:", comp_name))
        ggsave(file.path(gsea_dir, "GSEA_Hallmark_dotplot.pdf"), p_dot_h,
               width = 14, height = max(5, min(20, nrow(df_h)) * 0.4 + 3), dpi = 300)
        cat("  GSEA Hallmark:", nrow(df_h), "gene sets\n")
      }
    }, error = function(e) cat("  GSEA Hallmark error:", e$message, "\n"))
  } else {
    cat("  Not enough ranked genes for GSEA, skip\n")
  }
}

cat("\nEnrichment analyses complete. Results:", ENR_DIR, "\n")
