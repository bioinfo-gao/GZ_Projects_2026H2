#!/usr/bin/env Rscript
# 运行环境: conda activate DE_R45
# 运行方法: cd /home/gao/projects_2026H2/9_Lijian_Wu_Mouse/scripts && Rscript 5_run_enrichment.R
#
# 分析内容-Standard Pathway Enrichment:
#   1. GO enrichment analysis  (BP / MF / CC)
#   2. KEGG pathway analysis
#   3. GSEA (Gene Set Enrichment Analysis)  — GO BP + KEGG + MSigDB Hallmark
#
# 前置条件: 先运行 4_run_DE_PCA.R 生成 res_list.rds / meta.rds / counts_raw.rds

library(clusterProfiler) # clusterProfiler
library(enrichplot)
library(org.Mm.eg.db)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# ================= 1. 路径设置 =================
setwd("/home/gao/projects_2026H2/9_Lijian_Wu_Mouse/scripts/")

# 自动找最新的 Data_Analysis_YYYYMMDD 文件夹
data_dirs <- sort(
  list.dirs("..", full.names = TRUE, recursive = FALSE),
  decreasing = TRUE
)
data_dirs <- data_dirs[grepl("/Data_Analysis_[0-9]{8}$", data_dirs)]
if (length(data_dirs) == 0) {
  stop("No Data_Analysis_YYYYMMDD folder found. Run 4_run_DE_PCA.R first.")
}
DATA_DIR <- data_dirs[1]
cat("Using:", DATA_DIR, "\n")

DE_DIR <- "rds_cache"  # internal handoff objects (not the client-facing DE_PCA_Results/ folder)
ENR_DIR <- file.path(DATA_DIR, "Enrichment")
dir.create(ENR_DIR, showWarnings = FALSE, recursive = TRUE)

# ================= 2. 加载 DE 结果 =================
res_list <- readRDS(file.path(DE_DIR, "res_list.rds"))
meta <- readRDS(file.path(DE_DIR, "meta.rds"))
counts_raw <- readRDS(file.path(DE_DIR, "counts_raw.rds"))
summary(res_list)
summary(meta)
summary(counts_raw)

sig_col <- "sig (padj<=0.05 & |log2FC|>=0.263)"

cat(
  "✅ DE results loaded. Comparisons:",
  paste(names(res_list), collapse = ", "),
  "\n"
)

# ================= 3. 辅助函数 =================

# Ensembl → Entrez ID 转换（鼠: org.Mm.eg.db）
# 返回双列 data.frame: ENSEMBL | ENTREZID（一对多时会扩展行数）
ensembl_to_entrez <- function(ensembl_ids) {
  # nf-core salmon 输出的 gene_id 含版本号后缀，需先去掉才能匹配 OrgDb e.g. ENSMUSG00000051951.6 → ENSMUSG00000051951
  clean_ids <- sub("\\..*$", "", ensembl_ids)
  # bitr = Biological Id TRanslator；不能映射的 ID 会被静默丢弃
  map <- bitr(
    clean_ids,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  return(map)
}

# 保存富集结果 CSV + dotplot
# enr_res : enrichResult 或 gseaResult 对象
# prefix  : 输出文件名前缀，例如 "GO_BP_UP"
# out_dir : 输出目录（不存在会自动创建）
# n_show  : dotplot 最多显示的 term 数量
save_enrichment <- function(enr_res, prefix, out_dir, n_show = 20) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(enr_res) || nrow(as.data.frame(enr_res)) == 0) {
    cat("  ⚠️  No significant results for:", prefix, "\n")
    return(invisible(NULL)) # invisible: suppresses auto-print when caller assigns with <-   (让caller 用 <- 接返回值时不自动打印)
  }

  df <- as.data.frame(enr_res)
  write_csv(df, file.path(out_dir, paste0(prefix, ".csv")))

  p <- dotplot(enr_res, showCategory = min(n_show, nrow(df)), title = prefix) +
    theme(axis.text.y = element_text(size = 7))

  # 图高随 term 数量动态调整，避免文字堆叠或空白过多
  ggsave(
    file.path(out_dir, paste0(prefix, "_dotplot.pdf")),
    p,
    width = 10,
    height = max(5, min(nrow(df), n_show) * 0.35 + 3),
    dpi = 300
  )

  cat("  ✅ Saved:", prefix, "(", nrow(df), "terms )\n")
  return(invisible(df)) # returns df for optional downstream use without auto-printing  (返回 data.frame 供调用方按需使用，但不触发自动打印)
}

# ================= 4. 主循环：每个对比 =================
for (comp_name in names(res_list)) {
  cat("\n", strrep("=", 60), "\n")
  cat("Processing:", comp_name, "\n")
  cat(strrep("=", 60), "\n")

  comp_dir <- file.path(ENR_DIR, comp_name)
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

  df <- res_list[[comp_name]]

  # padj / log2FC 含 NA 的行（低表达或过滤掉的基因）不参与富集
  df_clean <- df %>% filter(!is.na(padj), !is.na(log2FoldChange))

  # ORA 支持方向性分析：UP/DOWN 单独跑可识别方向特异通路；ALL 做总体富集
  sig_up   <- df_clean %>% filter(.data[[sig_col]] == "Up")
  sig_down <- df_clean %>% filter(.data[[sig_col]] == "Down")
  sig_all  <- df_clean %>% filter(.data[[sig_col]] != "NS")

  # 背景集用 df_clean 全部基因（而非仅 DEG），这是 ORA 正确的 universe 设置
  # bitr 一对多时会产生重复 ENSEMBL，entrez_bg 保留全部以覆盖最大背景
  id_map <- ensembl_to_entrez(sub("\\..*$", "", df_clean$gene_id))

  # 闭包：捕获上方 id_map，从任意子 data.frame 中提取对应 Entrez ID
  get_entrez <- function(sub_df) {
    clean <- sub("\\..*$", "", sub_df$gene_id)
    id_map$ENTREZID[id_map$ENSEMBL %in% clean]
  }
  entrez_up   <- get_entrez(sig_up)
  entrez_down <- get_entrez(sig_down)
  entrez_all  <- get_entrez(sig_all)
  entrez_bg   <- id_map$ENTREZID # ORA universe：能映射到 Entrez 的所有检测基因

  cat(
    "  DEGs — Up:",
    length(entrez_up),
    " Down:",
    length(entrez_down),
    " Total sig:",
    length(entrez_all),
    "\n"
  )

  # ---- 4a. GO 富集 ----------------------------------------
  # 对 BP/MF/CC 三个本体分别做 ORA，每个本体再按 UP/DOWN/ALL 分方向跑
  cat("\n--- GO Enrichment ---\n")
  go_dir <- file.path(comp_dir, "GO")

  for (ont in c("BP", "MF", "CC")) {
    for (direction in list(
      list(ids = entrez_up,   name = "UP"),
      list(ids = entrez_down, name = "DOWN"),
      list(ids = entrez_all,  name = "ALL")
    )) {
      if (length(direction$ids) < 5) {
        next
      } # 基因太少时 ORA 统计不可靠，直接跳过

      tryCatch(
        {
          enr <- enrichGO(
            gene = direction$ids,
            universe = entrez_bg,
            OrgDb = org.Mm.eg.db,
            ont = ont,
            pAdjustMethod = "BH", # Benjamini-Hochberg FDR 校正
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2, # q-value 作为二级过滤（比 padj 宽松）
            readable = TRUE
          ) # 自动将结果中 Entrez ID 替换为基因 symbol
          save_enrichment(enr, paste0("GO_", ont, "_", direction$name), go_dir)
        },
        error = function(e) {
          cat("  ⚠️  GO", ont, direction$name, "error:", e$message, "\n")
        }
      )
    }
  }

  # ---- 4b. KEGG 富集 --------------------------------------
  # enrichKEGG 不支持 readable=TRUE 参数，需要在得到结果后手动调用 setReadable()
  cat("\n--- KEGG Enrichment ---\n")
  kegg_dir <- file.path(comp_dir, "KEGG")

  for (direction in list(
    list(ids = entrez_up,   name = "UP"),
    list(ids = entrez_down, name = "DOWN"),
    list(ids = entrez_all,  name = "ALL")
  )) {
    if (length(direction$ids) < 5) {
      next
    }

    tryCatch(
      {
        enr <- enrichKEGG(
          gene = direction$ids,
          organism      = "mmu", # Mus musculus KEGG 物种代码
          universe      = entrez_bg,
          pAdjustMethod = "BH",
          pvalueCutoff  = 0.05,
          qvalueCutoff  = 0.2
        )

        # enrichKEGG 结果中 geneID 列默认是 Entrez 数字；转为 symbol 方便报告阅读
        if (!is.null(enr) && nrow(as.data.frame(enr)) > 0) {
          enr <- setReadable(enr, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
        }
        save_enrichment(enr, paste0("KEGG_", direction$name), kegg_dir)
      },
      error = function(e) {
        cat("  ⚠️  KEGG", direction$name, "error:", e$message, "\n")
      }
    )
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
    tryCatch(
      {
        gsea_go <- gseGO(
          geneList  = gene_ranks,
          OrgDb     = org.Mm.eg.db,
          ont       = "BP",
          keyType   = "ENTREZID",
          minGSSize = 15,
          maxGSSize = 500,
          pvalueCutoff = 0.05,
          verbose   = FALSE
        )
        if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
          gsea_go <- setReadable(gsea_go, OrgDb = org.Mm.eg.db)
          df_g <- as.data.frame(gsea_go)
          write_csv(df_g, file.path(gsea_dir, "GSEA_GO_BP.csv"))

          # Ridge plot (top 20)
          p_ridge <- ridgeplot(gsea_go, showCategory = min(20, nrow(df_g))) +
            labs(title = paste("GSEA GO BP:", comp_name)) +
            theme(axis.text.y = element_text(size = 6))
          ggsave(
            file.path(gsea_dir, "GSEA_GO_BP_ridgeplot.pdf"),
            p_ridge,
            width  = 12,
            height = 8,
            dpi    = 300
          )

          # GSEA plot — top 3 enriched/depleted
          top_ids <- head(df_g$ID[df_g$NES > 0], 3)
          bot_ids <- head(df_g$ID[df_g$NES < 0], 3)
          for (id in c(top_ids, bot_ids)) {
            tryCatch(
              {
                p_gsea <- gseaplot2(
                  gsea_go,
                  geneSetID = id,
                  title = df_g$Description[df_g$ID == id][1]
                )
                ggsave(
                  file.path(gsea_dir, paste0("GSEA_GO_", id, ".pdf")),
                  p_gsea,
                  width = 8,
                  height = 5,
                  dpi = 300
                )
              },
              error = function(e) NULL
            )
          }
          cat("  ✅ GSEA GO BP:", nrow(df_g), "terms\n")
        }
      },
      error = function(e) cat("  ⚠️  GSEA GO BP error:", e$message, "\n")
    )

    # GSEA — KEGG
    tryCatch(
      {
        gsea_kegg <- gseKEGG(
          geneList = gene_ranks,
          organism = "mmu",
          minGSSize = 15,
          maxGSSize = 500,
          pvalueCutoff = 0.05,
          verbose = FALSE
        )
        if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
          gsea_kegg <- setReadable(
            gsea_kegg,
            OrgDb = org.Mm.eg.db,
            keyType = "ENTREZID"
          )
          df_k <- as.data.frame(gsea_kegg)
          write_csv(df_k, file.path(gsea_dir, "GSEA_KEGG.csv"))
          p_dot_k <- dotplot(
            gsea_kegg,
            showCategory = min(20, nrow(df_k)),
            split = ".sign"
          ) +
            facet_grid(. ~ .sign) +
            theme(axis.text.y = element_text(size = 7)) +
            labs(title = paste("GSEA KEGG:", comp_name))
          ggsave(
            file.path(gsea_dir, "GSEA_KEGG_dotplot.pdf"),
            p_dot_k,
            width = 14,
            height = max(5, min(20, nrow(df_k)) * 0.4 + 3),
            dpi = 300
          )
          cat("  ✅ GSEA KEGG:", nrow(df_k), "pathways\n")
        }
      },
      error = function(e) cat("  ⚠️  GSEA KEGG error:", e$message, "\n")
    )

    # GSEA — MSigDB Hallmark (H collection, mouse)
    tryCatch(
      {
        msig_h <- msigdbr(species = "Mus musculus", category = "H") %>%
          select(gs_name, entrez_gene) %>%
          rename(term = gs_name, gene = entrez_gene) %>%
          mutate(gene = as.character(gene))
        gene_ranks_char <- setNames(gene_ranks, names(gene_ranks))

        gsea_h <- GSEA(
          geneList = gene_ranks,
          TERM2GENE = msig_h,
          minGSSize = 10,
          maxGSSize = 500,
          pvalueCutoff = 0.05,
          verbose = FALSE
        )
        if (!is.null(gsea_h) && nrow(as.data.frame(gsea_h)) > 0) {
          df_h <- as.data.frame(gsea_h)
          write_csv(df_h, file.path(gsea_dir, "GSEA_Hallmark.csv"))
          p_dot_h <- dotplot(
            gsea_h,
            showCategory = min(20, nrow(df_h)),
            split = ".sign"
          ) +
            facet_grid(. ~ .sign) +
            theme(axis.text.y = element_text(size = 7)) +
            labs(title = paste("GSEA Hallmark:", comp_name))
          ggsave(
            file.path(gsea_dir, "GSEA_Hallmark_dotplot.pdf"),
            p_dot_h,
            width = 14,
            height = max(5, min(20, nrow(df_h)) * 0.4 + 3),
            dpi = 300
          )
          cat("  ✅ GSEA Hallmark:", nrow(df_h), "gene sets\n")
        }
      },
      error = function(e) cat("  ⚠️  GSEA Hallmark error:", e$message, "\n")
    )
  } else {
    cat("  ⚠️  Not enough ranked genes for GSEA, skip\n")
  }
}

cat("\nEnrichment analyses complete. Results:", ENR_DIR, "\n")
