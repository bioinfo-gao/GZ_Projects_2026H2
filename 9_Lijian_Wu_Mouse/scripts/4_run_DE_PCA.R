#!/usr/bin/env Rscript
# 运行环境: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript
# 运行方法: cd /home/gao/projects_2026H2/9_Lijian_Wu_Mouse/scripts && Rscript 4_run_DE_PCA.R
#
# 分析设计: Control = 对照组 (client note: "control vs other groups")
#   TNFa vs Control, Tumor vs Control, Tumor_TNF vs Control

library(DESeq2)
library(ashr)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(readr)
library(readxl)
library(tidyr)
library(ggrepel)

# ================= 1. 路径设置 =================
setwd("/home/gao/projects_2026H2/9_Lijian_Wu_Mouse/scripts/")

TODAY      <- format(Sys.Date(), "%Y%m%d")
META_FILE  <- "../Sample_Sheet_Lijian_Wu.xlsx"
COUNT_FILE <- "../output_results/star_salmon/salmon.merged.gene_counts.tsv"
TPM_FILE   <- "../output_results/star_salmon/salmon.merged.gene_tpm.tsv"
DATA_DIR   <- paste0("../Data_Analysis_", TODAY)
OUT_DIR    <- file.path(DATA_DIR, "DE_PCA_Results")
READS_DIR  <- file.path(DATA_DIR, "Reads")

dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(READS_DIR, showWarnings = FALSE, recursive = TRUE)

# ================= 2. 拷贝原始计数文件 =================
for (f in list(list(src=COUNT_FILE, dst="All_sample_gene_counts.tsv", req=TRUE),
               list(src=TPM_FILE,   dst="All_sample_gene_tpm.tsv",    req=FALSE))) {
  if (file.exists(f$src)) {
    file.copy(f$src, file.path(READS_DIR, f$dst), overwrite=TRUE)
    cat("Copied:", f$dst, "\n")
  } else {
    if (f$req) stop("Required file missing: ", f$src)
    cat("Optional file missing, skipped:", f$src, "\n")
  }
}

# ================= 2.5 拷贝基因注释 + QC 文件 =================
DATA_ANALYSIS_DIR <- DATA_DIR

annot_files <- sort(
  list.files("/Work_bio/references/Genes",
             pattern = "^mouse_Gene_annotation_.*\\.xlsx$", full.names = TRUE),
  decreasing = TRUE
)
if (length(annot_files) > 0) {
  ANNOT_SRC  <- annot_files[1]
  ANNOT_DEST <- file.path(DATA_ANALYSIS_DIR, basename(ANNOT_SRC))
  file.copy(ANNOT_SRC, ANNOT_DEST, overwrite = TRUE)
  cat("Gene annotation copied:", basename(ANNOT_SRC), "\n")
} else {
  cat("No mouse_Gene_annotation_*.xlsx found, skipping\n")
}

QC_DEST_BASE <- file.path(DATA_ANALYSIS_DIR, "QC")
dir.create(QC_DEST_BASE, showWarnings = FALSE, recursive = TRUE)
if (dir.exists("../output_results/multiqc")) {
  file.copy("../output_results/multiqc", QC_DEST_BASE, recursive = TRUE, overwrite = TRUE)
  cat("QC copied: multiqc ->", QC_DEST_BASE, "\n")
}

# ================= 3. 元数据 =================
# sample_id 直接取 "Name in File" 列 (e.g. cDC1_TNFa_1)，与 nf-core samplesheet 'sample' 列一致
meta_raw <- read_excel(META_FILE)
meta_raw <- meta_raw[!is.na(meta_raw$Group), ]

meta <- meta_raw %>%
  transmute(
    sample_id = `Name in File`,
    Group     = factor(as.character(Group), levels = c("Control", "TNFa", "Tumor", "Tumor_TNF"))  # Control = 对照组
  )

cat("Metadata loaded:", nrow(meta), "samples\n")
print(meta)

# ================= 4. 表达矩阵预处理 =================
counts_raw <- read_tsv(COUNT_FILE, col_types = cols())

valid_samples <- colnames(counts_raw)[3:ncol(counts_raw)]
meta <- meta[meta$sample_id %in% valid_samples, ]
counts_mat <- as.matrix(counts_raw[, meta$sample_id])
rownames(counts_mat) <- counts_raw$gene_id
counts_mat <- round(counts_mat)

# ================= 5. 基因过滤 (鼠适配 regex) =================
gene_info <- counts_raw[, c("gene_id", "gene_name")]

ribo_pattern <- "^Rpl|^Rps|^Mrpl|^Mrps|^Rplp|^Rpsa"
non_ribosomal <- !grepl(ribo_pattern, gene_info$gene_name, ignore.case = TRUE)

noncoding_pattern <- paste0(
  "^mt-",
  "|^Snord|^Snora|^Rnu",
  "|^Malat1|^Neat1|^Xist",
  "|^Mir[0-9]|^let-",
  "|^Linc|^Gm[0-9]",
  "|pseudogene|antisense"
)
protein_coding_like <- !grepl(noncoding_pattern, gene_info$gene_name, ignore.case = TRUE)

regex_filter <- non_ribosomal & protein_coding_like

n_samples <- ncol(counts_mat)
low_count_filter <- rowSums(counts_mat >= 10) >= (n_samples - 2)

final_filter <- regex_filter & low_count_filter
counts_mat_filtered <- counts_mat[final_filter, ]

cat("Original genes:", nrow(counts_mat), "\n")
cat("After regex filter:", sum(regex_filter), "\n")
cat("After low-count filter:", sum(final_filter), "\n")
cat("Final genes for DESeq2:", nrow(counts_mat_filtered), "\n")

# ================= 6. DESeq2 建模 =================
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat_filtered,
  colData   = meta,
  design    = ~Group
)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# ================= 7. PCA =================
pca_data   <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = Group, label = name)) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_text_repel(size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 12) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 10),
        axis.text       = element_text(size = 10),
        axis.title      = element_text(size = 11))

ggsave(file.path(OUT_DIR, "PCA.pdf"), p_pca, width = 8, height = 6, dpi = 300)
cat("PCA saved\n")

# ================= 8. 差异表达分析: TNFa/Tumor/Tumor_TNF vs Control =================
sig_col  <- "sig (padj<=0.05 & |log2FC|>=0.263)"
contrasts <- list(
  c("Group", "TNFa",      "Control"),
  c("Group", "Tumor",      "Control"),
  c("Group", "Tumor_TNF",  "Control")
)

res_list <- list()

for (comp in contrasts) {
  grp_trt  <- comp[2]
  grp_ctrl <- comp[3]
  comp_name <- paste(grp_trt, "vs", grp_ctrl, sep = "_")
  cat("\nAnalysing:", comp_name, "\n")

  res <- lfcShrink(dds, contrast = comp, type = "ashr")

  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)

  samples_trt  <- meta$sample_id[meta$Group == grp_trt]
  samples_ctrl <- meta$sample_id[meta$Group == grp_ctrl]
  raw_sub <- counts_raw[, c("gene_id", "gene_name", samples_trt, samples_ctrl), drop=FALSE]

  res_df <- left_join(res_df, raw_sub, by = "gene_id")

  final_cols <- c("gene_id", "gene_name", samples_trt, samples_ctrl,
                  "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
  res_df <- res_df[, final_cols] %>% arrange(padj)

  res_df[[sig_col]] <- case_when(
    res_df$padj <= 0.05 & res_df$log2FoldChange >=  0.263 ~ "Up",
    res_df$padj <= 0.05 & res_df$log2FoldChange <= -0.263 ~ "Down",
    TRUE ~ "NS"
  )

  write_csv(res_df, file.path(OUT_DIR, paste0("DEG_", comp_name, ".csv")))
  res_list[[comp_name]] <- res_df

  # 火山图
  res_df$negLogPadj <- -log10(res_df$padj)
  res_df$negLogPadj[!is.finite(res_df$negLogPadj)] <- NA

  top_labels <- res_df %>%
    filter(.data[[sig_col]] != "NS", !is.na(negLogPadj)) %>%
    arrange(padj) %>% head(10) %>%
    mutate(label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name))

  p_vol <- ggplot(res_df, aes(x=log2FoldChange, y=negLogPadj, color=.data[[sig_col]])) +
    geom_point(alpha=0.7, size=0.5) +
    scale_color_manual(values = c("Up"="#E41A1C","Down"="#377EB8","NS"="grey80"),
                       labels = c("Up"="Upregulated","Down"="Downregulated","NS"="Not Significant")) +
    theme_bw(base_size = 10) +
    labs(title = paste(grp_trt, "vs", grp_ctrl, "(log2FC>0 =", grp_trt, "upregulated)"),
         x = "log2 Fold Change", y = "-log10(adj. P-value)") +
    theme(plot.title=element_text(hjust=0.5, size=9), legend.position="bottom",
          legend.text=element_text(size=7), axis.text=element_text(size=8)) +
    geom_text_repel(data=top_labels, aes(label=label), size=2,
                    box.padding=0.3, max.overlaps=20, color="black", fontface="plain")

  ggsave(file.path(OUT_DIR, paste0("Volcano_", comp_name, ".png")),
         p_vol, width=8, height=6, dpi=300)

  cat(comp_name, "significant DEGs:", sum(res_df[[sig_col]] != "NS"), "\n")
}

# ================= 9. 热图 =================
for (comp_name in names(res_list)) {
  cat("\nHeatmap:", comp_name, "\n")

  top50 <- res_list[[comp_name]] %>%
    filter(.data[[sig_col]] != "NS") %>% arrange(padj) %>% head(50)

  if (nrow(top50) == 0) { cat("No significant DEGs, skip heatmap\n"); next }

  parts    <- strsplit(comp_name, "_vs_")[[1]]
  grp_trt  <- parts[1]; grp_ctrl <- parts[2]
  smp_comp <- c(meta$sample_id[meta$Group == grp_trt],
                meta$sample_id[meta$Group == grp_ctrl])

  mat <- assay(vsd)[top50$gene_id, smp_comp, drop=FALSE]
  mat <- t(scale(t(mat)))

  gene_labels <- ifelse(is.na(top50$gene_name) | top50$gene_name == "",
                        top50$gene_id, top50$gene_name)
  if (any(duplicated(gene_labels))) {
    dup_idx <- duplicated(gene_labels)
    gene_labels[dup_idx] <- paste0(gene_labels[dup_idx], "_", top50$gene_id[dup_idx])
  }
  rownames(mat) <- gene_labels

  ann_df <- data.frame(Group = as.character(meta$Group[meta$sample_id %in% smp_comp]),
                       row.names = smp_comp)

  pheatmap(mat, annotation_col=ann_df,
           filename=file.path(OUT_DIR, paste0("Heatmap_top50_", comp_name, ".pdf")),
           show_rownames=TRUE,
           main=paste("Top 50 DEGs:", comp_name),
           fontsize=7, fontsize_row=5, fontsize_col=7, fontfamily="sans")

  cat("Heatmap saved:", comp_name, "\n")
}

# ================= 10. 保存 DEG 汇总 + res_list (供 script 5 使用) =================
saveRDS(res_list,    file.path(OUT_DIR, "res_list.rds"))
saveRDS(vsd,         file.path(OUT_DIR, "vsd.rds"))
saveRDS(meta,        file.path(OUT_DIR, "meta.rds"))
saveRDS(counts_raw,  file.path(OUT_DIR, "counts_raw.rds"))
saveRDS(list(
  n_original    = nrow(counts_mat),
  n_after_regex = sum(regex_filter),
  n_final       = nrow(counts_mat_filtered),
  n_samples     = n_samples,
  low_count_min = 10,
  low_count_min_samples = n_samples - 2
), file.path(OUT_DIR, "filter_stats.rds"))
cat("\nRDS objects saved for enrichment analysis (script 5)\n")

cat("\nAll DE/PCA analyses complete. Results:", OUT_DIR, "\n")
