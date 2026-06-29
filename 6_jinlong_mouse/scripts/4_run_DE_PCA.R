#!/usr/bin/env Rscript
# 运行环境: conda activate DE_R45
# 运行方法: cd /home/gao/projects_2026H2/6_jinlong_mouse/scripts && Rscript 4_run_DE_PCA.R
#
# 分析设计: Group1/2/3 vs Group4 (对照)
#   G1 vs G4, G2 vs G4, G3 vs G4
# log2FC > 0 = 在处理组上调; log2FC < 0 = 在对照组(G4)上调

library(DESeq2)
library(ashr)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(readr)
library(tidyr)
library(ggrepel)

# ================= 1. 路径设置 =================
setwd("/home/gao/projects_2026H2/6_jinlong_mouse/scripts/")

TODAY     <- format(Sys.Date(), "%Y%m%d")
META_FILE  <- "../jinlong.csv"
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
    cat("✅ Copied:", f$dst, "\n")
  } else {
    if (f$req) stop("❌ Required file missing: ", f$src)
    cat("⚠️  Optional file missing, skipped:", f$src, "\n")
  }
}

# ================= 2.5 拷贝基因注释 + QC 文件 =================
DATA_ANALYSIS_DIR <- DATA_DIR

# 基因注释：自动选取 Genes/ 下日期最新的 mouse xlsx，拷贝到 Data_Analysis 顶层
annot_files <- sort(
  list.files("/home/gao/projects_2026H1/Genes",
             pattern = "^mouse_Gene_annotation_.*\\.xlsx$", full.names = TRUE),
  decreasing = TRUE
)
if (length(annot_files) > 0) {
  ANNOT_SRC  <- annot_files[1]
  ANNOT_DEST <- file.path(DATA_ANALYSIS_DIR, basename(ANNOT_SRC))
  file.copy(ANNOT_SRC, ANNOT_DEST, overwrite = TRUE)
  cat("✅ 基因注释已拷贝:", basename(ANNOT_SRC), "\n")
} else {
  cat("⚠️  未找到 mouse_Gene_annotation_*.xlsx，跳过\n")
}

# QC 文件夹：multiqc / fastqc / pipeline_info → Data_Analysis/QC/
QC_DEST_BASE <- file.path(DATA_ANALYSIS_DIR, "QC")
dir.create(QC_DEST_BASE, showWarnings = FALSE, recursive = TRUE)
for (qc_src in c("../output_results/multiqc",
                  # "../output_results/fastqc",
                  # "../output_results/pipeline_info"
                )) {
  if (dir.exists(qc_src)) {
    file.copy(qc_src, QC_DEST_BASE, recursive = TRUE, overwrite = TRUE)
    cat("✅ QC 已拷贝:", basename(qc_src), "->", QC_DEST_BASE, "\n")
  }
}

# ================= 3. 元数据 =================
# Name in File: 数字 → sample_id J_XXX; 字母 → 保持原样
meta_raw <- read_csv(META_FILE)

meta <- meta_raw %>%
  select(Group, `Name in File`) %>%
  filter(!is.na(Group)) %>%
  mutate(
    sample_id = ifelse(grepl("^[0-9]+$", as.character(`Name in File`)),
                       paste0("J_", `Name in File`),
                       as.character(`Name in File`)),
    Group = paste0("G", Group),
    Group = factor(Group, levels = c("G4", "G1", "G2", "G3"))  # G4 = 对照/参照组
  ) %>%
  select(sample_id, Group)

cat("✅ Metadata loaded:", nrow(meta), "samples\n")
print(meta)

# ================= 4. 表达矩阵预处理 =================
counts_raw <- read_tsv(COUNT_FILE, col_types = cols())
head(counts_raw)

valid_samples <- colnames(counts_raw)[3:ncol(counts_raw)]
meta <- meta[meta$sample_id %in% valid_samples, ]
counts_mat <- as.matrix(counts_raw[, meta$sample_id])
rownames(counts_mat) <- counts_raw$gene_id
counts_mat <- round(counts_mat)

# ================= 5. 基因过滤 (鼠适配 regex) =================
# 没有鼠专用基因注释 xlsx，直接用 regex 过滤
gene_info <- counts_raw[, c("gene_id", "gene_name")]

# 核糖体基因 (鼠命名: Rpl/Rps/Mrpl/Mrps, 首字母大写)
ribo_pattern <- "^Rpl|^Rps|^Mrpl|^Mrps|^Rplp|^Rpsa"
non_ribosomal <- !grepl(ribo_pattern, gene_info$gene_name, ignore.case = TRUE)

# 非编码 RNA、线粒体基因 (鼠 mt 基因命名: mt-Co1 等, 小写 mt)
noncoding_pattern <- paste0(
  "^mt-",                         # 线粒体基因
  "|^Snord|^Snora|^Rnu",          # snoRNA / snRNA
  "|^Malat1|^Neat1|^Xist",        # 著名 lncRNA
  "|^Mir[0-9]|^let-",             # miRNA
  "|^Linc|^Gm[0-9]",             # lincRNA / 预测基因
  "|pseudogene|antisense"
)
protein_coding_like <- !grepl(noncoding_pattern, gene_info$gene_name, ignore.case = TRUE)

regex_filter <- non_ribosomal & protein_coding_like

# 低表达过滤：至少 (n-2) 个样本中 counts >= 10
n_samples <- ncol(counts_mat)
low_count_filter <- rowSums(counts_mat >= 10) >= (n_samples - 2)

final_filter <- regex_filter & low_count_filter
counts_mat_filtered <- counts_mat[final_filter, ]

cat("✅ Original genes:", nrow(counts_mat), "\n")
cat("✅ After regex filter:", sum(regex_filter), "\n")
cat("✅ After low-count filter:", sum(final_filter), "\n")
cat("✅ Final genes for DESeq2:", nrow(counts_mat_filtered), "\n")

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
cat("✅ PCA saved\n")

# ================= 8. 差异表达分析 =================
# G4 = 对照 (DESeq2 参照组 = 因子第一水平 = G4)
sig_col  <- "sig (padj<=0.05 & |log2FC|>=0.263)"
contrasts <- list(
  c("Group", "G1", "G4"),
  c("Group", "G2", "G4"),
  c("Group", "G3", "G4")
)

res_list <- list()

for (comp in contrasts) {
  grp_trt  <- comp[2]
  grp_ctrl <- comp[3]
  comp_name <- paste(grp_trt, "vs", grp_ctrl, sep = "_")
  cat("\n🔍 Analysing:", comp_name, "\n")

  res <- results(dds, contrast = comp)
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

  cat("✅", comp_name, "— significant DEGs:", sum(res_df[[sig_col]] != "NS"), "\n")
}

# ================= 9. 热图 =================
for (comp_name in names(res_list)) {
  cat("\n🔍 Heatmap:", comp_name, "\n")

  top50 <- res_list[[comp_name]] %>%
    filter(.data[[sig_col]] != "NS") %>% arrange(padj) %>% head(50)

  if (nrow(top50) == 0) { cat("⚠️  No significant DEGs, skip heatmap\n"); next }

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

  cat("✅ Heatmap saved:", comp_name, "\n")
}

# ================= 10. 保存 DEG 汇总 + res_list (供 script 5 使用) =================
saveRDS(res_list,    file.path(OUT_DIR, "res_list.rds"))
saveRDS(vsd,         file.path(OUT_DIR, "vsd.rds"))
saveRDS(meta,        file.path(OUT_DIR, "meta.rds"))
saveRDS(counts_raw,  file.path(OUT_DIR, "counts_raw.rds"))
cat("\n✅ RDS objects saved for enrichment analysis (script 5)\n")

cat("\n🎉 All DE/PCA analyses complete. Results:", OUT_DIR, "\n")
