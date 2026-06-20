#!/usr/bin/env Rscript
# ==========================================================
# 单独核查 padj<0.05 的显著基因 (pi5 vs NC)
# 背景：4_run_DE_PCA_final.R 里 ashr shrinkage 把效应量压得太小，
# 导致"padj<0.05 且 |log2FC|>=0.585"的双重过滤后 0 个基因入选；
# 但 padj 本身不受 shrinkage 影响，仍有 16 个基因 padj<0.05。
# 这里把这16个基因单独抽出来：
#   1. 输出每个样本的原始count + DESeq2标准化count + 收缩前/收缩后log2FC，方便人工核查
#      是否由单个样本的离群值驱动（尤其NC/pi5两组文库大小相差约27%，需要排除是否是
#      标准化没扣干净的残余效应）
#   2. 画热图 (z-score, VST标准化)
#   3. 画逐基因的样本点图 (per-sample dot plot)，比热图更容易看出是否单样本驱动
# ==========================================================

library(DESeq2)
library(ashr)
library(dplyr)
library(readr)
library(pheatmap)
library(ggplot2)
library(tidyr)

setwd("/home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/")

META_FILE <- "op.csv"
COUNT_FILE <- "../output_results/star_salmon/salmon.merged.gene_counts.tsv"
OUT_DIR <- "../Data_Analysis/DE_PCA_Results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. 复现 4_run_DE_PCA_final.R 里同样的建模流程，保证dds和之前完全一致 ----
meta_raw <- read_csv(META_FILE, show_col_types = FALSE)
meta <- meta_raw %>%
  select(Group, `Name in File`) %>%
  rename(sample_id = `Name in File`) %>%
  filter(!is.na(Group)) %>%
  mutate(Group = factor(Group, levels = c("NC", "pi5")))

counts_raw <- read_tsv(COUNT_FILE, col_types = cols())
all_sample_cols <- colnames(counts_raw)[3:ncol(counts_raw)]
meta <- meta[meta$sample_id %in% all_sample_cols, ]

counts_mat <- as.matrix(counts_raw[, all_sample_cols])
rownames(counts_mat) <- counts_raw$gene_id
counts_mat <- counts_mat[, meta$sample_id]
counts_mat <- round(counts_mat)
keep <- rowSums(counts_mat >= 10) >= 4
counts_mat <- counts_mat[keep, ]

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = meta, design = ~Group)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

PADJ_THRESHOLD <- 0.05 # 与 4_run_DE_PCA_final.R 保持一致

# ---- 2. 取未收缩(MLE)和ashr收缩后两版结果，分别记录 ----
res_unshrunk <- results(dds, contrast = c("Group", "pi5", "NC"))
res_shrunk <- lfcShrink(dds, contrast = c("Group", "pi5", "NC"), type = "ashr")

unshrunk_df <- as.data.frame(res_unshrunk) %>%
  mutate(gene_id = rownames(res_unshrunk)) %>%
  select(gene_id, log2FC_unshrunk = log2FoldChange, lfcSE_unshrunk = lfcSE)

shrunk_df <- as.data.frame(res_shrunk) %>%
  mutate(gene_id = rownames(res_shrunk)) %>%
  select(gene_id, log2FC_shrunk_ashr = log2FoldChange, lfcSE_shrunk = lfcSE, baseMean, pvalue, padj)

# ---- 3. 只保留 padj<0.05 的基因 (不论 |log2FC| 是否达到 0.585) ----
sig_genes <- shrunk_df %>%
  filter(!is.na(padj), padj < PADJ_THRESHOLD) %>%
  left_join(unshrunk_df, by = "gene_id") %>%
  left_join(counts_raw[, c("gene_id", "gene_name")], by = "gene_id") %>%
  arrange(padj)

cat("✅ padj <", PADJ_THRESHOLD, "的基因数:", nrow(sig_genes), "\n")

# ---- 4. 拼上每个样本的原始count + DESeq2标准化count，供人工核查 ----
raw_counts_sub <- as.data.frame(counts_mat[sig_genes$gene_id, , drop = FALSE])
colnames(raw_counts_sub) <- paste0("raw_", colnames(raw_counts_sub))
raw_counts_sub$gene_id <- rownames(raw_counts_sub)

norm_counts <- counts(dds, normalized = TRUE)
norm_counts_sub <- as.data.frame(norm_counts[sig_genes$gene_id, , drop = FALSE])
colnames(norm_counts_sub) <- paste0("norm_", colnames(norm_counts_sub))
norm_counts_sub$gene_id <- rownames(norm_counts_sub)

check_table <- sig_genes %>%
  left_join(raw_counts_sub, by = "gene_id") %>%
  left_join(norm_counts_sub, by = "gene_id") %>%
  select(
    gene_id, gene_name, baseMean,
    log2FC_unshrunk, lfcSE_unshrunk, log2FC_shrunk_ashr, lfcSE_shrunk,
    pvalue, padj,
    starts_with("raw_"), starts_with("norm_")
  )

write_csv(check_table, file.path(OUT_DIR, "Sig_padj_genes_manual_check.csv"))
cat("✅ 人工核查表已保存:", file.path(OUT_DIR, "Sig_padj_genes_manual_check.csv"), "\n")

# ---- 5. 热图 (VST标准化, z-score by gene) ----
mat <- assay(vsd)[sig_genes$gene_id, , drop = FALSE]
mat <- t(scale(t(mat)))

gene_labels <- ifelse(is.na(sig_genes$gene_name) | sig_genes$gene_name == "",
                       sig_genes$gene_id, sig_genes$gene_name)
if (any(duplicated(gene_labels))) {
  dup_idx <- duplicated(gene_labels)
  gene_labels[dup_idx] <- paste0(gene_labels[dup_idx], "_", sig_genes$gene_id[dup_idx])
}
rownames(mat) <- gene_labels

annotation_df <- data.frame(Group = as.character(meta$Group), row.names = meta$sample_id)

pheatmap(
  mat,
  annotation_col = annotation_df,
  filename = file.path(OUT_DIR, "Heatmap_padj_sig_genes_pi5_vs_NC.pdf"),
  show_rownames = TRUE,
  main = paste0("padj<", PADJ_THRESHOLD, " genes (n=", nrow(sig_genes), "): pi5 vs NC"),
  fontsize = 9,
  fontsize_row = 8,
  fontsize_col = 9,
  fontfamily = "sans"
)
cat("✅ 热图已保存:", file.path(OUT_DIR, "Heatmap_padj_sig_genes_pi5_vs_NC.pdf"), "\n")

# ---- 6. 逐基因样本点图：比热图更容易看出是否被单个样本驱动 ----
norm_long <- as.data.frame(norm_counts[sig_genes$gene_id, , drop = FALSE]) %>%
  mutate(gene_id = rownames(.)) %>%
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "norm_count") %>%
  left_join(data.frame(sample_id = meta$sample_id, Group = meta$Group), by = "sample_id") %>%
  left_join(sig_genes[, c("gene_id", "gene_name", "padj")], by = "gene_id") %>%
  mutate(
    gene_label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name),
    gene_label = factor(gene_label, levels = unique(gene_label[order(padj)]))
  )

p_check <- ggplot(norm_long, aes(x = Group, y = norm_count, color = Group)) +
  geom_point(size = 2, position = position_jitter(width = 0.1, height = 0)) +
  facet_wrap(~gene_label, scales = "free_y", ncol = 4) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 9) +
  labs(
    title = paste0("Per-sample normalized counts, padj<", PADJ_THRESHOLD, " genes (pi5 vs NC)"),
    subtitle = "每个点=一个样本；若某组4个点里有1个明显偏离其余3个，说明该基因的显著性可能是被单样本驱动的",
    y = "DESeq2 normalized count",
    x = NULL
  ) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7))

ggsave(
  file.path(OUT_DIR, "Check_padj_sig_genes_per_sample_dotplot.png"),
  p_check, width = 12, height = 10, dpi = 300
)
cat("✅ 逐基因样本点图已保存:", file.path(OUT_DIR, "Check_padj_sig_genes_per_sample_dotplot.png"), "\n")

cat("\n🎉 核查完成，请人工查看以上3个输出文件\n")
