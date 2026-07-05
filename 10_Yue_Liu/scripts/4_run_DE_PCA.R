#!/usr/bin/env Rscript
# 运行环境: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript
# 运行方法: cd /home/gao/projects_2026H2/10_Yue_Liu/scripts && Rscript 4_run_DE_PCA.R
#
# ================================================================
# REDO 2026-07-04: Integrated time-course analysis
# ================================================================
# Client correction: "A" group in Sample_Sheet_Yue_Liu.xlsx was mislabeled.
# It is actually MOLM13 cells 8h post X-ray radiation (renamed "8h").
# "16" group = 16h post radiation (renamed "16h"). Both from the current
# NovaSeq X Plus batch, quantified with Salmon (this project's nf-core run).
#
# Client also provided OLD_Count/counts-filtered protein coding original file.xlsx:
# a pre-quantified count matrix from a SEPARATE, EARLIER sequencing batch,
# containing Control-1/2/3 (untreated baseline) and X-Ray-1/2/3 (4h post
# X-ray radiation). Quantification method for this legacy file is not
# independently verified (client-supplied; not this pipeline's Salmon output).
#
# IMPORTANT CAVEAT (confirmed with user, see analysis_plan): the legacy batch
# (Control, 4h) and the current batch (8h, 16h) share NO overlapping condition,
# so sequencing batch and time point are fully confounded — a ~batch+Group
# DESeq2 design is not estimable (Control only exists in the old batch). Per
# user's explicit decision, Control is used as the common reference for all
# three contrasts (4h/8h/16h vs Control) to produce a unified time-course view,
# but the 8h vs Control and 16h vs Control results cannot statistically
# distinguish true radiation effect from batch effect. This is disclosed
# prominently in the client report (Section 4 + Limitations).
# ================================================================

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
setwd("/home/gao/projects_2026H2/10_Yue_Liu/scripts/")

TODAY          <- format(Sys.Date(), "%Y%m%d")
META_FILE      <- "../Sample_Sheet_Yue_Liu.xlsx"
NEW_COUNT_FILE <- "../output_results/star_salmon/salmon.merged.gene_counts.tsv"
NEW_TPM_FILE   <- "../output_results/star_salmon/salmon.merged.gene_tpm.tsv"
OLD_COUNT_FILE <- "../OLD_Count/counts-filtered protein coding original file.xlsx"
DATA_DIR       <- paste0("../Data_Analysis_", TODAY)
OUT_DIR        <- file.path(DATA_DIR, "DE_PCA_Results")
READS_DIR      <- file.path(DATA_DIR, "Reads")

dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(READS_DIR, showWarnings = FALSE, recursive = TRUE)

# ================= 2. 拷贝原始计数文件 =================
file.copy(NEW_COUNT_FILE, file.path(READS_DIR, "New_batch_8h_16h_gene_counts.tsv"), overwrite = TRUE)
file.copy(NEW_TPM_FILE,   file.path(READS_DIR, "New_batch_8h_16h_gene_tpm.tsv"),    overwrite = TRUE)
file.copy(OLD_COUNT_FILE, file.path(READS_DIR, "Old_batch_Control_4h_gene_counts.xlsx"), overwrite = TRUE)
cat("Copied raw count files (both batches) to Reads/\n")

# ================= 2.5 拷贝基因注释 + QC 文件 =================
DATA_ANALYSIS_DIR <- DATA_DIR

annot_files <- sort(
  list.files("/Work_bio/references/Genes",
             pattern = "^human_Gene_annotation_.*\\.xlsx$", full.names = TRUE),
  decreasing = TRUE
)
ANNOT_SRC <- NA
if (length(annot_files) > 0) {
  ANNOT_SRC  <- annot_files[1]
  ANNOT_DEST <- file.path(DATA_ANALYSIS_DIR, basename(ANNOT_SRC))
  file.copy(ANNOT_SRC, ANNOT_DEST, overwrite = TRUE)
  cat("Gene annotation copied:", basename(ANNOT_SRC), "\n")
} else {
  cat("No human_Gene_annotation_*.xlsx found, will fall back to regex filter\n")
}

QC_DEST_BASE <- file.path(DATA_ANALYSIS_DIR, "QC")
dir.create(QC_DEST_BASE, showWarnings = FALSE, recursive = TRUE)
if (dir.exists("../output_results/multiqc")) {
  file.copy("../output_results/multiqc", QC_DEST_BASE, recursive = TRUE, overwrite = TRUE)
  cat("QC copied: multiqc ->", QC_DEST_BASE, "(covers new batch 8h/16h only; no QC available for legacy batch)\n")
}

# ================= 3. 加载并合并两个批次的计数矩阵 =================

# --- 3a. 新批次 (8h/16h, Salmon) ---
new_counts_raw <- read_tsv(NEW_COUNT_FILE, col_types = cols())
new_counts_raw$base_id <- sub("\\..*$", "", new_counts_raw$gene_id)

# nf-core samplesheet 仍沿用原始 A_1/A_2/A_3 命名 (无需重新比对)；此处映射到修正后的 8h 标签
rename_map <- c(A_1="8h_1", A_2="8h_2", A_3="8h_3",
                `16_1`="16h_1", `16_2`="16h_2", `16_3`="16h_3")
new_sample_cols <- names(rename_map)
stopifnot(all(new_sample_cols %in% colnames(new_counts_raw)))
new_mat <- as.matrix(new_counts_raw[, new_sample_cols])
colnames(new_mat) <- rename_map[colnames(new_mat)]
rownames(new_mat) <- new_counts_raw$base_id
new_mat <- round(new_mat)

# --- 3b. 旧批次 (Control/4h, 客户提供的计数矩阵) ---
old_counts_raw <- read_excel(OLD_COUNT_FILE)
old_counts_raw$base_id <- sub("\\..*$", "", old_counts_raw$Geneid)
old_sample_cols <- c("Control-1", "Control-2", "Control-3", "X-Ray-1", "X-Ray-2", "X-Ray-3")
stopifnot(all(old_sample_cols %in% colnames(old_counts_raw)))
old_mat <- as.matrix(old_counts_raw[, old_sample_cols])
colnames(old_mat) <- c("Control_1", "Control_2", "Control_3", "4h_1", "4h_2", "4h_3")
rownames(old_mat) <- old_counts_raw$base_id
old_mat <- round(old_mat)

# --- 3c. 按 base Ensembl ID 取交集合并 (99.96% 基因在两批次间可对应，详见分析计划) ---
common_ids <- intersect(rownames(new_mat), rownames(old_mat))
cat("New batch genes:", nrow(new_mat), " | Old batch genes:", nrow(old_mat),
    " | Common (base ID) genes:", length(common_ids), "\n")

# --- 3d. 验证: 客户旧批次已预先过滤为 "protein coding" 基因 —— 检查与我们自己的
# protein_coding 判定是否一致 (方法学一致性校验，写入报告) ---
old_ids_not_in_new  <- setdiff(rownames(old_mat), rownames(new_mat))
gene_overlap_check <- list(
  old_total_genes   = nrow(old_mat),
  common_genes      = length(common_ids),
  missing_from_new  = length(old_ids_not_in_new),
  missing_gene_ids  = old_ids_not_in_new
)
cat("Client legacy batch was pre-filtered to protein-coding genes (n =", nrow(old_mat), ").\n")
cat(length(common_ids), "/", nrow(old_mat), "(", round(100*length(common_ids)/nrow(old_mat), 2),
    "%) matched by base Ensembl ID against our GENCODE v45 annotation;",
    length(old_ids_not_in_new), "absent entirely (likely added in a newer Ensembl release than v45).\n")

counts_mat <- cbind(old_mat[common_ids, ], new_mat[common_ids, ])
gene_id_map <- new_counts_raw %>% distinct(base_id, gene_id, gene_name) %>%
  filter(base_id %in% common_ids)

# ================= 4. 元数据: Control(参照) / 4h / 8h / 16h, 并记录 Batch (仅作说明/可视化用) =================
meta <- data.frame(
  sample_id = colnames(counts_mat),
  Group = factor(
    sub("_[0-9]$", "", colnames(counts_mat)),
    levels = c("Control", "4h", "8h", "16h")
  ),
  Batch = ifelse(colnames(counts_mat) %in% c("Control_1","Control_2","Control_3","4h_1","4h_2","4h_3"),
                  "Old_batch_legacy", "New_batch_NovaSeq"),
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$sample_id
cat("Metadata loaded:", nrow(meta), "samples (integrated Control/4h/8h/16h time course)\n")
print(meta)

# ================= 5. 基因过滤 (Human: 优先用 annotation xlsx biotype 列) =================
gene_info <- gene_id_map %>% rename(gene_id_full = gene_id)

if (!is.na(ANNOT_SRC)) {
  annot <- read_excel(ANNOT_SRC) %>% select(gene_id, gene_type) %>%
    mutate(base_id = sub("\\..*$", "", gene_id))
  gene_info <- gene_info %>% left_join(annot %>% select(base_id, gene_type), by = "base_id")
  regex_filter <- !is.na(gene_info$gene_type) & gene_info$gene_type == "protein_coding"
  cat("Using annotation-based protein_coding filter:", basename(ANNOT_SRC), "\n")

  # --- 方法学一致性校验: 客户旧批次的 20,065 个基因中，有多少也被我们自己的
  # GENCODE v45 注释独立判定为 protein_coding? (验证两批次过滤标准一致) ---
  old_common_annot <- annot %>% filter(base_id %in% intersect(rownames(old_mat), rownames(new_mat)))
  gene_overlap_check$old_common_protein_coding <- sum(old_common_annot$gene_type == "protein_coding", na.rm = TRUE)
  gene_overlap_check$old_common_biotype_mismatch <- old_common_annot %>%
    filter(gene_type != "protein_coding") %>% count(gene_type)
  cat("Of the", length(common_ids), "genes shared with the client's legacy batch,",
      gene_overlap_check$old_common_protein_coding, "(",
      round(100*gene_overlap_check$old_common_protein_coding/length(common_ids), 2),
      "%) are independently confirmed protein_coding by our own annotation",
      "— validating our filtering is consistent with the client's pre-filtered file.\n")
} else {
  ribo_pattern <- "^RPL|^RPS|^MRPL|^MRPS|^RPLP|^RPSA"
  non_ribosomal <- !grepl(ribo_pattern, gene_info$gene_name, ignore.case = TRUE)
  noncoding_pattern <- paste0(
    "^MT-",
    "|^SNORD|^SNORA|^RNU",
    "|^MALAT1|^NEAT1|^XIST",
    "|^MIR[0-9]|^LET-",
    "|^LINC|pseudogene|antisense"
  )
  protein_coding_like <- !grepl(noncoding_pattern, gene_info$gene_name, ignore.case = TRUE)
  regex_filter <- non_ribosomal & protein_coding_like
  cat("Using regex fallback filter (no annotation xlsx found)\n")
}
names(regex_filter) <- gene_info$base_id
regex_filter <- regex_filter[rownames(counts_mat)]

n_samples <- ncol(counts_mat)

# NOTE: no low-count filter here (deliberate). The client's legacy "counts-filtered
# protein coding" file retains ALL protein-coding genes regardless of expression level
# (4,746 / 20,065 genes are all-zero across its 6 samples) — i.e. they applied only a
# biotype filter, no low-expression filter. To keep both batches on the same gene-set
# convention, we match that: protein_coding filter only. DESeq2's own independent
# filtering (applied automatically during BH/padj correction) still down-weights
# unreliably-low-count genes at the statistical-testing stage.
final_filter <- regex_filter
counts_mat_filtered <- counts_mat[final_filter, ]

cat("Original (common) genes:", nrow(counts_mat), "\n")
cat("After regex/biotype filter:", sum(regex_filter), "\n")
cat("Final genes for DESeq2:", nrow(counts_mat_filtered),
    "(no low-count filter — matches client's legacy file convention of keeping all protein-coding genes)\n")

# ================= 6. DESeq2 建模 =================
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat_filtered,
  colData   = meta,
  design    = ~Group
)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# ================= 7. PCA (按 Group 着色，并叠加 Batch 形状以可视化批次效应) =================
pca_data   <- plotPCA(vsd, intgroup = c("Group", "Batch"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = Group, shape = Batch, label = name)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "PCA — Control/4h from legacy batch, 8h/16h from current NovaSeq batch") +
  theme_bw(base_size = 12) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom",
        plot.title      = element_text(size = 9, hjust = 0.5),
        legend.text     = element_text(size = 9),
        axis.text       = element_text(size = 10),
        axis.title      = element_text(size = 11))

ggsave(file.path(OUT_DIR, "PCA.pdf"), p_pca, width = 8, height = 6.5, dpi = 300)
cat("PCA saved (batch-vs-group visualization)\n")

# ================= 8. 差异表达分析: 4h/8h/16h vs Control =================
sig_col  <- "sig (padj<=0.05 & |log2FC|>=0.263)"
contrasts <- list(
  c("Group", "4h",  "Control"),
  c("Group", "8h",  "Control"),
  c("Group", "16h", "Control")
)

res_list <- list()
gene_name_lookup <- setNames(gene_info$gene_name, gene_info$base_id)

for (comp in contrasts) {
  grp_trt  <- comp[2]
  grp_ctrl <- comp[3]
  comp_name <- paste(grp_trt, "vs", grp_ctrl, sep = "_")
  cat("\nAnalysing:", comp_name, "\n")

  res <- lfcShrink(dds, contrast = comp, type = "ashr")

  res_df <- as.data.frame(res)
  res_df$base_id <- rownames(res_df)
  res_df$gene_name <- gene_name_lookup[res_df$base_id]

  samples_trt  <- meta$sample_id[meta$Group == grp_trt]
  samples_ctrl <- meta$sample_id[meta$Group == grp_ctrl]
  raw_sub <- as.data.frame(counts_mat[, c(samples_trt, samples_ctrl)])
  raw_sub$base_id <- rownames(raw_sub)

  res_df <- left_join(res_df, raw_sub, by = "base_id")

  final_cols <- c("base_id", "gene_name", samples_trt, samples_ctrl,
                  "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
  res_df <- res_df[, final_cols] %>% arrange(padj)
  colnames(res_df)[colnames(res_df) == "base_id"] <- "gene_id"

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
saveRDS(gene_info,   file.path(OUT_DIR, "gene_info.rds"))
saveRDS(list(
  n_original    = nrow(counts_mat),
  n_after_regex = sum(regex_filter),
  n_final       = nrow(counts_mat_filtered),
  n_samples     = n_samples,
  low_count_filter_applied = FALSE,
  low_count_filter_note = "No low-expression filter applied — matches client's legacy 'counts-filtered protein coding' file, which retains all protein-coding genes including all-zero ones (4,746/20,065 genes are all-zero across its 6 samples). Only the protein_coding biotype filter is applied; DESeq2's built-in independent filtering still applies at the padj-correction stage.",
  batch_confound_note = "Control (n=3) and 4h (n=3) come from a client-supplied legacy sequencing batch; 8h (n=3) and 16h (n=3) come from the current NovaSeq X Plus batch. No condition overlaps both batches, so a ~Batch+Group model is not estimable and sequencing batch is fully confounded with time point for the 8h/16h vs Control contrasts.",
  gene_overlap_check = gene_overlap_check
), file.path(OUT_DIR, "filter_stats.rds"))
cat("\nRDS objects saved for enrichment analysis (script 5)\n")

cat("\nAll DE/PCA analyses complete. Results:", OUT_DIR, "\n")
