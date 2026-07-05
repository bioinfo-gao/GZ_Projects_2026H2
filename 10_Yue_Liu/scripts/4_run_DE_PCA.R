#!/usr/bin/env Rscript
# 运行环境: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript
# 运行方法: cd /home/gao/projects_2026H2/10_Yue_Liu/scripts && Rscript 4_run_DE_PCA.R
#
# ================================================================
# REDO #2 (2026-07-04): split into RELIABLE (same-batch) vs UNRELIABLE (cross-batch)
# ================================================================
# Background: Control/4h come from a client-supplied LEGACY sequencing batch;
# 8h/16h come from the CURRENT NovaSeq batch (this project's own nf-core run).
# No condition was sequenced in both batches.
#
# CRITICAL FINDING that triggered this redo: fitting ONE DESeq2 model across
# all 12 samples (~Group, 4 levels) contaminates the dispersion/mean-trend
# estimation so badly that even the *same-batch* 4h vs Control contrast
# produced false positives. Diagnostic case: FSCN1 (ENSG00000075618) —
# raw counts Control=(528,618,415) vs 4h=(593,512,434), essentially flat —
# was flagged padj=0.0017 "Up" in the merged 4-group model, but padj=0.9996
# (no effect) when Control+4h are modeled in ISOLATION. Root cause: baseMean
# for this gene jumps from 521 (Control+4h only) to 3020 when 8h/16h samples
# (which express this gene 10x higher) are included in the same fit, which
# drags the mean-dispersion trend curve and corrupts the shrinkage for every
# gene, not just the ones genuinely affected by 8h/16h.
#
# FIX: run TWO fully independent DESeq2 models (no shared dispersion estimation):
#   Model A (RELIABLE) — legacy batch only:  Control vs 4h        (6 samples)
#   Model B (RELIABLE) — new batch only:     8h vs 16h            (6 samples)
#   Model C (UNRELIABLE, for reference only) — all 12 samples merged, used
#            only to produce the "4-timepoint" cross-batch view the client
#            requested for visual reference: PCA (all 4 groups) + the two
#            cross-batch contrasts (8h vs Control, 16h vs Control). These
#            numbers are flagged unreliable and are NOT used for any
#            conclusion, and are NOT run through enrichment analysis.
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
NEW_COUNT_FILE <- "../output_results/star_salmon/salmon.merged.gene_counts.tsv"
NEW_TPM_FILE   <- "../output_results/star_salmon/salmon.merged.gene_tpm.tsv"
OLD_COUNT_FILE <- "../OLD_Count/counts-filtered protein coding original file.xlsx"
DATA_DIR       <- paste0("../Data_Analysis_", TODAY)
REL_DIR        <- file.path(DATA_DIR, "DE_PCA_Results_Reliable")
UNREL_DIR      <- file.path(DATA_DIR, "DE_PCA_Results_Unreliable_CrossBatch")
READS_DIR      <- file.path(DATA_DIR, "Reads")

dir.create(REL_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(UNREL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(READS_DIR, showWarnings = FALSE, recursive = TRUE)

# ================= 2. 拷贝原始计数文件 =================
file.copy(NEW_COUNT_FILE, file.path(READS_DIR, "New_batch_8h_16h_gene_counts.tsv"), overwrite = TRUE)
file.copy(NEW_TPM_FILE,   file.path(READS_DIR, "New_batch_8h_16h_gene_tpm.tsv"),    overwrite = TRUE)
file.copy(OLD_COUNT_FILE, file.path(READS_DIR, "Old_batch_Control_4h_gene_counts.xlsx"), overwrite = TRUE)
cat("Copied raw count files (both batches) to Reads/\n")

# ================= 2.5 拷贝基因注释 + QC 文件 =================
annot_files <- sort(
  list.files("/Work_bio/references/Genes",
             pattern = "^human_Gene_annotation_.*\\.xlsx$", full.names = TRUE),
  decreasing = TRUE
)
ANNOT_SRC <- NA
if (length(annot_files) > 0) {
  ANNOT_SRC  <- annot_files[1]
  file.copy(ANNOT_SRC, file.path(DATA_DIR, basename(ANNOT_SRC)), overwrite = TRUE)
  cat("Gene annotation copied:", basename(ANNOT_SRC), "\n")
} else {
  cat("No human_Gene_annotation_*.xlsx found, will fall back to regex filter\n")
}

QC_DEST_BASE <- file.path(DATA_DIR, "QC")
dir.create(QC_DEST_BASE, showWarnings = FALSE, recursive = TRUE)
if (dir.exists("../output_results/multiqc")) {
  file.copy("../output_results/multiqc", QC_DEST_BASE, recursive = TRUE, overwrite = TRUE)
  cat("QC copied: multiqc ->", QC_DEST_BASE, "(covers new batch 8h/16h only; no QC available for legacy batch)\n")
}

# ================= 3. 加载并合并两个批次的计数矩阵 =================
new_counts_raw <- read_tsv(NEW_COUNT_FILE, col_types = cols())
new_counts_raw$base_id <- sub("\\..*$", "", new_counts_raw$gene_id)

rename_map <- c(A_1="8h_1", A_2="8h_2", A_3="8h_3",
                `16_1`="16h_1", `16_2`="16h_2", `16_3`="16h_3")
new_sample_cols <- names(rename_map)
stopifnot(all(new_sample_cols %in% colnames(new_counts_raw)))
new_mat <- as.matrix(new_counts_raw[, new_sample_cols])
colnames(new_mat) <- rename_map[colnames(new_mat)]
rownames(new_mat) <- new_counts_raw$base_id
new_mat <- round(new_mat)

old_counts_raw <- read_excel(OLD_COUNT_FILE)
old_counts_raw$base_id <- sub("\\..*$", "", old_counts_raw$Geneid)
old_sample_cols <- c("Control-1", "Control-2", "Control-3", "X-Ray-1", "X-Ray-2", "X-Ray-3")
stopifnot(all(old_sample_cols %in% colnames(old_counts_raw)))
old_mat <- as.matrix(old_counts_raw[, old_sample_cols])
colnames(old_mat) <- c("Control_1", "Control_2", "Control_3", "4h_1", "4h_2", "4h_3")
rownames(old_mat) <- old_counts_raw$base_id
old_mat <- round(old_mat)

common_ids <- intersect(rownames(new_mat), rownames(old_mat))
cat("New batch genes:", nrow(new_mat), " | Old batch genes:", nrow(old_mat),
    " | Common (base ID) genes:", length(common_ids), "\n")

old_ids_not_in_new <- setdiff(rownames(old_mat), rownames(new_mat))
gene_overlap_check <- list(
  old_total_genes  = nrow(old_mat),
  common_genes     = length(common_ids),
  missing_from_new = length(old_ids_not_in_new),
  missing_gene_ids = old_ids_not_in_new
)

counts_mat  <- cbind(old_mat[common_ids, ], new_mat[common_ids, ])
gene_id_map <- new_counts_raw %>% distinct(base_id, gene_id, gene_name) %>%
  filter(base_id %in% common_ids)
gene_info <- gene_id_map %>% rename(gene_id_full = gene_id)

# ================= 4. 基因过滤 (protein_coding only, 无低表达过滤 — 与客户旧文件一致) =================
if (!is.na(ANNOT_SRC)) {
  annot <- read_excel(ANNOT_SRC) %>% select(gene_id, gene_type) %>%
    mutate(base_id = sub("\\..*$", "", gene_id))
  gene_info <- gene_info %>% left_join(annot %>% select(base_id, gene_type), by = "base_id")
  regex_filter <- !is.na(gene_info$gene_type) & gene_info$gene_type == "protein_coding"
  cat("Using annotation-based protein_coding filter:", basename(ANNOT_SRC), "\n")

  old_common_annot <- annot %>% filter(base_id %in% common_ids)
  gene_overlap_check$old_common_protein_coding <- sum(old_common_annot$gene_type == "protein_coding", na.rm = TRUE)
  gene_overlap_check$old_common_biotype_mismatch <- old_common_annot %>%
    filter(gene_type != "protein_coding") %>% count(gene_type)
} else {
  ribo_pattern <- "^RPL|^RPS|^MRPL|^MRPS|^RPLP|^RPSA"
  non_ribosomal <- !grepl(ribo_pattern, gene_info$gene_name, ignore.case = TRUE)
  noncoding_pattern <- paste0(
    "^MT-", "|^SNORD|^SNORA|^RNU", "|^MALAT1|^NEAT1|^XIST",
    "|^MIR[0-9]|^LET-", "|^LINC|pseudogene|antisense"
  )
  protein_coding_like <- !grepl(noncoding_pattern, gene_info$gene_name, ignore.case = TRUE)
  regex_filter <- non_ribosomal & protein_coding_like
  cat("Using regex fallback filter (no annotation xlsx found)\n")
}
names(regex_filter) <- gene_info$base_id
regex_filter <- regex_filter[rownames(counts_mat)]

# NOTE: deliberately no low-count filter (client's legacy file keeps all-zero genes too).
counts_mat_filtered <- counts_mat[regex_filter, ]
gene_name_lookup <- setNames(gene_info$gene_name, gene_info$base_id)

cat("Original (common) genes:", nrow(counts_mat), "\n")
cat("After protein_coding filter:", nrow(counts_mat_filtered), "\n")

write_tsv(as.data.frame(counts_mat_filtered) %>%
            mutate(gene_id = rownames(counts_mat_filtered), gene_name = gene_name_lookup[gene_id], .before = 1),
          file.path(UNREL_DIR, "All_12samples_gene_counts.tsv"))

# ================= 辅助函数: 跑一个独立 DESeq2 模型 (同批次, 无跨批次污染) =================
run_one_contrast <- function(count_subset, meta_sub, grp_trt, grp_ctrl, out_dir, pca_title) {
  comp_name <- paste(grp_trt, "vs", grp_ctrl, sep = "_")
  cat("\n====", comp_name, "(", pca_title, ") ====\n")

  dds <- DESeqDataSetFromMatrix(countData = count_subset, colData = meta_sub, design = ~Group)
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = FALSE)

  # PCA
  pca_data   <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  p_pca <- ggplot(pca_data, aes(PC1, PC2, color = Group, label = name)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 20) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(title = pca_title) +
    theme_bw(base_size = 12) +
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = "bottom", plot.title = element_text(size = 9, hjust = 0.5))
  ggsave(file.path(out_dir, paste0("PCA_", comp_name, ".pdf")), p_pca, width = 8, height = 6.5, dpi = 300)

  # DE
  sig_col <- "sig (padj<=0.05 & |log2FC|>=0.263)"
  res <- lfcShrink(dds, contrast = c("Group", grp_trt, grp_ctrl), type = "ashr")
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_name <- gene_name_lookup[res_df$gene_id]

  samples_trt  <- meta_sub$sample_id[meta_sub$Group == grp_trt]
  samples_ctrl <- meta_sub$sample_id[meta_sub$Group == grp_ctrl]
  raw_sub <- as.data.frame(count_subset[, c(samples_trt, samples_ctrl)])
  raw_sub$gene_id <- rownames(raw_sub)
  res_df <- left_join(res_df, raw_sub, by = "gene_id")

  final_cols <- c("gene_id", "gene_name", samples_trt, samples_ctrl,
                  "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
  res_df <- res_df[, final_cols] %>% arrange(padj)
  res_df[[sig_col]] <- case_when(
    res_df$padj <= 0.05 & res_df$log2FoldChange >=  0.263 ~ "Up",
    res_df$padj <= 0.05 & res_df$log2FoldChange <= -0.263 ~ "Down",
    TRUE ~ "NS"
  )
  write_csv(res_df, file.path(out_dir, paste0("DEG_", comp_name, ".csv")))

  # Volcano
  res_df$negLogPadj <- -log10(res_df$padj)
  res_df$negLogPadj[!is.finite(res_df$negLogPadj)] <- NA
  top_labels <- res_df %>% filter(.data[[sig_col]] != "NS", !is.na(negLogPadj)) %>%
    arrange(padj) %>% head(10) %>%
    mutate(label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name))
  p_vol <- ggplot(res_df, aes(x=log2FoldChange, y=negLogPadj, color=.data[[sig_col]])) +
    geom_point(alpha=0.7, size=0.5) +
    scale_color_manual(values = c("Up"="#E41A1C","Down"="#377EB8","NS"="grey80"),
                       labels = c("Up"="Upregulated","Down"="Downregulated","NS"="Not Significant")) +
    theme_bw(base_size = 10) +
    labs(title = paste(grp_trt, "vs", grp_ctrl, "(log2FC>0 =", grp_trt, "upregulated)"),
         x = "log2 Fold Change", y = "-log10(adj. P-value)") +
    theme(plot.title=element_text(hjust=0.5, size=9), legend.position="bottom") +
    geom_text_repel(data=top_labels, aes(label=label), size=2, box.padding=0.3, max.overlaps=20, color="black")
  ggsave(file.path(out_dir, paste0("Volcano_", comp_name, ".png")), p_vol, width=8, height=6, dpi=300)

  # Heatmap
  top50 <- res_df %>% filter(.data[[sig_col]] != "NS") %>% arrange(padj) %>% head(50)
  if (nrow(top50) > 0) {
    mat <- assay(vsd)[top50$gene_id, , drop = FALSE]
    mat <- t(scale(t(mat)))
    gene_labels <- ifelse(is.na(top50$gene_name) | top50$gene_name == "", top50$gene_id, top50$gene_name)
    if (any(duplicated(gene_labels))) {
      dup_idx <- duplicated(gene_labels)
      gene_labels[dup_idx] <- paste0(gene_labels[dup_idx], "_", top50$gene_id[dup_idx])
    }
    rownames(mat) <- gene_labels
    ann_df <- data.frame(Group = as.character(meta_sub$Group), row.names = meta_sub$sample_id)
    pheatmap(mat, annotation_col = ann_df,
             filename = file.path(out_dir, paste0("Heatmap_top50_", comp_name, ".pdf")),
             show_rownames = TRUE, main = paste("Top 50 DEGs:", comp_name),
             fontsize = 7, fontsize_row = 5, fontsize_col = 7, fontfamily = "sans")
  } else {
    cat("No significant DEGs, skip heatmap\n")
  }

  cat(comp_name, "— significant DEGs:", sum(res_df[[sig_col]] != "NS"), "(baseMean-sanity, isolated model)\n")
  list(res_df = res_df, vsd = vsd)
}

# ================= 5. MODEL A (RELIABLE): 旧批次单独 — Control vs 4h =================
meta_A <- data.frame(
  sample_id = c("Control_1","Control_2","Control_3","4h_1","4h_2","4h_3"),
  Group = factor(c("Control","Control","Control","4h","4h","4h"), levels = c("Control","4h")),
  row.names = c("Control_1","Control_2","Control_3","4h_1","4h_2","4h_3")
)
out_A <- run_one_contrast(counts_mat_filtered[, meta_A$sample_id], meta_A, "4h", "Control",
                           REL_DIR, "PCA — Control vs 4h (legacy batch only, same-batch, reliable)")

# ================= 6. MODEL B (RELIABLE): 新批次单独 — 16h vs 8h =================
meta_B <- data.frame(
  sample_id = c("8h_1","8h_2","8h_3","16h_1","16h_2","16h_3"),
  Group = factor(c("8h","8h","8h","16h","16h","16h"), levels = c("8h","16h")),
  row.names = c("8h_1","8h_2","8h_3","16h_1","16h_2","16h_3")
)
out_B <- run_one_contrast(counts_mat_filtered[, meta_B$sample_id], meta_B, "16h", "8h",
                           REL_DIR, "PCA — 8h vs 16h (current NovaSeq batch only, same-batch, reliable)")

res_list_reliable <- list(
  "4h_vs_Control" = out_A$res_df,
  "16h_vs_8h"      = out_B$res_df
)

# ================= 7. MODEL C (UNRELIABLE, reference only): 全部12样本合并 =================
meta_full <- data.frame(
  sample_id = colnames(counts_mat_filtered),
  Group = factor(sub("_[0-9]$", "", colnames(counts_mat_filtered)), levels = c("Control","4h","8h","16h")),
  Batch = ifelse(colnames(counts_mat_filtered) %in% c("Control_1","Control_2","Control_3","4h_1","4h_2","4h_3"),
                 "Old_batch_legacy", "New_batch_NovaSeq"),
  stringsAsFactors = FALSE
)
rownames(meta_full) <- meta_full$sample_id

dds_full <- DESeqDataSetFromMatrix(counts_mat_filtered, meta_full, design = ~Group)
dds_full <- DESeq(dds_full)
vsd_full <- vst(dds_full, blind = FALSE)

pca_data   <- plotPCA(vsd_full, intgroup = c("Group","Batch"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p_pca_full <- ggplot(pca_data, aes(PC1, PC2, color = Group, shape = Batch, label = name)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(size = 3, box.padding = 0.3, point.padding = 0.3, max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "UNRELIABLE — all 4 timepoints merged (cross-batch; do not use for conclusions)") +
  theme_bw(base_size = 12) + scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom", plot.title = element_text(size = 9, hjust = 0.5, color = "firebrick"))
ggsave(file.path(UNREL_DIR, "PCA_all_4_groups.pdf"), p_pca_full, width = 8, height = 6.5, dpi = 300)

sig_col <- "sig (padj<=0.05 & |log2FC|>=0.263)"
res_list_unreliable <- list()
for (grp_trt in c("8h", "16h")) {
  comp_name <- paste(grp_trt, "vs_Control", sep = "_")
  res <- lfcShrink(dds_full, contrast = c("Group", grp_trt, "Control"), type = "ashr")
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_name <- gene_name_lookup[res_df$gene_id]
  samples_trt  <- meta_full$sample_id[meta_full$Group == grp_trt]
  samples_ctrl <- meta_full$sample_id[meta_full$Group == "Control"]
  raw_sub <- as.data.frame(counts_mat_filtered[, c(samples_trt, samples_ctrl)])
  raw_sub$gene_id <- rownames(raw_sub)
  res_df <- left_join(res_df, raw_sub, by = "gene_id")
  res_df <- res_df[, c("gene_id","gene_name", samples_trt, samples_ctrl, "baseMean","log2FoldChange","lfcSE","pvalue","padj")] %>% arrange(padj)
  res_df[[sig_col]] <- case_when(
    res_df$padj <= 0.05 & res_df$log2FoldChange >=  0.263 ~ "Up",
    res_df$padj <= 0.05 & res_df$log2FoldChange <= -0.263 ~ "Down",
    TRUE ~ "NS"
  )
  write_csv(res_df, file.path(UNREL_DIR, paste0("DEG_", comp_name, ".csv")))
  res_list_unreliable[[comp_name]] <- res_df

  res_df$negLogPadj <- -log10(res_df$padj); res_df$negLogPadj[!is.finite(res_df$negLogPadj)] <- NA
  top_labels <- res_df %>% filter(.data[[sig_col]] != "NS", !is.na(negLogPadj)) %>% arrange(padj) %>% head(10) %>%
    mutate(label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name))
  p_vol <- ggplot(res_df, aes(x=log2FoldChange, y=negLogPadj, color=.data[[sig_col]])) +
    geom_point(alpha=0.7, size=0.5) +
    scale_color_manual(values = c("Up"="#E41A1C","Down"="#377EB8","NS"="grey80")) +
    theme_bw(base_size = 10) +
    labs(title = paste0("UNRELIABLE (cross-batch): ", grp_trt, " vs Control"),
         x = "log2 Fold Change", y = "-log10(adj. P-value)") +
    theme(plot.title=element_text(hjust=0.5, size=9, color="firebrick"), legend.position="bottom") +
    geom_text_repel(data=top_labels, aes(label=label), size=2, box.padding=0.3, max.overlaps=20, color="black")
  ggsave(file.path(UNREL_DIR, paste0("Volcano_", comp_name, ".png")), p_vol, width=8, height=6, dpi=300)

  top50 <- res_df %>% filter(.data[[sig_col]] != "NS") %>% arrange(padj) %>% head(50)
  if (nrow(top50) > 0) {
    mat <- assay(vsd_full)[top50$gene_id, c(samples_trt, samples_ctrl), drop = FALSE]
    mat <- t(scale(t(mat)))
    gene_labels <- ifelse(is.na(top50$gene_name) | top50$gene_name == "", top50$gene_id, top50$gene_name)
    if (any(duplicated(gene_labels))) {
      dup_idx <- duplicated(gene_labels)
      gene_labels[dup_idx] <- paste0(gene_labels[dup_idx], "_", top50$gene_id[dup_idx])
    }
    rownames(mat) <- gene_labels
    ann_df <- data.frame(Group = as.character(meta_full$Group[meta_full$sample_id %in% c(samples_trt, samples_ctrl)]),
                         row.names = c(samples_trt, samples_ctrl))
    pheatmap(mat, annotation_col = ann_df,
             filename = file.path(UNREL_DIR, paste0("Heatmap_top50_", comp_name, ".pdf")),
             show_rownames = TRUE, main = paste0("UNRELIABLE — Top 50: ", comp_name),
             fontsize = 7, fontsize_row = 5, fontsize_col = 7, fontfamily = "sans")
  }
  cat(comp_name, "(UNRELIABLE, cross-batch) significant DEGs:", sum(res_df[[sig_col]] != "NS"), "\n")
}

# ================= 8. 保存 RDS =================
saveRDS(res_list_reliable, file.path(REL_DIR, "res_list.rds"))
saveRDS(meta_A, file.path(REL_DIR, "meta_A.rds"))
saveRDS(meta_B, file.path(REL_DIR, "meta_B.rds"))
saveRDS(out_A$vsd, file.path(REL_DIR, "vsd_A.rds"))
saveRDS(out_B$vsd, file.path(REL_DIR, "vsd_B.rds"))
saveRDS(gene_info, file.path(REL_DIR, "gene_info.rds"))

saveRDS(res_list_unreliable, file.path(UNREL_DIR, "res_list.rds"))
saveRDS(meta_full, file.path(UNREL_DIR, "meta.rds"))

saveRDS(list(
  n_original    = nrow(counts_mat),
  n_after_regex = sum(regex_filter),
  n_final       = nrow(counts_mat_filtered),
  low_count_filter_applied = FALSE,
  low_count_filter_note = "No low-expression filter applied — matches client's legacy 'counts-filtered protein coding' file, which retains all protein-coding genes including all-zero ones (4,746/20,065 genes are all-zero across its 6 samples). Only the protein_coding biotype filter is applied; DESeq2's built-in independent filtering still applies at the padj-correction stage.",
  batch_confound_note = "Control (n=3) and 4h (n=3) come from a client-supplied legacy sequencing batch; 8h (n=3) and 16h (n=3) come from the current NovaSeq X Plus batch. No condition overlaps both batches.",
  contamination_finding = "A merged 4-group (~Group, 12-sample) DESeq2 model was found to CONTAMINATE the dispersion/mean-trend estimation badly enough that even the same-batch 4h vs Control contrast produced a false positive (FSCN1/ENSG00000075618: padj=0.0017 in the merged model vs padj=0.9996 when Control+4h are modeled in isolation). Fix: Control-vs-4h and 8h-vs-16h are now each fit in their own independent 2-group DESeq2 model (no shared dispersion estimation). The merged 12-sample model is retained ONLY to produce the cross-batch reference PCA/DEG lists in DE_PCA_Results_Unreliable_CrossBatch/, which are explicitly flagged as not usable for conclusions.",
  gene_overlap_check = gene_overlap_check
), file.path(REL_DIR, "filter_stats.rds"))

cat("\nAll DE/PCA analyses complete.\n")
cat("RELIABLE (same-batch) results:", REL_DIR, "\n")
cat("UNRELIABLE (cross-batch, reference only) results:", UNREL_DIR, "\n")
