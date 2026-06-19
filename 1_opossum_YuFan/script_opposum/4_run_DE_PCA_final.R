#!/usr/bin/env Rscript
# env in bash 各种包装在了 DE_R45 环境  
# mamba activate DE_R45               # mamba activate regular_bioinfo # === Regular_bioinfo lacks ggrepel and ashr ， by CC
# [MODIFIED-0] 补充R解释器路径注释（运行环境用 mamba 创建，非 conda）
# R interpreter: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/R

# [MODIFIED-1] setwd 改为本项目(opossum)的脚本目录，原为 2026_Item16_ZhenYan/scripts/
setwd("/home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/")
getwd()
# 跑完的输出文件：
# DEG_Test_1vsControl.csv
# DEG_Test_2vsControl.csv
# Volcano_Test_1vsControl.png
# Volcano_Test_2vsControl.png
# PCA.pdf
# Heatmap_top50_Test_1_vs_Control.pdf
# All_sample_gene_counts.tsv (拷贝的原始文件)

#!/usr/bin/env Rscript
# ==========================================================
# nf-core RNA-seq 下游分析：PCA + 差异表达 + 可视化 (针对LZJ项目修正版)
# 对比顺序: c("分组变量", "处理组/分子", "对照组/分母")
# log2FC = log2(处理组 / 对照组)
# 正数 = 在处理组中上调；负数 = 在对照组中上调
# ==========================================================

# ================= 0. 依赖加载 =================
library(DESeq2)
library(ashr)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(readr)
library(tidyr)
library(ggrepel)

# ================= 1. 路径设置 =================
META_FILE  <- "op.csv"                 
COUNT_FILE <- "../output_results/star_salmon/salmon.merged.gene_counts.tsv"
TPM_FILE   <- "../output_results/star_salmon/salmon.merged.gene_tpm.tsv"
OUT_DIR    <- "../Data_Analysis/DE_PCA_Results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

READS_DIR <- "../Data_Analysis/Reads"
dir.create(READS_DIR, showWarnings = FALSE, recursive = TRUE)

# ================= 2. 拷贝原始计数文件和TPM文件到目标目录 =================
# Create subdirectory for raw/normalized data

READS_DIR

# 定义需要拷贝的文件列表及其目标文件名
files_to_copy <- list(
  list(
    source = COUNT_FILE,
    dest = "All_sample_gene_counts.tsv",
    required = TRUE
  ),
  list(source = TPM_FILE, dest = "All_sample_gene_tpm.tsv", required = FALSE)
)

# 统一执行拷贝逻辑
for (file_info in files_to_copy) {
  src <- file_info$source
  dst_name <- file_info$dest
  is_required <- file_info$required

  if (file.exists(src)) {
    dst_path <- file.path(READS_DIR, dst_name)
    file.copy(src, dst_path, overwrite = TRUE)
    cat("✅ 文件已拷贝:", dst_name, "\n")
  } else {
    if (is_required) {
      stop("❌ 错误: 找不到必需文件: ", src)
    } else {
      cat("⚠️  警告: 可选文件不存在，已跳过: ", src, "\n")
    }
  }
}

# ================= 2.5 拷贝 QC 文件夹和生成分析报告 =================
# [MODIFIED-3] 删除了原来拷贝 human_Gene_annotation_20260202.xlsx 的代码块。
# 原因：本项目物种为 Didelphis virginiana（opossum），人类基因注释表完全不适用；
# 基因注释改由同目录下的 5B_annotation.R（GTF转注释表）单独处理，此处跳过。

# 定义 QC 文件夹源路径 (请根据实际路径调整，例如 nf-core 的 multiqc 输出或自定义 QC 文件夹)
# 假设 QC 文件夹位于项目根目录下的 output_results/multiqc 或当前目录下的 QC 文件夹
QC_SOURCE_CANDIDATES <- c(
  "../output_results/multiqc",
  "../output_results/pipeline_info",
  "QC"
)

QC_SRC <- NULL
for (path in QC_SOURCE_CANDIDATES) {
  if (dir.exists(path)) {
    QC_SRC <- path
    break
  }
}

QC_SRC


if (!is.null(QC_SRC)) {
  # Define the specific destination path for QC files
  # [MODIFIED-4] QC_DEST_DIR 改为本项目的 Data_Analysis/QC，原为 2026_Item16_ZhenYan/Data_Analysis/QC
  QC_DEST_DIR <- "/home/gao/projects_2026H2/1_opossum_YuFan/Data_Analysis/QC"

  # Create the destination directory if it doesn't exist
  dir.create(QC_DEST_DIR, showWarnings = FALSE, recursive = TRUE)

  # Copy contents of QC_SRC into QC_DEST_DIR
  # file.copy with recursive=TRUE copies the source folder INTO the destination if destination exists
  success <- file.copy(
    from = QC_SRC,
    to = QC_DEST_DIR,
    recursive = TRUE,
    overwrite = TRUE
  )

  if (any(success)) {
    cat("✅ QC 文件夹已拷贝:", QC_SRC, "->", QC_DEST_DIR, "\n")
  } else {
    cat("⚠️  警告: QC 文件夹拷贝失败\n")
  }
} else {
  cat("⚠️  警告: 未找到常见的 QC 文件夹路径，跳过拷贝\n")
}

# ================= 3. 读取并清洗元数据 =================
# ZG 需要按照真实组名修改 <<====================   # mutate(Group = factor(Group, levels = c("Control", "Test_1", "Test_2"))) ！！！
meta_raw <- read_csv(META_FILE)
meta_raw
# meta <- meta_raw %>%
#   select(Group, `Name in File`) 

# meta <- meta_raw %>%
#   select(Group, `Name in File`) %>%
#   rename(sample_id = `Name in File`) %>%
#   filter(!is.na(Group)) 
meta <- meta_raw %>%
  select(Group, `Name in File`) %>%       # 选择两列 “Group` 和 “Name in File”
  rename(sample_id = `Name in File`) %>%  # 重命名“Name in File”列名为 sample_id
  filter(!is.na(Group))           %>%       # 删除 Group 为 NA 的行
  # [MODIFIED-5] Group 水平改为本项目实际的两组 "NC","pi5"，原为 "CTRL","SMA4","SMC2","ME13" 四组
  mutate(Group = factor(Group, levels = c("NC", "pi5"))) # Group = factor(...): 将 Group 列转换为因子（factor）类型，并指定因子的水平（levels）顺序。

# levels = c("NC", "pi5"): 明确定义因子的水平顺序为："NC" (对照组) "pi5" (处理组)
# 为什么要这样做？
# 在 R 中，因子的水平顺序非常重要，特别是在进行统计分析时：
# DESeq2 差异表达分析会将第一个水平（"CTRL"）作为参照组（baseline/reference）
# 其他组会与参照组进行比较
# 因子水平的顺序会影响结果的解释和可视化（如 PCA 图中组别的颜色顺序）
# 通过显式设置 levels，可以确保分析的一致性和可重复性，避免 R 自动按字母顺序排列因子水平（那样可能会把 "CTRL" 排在后面）。

meta

cat("✅ 元数据加载完成，有效样本数:", nrow(meta), "\n")

# ================= 4. 读取表达矩阵 & 预处理 =================
counts_raw <- read_tsv(COUNT_FILE, col_types = cols())
head(counts_raw)

all_sample_cols <- colnames(counts_raw)[3:ncol(counts_raw)]
all_sample_cols

# valid_samples <- all_sample_cols[!all_sample_cols %in% c("TperMix", "TtriMix")]
valid_samples <- all_sample_cols

meta <- meta[meta$sample_id %in% valid_samples, ]
counts_mat <- as.matrix(counts_raw[, valid_samples])
rownames(counts_mat) <- counts_raw$gene_id
counts_mat <- counts_mat[, meta$sample_id]
head(counts_mat)

counts_mat <- round(counts_mat)                                   # to integer
# [MODIFIED-6] 过滤阈值改为 >=4，原为 >=6（原项目12样本的一半=6；本项目共8样本，一半=4）
keep <- rowSums(counts_mat >= 10) >= 4                            # 删除极低表达基因，至少在8/2 = 4个样本中至少有 10 个 reads才保留
counts_mat <- counts_mat[keep, ]


# [未修改/原项目遗留注释，与本项目(opossum)无关，仅保留作历史参考]
# gene_id	gene_name	ME13-1	ME13-2	ME13-3	ME13-4	CTRL-1	CTRL-2	CTRL-3	CTRL-4	baseMean	log2FoldChange	lfcSE	pvalue	padj	sig (padj<0.05 & |log2FC|>=0.585)
# ENSG00000278233.1	RNA5-8SN2	77.0	171.0	14.0	15.0	0.0	143.0	0.0	0.0	21.963509482370625	0.00668096140267046	0.35611598718401	0.851120448689936	0.9072712331508488	NS
# 本实验出现一个意外结果，有一个基因在7个样本中均没有 reads，5个样本数中还不少。 导致 log2FC 很成问题计算，因此把该基因删除
# 我的整体思想是最少一半的样本表达， 所以伤寒改成6


head(counts_mat)
cat("✅ 表达矩阵加载完成，过滤后保留基因数:", nrow(counts_mat), "\n")


# ================= 5. DESeq2 建模 =================
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = meta,
  design = ~Group
)

dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)


# ================= 5. PCA 分析 (🔻 字体已调整) =================
pca_data <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = Group, label = name)) +
  geom_point(size = 2.5, alpha = 0.9) + # 🔻 增大点大小和透明度以提高可见性
  geom_text_repel(
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = 20
  ) + # 🔻 使用 geom_text_repel 避免标签与点重合
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 12) + # 🔻 恢复全局字体大小以提高可读性
  scale_color_brewer(palette = "Set1") + # 🔻 使用对比度更高的颜色 palette
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )
ggsave(file.path(OUT_DIR, "PCA.pdf"), p_pca, width = 8, height = 6, dpi = 300)
cat("✅ PCA plot saved (optimized display)\n")


# ================= 7. 差异表达分析 (2 组对比) =================1
# ✅ 对比组定义：c("分组列名", "处理组/分子", "对照组/分母")
# log2FC = log2(处理组均值 / 对照组均值)
# 正数 = 在处理组(第一个)中上调；负数 = 在对照组(第二个)中上调

# [MODIFIED-7] contrasts 改为本项目唯一的一组对比 pi5 vs NC，原为 3组对比 (SMA4/SMC2/ME13 各 vs CTRL)
# c("NC", "pi5")
contrasts <- list(
  c("Group", "pi5", "NC")               # ✅ pi5 vs NC → log2FC>0 = pi5 上调
)

res_list <- list()
#sig_col_name <- "sig (padj<0.05 & |log2FC|>=1)" # log2(1.5) =0.585
sig_col_name <- "sig (padj<0.05 & |log2FC|>=0.585)" # log2(1.5) =0.585

for (comp in contrasts) {
  grp_treatment <- comp[2] # 处理组（分子）
  grp_control <- comp[3] # 对照组（分母）
  comp_name <- paste(grp_treatment, "vs", grp_control, sep = "_")
  cat(paste0(
    "\n🔍 正在分析: ",
    comp_name,
    " (处理:",
    grp_treatment,
    " vs 对照:",
    grp_control,
    ")\n"
  ))

  # 提取 DESeq2 结果并收缩
  res <- results(dds, contrast = comp)
  res <- lfcShrink(dds, contrast = comp, type = "ashr")

  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)

  # 提取当前对比两组的原始 reads
  samples_treatment <- meta$sample_id[meta$Group == grp_treatment]
  samples_control <- meta$sample_id[meta$Group == grp_control]
  raw_sub <- counts_raw[,
    c("gene_id", "gene_name", samples_treatment, samples_control),
    drop = FALSE
  ]

  # 合并结果与原始计数
  res_df <- left_join(res_df, raw_sub, by = "gene_id")

  # ✅ 严格按指定顺序排列列（处理组在前，对照组在后）
  final_cols <- c(
    "gene_id",
    "gene_name",
    samples_treatment,
    samples_control,
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "pvalue",
    "padj"
  )
  res_df <- res_df[, final_cols]

  # ✅ 按 padj 从小到大排序
  res_df <- res_df %>% arrange(padj)

  # ✅ 添加 sig 列（列名已含标准）
  res_df[[sig_col_name]] <- case_when(
    res_df$padj < 0.05 & res_df$log2FoldChange >= 1 ~ "Up",
    res_df$padj < 0.05 & res_df$log2FoldChange <= -1 ~ "Down",
    TRUE ~ "NS"
  )

  # 保存 CSV
  write_csv(res_df, file.path(OUT_DIR, paste0("DEG_", comp_name, ".csv")))
  res_list[[comp_name]] <- res_df

  # 🌋 火山图 (🔻 字体已调整，去除黑体)
  res_df$negLog10Padj <- -log10(res_df$padj)
  res_df$negLog10Padj[!is.finite(res_df$negLog10Padj)] <- NA

  # 选取 Top 10 显著基因用于标注
  top_labels <- res_df %>%
    filter(.data[[sig_col_name]] != "NS", !is.na(negLog10Padj)) %>%
    arrange(padj) %>%
    head(10) %>%
    mutate(
      label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)
    )

  p_vol <- ggplot(
    res_df,
    aes(x = log2FoldChange, y = negLog10Padj, color = .data[[sig_col_name]])
  ) +
    geom_point(alpha = 0.7, size = 0.5) + # 🔻 点大小：0.6 → 0.5
    scale_color_manual(
      values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey80"),
      labels = c(
        "Up" = "Upregulated",
        "Down" = "Downregulated",
        "NS" = "Not Significant"
      )
    ) +
    theme_bw(base_size = 10) + # 🔻 全局字体：12 → 10
    labs(
      title = paste(
        grp_treatment,
        "vs",
        grp_control,
        "(log2FC > 0 =",
        grp_treatment,
        "upregulated)"
      ),
      x = "log2 Fold Change",
      y = "-log10(adj. P-value)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    ) +
    # ✅ 去除黑体，使用普通字体
    geom_text_repel(
      data = top_labels,
      aes(label = label),
      size = 2, # 🔻 标签大小：2.5 → 2
      box.padding = 0.3,
      max.overlaps = 20,
      color = "black",
      fontface = "plain"
    ) # ✅ 去除 bold

  ggsave(
    file.path(OUT_DIR, paste0("Volcano_", comp_name, ".png")),
    p_vol,
    width = 8,
    height = 6,
    dpi = 300
  )
  cat(paste0(
    "✅ ",
    comp_name,
    " completed, significant DEGs: ",
    sum(res_df[[sig_col_name]] != "NS"),
    "\n"
  ))
}

# ================= 8. 热图 (🔻 字体已调整，去除黑体) =================
# Generate separate heatmap for each comparison using their own top 50 DEGs
heatmap_files <- c()

for (comp_name in names(res_list)) {
  cat(paste0("\n🔍 正在生成热图: ", comp_name, "\n"))
  
  # Get top 50 significant genes for this specific comparison
  top_genes_df <- res_list[[comp_name]] %>%
    filter(.data[[sig_col_name]] != "NS") %>%
    arrange(padj) %>%
    head(50)
  
  if (nrow(top_genes_df) == 0) {
    cat(paste0("⚠️  警告: ", comp_name, " 没有显著差异基因，跳过热图生成\n"))
    next
  }
  
  top_genes_ids <- top_genes_df$gene_id
  
  # Extract the specific treatment and control groups from comp_name
  # comp_name format is "treatment_vs_control"
  comp_parts <- strsplit(comp_name, "_vs_")[[1]]
  grp_treatment <- comp_parts[1]
  grp_control <- comp_parts[2]
  
  # Get only the samples for this specific comparison
  samples_treatment <- meta$sample_id[meta$Group == grp_treatment]
  samples_control <- meta$sample_id[meta$Group == grp_control]
  comp_samples <- c(samples_treatment, samples_control)
  
  # Extract VST normalized expression for these genes and only these samples
  mat <- assay(vsd)[top_genes_ids, comp_samples, drop = FALSE]
  mat <- t(scale(t(mat)))
  
  # ✅ 将行名替换为 gene_name，如果 gene_name 为空则使用 gene_id
  gene_names_for_plot <- top_genes_df$gene_name
  # 处理缺失或空的 gene_name
  gene_names_for_plot[
    is.na(gene_names_for_plot) | gene_names_for_plot == ""
  ] <- top_genes_df$gene_id[
    is.na(gene_names_for_plot) | gene_names_for_plot == ""
  ]
  # 确保行名唯一，如果有重复，可以添加后缀或保留原ID，这里简单处理：
  # 如果存在重复的 gene_name，pheatmap 可能会报错或显示不全。
  # 为了安全，我们检查重复并必要时回退到 ID
  if (any(duplicated(gene_names_for_plot))) {
    # 简单的去重策略：如果重复，追加 gene_id 以确保唯一性
    duplicated_idx <- duplicated(gene_names_for_plot)
    gene_names_for_plot[duplicated_idx] <- paste0(
      gene_names_for_plot[duplicated_idx],
      "_",
      top_genes_df$gene_id[duplicated_idx]
    )
  }
  
  rownames(mat) <- gene_names_for_plot
  
  # 创建正确的annotation数据框，只包含当前比较的样品
  annotation_df <- data.frame(
    Group = as.character(meta$Group[meta$sample_id %in% comp_samples]),
    row.names = comp_samples
  )
  
  # Generate heatmap filename
  heatmap_filename <- paste0("Heatmap_top50_", gsub("_vs_", "_", comp_name), ".pdf")
  heatmap_files <- c(heatmap_files, heatmap_filename)
  
  pheatmap(
    mat,
    annotation_col = annotation_df,
    filename = file.path(OUT_DIR, heatmap_filename),
    show_rownames = TRUE,
    main = paste("Top 50 DEGs:", comp_name),
    fontsize = 7,
    fontsize_row = 5,
    fontsize_col = 7,
    fontfamily = "sans",
    legend_labels = c("Low", "High")
  ) # Ensure legend labels are English if auto-generated
  
  cat(paste0("✅ ", comp_name, " 热图已生成: ", heatmap_filename, "\n"))
}

cat("✅ 所有热图生成完成\n")

# # ================= 9. 保存标准化计数矩阵 =================
# norm_counts <- counts(dds, normalized = TRUE)
# norm_counts_df <- as.data.frame(norm_counts)
# norm_counts_df$gene_id <- rownames(norm_counts_df)
# # 获取基因名称映射
# gene_names <- counts_raw[, c("gene_id", "gene_name")]
# norm_counts_df <- left_join(norm_counts_df, gene_names, by = "gene_id")
# # 重新排列列
# final_norm_cols <- c("gene_id", "gene_name", colnames(norm_counts))
# norm_counts_df <- norm_counts_df[, final_norm_cols]
# write_csv(norm_counts_df, file.path(OUT_DIR, "Normalized_Counts.csv"))
# cat("✅ 标准化计数矩阵已保存\n")


# ================= 10. 生成分析报告 =================
# Generate Bioinformatics_Analysis_Report.md with detailed content
report_file <- file.path(OUT_DIR, "Bioinformatics_Analysis_Report.md")


# Calculate summary statistics for the report
deg_summary <- lapply(res_list, function(df) {
  up <- sum(df[[sig_col_name]] == "Up", na.rm = TRUE)
  down <- sum(df[[sig_col_name]] == "Down", na.rm = TRUE)
  total_sig <- up + down
  list(up = up, down = down, total = total_sig)
})

# Build report content
report_content <- c(
  "# Bioinformatics Analysis Report",
  "",
  paste("Date:", Sys.Date()),
  # [MODIFIED-8] Project 名称改为本项目标识，原为 "2026_Item16_ZhenYan"
  paste("Project:", "1_opossum_YuFan (Didelphis virginiana, NC vs pi5)"),
  "",
  "## 1. Overview",
  "This report summarizes the differential expression analysis and quality control metrics for the RNA-seq dataset.",
  "- **Analysis Tool**: DESeq2",
  "- **Normalization**: VST (Variance Stabilizing Transformation) for PCA/Heatmap, Median-of-ratios for DE",
  "- **Significance Thresholds**: padj < 0.05, |log2FoldChange| >= 0.585",
  "",
  "## 2. Quality Control (QC)",
  ifelse(
    !is.null(QC_SRC),
    c(
      "- QC reports were generated using MultiQC.",
      "- Raw data quality and alignment metrics are available in the `QC/` directory.",
      "- Please refer to `QC/multiqc_report.html` for detailed interactive plots."
    ),
    c(
      "- ⚠️ QC folder not found or not copied.",
      "- Ensure raw data QC was performed prior to this step."
    )
  ),
  "",
  "## 3. Differential Expression Analysis Results",
  ""
)

# Add DEG statistics for each contrast
# heatmap ARE changed !, fold are log2(1.5)= 0.585
for (name in names(deg_summary)) {
  stats <- deg_summary[[name]]
  report_content <- c(
    report_content,
    paste0("### Contrast: ", name),
    paste0("- Total Significant Genes: ", stats$total),
    paste0("  - Upregulated: ", stats$up),
    paste0("  - Downregulated: ", stats$down),
    paste0("- Output File: `DEG_", name, ".csv`"),
    ""
  )
}

report_content <- c(
  report_content,
  "## 4. Visualizations",
  "",
  "### Principal Component Analysis (PCA)",
  "- **File**: `PCA.pdf`",
  "- **Description**: Shows sample clustering based on the top 500 most variable genes. Samples should cluster by biological group if the treatment effect is strong.",
  "",
  "### Volcano Plots",
  "- **Files**: `Volcano_*.png`",
  "- **Description**: Displays the relationship between statistical significance (-log10 padj) and magnitude of change (log2FC). Red points indicate upregulated genes, blue points indicate downregulated genes.",
  "",
  "### Heatmaps",
  "- **Files**: `Heatmap_top50_*.pdf`",
  "- **Description**: Hierarchical clustering of the top 50 differentially expressed genes for each contrast separately. Each heatmap shows expression patterns (Z-score normalized) for the most significant genes in that specific comparison across all samples.",
  "",
  "## 5. Generated Data Files",
  "",
  "| File Name | Description |",
  "| :--- | :--- |",
  "| `All_sample_gene_counts.tsv` | Raw count matrix for all samples. |",
  "| `All_sample_gene_tpm.tsv` | TPM (Transcripts Per Million) matrix, if available. |",
  "| `Normalized_Counts.csv` | DESeq2 normalized counts for downstream analysis. |",
  "| `DEG_*.csv` | Detailed differential expression results including log2FC, p-values, and base means. |",
  "| `PCA.pdf` | PCA plot showing sample relationships. |",
  "| `Volcano_*.png` | Volcano plots for each contrast. |",
  "| `Heatmap_top50_*.pdf` | Heatmaps of top 50 DEGs for each contrast separately. |",
  ifelse(
    !is.null(QC_SRC),
    "| `QC/` | Directory containing MultiQC and other QC reports. |",
    "| `QC/` | Not available. |"
  )
)

# report_file
# [1] "../Data_Analysis/DE_PCA_Results/Bioinformatics_Analysis_Report.md"

writeLines(report_content, con = report_file)
cat("✅ 详细分析报告已生成:", report_file, "\n")

cat("\n🎉 全部分析完成！结果已保存至:", OUT_DIR, "\n")
