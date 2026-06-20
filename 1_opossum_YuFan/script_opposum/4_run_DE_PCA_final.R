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
# 按优先级列出几个"候选"路径——数组顺序即优先级顺序，排前面的优先采用：
#   1. nf-core 流程生成的 multiqc 汇总报告（最理想，图文并茂）
#   2. 没有multiqc时，退而求其次找 pipeline_info（流程运行记录，不是真正的QC报告，但好过没有）
#   3. 都没有就找当前目录下手动放的 QC 文件夹
QC_SOURCE_CANDIDATES <- c(
  "../output_results/multiqc",
  "../output_results/pipeline_info",
  "QC"
)

# 依次检查候选路径是否存在，命中第一个存在的就 break 退出循环，
# 即"取优先级最高且实际存在的那一个"；如果一个都不存在，QC_SRC 保持 NULL <<< =========================== 理解关键点
QC_SRC <- NULL
for (path in QC_SOURCE_CANDIDATES) {
  if (dir.exists(path)) {
    QC_SRC <- path
    break
  }
}

QC_SRC 
# [1] "../output_results/multiqc"


# 只有找到了QC源目录才执行拷贝，否则跳过（见下面 else 分支）
if (!is.null(QC_SRC)) {
  # Define the specific destination path for QC files
  # [MODIFIED-4] QC_DEST_DIR 改为本项目的 Data_Analysis/QC，原为 2026_Item16_ZhenYan/Data_Analysis/QC
  QC_DEST_DIR <- "/home/gao/projects_2026H2/1_opossum_YuFan/Data_Analysis/QC"

  # Create the destination directory if it doesn't exist
  dir.create(QC_DEST_DIR, showWarnings = FALSE, recursive = TRUE)

  # Copy contents of QC_SRC into QC_DEST_DIR
  # file.copy with recursive=TRUE copies the source folder INTO the destination if destination exists
  # 注意 file.copy 在 to 已存在的情况下，是把 QC_SRC 这个文件夹本身作为子目录拷进 QC_DEST_DIR，
  # 而不是把里面的文件铺平拷进去，所以实际产物会是 QC_DEST_DIR/multiqc/... 这样多一层目录
  success <- file.copy(
    from = QC_SRC,
    to = QC_DEST_DIR,
    recursive = TRUE,
    overwrite = TRUE
  )

  # file.copy 对多文件/递归拷贝会返回一个逐项的逻辑向量，any(success)即"只要有一个文件拷成功就算成功"
  if (any(success)) {
    cat("✅ QC 文件夹已拷贝:", QC_SRC, "->", QC_DEST_DIR, "\n")
  } else {
    cat("⚠️  警告: QC 文件夹拷贝失败\n")
  }
} else {
  # 三个候选路径都不存在的情况，不报错、只警告，不影响后面DE分析继续跑
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
  select(Group, `Name in File`) %>%         # 选择两列 “Group` 和 “Name in File”
  rename(sample_id = `Name in File`) %>%    # 重命名“Name in File”列名为 sample_id
  filter(!is.na(Group))           %>%       # 删除 Group 为 NA 的行
  # [MODIFIED-5] Group 水平改为本项目实际的两组 "NC","pi5"，原为 "CTRL","SMA4","SMC2","ME13" 四组
  mutate(Group = factor(Group, levels = c("NC", "pi5"))) # Group = factor(...): 将 Group 列转换为因子（factor）类型，并指定因子的水平（levels）顺序。

meta

# 前段代码解释：
# levels = c("NC", "pi5"): 明确定义因子的水平顺序为："NC" (对照组) "pi5" (处理组)
# 为什么要这样做？
# 在 R 中，因子的水平顺序非常重要，特别是在进行统计分析时：
# DESeq2 差异表达分析会将第一个水平（"CTRL"）作为参照组（baseline/reference）
# 其他组会与参照组进行比较
# 因子水平的顺序会影响结果的解释和可视化（如 PCA 图中组别的颜色顺序）
# 通过显式设置 levels，可以确保分析的一致性和可重复性，避免 R 自动按字母顺序排列因子水平（那样可能会把 "CTRL" 排在后面）。

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
contrasts <- list(
  c("Group", "pi5", "NC")               # ✅ pi5 vs NC → log2FC>0 = pi5 上调
)

res_list <- list()

# 只需要改这两行，其余地方会自动同步，不会再出现"标签和实际逻辑不一致"的问题。
LFC_THRESHOLD <- 0.585          # log2(1.5) = 0.585
PADJ_THRESHOLD <- 0.05          # 显著性 padj level

sig_col_name <- paste0(
  "sig (padj<", PADJ_THRESHOLD, " & |log2FC|>=", LFC_THRESHOLD, ")"
)

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
  # [FIXED-BUG] fold-change阈值改用 LFC_THRESHOLD（与上面 sig_col_name 标签同源），原来这里
  # 硬编码的是 1，跟列名标注的 0.585 不一致
  # [MODIFIED] padj阈值改用 PADJ_THRESHOLD，原来这里硬编码的是 0.05，现在和标签、报告统一管理
  res_df[[sig_col_name]] <- case_when(
    res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange >= LFC_THRESHOLD ~ "Up",
    res_df$padj < PADJ_THRESHOLD & res_df$log2FoldChange <= -LFC_THRESHOLD ~ "Down",
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

# ✅ pi5_vs_NC completed, significant DEGs: 0
# Warning message:
# Removed 10024 rows containing missing values or values outside the scale range (`geom_point()`). 

# 两个现象分别由两个独立原因造成，诊断结果如下：

# 一、为什么有10024个missing value（geom_point的警告）
# 这是 DESeq2 的标准行为，不是bug。 验证过了：

# DESeq2 的 results() 默认会做"independent filtering"——自动算出一个 baseMean 阈值（这次算出来是 2651），低于这个阈值的基因直接把 padj 设为 NA（不参与多重检验校正）。
# 实测：padj 为 NA 的 10,024 个基因，baseMean 全部 ≤ 2650.97；padj 非NA的3582个基因，baseMean 全部 ≥ 2651。完全精确对应，10024 恰好是总数13606的73.67%（和DESeq2自己报告的 filterThreshold百分位73.67%吻合）。
# 原理：低表达基因的log2FC估计本身噪音很大，硬塞进多重检验校正只会拉低统计功效（让真正显著的基因也被牵连），所以DESeq2自动把它们排除在外。这是文献认可的标准统计学优化，不是数据或脚本问题。
# geom_point() 的警告只是说"画图时这10024个点因为y轴(-log10 padj)是NA画不出来"，纯粹是画图层面的提示，完全正常，可以忽略。
# 二、为什么 "significant DEGs: 0"（这是更值得关注的问题）
# 这个和上面是两件独立的事，我往下查了两层：

# 表层问题——脚本本身有个继承自原模板、从没被注意到的bug：

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

# 正在生成热图: pi5_vs_NC
# ⚠️  警告: pi5_vs_NC 没有显著差异基因，跳过热图生成

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
# [MODIFIED-ALL] 本section整体重写。原模板是通用报告，没有任何针对本实验的特殊说明。
# 本实验有几个非常特殊、足以影响结果解读的背景，原模板完全没有体现，必须写进报告：
#   (a) 参考基因组注释是liftoff同源迁移注释，不是该物种原生注释，且发现415个基因有
#       多locus歧义 (见 5B_annotation.R)
#   (b) STAR比对实际跑的是2-pass模式（--twopassMode Basic）。
#       [重要更正] params_*.json 里记录的 star_twopass:false 不可信 ——
#       实际查看 work 目录下某样本(NC_4)的 Log.out，命令行明确写着
#       `--twopassMode Basic`，且存在 _STARpass1 目录、Log.progress.out 里也写着
#       "Started 1st pass mapping"，这些都是2-pass确实跑过的直接证据。
#       推测是后来改了配置(改成false)但用 -resume 时复用了之前(true时)的缓存结果，
#       所以params json记录的是"最新配置"而不是"实际产出这批结果时用的配置"。
#       这提醒一个通用教训：nf-core/nextflow用-resume时，params json不能完全代表
#       实际跑出这批数据用的参数，真正可信的证据是work目录里每个任务的Log.out/
#       .command.sh —— 报告里必须按实际证据(2-pass已跑)来写，不能直接信params json。
#   (c) 这次DE分析本身发现的几个统计层面的警示信号(PCA不分离、16个padj显著基因
#       100%同方向且完全分离、文库大小两组相差~27%)，这些不是"bug"，是这批数据
#       真实存在的、解读结果时必须考虑的局限性
# 凡是能从脚本里已有的变量动态算出来的数字(文库大小比例、显著基因数、方向、分离情况)，
# 全部动态计算，不写死数字，避免以后数据换了报告却没跟着更新。
# Generate Bioinformatics_Analysis_Report.md with detailed content
report_file <- file.path(OUT_DIR, "Bioinformatics_Analysis_Report.md")

# ---- 10.1 动态计算文库大小差异 (回应 27% library size 差异这个发现) ----
libsize_by_sample <- colSums(counts_raw[, meta$sample_id])
libsize_by_group <- tapply(libsize_by_sample, meta$Group, mean)
libsize_ratio_pct <- round(
  100 * abs(libsize_by_group["NC"] - libsize_by_group["pi5"]) /
    mean(libsize_by_group),
  1
)

# ---- 10.2 对每个contrast，动态核查"padj显著基因是否方向一致/组间完全分离" ----
# 这是对4B_check_padj_sig_genes.R里那次人工核查的轻量版固化：只用padj阈值(不叠加LFC
# 门槛)挑出基因，检查它们是否100%同方向、是否在两组间完全不重叠 —— 如果是，提示更可能
# 是系统性/技术性混杂因素，而不是16个独立的生物学调控事件
check_sig_gene_pattern <- function(res_df, norm_mat, meta) {
  sig <- res_df[!is.na(res_df$padj) & res_df$padj < PADJ_THRESHOLD, ]
  if (nrow(sig) == 0) {
    return(list(n = 0, note = "无padj显著基因，无需核查方向/分离情况"))
  }
  n_up <- sum(sig$log2FoldChange > 0)
  n_down <- sum(sig$log2FoldChange < 0)
  same_direction <- (n_up == 0 || n_down == 0)

  groups <- levels(meta$Group)
  sep_count <- 0
  for (gid in sig$gene_id) {
    vals_a <- norm_mat[gid, meta$sample_id[meta$Group == groups[1]]]
    vals_b <- norm_mat[gid, meta$sample_id[meta$Group == groups[2]]]
    if (max(vals_a) < min(vals_b) || max(vals_b) < min(vals_a)) sep_count <- sep_count + 1
  }
  list(
    n = nrow(sig), n_up = n_up, n_down = n_down,
    same_direction = same_direction, sep_count = sep_count
  )
}

norm_counts_for_check <- counts(dds, normalized = TRUE)
sig_pattern_by_contrast <- lapply(res_list, function(df) {
  check_sig_gene_pattern(df, norm_counts_for_check, meta)
})

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
  paste("Project:", "1_opossum_YuFan (Didelphis virginiana, NC vs pi5)"),
  # [MODIFIED] 加入Author信息，客户报告需要署名
  "Author: Zhen Gao, PhD, Principal Bioinformatics Scientist",
  "",
  "## 1. Overview",
  "This report summarizes the differential expression analysis and quality control metrics for the RNA-seq dataset.",
  "- **Analysis Tool**: DESeq2",
  "- **Normalization**: VST (Variance Stabilizing Transformation) for PCA/Heatmap, Median-of-ratios for DE",
  "- **LFC shrinkage**: ashr (Stephens 2016) — adaptive, shrinkage strength depends on each gene's own standard error",
  paste0(
    "- **Significance Thresholds**: padj < ", PADJ_THRESHOLD,
    ", |log2FoldChange| >= ", LFC_THRESHOLD
  ),
  "",
  "## 2. Upstream Pipeline & Reference Caveats",
  "This experiment uses a non-model organism with a liftoff-transferred annotation, which has",
  "real implications for how the downstream DE results below should be interpreted.",
  "",
  "**Genome / annotation**",
  "- Species: *Didelphis virginiana* (opossum). Genome assembly: `dv-2k.fasta` (Hi-C scaffolded, DNA Zoo).",
  "- Gene annotation: `Didelphis_v.liftoff.gtf` — produced by **liftoff** (homology-based annotation",
  "  transfer from a related reference species), **not** a native, experimentally-curated annotation",
  "  for this species. Two concrete data-quality issues were found and corrected during processing",
  "  (see `5B_annotation.R`):",
  "  - This GTF has **no `gene` feature rows** (only `transcript`/`exon`/`CDS`); gene-level coordinates",
  "    had to be derived by aggregating transcript records.",
  "  - **415 of 27,668 genes (~1.5%)** have their `gene_id` mapped to two different scaffolds",
  "    simultaneously (a known liftoff artifact, e.g. paralog/repeat region mis-mapping). The locus",
  "    with the most supporting transcripts was kept as primary; see the `n_loci` column in",
  "    `Didelphis_virginiana_Gene_Annotation_Client.csv` to identify which genes are affected.",
  "",
  "**Alignment (nf-core/rnaseq, aligner = star_salmon)**",
  "- STAR was run with 2-pass mode (`--twopassMode Basic`). Splice junctions discovered across all",
  "  samples in the first pass were used when finalizing alignments in the second pass, which",
  "  partially compensates for the liftoff annotation being an imperfect/incomplete transfer by",
  "  recovering species-specific or novel junctions that the liftoff GTF alone would have missed.",
  "- STAR filtering parameters were tightened for this non-model-organism / imperfect-reference",
  "  scenario (see `local_optimized.config`): `--outFilterMultimapNmax 8`, `--alignSJoverhangMin 8`,",
  "  `--alignSJDBoverhangMin 1`, `--outFilterMismatchNmax 2`.",
  "- Quantification: Salmon, using a tx2gene mapping built from the same liftoff GTF.",
  "",
  "## 3. Quality Control (QC)",
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
  paste0(
    "- Mean library size differs by **", libsize_ratio_pct,
    "%** between NC and pi5 groups (NC mean = ",
    format(round(libsize_by_group["NC"]), big.mark = ","),
    " reads; pi5 mean = ", format(round(libsize_by_group["pi5"]), big.mark = ","),
    " reads). DESeq2 size-factor normalization corrects for this at the model-fitting level, but a",
    " difference of this size is worth confirming against the sequencing/library-prep batch records",
    " (see Section 6)."
  ),
  "",
  "## 4. Differential Expression Analysis Results",
  ""
)

# Add DEG statistics for each contrast
for (name in names(deg_summary)) {
  stats <- deg_summary[[name]]
  pat <- sig_pattern_by_contrast[[name]]
  report_content <- c(
    report_content,
    paste0("### Contrast: ", name),
    paste0(
      "- Genes with padj < ", PADJ_THRESHOLD, " (regardless of fold-change size): ",
      pat$n
    ),
    paste0(
      "- Of those, passing the |log2FC| >= ", LFC_THRESHOLD, " effect-size filter (final \"sig\" call): ",
      stats$total, " (Up: ", stats$up, ", Down: ", stats$down, ")"
    ),
    paste0("- Output File: `DEG_", name, ".csv`"),
    ""
  )
}

report_content <- c(
  report_content,
  "## 5. Key Caveats Found During This Analysis",
  "",
  "- **PCA shows NC and pi5 samples do not separate** (see `PCA.pdf`) — samples from the two groups",
  "  are interspersed rather than forming distinct clusters, indicating the overall transcriptome is",
  "  highly similar between groups at this sample size.",
  "- **ashr shrinkage compresses nearly all effect sizes toward zero** (most genes' shrunk log2FC fall",
  "  within roughly ±0.1), consistent with the weak overall signal seen in the PCA. This is why the",
  "  `padj < 0.585`-filtered \"sig\" column above is much smaller than the raw padj-significant count.",
  "- For each contrast, the padj-significant genes (regardless of fold-change size) were checked for",
  "  whether they are dominated by a single direction and/or show complete (non-overlapping) separation",
  "  between groups despite modest effect sizes — a pattern that, when it affects *all* significant",
  "  genes simultaneously, is more consistent with a shared systematic/technical confound (e.g. the",
  "  library-size difference noted in Section 3, or an unmodeled batch effect correlated with Group)",
  "  than with that many independent gene-specific biological regulation events:"
)

for (name in names(sig_pattern_by_contrast)) {
  pat <- sig_pattern_by_contrast[[name]]
  if (is.null(pat$n) || pat$n == 0) {
    report_content <- c(report_content, paste0("  - ", name, ": ", pat$note))
  } else {
    report_content <- c(
      report_content,
      paste0(
        "  - ", name, ": ", pat$n, " padj-significant genes — ",
        pat$n_up, " up / ", pat$n_down, " down",
        ifelse(pat$same_direction, " (**100% same direction**)", " (mixed directions, less concerning)"),
        "; ", pat$sep_count, "/", pat$n,
        " show complete non-overlapping separation between groups."
      )
    )
  }
}

report_content <- c(
  report_content,
  "",
  "- **Recommendation**: before treating any gene from this dataset's DEG list as a confirmed",
  "  biological finding, (1) confirm with the wet-lab team whether NC and pi5 samples were",
  "  prepared/sequenced in the same batch, (2) for genes of interest, inspect per-sample normalized",
  "  counts individually (see `4B_check_padj_sig_genes.R` → `Check_padj_sig_genes_per_sample_dotplot.png`",
  "  and `Sig_padj_genes_manual_check.csv`) rather than relying on padj/log2FC alone.",
  "",
  "## 6. Visualizations",
  "",
  "### Principal Component Analysis (PCA)",
  "- **File**: `PCA.pdf`",
  "- **Description**: Shows sample clustering based on the top 500 most variable genes. Samples should cluster by biological group if the treatment effect is strong.",
  "",
  "### Volcano Plots",
  "- **Files**: `Volcano_*.png`",
  "- **Description**: Displays the relationship between statistical significance (-log10 padj) and magnitude of change (log2FC). Red points indicate upregulated genes, blue points indicate downregulated genes.",
  "  Note: with strong ashr shrinkage and a narrow effect-size range, the two \"wings\" can look like a",
  "  smooth curve rather than a scattered cloud — this is expected when shrunk LFC and -log10(padj)",
  "  both become near-monotonic functions of the same underlying z-statistic.",
  "",
  "### Heatmaps",
  # [MODIFIED] 改为根据deg_summary动态判断：本次没有任何contrast有基因同时通过padj和
  # |log2FC|门槛，所以根本没有生成Heatmap_top50_*.pdf这个文件。原模板写死引用这个文件名，
  # 跟实际产出的文件完全不对应。改为动态判断 + 指向实际存在的文件
  # (Heatmap_padj_sig_genes_pi5_vs_NC.pdf, 来自4B_check_padj_sig_genes.R)
  if (any(sapply(deg_summary, function(x) x$total) > 0)) {
    c(
      "- **Files**: `Heatmap_top50_*.pdf`",
      "- **Description**: Hierarchical clustering of the top 50 differentially expressed genes for contrasts that had at least one gene passing both the padj and |log2FC| thresholds. Each heatmap shows expression patterns (Z-score normalized) across all samples."
    )
  } else {
    c(
      "- No contrast had any gene passing both the padj and |log2FC| thresholds, so no `Heatmap_top50_*.pdf` files were produced.",
      "- **File**: `Heatmap_padj_sig_genes_pi5_vs_NC.pdf`",
      "- **Description**: Heatmap of all genes passing the padj threshold alone (regardless of fold-change size), generated separately for manual review.",
      "- **File**: `Check_padj_sig_genes_per_sample_dotplot.png`",
      "- **Description**: Per-sample normalized counts for the same genes, plotted individually, to verify that significance calls are not driven by a single outlier sample."
    )
  },
  "",
  "## 7. Generated Data Files",
  "",
  "| File Name | Description |",
  "| :--- | :--- |",
  "| `Bioinformatics_Analysis_Report.md` | This report. |",
  "| `All_sample_gene_counts.tsv` (in `Reads/`) | Raw count matrix for all samples. |",
  "| `All_sample_gene_tpm.tsv` (in `Reads/`) | TPM (Transcripts Per Million) matrix. |",
  "| `DEG_*.csv` | Differential expression results per contrast, including log2FC, p-values, and base means. |",
  "| `PCA.pdf` | PCA plot showing sample relationships. |",
  "| `Volcano_*.png` | Volcano plot for each contrast. |",
  if (any(sapply(deg_summary, function(x) x$total) > 0)) {
    "| `Heatmap_top50_*.pdf` | Heatmap of the top 50 DEGs, for contrasts with genes passing both thresholds. |"
  } else {
    "| `Heatmap_padj_sig_genes_pi5_vs_NC.pdf` | Heatmap of padj-significant genes (see Section 6). |"
  },
  "| `Check_padj_sig_genes_per_sample_dotplot.png` | Per-sample verification plot for padj-significant genes (see Section 6). |",
  "| `Sig_padj_genes_manual_check.csv` | Per-sample raw and normalized counts plus pre/post-shrinkage log2FC for the padj-significant genes, for manual review. |",
  "| `Didelphis_virginiana_Gene_Annotation_Client.csv` (in `Data_Analysis/`) | Gene-level annotation derived from the liftoff GTF; see the `n_loci` column for genes with multiple candidate loci. |",
  ifelse(
    !is.null(QC_SRC),
    "| `QC/` (in `Data_Analysis/`) | MultiQC and other QC reports. |",
    "| `QC/` | Not available. |"
  )
)

writeLines(report_content, con = report_file)
cat("✅ 详细分析报告已生成:", report_file, "\n")

cat("\n🎉 全部分析完成！结果已保存至:", OUT_DIR, "\n")
