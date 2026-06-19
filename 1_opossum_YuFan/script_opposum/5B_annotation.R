#!/usr/bin/env Rscript
# ==========================================================
# [MODIFIED-ALL] 整个文件已重写。原脚本是植物(T_majus)模板：
#   - GTF_FILE 指向不存在的占位文件 "T_majus.final.gtf"
#   - 核心聚合逻辑被注释掉了，且最后一行 cat() 引用了从未生成的 client_df
#     变量，原脚本若直接运行必定报错
# 现改为真正适配本项目物种 Didelphis virginiana (opossum) 的版本。
#
# Liftoff GTF 转客户友好型基因注释表
# 物种: Didelphis virginiana (opossum)
# 输入: Didelphis_v.liftoff.gtf (来自 liftoff 同源迁移注释，无 "gene" feature，
#       只有 transcript/exon/CDS，需按 gene_id 从 transcript 行聚合出基因坐标)
# 输出: Didelphis_virginiana_Gene_Annotation_Client.csv
# ==========================================================

library(rtracklayer)
library(dplyr)
library(readr)

# [MODIFIED-1] 新增 setwd，原脚本没有设置工作目录
setwd("/home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/")

# [MODIFIED-2] GTF_FILE 改为本项目实际使用的 liftoff 注释（与 nf-core/rnaseq 跑流程
# 时 --gtf 参数完全一致，见 output_results/pipeline_info/params_*.json）
# 原值: "T_majus.final.gtf"（占位符，文件不存在）
GTF_FILE <- "/Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/Didelphis_v.liftoff.gtf"
# [MODIFIED-3] 输出目录改为与 4_run_DE_PCA_final.R 共用的 Data_Analysis 目录
OUT_DIR  <- "../Data_Analysis"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
# [MODIFIED-4] 输出文件名改为本物种命名，原为 "T_majus_Gene_Annotation_Client.csv"
OUT_CSV  <- file.path(OUT_DIR, "Didelphis_virginiana_Gene_Annotation_Client.csv")

cat("📖 正在读取 GTF 注释文件...\n")
gtf <- rtracklayer::import(GTF_FILE)
df  <- as.data.frame(gtf)

head(df)
table(df$type)

# [MODIFIED-5] 核心聚合逻辑完全重写。原脚本这部分代码全部被注释掉(filter(type=="gene")
# 那段)，而且这个GTF本身就没有 "gene" feature(只有transcript/exon/CDS)，按原逻辑
# 跑会得到0行。现改为按 transcript 行聚合到基因层级。
# 另外发现并修正了一个数据陷阱：liftoff 同源迁移注释里，约1.5%的 gene_id 会同时出现
# 在两个不同的scaffold上（如 KCTD1 同时在 HiC_scaffold_1 和 HiC_scaffold_22483），
# 若直接按 gene_id 取全局 min(start)/max(end)会把两个不同scaffold的坐标
# 错误合并成跨越上亿bp的"假基因"。这里先按 (gene_id, chromosome) 聚合出
# 每个locus，再为每个 gene_id 选转录本数最多的 locus 作为主要位置，
# 同时用 n_loci 标记该基因是否存在多个候选locus，方便后续核查。
loci_df <- df %>%
  filter(type == "transcript") %>%
  group_by(gene_id, chromosome = as.character(seqnames)) %>%
  summarise(
    gene_name     = first(gene_name),
    start         = min(start),
    end           = max(end),
    strand        = first(as.character(strand)),
    n_transcripts = n(),
    .groups = "drop"
  )

client_df <- loci_df %>%
  group_by(gene_id) %>%
  mutate(n_loci = n()) %>%
  slice_max(n_transcripts, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    # liftoff 注释里 gene_name 基本等于 gene_id；缺失时标记为预测基因
    gene_name = ifelse(is.na(gene_name) | gene_name == "",
                        paste0("Predicted_", gene_id), gene_name)
  ) %>%
  arrange(chromosome, start)

# 导出
write_csv(client_df, OUT_CSV)

cat("✅ 完成！客户友好型注释表已保存至:", OUT_CSV, "\n")
cat("📊 共提取", nrow(client_df), "个基因\n")
# [MODIFIED-6] 新增多locus统计提示，原脚本最后一行 cat() 引用了不存在的 client_df$多余统计，会直接报错
cat("⚠️  其中", sum(client_df$n_loci > 1), "个基因在GTF中有多个候选locus，已选转录本数最多的一个，n_loci列标记了候选数\n")
cat("💡 提示: 可直接用 Excel 打开，或作为差异分析表的基因注释底表\n")
