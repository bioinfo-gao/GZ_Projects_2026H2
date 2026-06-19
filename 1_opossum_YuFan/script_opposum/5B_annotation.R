#!/usr/bin/env Rscript
# ==========================================================
# Liftoff GTF 转客户友好型基因注释表
# 物种: Didelphis virginiana (opossum)
# 输入: Didelphis_v.liftoff.gtf (来自 liftoff 同源迁移注释，无 "gene" feature，
#       只有 transcript/exon/CDS，需按 gene_id 从 transcript 行聚合出基因坐标)
# 输出: Didelphis_virginiana_Gene_Annotation_Client.csv
# ==========================================================

library(rtracklayer)
library(dplyr)
library(readr)

setwd("/home/gao/projects_2026H2/1_opossum_YuFan/script_opposum/")

GTF_FILE <- "/Work_bio/references/Didelphis_virginiana/mDidVir1/DNA_Zoo/Didelphis_v.liftoff.gtf"
OUT_DIR  <- "../Data_Analysis"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_CSV  <- file.path(OUT_DIR, "Didelphis_virginiana_Gene_Annotation_Client.csv")

cat("📖 正在读取 GTF 注释文件...\n")
gtf <- rtracklayer::import(GTF_FILE)
df  <- as.data.frame(gtf)

head(df)
table(df$type)

# 该 GTF 没有 "gene" feature，按 transcript 行聚合到基因层级
# （同一 gene_id 取最小 start / 最大 end，覆盖该基因全部转录本范围）
client_df <- df %>%
  filter(type == "transcript") %>%
  group_by(gene_id) %>%
  summarise(
    gene_name   = first(gene_name),
    chromosome  = first(as.character(seqnames)),
    start       = min(start),
    end         = max(end),
    strand      = first(as.character(strand)),
    n_transcripts = n()
  ) %>%
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
cat("💡 提示: 可直接用 Excel 打开，或作为差异分析表的基因注释底表\n")
