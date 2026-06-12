# ==========================================================
# 🌱 草图基因组 GTF 转客户友好型注释表 (一键运行)
# 输入: T_majus.final.gtf
# 输出: T_majus_Gene_Annotation_Client.csv
# ==========================================================

library(rtracklayer)
library(dplyr)
library(readr)

GTF_FILE <- "T_majus.final.gtf"
OUT_CSV  <- "T_majus_Gene_Annotation_Client.csv"

cat("📖 正在读取 GTF 注释文件...\n")
gtf <- rtracklayer::import(GTF_FILE)
df  <- as.data.frame(gtf)

head(df)

# # 仅保留基因级别注释
# gene_df <- df %>% filter(type == "gene")

# # 安全提取函数（兼容不同版本/质量的 GTF，缺失列自动补 NA 不报错）
# safe_col <- function(df, name) {
#   if (name %in% colnames(df)) df[[name]] else rep(NA_character_, nrow(df))
# }

# # 构建客户可读表格
# client_df <- tibble(
#   gene_id      = safe_col(gene_df, "gene_id"),
#   gene_name    = safe_col(gene_df, "gene_name"),
#   chromosome   = as.character(safe_col(gene_df, "seqnames")),
#   start        = safe_col(gene_df, "start"),
#   end          = safe_col(gene_df, "end"),
#   strand       = as.character(safe_col(gene_df, "strand")),
#   biotype      = safe_col(gene_df, "gene_biotype")
# ) %>%
#   mutate(
#     # 若无 gene_name 或与 ID 完全一致，自动标记为预测基因
#     gene_name = ifelse(is.na(gene_name) | gene_name == "" | gene_name == gene_id,
#                        paste0("Predicted_", gene_id), gene_name),
#     # 清理空白的生物类型
#     biotype   = ifelse(is.na(biotype) | biotype == "", "unknown", biotype)
#   ) %>%
#   arrange(chromosome, start)

# 导出
#write_csv(client_df, OUT_CSV)
write_csv(df, OUT_CSV)

cat("✅ 完成！客户友好型注释表已保存至:", OUT_CSV, "\n")
cat("📊 共提取", nrow(client_df), "个基因位点\n")
cat("💡 提示: 可直接用 Excel 打开，或作为差异分析表的基因注释底表\n")
