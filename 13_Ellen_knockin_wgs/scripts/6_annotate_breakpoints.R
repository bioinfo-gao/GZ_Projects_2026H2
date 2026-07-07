#!/usr/bin/env Rscript
# ============================================================================
# Step 6 — 整合断点注释到小鼠基因（GENCODE vM35）
#   读取 Step 4 的候选整合位点 TSV，用 GRCm39 基因模型注释每个断点落在
#   哪个小鼠基因（外显子/内含子/基因间），判断是否破坏内源基因。
#
#   运行: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript 6_annotate_breakpoints.R RAGH_153
#   依赖: GenomicRanges, rtracklayer（DE_R45 环境）
# ============================================================================
suppressMessages({ library(GenomicRanges); library(rtracklayer) })

args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else stop("用法: Rscript 6_annotate_breakpoints.R <sample>")
PROJ   <- "/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
GTF    <- "/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/gencode.vM35.annotation.gtf"
INDIR  <- file.path(PROJ, "analysis/integration", SAMPLE)
OUTDIR <- file.path(PROJ, "analysis/annotation", SAMPLE); dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# 载入基因模型（gene + exon）
message("载入 GTF ...")
gtf   <- import(GTF)
genes <- gtf[gtf$type == "gene"]
exons <- gtf[gtf$type == "exon"]

# 汇总所有构建体的候选整合位点
sites_files <- list.files(INDIR, pattern = "candidate_integration_sites\\.tsv$", full.names = TRUE)
if (length(sites_files) == 0) stop("未找到候选整合位点文件，先跑 4_integration_analysis.sh")

all_out <- list()
for (f in sites_files) {
  tg <- sub("\\.candidate_integration_sites\\.tsv$", "", basename(f))
  df <- tryCatch(read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) next
  gr <- GRanges(df$chrom, IRanges(df$start, df$end))

  hit_gene <- findOverlaps(gr, genes)
  hit_exon <- findOverlaps(gr, exons)
  gene_name <- rep("intergenic", length(gr))
  gene_name[queryHits(hit_gene)] <- genes$gene_name[subjectHits(hit_gene)]
  in_exon <- rep("intron/intergenic", length(gr))
  in_exon[queryHits(hit_exon)] <- "EXON"

  # 最近基因（当落在基因间时）
  nearest_idx <- nearest(gr, genes)
  nearest_gene <- genes$gene_name[nearest_idx]
  dist_near <- distance(gr, genes[nearest_idx])

  out <- data.frame(construct = tg, df,
                    overlap_gene = gene_name, feature = in_exon,
                    nearest_gene = nearest_gene, dist_to_nearest = dist_near,
                    stringsAsFactors = FALSE)
  all_out[[tg]] <- out
}

res <- do.call(rbind, all_out)
outfile <- file.path(OUTDIR, paste0(SAMPLE, "_integration_annotated.tsv"))
write.table(res, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
message("写入: ", outfile)
print(res)

# 定点打靶预期位点提示
cat("\n预期靶位点（定点打靶）：RAGH→Rag2(chr2:101.45Mb) | MTTH→Htt(chr5:34.9Mb)\n")
cat("落在预期基因上的位点=on-target 验证；落在别处的高支持位点=脱靶随机整合，需重点核查。\n")
