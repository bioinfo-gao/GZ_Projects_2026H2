#!/usr/bin/env Rscript
# ============================================================================
# Step 7 — 客户报告生成（英文，CLAUDE.md 规范）—— SCAFFOLD
#   汇总 Step 4/5/6 各样本输出 → English client report。
#   注意：报告的 Objectives / Key Findings / Conclusions 依赖客户确认的分析目标，
#   最终版在 (a) 客户回复分析目标、(b) 全部 6 样(含 CD1A)跑完 后定稿。
#   本脚本先把"结果表"部分自动填好，叙述部分留占位待补。
#
#   运行: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript 7_generate_report.R
# ============================================================================
suppressMessages(library(tools))
PROJ <- "/home/gao/projects_2026H2/13_Ellen_knockin_wgs"
TODAY <- format(Sys.Date(), "%m%d")
DELIV <- file.path(PROJ, paste0("custom_research_report_", format(Sys.Date(), "%Y%m%d")))
dir.create(DELIV, showWarnings = FALSE, recursive = TRUE)
REPORT <- file.path(DELIV, paste0("Ellen_KnockIn_WGS_", TODAY, ".md"))

samples <- list.dirs(file.path(PROJ, "analysis/copy_number"), recursive = FALSE, full.names = FALSE)
con <- file(REPORT, "w")
w <- function(...) cat(..., "\n", file = con, sep = "")

w("# Ellen Humanized-Mouse WGS — Integration, Copy Number & Structural Analysis\n")
w("**Report Date:** ", format(Sys.Date(), "%Y-%m-%d"))
w("**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics")
w("**Analysis Platform:** Linux HPC server\n")

w("## 1. Objectives")
w("_[TODO: 依据客户确认的分析目标填写 — 确认定点整合 / 拷贝数 / 脱靶筛查 / 合子型 / Neo 删除验证]_\n")

w("## 2. Key Findings")
w("_[TODO: 3–5 条大白话结论，全部样本跑完后填写]_\n")

w("## 3. Sample Information\n")
w("| Line | Sample | Targeted locus (GRCm39) | Insert |")
w("| :--- | :---: | :---: | :---: |")
w("| RAGH | RAGH_153, RAGH_273 | Rag2 (chr2:101.45Mb) | human G-CSF/M-CSF/IL-6/IL-1b/IL-7/IL-15 |")
w("| MTTH | MTTH_284/412/524 | Htt (chr5:34.9Mb) | human HTT full-genomic (~170kb) |")
w("| CD1A | CD1A_B125 | TBD | TBD (sequence pending) |")
w("")

w("## 4. Analysis Rationale and Decision Criteria")
w("These are **targeted (homologous-recombination) knock-ins**, so the intended integration site is known by design; the analysis validates on-target integration and screens for off-target events, copy number, zygosity, and Neo removal. Reads were aligned to a **hybrid reference** (GRCm39 + construct contigs); each sample covers only its own construct, giving built-in cross-line negative controls.\n")

w("## 5. Methods")
w("| Step | Tool | Key parameters |")
w("| :--- | :---: | :---: |")
w("| Alignment/QC/dedup/SV | nf-core/sarek 3.8.1 | bwa-mem2; skip BQSR; --tools manta,tiddit; hybrid ref |")
w("| Integration junctions | samtools + Manta/DELLY | discordant + split reads; BND on TG_ contigs |")
w("| Copy number | mosdepth / samtools coverage | construct depth / autosomal single-copy baseline, MAPQ>=20 |")
w("| Breakpoint annotation | GENCODE vM35 | GenomicRanges overlap |")
w("")

w("## 6. Results\n")
w("### 6.1 Copy number (per construct, per sample)\n")
w("| Sample | Construct | Mean depth | Breadth % | Baseline | CN ratio |")
w("| :--- | :---: | :---: | :---: | :---: | :---: |")
for (s in samples) {
  cnf <- file.path(PROJ, "analysis/copy_number", s, "copy_number.tsv")
  if (file.exists(cnf)) {
    d <- read.table(cnf, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    for (i in seq_len(nrow(d)))
      w("| ", s, " | ", d$construct[i], " | ", d$mean_depth[i], " | ",
        d$breadth_pct[i], " | ", d$baseline[i], " | ", d$copy_number_ratio[i], " |")
  }
}
w("")
w("### 6.2 Integration sites\n")
w("_[TODO: 汇总 analysis/integration/<sample>/*.candidate_integration_sites.tsv 与 annotation 结果]_\n")

w("## 7. Conclusions")
w("_[TODO: 每样本 on-target 验证 / 脱靶 / 拷贝数 / 合子型 / Neo 结论表]_\n")

w("## 8. Deliverable Files")
w("_[TODO: 交付文件清单]_\n")
w("---\n")
w("*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*")
close(con)
cat("报告 scaffold written:", REPORT, "\n")
cat("注意：叙述性小节(Objectives/Key Findings/Conclusions)待客户确认目标 + 全样本完成后定稿。\n")
