#!/usr/bin/env Rscript
# ============================================================================
# Step 10 — 客户报告生成（英文，CLAUDE.md 规范）—— SCAFFOLD
#   汇总 Step 4/5/6/7/8/9 各样本输出 → English client report。
#   注意：Key Findings / Conclusions 的具体数值结论仍需在全部 6 样跑完后人工核定
#   （交付前自查：PCA 类比不适用于此项目，但拷贝数/合子型/脱靶的跨样本一致性
#   仍需人工核对，不能只信自动填表）。
#
#   运行: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript 10_generate_report.R
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
w("Per client confirmation (2026-07-07/08): (1) confirm the targeting vector is correctly integrated at the intended gene locus, replacing the wild-type allele with the human sequence; (2) verify the targeting-vector copy number; (3) confirm the complete integrity of the human knock-in sequence. Zygosity, off-target screening, and Neo-cassette status are natural extensions of these objectives and are reported alongside.\n")

w("## 2. Key Findings")
w("_[TODO: 3–5 条大白话结论，全部样本跑完并人工核对后填写——不要只信自动填表，尤其拷贝数/合子型的跨样本一致性]_\n")

w("## 3. Sample Information\n")
w("| Line | Sample(s) | Targeted locus (GRCm39) | Insert |")
w("| :--- | :---: | :---: | :---: |")
w("| RAGH | RAGH_153, RAGH_273 | Rag2 (chr2:101.45Mb) | human G-CSF/M-CSF/IL-6/IL-1b/IL-7/IL-15 cytokine cassette |")
w("| MTTH | MTTH_284/412/524 | Htt (chr5:34.9Mb) | full-genomic human HTT (~170kb, 67 exons) |")
w("| CD1A | CD1A_B125 | Cd1d1+Cd1d2 (chr3:86.9Mb) | **entire human CD1 gene cluster**: CD1D+CD1A+CD1C+CD1B+CD1E (~127kb, 5 genes) |")
w("")
w("_Note: the CD1A line's insert scope (whole CD1 gene cluster, confirmed by direct sequence alignment) is broader than the client's one-line description (\"CD1D1/D2 humanization\"); confirmed with client as an added detail, not a correction._\n")

w("## 4. Analysis Rationale and Decision Criteria")
w("These are **targeted (homologous-recombination) knock-ins** (mouse homology arms flank each human insert), so the intended integration site is known by design. Reads were aligned to a **hybrid reference** (GRCm39 + 3 construct contigs); each sample covers only its own construct, giving built-in cross-line negative controls. Because homology arms are sequence-identical to the endogenous locus, on-target validation uses arm-assisted bridging reads (no MAPQ filter); off-target screening and copy number use MAPQ>=20 unique reads restricted to the human-specific region only (whole-contig averaging underestimates copy number via arm dilution — validated in the RAGH_153 pilot run). Client-supplied wild-type (WT) allele sequences for all three lines enabled base-precise arm/insert boundaries via WT-vs-KI alignment, and a zygosity method based on coverage over the excised mouse span between the arms.\n")

w("## 5. Methods")
w("| Step | Tool | Key parameters |")
w("| :--- | :---: | :---: |")
w("| Alignment/QC/dedup/SV | nf-core/sarek 3.8.1 | bwa-mem2; skip BQSR; --tools tiddit (Manta dropped: pilot run showed >7h/sample with zero true junction calls after MAPQ filtering); hybrid ref |")
w("| Integration junctions | samtools (chimeric reads) | on-target arm-assisted bridging (no MAPQ) + off-target screen (MAPQ>=20, unique) + high-depth artifact blacklist (>5x baseline) |")
w("| Copy number | mosdepth / samtools coverage | human-specific-region depth / autosomal baseline, MAPQ>=20 |")
w("| KI integrity | mosdepth (windowed) | 500bp sliding-window depth uniformity across human insert |")
w("| Zygosity | minimap2 + samtools coverage | depth over mouse span excised between homology arms (arm coordinates re-mapped to GRCm39 at runtime) |")
w("| CD1A Neo status | samtools coverage | depth at NeoR/KanR cassette (coordinates parsed at runtime from client's construct file) |")
w("| Breakpoint annotation | GENCODE vM35 | GenomicRanges overlap |")
w("")

w("## 6. Results\n")
w("### 6.1 Copy number (per construct, per sample; human-specific region only)\n")
w("| Sample | Construct | Human region | Mean depth | Breadth % | Baseline | CN ratio | Interpretation |")
w("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
for (s in samples) {
  cnf <- file.path(PROJ, "analysis/copy_number", s, "copy_number.tsv")
  if (file.exists(cnf)) {
    d <- read.table(cnf, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    for (i in seq_len(nrow(d)))
      w("| ", s, " | ", d$construct[i], " | ", d$human_region[i], " | ", d$human_meandepth[i], " | ",
        d$human_breadth_pct[i], " | ", d$baseline[i], " | ", d$copy_number_ratio[i], " | ", d$interpretation[i], " |")
  }
}
w("")

w("### 6.2 Integration sites (on-target / off-target)\n")
w("_[TODO: 汇总 analysis/integration/<sample>/integration_summary.tsv]_\n")

w("### 6.3 KI sequence integrity\n")
w("_[TODO: 汇总 analysis/ki_integrity/<sample>/*.integrity_flags.tsv — 有异常窗口才列出，否则写\"no anomalies detected\"]_\n")

w("### 6.4 Zygosity\n")
w("_[TODO: 汇总 analysis/zygosity/<sample>/zygosity_summary.tsv]_\n")

w("### 6.5 CD1A Neo-cassette status\n")
w("_[TODO: analysis/cd1a_neo_status/CD1A_B125/verdict.txt —— 明确写出 intact 还是 deleted，不臆测]_\n")

w("## 7. Conclusions")
w("_[TODO: 每样本 on-target / 脱靶 / 拷贝数 / KI完整性 / 合子型 结论表；CD1A 额外含 Neo 状态]_\n")

w("## 8. Deliverable Files")
w("_[TODO: 交付文件清单，按 qc/ integration/ copy_number/ ki_integrity/ zygosity/ cd1a_neo_status/ 分文件夹]_\n")
w("---\n")
w("*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*")
close(con)
cat("报告 scaffold written:", REPORT, "\n")
cat("注意：叙述性小节(Key Findings/6.2-6.5/Conclusions/Deliverables)待全样本跑完+人工核对后定稿。\n")
