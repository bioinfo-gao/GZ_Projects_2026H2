#!/usr/bin/env Rscript
# P17 — Phase 1 下游分析与出图（assembly-free taxonomy + diversity）
# 输入：output_results/  (MetaPhlAn combined + Bracken combined)
# 输出：analysis/figures/ (PDF+PNG 300dpi), analysis/tables/ (tsv), analysis/rds_cache/
# 图风格：colorblind-safe 固定顺序调色板、thin marks、direct labels、clean theme。
# env: regular_bioinfo (全局 R_LIBS 统一指向该 env 的 R library；vegan/ggrepel/patchwork/ape 装于此)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(vegan); library(ggrepel); library(patchwork); library(tibble); library(stringr)
})
set.seed(1)

PROJ <- "/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics"
OUT  <- file.path(PROJ, "output_results")
FIG  <- file.path(PROJ, "analysis", "figures"); TAB <- file.path(PROJ, "analysis", "tables")
RDS  <- file.path(PROJ, "analysis", "rds_cache")
for (d in c(FIG, TAB, RDS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---------- aesthetics ----------
# 组别配色：Okabe-Ito CVD-safe pair
GRP_COL <- c(AL = "#0072B2", IF = "#D55E00")
# 分类堆叠 Top-N 配色（12 色，尽量可区分）+ Other 灰
TAXA_PAL <- c("#4E79A7","#F28E2B","#59A14F","#E15759","#B07AA1","#76B7B2",
              "#EDC948","#FF9DA7","#9C755F","#BAB0AC","#86BCB6","#D37295")
OTHER_COL <- "#CCCCCC"

theme_pub <- function(base=12){
  theme_minimal(base_size=base, base_family="sans") +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          axis.line=element_line(color="grey30", linewidth=0.3),
          axis.ticks=element_line(color="grey30", linewidth=0.3),
          plot.title=element_text(face="bold", size=base+2),
          plot.subtitle=element_text(color="grey30", size=base-1),
          legend.key.size=unit(11,"pt"),
          strip.text=element_text(face="bold"))
}
save_fig <- function(p, name, w=8, h=5.5){
  ggsave(file.path(FIG, paste0(name,".pdf")), p, width=w, height=h, device=cairo_pdf)
  ggsave(file.path(FIG, paste0(name,".png")), p, width=w, height=h, dpi=300)
  message("saved: ", name)
}
clean_id <- function(x){  # "HFD_AL_6_05_12_HFD_AL_6_05_12_L4_..." -> "HFD_AL_6_05_12"
  m <- str_match(x, "HFD_(AL|IF)_[0-9]+_[0-9]+_[0-9]+")[,1]
  ifelse(is.na(m), x, m)
}
grp_of <- function(id) ifelse(grepl("_AL_", id), "AL", "IF")

## ---------- load MetaPhlAn combined (relative abundance %, with lineage names) ----------
mp_file <- Sys.glob(file.path(OUT, "metaphlan", "*combined_reports.txt"))[1]
stopifnot(!is.na(mp_file))
# 文件首行是 "#mpa_vJan25..." 注释，真正表头在 "clade_name" 行 —— 定位后再读
mp_lines <- readLines(mp_file)
hdr <- grep("^clade_name", mp_lines)[1]
mp <- read.delim(text=mp_lines[hdr:length(mp_lines)], check.names=FALSE)
names(mp)[1] <- "clade_name"
# 样品列改名
samp_cols <- setdiff(names(mp), c("clade_name","NCBI_tax_id","clade_taxid"))
clean_map <- setNames(clean_id(samp_cols), samp_cols)
mp2 <- mp %>% select(clade_name, all_of(samp_cols))
names(mp2)[-1] <- clean_map[names(mp2)[-1]]
long_level <- function(df, tag){  # tag e.g. "g__" genus, "s__" species (deepest rank)
  df %>%
    filter(str_detect(clade_name, fixed(tag)),
           !str_detect(clade_name, fixed(ifelse(tag=="g__","s__","t__")))) %>%
    mutate(taxon = str_extract(clade_name, paste0(tag, "[^|]+")) %>% str_remove(tag)) %>%
    select(-clade_name) %>% group_by(taxon) %>% summarise(across(everything(), sum), .groups="drop")
}
mp_g <- long_level(mp2, "g__"); mp_s <- long_level(mp2, "s__")

## ---------- load Bracken combined (counts + fractions, with names) ----------
br_file <- Sys.glob(file.path(OUT, "bracken", "*combined*.txt"))
if (length(br_file)==0) br_file <- Sys.glob(file.path(OUT, "bracken", "**", "*combined*"))
br_file <- br_file[1]
br_counts <- NULL
if (!is.na(br_file)){
  br <- read.delim(br_file, check.names=FALSE)
  # 典型列: name, taxonomy_id, taxonomy_lvl, <sample>.bracken_num, <sample>.bracken_frac ...
  num_cols <- grep("_num$|\\.num$|bracken_num", names(br), value=TRUE)
  if (length(num_cols)==0) num_cols <- grep("num", names(br), value=TRUE)
  br_counts <- br %>% select(name, all_of(num_cols))
  names(br_counts)[-1] <- clean_id(names(br_counts)[-1])
  br_counts <- br_counts %>% group_by(name) %>% summarise(across(everything(), sum), .groups="drop")
  # 去掉宿主/人源污染 taxa（Kraken2 std 库含 human/mouse decoy；去宿主只除了小鼠比对上的 read，
  # handling/kit 带入的 human read 会被判为 Homo sapiens）——不能算进微生物相对丰度。
  HOST_TAXA <- c("Homo sapiens", "Mus musculus")
  br_counts <- br_counts %>% filter(!name %in% HOST_TAXA)
}
saveRDS(list(mp_g=mp_g, mp_s=mp_s, br_counts=br_counts), file.path(RDS, "phase1_tables.rds"))

## ---------- metadata ----------
ids <- setdiff(names(mp_g), "taxon")
meta <- tibble(sample=ids, group=factor(grp_of(ids), levels=c("AL","IF")))

## ---------- FIG1: genus-level stacked relative abundance (MetaPhlAn) ----------
stacked_plot <- function(mat, level_name, topn=11){
  m <- mat %>% column_to_rownames("taxon") %>% as.matrix()
  ord <- names(sort(rowMeans(m), decreasing=TRUE))
  top <- head(ord, topn)
  df <- as.data.frame(m) %>% rownames_to_column("taxon") %>%
    mutate(taxon=ifelse(taxon %in% top, taxon, "Other")) %>%
    group_by(taxon) %>% summarise(across(everything(), sum), .groups="drop") %>%
    pivot_longer(-taxon, names_to="sample", values_to="ra") %>%
    left_join(meta, by="sample")
  lev <- c(top, "Other")
  df$taxon <- factor(df$taxon, levels=lev)
  pal <- c(setNames(TAXA_PAL[seq_along(top)], top), Other=OTHER_COL)
  ggplot(df, aes(sample, ra, fill=taxon)) +
    geom_col(width=0.82, color="white", linewidth=0.18) +
    facet_grid(~group, scales="free_x", space="free_x") +
    scale_fill_manual(values=pal, name=level_name) +
    scale_y_continuous(expand=expansion(mult=c(0,0.02)), labels=label_number(suffix="%")) +
    labs(title=paste0("Gut microbiome composition — ", level_name, " level (MetaPhlAn)"),
         subtitle="Relative abundance per sample, faceted by diet arm", x=NULL, y="Relative abundance") +
    theme_pub() + theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
                        legend.text=element_text(size=8, face="italic"))
}
save_fig(stacked_plot(mp_g, "Genus"),  "fig1_composition_genus",  w=9, h=6)
save_fig(stacked_plot(mp_s, "Species"),"fig1b_composition_species",w=9.5,h=6.5)

## ---------- FIG2: alpha diversity (Bracken counts) ----------
if (!is.null(br_counts)){
  cm <- br_counts %>% column_to_rownames("name") %>% as.matrix()
  cm <- cm[, meta$sample, drop=FALSE]; cmt <- t(cm)
  alpha <- tibble(sample=rownames(cmt),
                  Shannon=diversity(cmt, "shannon"),
                  Simpson=diversity(cmt, "simpson"),
                  Observed=rowSums(cmt>0)) %>% left_join(meta, by="sample") %>%
           pivot_longer(c(Shannon,Simpson,Observed), names_to="metric", values_to="value")
  alpha$metric <- factor(alpha$metric, levels=c("Observed","Shannon","Simpson"))
  pvals <- alpha %>% group_by(metric) %>%
    summarise(p=wilcox.test(value~group)$p.value, y=max(value)*1.06, .groups="drop")
  p2 <- ggplot(alpha, aes(group, value, fill=group)) +
    geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.85, color="grey25", linewidth=0.3) +
    geom_jitter(width=0.12, size=1.8, alpha=0.8, color="grey20") +
    geom_text(data=pvals, aes(x=1.5, y=y, label=sprintf("Wilcoxon p = %.3f", p)),
              inherit.aes=FALSE, size=3.2, color="grey30") +
    facet_wrap(~metric, scales="free_y") +
    scale_fill_manual(values=GRP_COL, guide="none") +
    labs(title="Alpha diversity by diet arm", subtitle="Bracken species-level counts; points = samples",
         x=NULL, y="Diversity index") + theme_pub()
  save_fig(p2, "fig2_alpha_diversity", w=8.5, h=4.2)

  ## ---------- FIG3: beta diversity PCoA (Bray-Curtis) + PERMANOVA ----------
  ra <- sweep(cmt, 1, rowSums(cmt), "/")
  bc <- vegdist(ra, "bray")
  pco <- cmdscale(bc, k=2, eig=TRUE)
  ve <- round(100*pco$eig[1:2]/sum(pco$eig[pco$eig>0]), 1)
  perm <- adonis2(bc ~ group, data=meta, permutations=999)
  pc <- as.data.frame(pco$points); names(pc) <- c("PCo1","PCo2")
  pc <- pc %>% rownames_to_column("sample") %>% left_join(meta, by="sample")
  p3 <- ggplot(pc, aes(PCo1, PCo2, color=group, fill=group)) +
    stat_ellipse(geom="polygon", alpha=0.12, color=NA, level=0.8) +
    geom_point(size=3.2, alpha=0.9) +
    geom_text_repel(aes(label=sample), size=2.7, color="grey30", max.overlaps=20, seg.color="grey70") +
    scale_color_manual(values=GRP_COL, name="Diet") + scale_fill_manual(values=GRP_COL, guide="none") +
    labs(title="Beta diversity — Bray-Curtis PCoA",
         subtitle=sprintf("PERMANOVA: R² = %.3f, p = %.3f (999 perms)",
                          perm$R2[1], perm$`Pr(>F)`[1]),
         x=sprintf("PCo1 (%.1f%%)", ve[1]), y=sprintf("PCo2 (%.1f%%)", ve[2])) +
    theme_pub()
  save_fig(p3, "fig3_beta_pcoa", w=7, h=5.5)
  write.table(as.data.frame(as.matrix(bc)), file.path(TAB,"braycurtis_distance.tsv"), sep="\t", quote=FALSE)

  ## ---------- FIG4: diet-associated trends among CORE abundant species ----------
  # Kraken2/Bracken 在 species 层产生大量近零丰度假阳性(数千 taxa)；直接按 |log2FC| 挑 top
  # 会全是环境/污染噪声菌。先按 丰度+流行度 过滤到核心菌，再检验，诚实标注显著性。
  # 全体 taxa 也做一版检验存表（供追溯），但出图/结论只用 core。
  g <- meta$group
  test_set <- function(taxa){
    lapply(taxa, function(t){
      v <- ra[,t]; if (sum(v>0) < 3) return(NULL)
      data.frame(taxon=t, meanAL=mean(v[g=="AL"]), meanIF=mean(v[g=="IF"]),
                 meanRA=mean(v), prev=sum(v>0),
                 log2FC=log2((mean(v[g=="IF"])+1e-6)/(mean(v[g=="AL"])+1e-6)),
                 p=suppressWarnings(wilcox.test(v~g)$p.value))
    }) %>% bind_rows()
  }
  diff_all <- test_set(colnames(ra)); diff_all$padj <- p.adjust(diff_all$p,"BH")
  write.table(diff_all %>% arrange(p), file.path(TAB,"differential_abundance_ALL_bracken.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  # 核心菌：平均相对丰度 ≥ 0.1% 且 ≥5/10 样本检出
  core <- colnames(ra)[colMeans(ra) >= 1e-3 & colSums(ra>0) >= 5]
  diffc <- test_set(core); diffc$padj <- p.adjust(diffc$p,"BH")
  diffc <- diffc %>% arrange(p)
  write.table(diffc, file.path(TAB,"differential_abundance_CORE_bracken.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  nsig <- sum(diffc$padj < 0.05, na.rm=TRUE)
  topd <- diffc %>% slice_max(meanRA, n=20) %>%
    mutate(taxon=factor(taxon, levels=taxon[order(log2FC)]),
           dir=ifelse(log2FC>0,"Higher in IF","Higher in AL"),
           sig=ifelse(padj<0.05,"*",""))
  sub4 <- if (nsig==0)
      sprintf("Core species only (mean RA ≥0.1%%, %d spp.); none significant at FDR<0.05 — trends", nrow(diffc))
    else sprintf("%d of %d core species significant at FDR<0.05", nsig, nrow(diffc))
  p4 <- ggplot(topd, aes(log2FC, taxon, color=dir)) +
    geom_segment(aes(x=0, xend=log2FC, yend=taxon), linewidth=0.5) +
    geom_point(aes(size=meanRA*100)) +
    geom_text(aes(label=sig), color="black", size=5, vjust=0.75, show.legend=FALSE) +
    geom_vline(xintercept=0, color="grey50", linewidth=0.3) +
    scale_color_manual(values=c("Higher in IF"=GRP_COL[["IF"]], "Higher in AL"=GRP_COL[["AL"]]), name=NULL) +
    scale_size_continuous(name="Mean RA (%)", range=c(1.5,7)) +
    labs(title="Diet-associated trends among core abundant species (Bracken)",
         subtitle=sub4, x=expression(log[2]~fold~change~(IF/AL)), y=NULL) +
    theme_pub() + theme(axis.text.y=element_text(face="italic", size=8))
  save_fig(p4, "fig4_differential_abundance", w=8.5, h=6)
}

## ---------- FIG5: cross-tool concordance (Bracken vs MetaPhlAn, species) ----------
if (!is.null(br_counts)){
  br_ra <- br_counts %>% column_to_rownames("name") %>% as.matrix()
  br_ra <- sweep(br_ra, 2, colSums(br_ra), "/")*100
  br_mean <- tibble(taxon=rownames(br_ra), bracken=rowMeans(br_ra[, meta$sample, drop=FALSE]))
  mp_mean <- mp_s %>% mutate(metaphlan=rowMeans(across(all_of(meta$sample)))) %>% select(taxon, metaphlan)
  # 物种名对齐（MetaPhlAn 用下划线，Bracken 用空格）
  mp_mean$key <- tolower(gsub("[ _]", "", mp_mean$taxon))
  br_mean$key <- tolower(gsub("[ _]", "", br_mean$taxon))
  cc <- inner_join(mp_mean, br_mean, by="key") %>% filter(metaphlan>0, bracken>0)
  if (nrow(cc) > 3){
    rho <- suppressWarnings(cor(log10(cc$metaphlan), log10(cc$bracken), method="spearman"))
    lab <- cc %>% slice_max(metaphlan+bracken, n=8)
    p5 <- ggplot(cc, aes(metaphlan, bracken)) +
      geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
      geom_point(color="#4E79A7", size=2.4, alpha=0.8) +
      geom_text_repel(data=lab, aes(label=taxon.x), size=2.6, fontface="italic", color="grey30", max.overlaps=15) +
      scale_x_log10() + scale_y_log10() +
      labs(title="Cross-tool concordance (species)",
           subtitle=sprintf("Mean relative abundance; Spearman ρ = %.2f", rho),
           x="MetaPhlAn (%)", y="Bracken (%)") + theme_pub()
    save_fig(p5, "fig5_crosstool_concordance", w=6.5, h=5.5)
  }
}

## ---------- export composition tables ----------
write.table(mp_g, file.path(TAB,"composition_genus_metaphlan.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(mp_s, file.path(TAB,"composition_species_metaphlan.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
message("=== downstream analysis DONE ===")
