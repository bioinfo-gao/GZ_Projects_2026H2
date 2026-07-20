#!/usr/bin/env Rscript
# P17 — Aitchison PCA (CLR-transformed species abundances)：成分数据正确版的 PCA 排序，
# 补充 Bray-Curtis PCoA(fig3)，给习惯看 "PCA" 的读者。风格 Okabe-Ito(AL=#0072B2 IF=#D55E00), env regular_bioinfo。
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })
PROJ <- "/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics"
TAX  <- file.path(PROJ,"custom_research_report_20260718","taxonomy")
DIV  <- file.path(PROJ,"custom_research_report_20260718","diversity")
GRP  <- c(AL="#0072B2", IF="#D55E00")

d <- read.table(file.path(TAX,"composition_species_metaphlan.tsv"), header=TRUE, sep="\t",
                row.names=1, check.names=FALSE, quote="")
X <- t(as.matrix(d)); X <- X[, colSums(X) > 0, drop=FALSE]     # samples x taxa, 去全零
ps <- min(X[X>0]) * 0.5                                         # 零替换(半最小非零) + 重闭合 + CLR
Xr <- X; Xr[Xr==0] <- ps; Xr <- Xr / rowSums(Xr)
clr <- t(apply(Xr, 1, function(v){ lv <- log(v); lv - mean(lv) }))
pc  <- prcomp(clr, center=TRUE, scale.=FALSE)
ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
df  <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2],
                  Diet=ifelse(grepl("IF",rownames(clr)),"IF","AL"),
                  lab=sub("HFD_","",rownames(clr)))

theme_pub <- function(b=12) theme_minimal(base_size=b) +
  theme(panel.grid.minor=element_blank(),
        axis.line=element_line(color="grey30",linewidth=0.3),
        axis.ticks=element_line(color="grey30",linewidth=0.3),
        plot.title=element_text(face="bold",size=b+2),
        plot.subtitle=element_text(color="grey30",size=b-1))
p <- ggplot(df, aes(PC1,PC2,color=Diet)) +
  geom_hline(yintercept=0,color="grey85") + geom_vline(xintercept=0,color="grey85") +
  stat_ellipse(aes(group=Diet,fill=Diet),geom="polygon",alpha=0.08,color=NA,type="norm",level=0.8) +
  geom_point(size=3.4) +
  geom_text_repel(aes(label=lab),size=3,show.legend=FALSE,max.overlaps=20,seed=1) +
  scale_color_manual(values=GRP) + scale_fill_manual(values=GRP) +
  labs(title="Aitchison PCA (CLR-transformed species abundances)",
       subtitle="Compositional ordination; complements Bray-Curtis PCoA (Fig 3). Ellipses = 80% normal.",
       x=sprintf("PC1 (%.1f%%)",ve[1]), y=sprintf("PC2 (%.1f%%)",ve[2])) +
  theme_pub()
ggsave(file.path(DIV,"fig10_aitchison_pca.png"), p, width=7.6, height=5.8, dpi=300)
ggsave(file.path(DIV,"fig10_aitchison_pca.pdf"), p, width=7.6, height=5.8)
cat("DONE fig10_aitchison_pca; PC1=",ve[1],"% PC2=",ve[2],"%\n",sep="")
