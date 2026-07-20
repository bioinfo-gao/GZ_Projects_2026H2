#!/usr/bin/env Rscript
# P17 — Phase 1b 功能通路下游分析与出图 (HUMAnN pathway relab)
# 输入: humann_merged/all_pathabundance_relab.tsv  (社区级 unstratified pathways)
# 输出: analysis/figures/fig6-8*  +  analysis/tables/
# env: regular_bioinfo。风格与 taxonomy 图一致 (Okabe-Ito, clean theme, honest significance)。
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(vegan); library(ggrepel); library(tibble); library(stringr)
})
set.seed(1)
PROJ <- "/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics"
FIG <- file.path(PROJ,"analysis","figures"); TAB <- file.path(PROJ,"analysis","tables")
dir.create(FIG,showWarnings=FALSE,recursive=TRUE); dir.create(TAB,showWarnings=FALSE,recursive=TRUE)

GRP_COL <- c(AL="#0072B2", IF="#D55E00")
theme_pub <- function(base=12){
  theme_minimal(base_size=base, base_family="sans") +
    theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
          axis.line=element_line(color="grey30",linewidth=0.3),
          axis.ticks=element_line(color="grey30",linewidth=0.3),
          plot.title=element_text(face="bold",size=base+2),
          plot.subtitle=element_text(color="grey30",size=base-1),
          strip.text=element_text(face="bold"))
}
save_fig <- function(p,name,w=8,h=5.5){
  ggsave(file.path(FIG,paste0(name,".pdf")),p,width=w,height=h,device=cairo_pdf)
  ggsave(file.path(FIG,paste0(name,".png")),p,width=w,height=h,dpi=300); message("saved: ",name)
}

## ---- load pathabundance (community-level, unstratified) ----
pa <- read.delim(file.path(PROJ,"humann_merged","all_pathabundance_relab.tsv"), check.names=FALSE)
names(pa)[1] <- "pathway"
names(pa)[-1] <- str_replace(names(pa)[-1], "_Abundance-RELAB$", "")
pa <- pa %>% filter(!str_detect(pathway, "\\|"),                       # unstratified only
                    !pathway %in% c("UNMAPPED","UNINTEGRATED"))
mat <- pa %>% column_to_rownames("pathway") %>% as.matrix()            # pathways x samples
ids <- colnames(mat)
meta <- tibble(sample=ids, group=factor(ifelse(grepl("_AL_",ids),"AL","IF"),levels=c("AL","IF")))
# renormalize to relative (rows sum vary; use column-relative within retained pathways)
ra <- t(sweep(mat, 2, colSums(mat), "/"))                              # samples x pathways
saveRDS(list(mat=mat, ra=ra, meta=meta), file.path(PROJ,"analysis","rds_cache","functional.rds"))

## ---- FIG6: functional beta diversity (Bray-Curtis PCoA + PERMANOVA) ----
bc <- vegdist(ra, "bray"); pco <- cmdscale(bc,k=2,eig=TRUE)
ve <- round(100*pco$eig[1:2]/sum(pco$eig[pco$eig>0]),1)
perm <- adonis2(bc~group, data=meta, permutations=999)
pc <- as.data.frame(pco$points); names(pc) <- c("PCo1","PCo2")
pc <- pc %>% rownames_to_column("sample") %>% left_join(meta,by="sample")
p6 <- ggplot(pc, aes(PCo1,PCo2,color=group,fill=group)) +
  stat_ellipse(geom="polygon",alpha=0.12,color=NA,level=0.8) +
  geom_point(size=3.2,alpha=0.9) +
  geom_text_repel(aes(label=sample),size=2.7,color="grey30",max.overlaps=20,seg.color="grey70") +
  scale_color_manual(values=GRP_COL,name="Diet") + scale_fill_manual(values=GRP_COL,guide="none") +
  labs(title="Functional beta diversity — pathway Bray-Curtis PCoA",
       subtitle=sprintf("PERMANOVA: R² = %.3f, p = %.3f (999 perms); HUMAnN MetaCyc pathways",
                        perm$R2[1], perm$`Pr(>F)`[1]),
       x=sprintf("PCo1 (%.1f%%)",ve[1]), y=sprintf("PCo2 (%.1f%%)",ve[2])) + theme_pub()
save_fig(p6,"fig6_functional_pcoa",w=7,h=5.5)

## ---- FIG7: differential pathways (core-filtered, Wilcoxon, BH, honest) ----
g <- meta$group
core <- colnames(ra)[colMeans(ra) >= 1e-3 & colSums(ra>0) >= 5]
dl <- lapply(core, function(t){
  v <- ra[,t]
  data.frame(pathway=t, meanAL=mean(v[g=="AL"]), meanIF=mean(v[g=="IF"]), meanRA=mean(v),
             log2FC=log2((mean(v[g=="IF"])+1e-6)/(mean(v[g=="AL"])+1e-6)),
             p=suppressWarnings(wilcox.test(v~g)$p.value))
}) %>% bind_rows()
dl$padj <- p.adjust(dl$p,"BH"); dl <- dl %>% arrange(p)
write.table(dl, file.path(TAB,"differential_pathways_humann.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
nsig <- sum(dl$padj<0.05,na.rm=TRUE)
shortname <- function(x) str_trunc(str_replace(x,"^[A-Za-z0-9+-]+:? ?",""), 42)
topd <- dl %>% slice_max(meanRA,n=20) %>%
  mutate(lab=str_trunc(pathway,52),
         lab=factor(lab, levels=lab[order(log2FC)]),
         dir=ifelse(log2FC>0,"Higher in IF","Higher in AL"))
sub7 <- if(nsig==0) sprintf("Core pathways (mean RA ≥0.1%%, %d); none significant at FDR<0.05 — trends",nrow(dl)) else sprintf("%d of %d core pathways significant at FDR<0.05",nsig,nrow(dl))
p7 <- ggplot(topd, aes(log2FC,lab,color=dir)) +
  geom_segment(aes(x=0,xend=log2FC,yend=lab),linewidth=0.5) +
  geom_point(aes(size=meanRA*100)) + geom_vline(xintercept=0,color="grey50",linewidth=0.3) +
  scale_color_manual(values=c("Higher in IF"=GRP_COL[["IF"]],"Higher in AL"=GRP_COL[["AL"]]),name=NULL) +
  scale_size_continuous(name="Mean RA (%)",range=c(1.5,6)) +
  labs(title="Diet-associated trends in metabolic pathways (HUMAnN)",
       subtitle=sub7, x=expression(log[2]~fold~change~(IF/AL)), y=NULL) +
  theme_pub() + theme(axis.text.y=element_text(size=7))
save_fig(p7,"fig7_pathway_differential",w=9.5,h=6)

## ---- FIG8: top pathways abundance heatmap (AL vs IF) ----
top <- names(sort(colMeans(ra),decreasing=TRUE))[1:25]
hd <- as.data.frame(ra[,top]) %>% rownames_to_column("sample") %>%
  pivot_longer(-sample,names_to="pathway",values_to="ra") %>% left_join(meta,by="sample") %>%
  mutate(pathway=str_trunc(pathway,50),
         sample=factor(sample, levels=meta$sample[order(meta$group)]))
hd$pathway <- factor(hd$pathway, levels=rev(unique(hd$pathway)))
p8 <- ggplot(hd, aes(sample,pathway,fill=ra*100)) +
  geom_tile(color="white",linewidth=0.3) +
  scale_fill_gradient(low="#f7fbff",high="#08519c",name="RA (%)") +
  facet_grid(~group,scales="free_x",space="free_x") +
  labs(title="Top 25 metabolic pathways by relative abundance", subtitle="HUMAnN MetaCyc, per sample", x=NULL,y=NULL) +
  theme_pub(11) + theme(axis.text.x=element_text(angle=45,hjust=1,size=7), axis.text.y=element_text(size=6.5))
save_fig(p8,"fig8_top_pathways_heatmap",w=10,h=8)

message("=== functional analysis DONE; sig pathways FDR<0.05: ", nsig, " ===")
