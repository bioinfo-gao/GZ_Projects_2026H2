#!/usr/bin/env Rscript
# P17 MAG — fig11 (bin quality scatter) + fig12 (recovered-MAG phylum composition), AL vs IF.
# 风格延续 fig1-10：Okabe-Ito(AL=#0072B2 IF=#D55E00), theme_pub, env regular_bioinfo.
suppressPackageStartupMessages({ library(ggplot2) })
PROJ <- "/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics"
OUT  <- file.path(PROJ, "custom_research_report_20260720", "assembly_binning")
GRP  <- c(AL="#0072B2", IF="#D55E00")

d <- read.table(file.path(OUT, "final_MAG_catalog.tsv"), header=TRUE, sep="\t",
                quote="", check.names=FALSE, na.strings=c("", "NA"))
d$completeness_pct  <- as.numeric(d$completeness_pct)
d$contamination_pct <- as.numeric(d$contamination_pct)
d$quality_tier_MIMAG <- factor(d$quality_tier_MIMAG, levels=c("High","Medium","Low"))

theme_pub <- function(b=12) theme_minimal(base_size=b) +
  theme(panel.grid.minor=element_blank(),
        axis.line=element_line(color="grey30",linewidth=0.3),
        axis.ticks=element_line(color="grey30",linewidth=0.3),
        plot.title=element_text(face="bold",size=b+2),
        plot.subtitle=element_text(color="grey30",size=b-1))

## Fig11: completeness vs contamination, MIMAG thresholds as reference lines
p11 <- ggplot(d, aes(completeness_pct, contamination_pct, color=group, shape=quality_tier_MIMAG)) +
  geom_vline(xintercept=c(50,90), linetype="dashed", color="grey70") +
  geom_hline(yintercept=c(5,10), linetype="dashed", color="grey70") +
  geom_point(size=2.6, alpha=0.85) +
  scale_color_manual(values=GRP, name="Group") +
  scale_shape_manual(values=c(High=16, Medium=17, Low=4), name="MIMAG tier") +
  labs(title="Final MAG catalog: completeness vs. contamination (n=184)",
       subtitle="Dashed lines = MIMAG thresholds (comp 50/90%, cont 5/10%); one point per non-redundant DASTool-refined genome",
       x="CheckM2 completeness (%)", y="CheckM2 contamination (%)") +
  theme_pub()
ggsave(file.path(OUT,"fig11_mag_quality_scatter.png"), p11, width=7.6, height=5.8, dpi=300)
ggsave(file.path(OUT,"fig11_mag_quality_scatter.pdf"), p11, width=7.6, height=5.8)

## Fig12: phylum-level composition of recovered MAGs, AL vs IF (classified bins only)
d$phylum <- ifelse(is.na(d$gtdb_classification) | d$gtdb_classification=="", NA,
                    sub(".*p__([^;]*);.*", "\\1", d$gtdb_classification))
dc <- d[!is.na(d$phylum) & d$phylum!="", ]
tab <- as.data.frame(table(group=dc$group, phylum=dc$phylum))
tab <- tab[tab$Freq>0,]
# collapse rare phyla (<3 total MAGs) into "Other" to keep legend readable
tot <- aggregate(Freq~phylum, tab, sum)
rare <- tot$phylum[tot$Freq<3]
tab$phylum <- ifelse(tab$phylum %in% rare, "Other", as.character(tab$phylum))
tab <- aggregate(Freq~group+phylum, tab, sum)

p12 <- ggplot(tab, aes(x=group, y=Freq, fill=phylum)) +
  geom_col(position="stack", width=0.6, color="white", linewidth=0.3) +
  scale_fill_brewer(palette="Set2", name="Phylum (GTDB)") +
  labs(title="Phylum-level composition of recovered MAGs",
       subtitle=sprintf("GTDB-Tk classified bins only (AL: %d/%d, IF: %d/%d); rare phyla (<3 MAGs) collapsed to 'Other'",
                         sum(dc$group=="AL"), sum(d$group=="AL"), sum(dc$group=="IF"), sum(d$group=="IF")),
       x=NULL, y="Number of MAGs") +
  theme_pub()
ggsave(file.path(OUT,"fig12_mag_taxonomy_composition.png"), p12, width=7.2, height=5.8, dpi=300)
ggsave(file.path(OUT,"fig12_mag_taxonomy_composition.pdf"), p12, width=7.2, height=5.8)

cat("DONE fig11 (n=",nrow(d),") fig12 (classified n=",nrow(dc),")\n", sep="")
