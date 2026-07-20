#!/usr/bin/env Rscript
# P17 — 样本质量/相似度 QC: Bray-Curtis 样本距离热图 + 层次聚类 + 离群/敏感性量化。
# 承接离群分析(IF_4_03_11 / AL_4_02_25)。风格: AL=#0072B2 IF=#D55E00 (Okabe-Ito), env regular_bioinfo。
suppressPackageStartupMessages({ library(pheatmap) })
PROJ  <- "/home/gao/projects_2026H2/17_Daniel_Mendes_gut_metagenomics"
DIV   <- file.path(PROJ, "custom_research_report_20260718", "diversity")

# --- 稳健读取(矩阵表头无角标: 首行10列名, 数据行=名+10值) ---
L    <- readLines(file.path(DIV,"braycurtis_distance.tsv"))
cn   <- strsplit(L[1], "\t")[[1]]                       # 10 列名
body <- do.call(rbind, strsplit(L[-1], "\t"))           # 10 行 x 11 列(名+值)
rn   <- body[,1]
M    <- apply(body[,-1,drop=FALSE], 2, as.numeric); rownames(M) <- rn; colnames(M) <- cn
M    <- M[order(rownames(M)), order(colnames(M))]
M    <- (M + t(M))/2; diag(M) <- 0                      # 对称 + 零对角

arm    <- ifelse(grepl("IF", rownames(M)), "IF", "AL")
ann    <- data.frame(Diet=arm); rownames(ann) <- rownames(M)
labels <- sub("HFD_", "", rownames(M))
hc     <- hclust(as.dist(M), method="average")

for (ext in c("png","pdf")) {
  pheatmap(M,
    color = colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(100),
    cluster_rows=hc, cluster_cols=hc,
    annotation_row=ann, annotation_col=ann,
    annotation_colors=list(Diet=c(AL="#0072B2", IF="#D55E00")),
    labels_row=labels, labels_col=labels,
    display_numbers=TRUE, number_format="%.2f", fontsize_number=7, number_color="grey25",
    main="Sample-to-sample Bray-Curtis distance (species) + hierarchical clustering",
    fontsize=10, border_color="white", treeheight_row=32, treeheight_col=32,
    width=8.6, height=7.2, filename=file.path(DIV, paste0("fig9_sample_distance_heatmap.",ext)))
}

# --- 离群量化: 每样本到同组其他的平均 Bray-Curtis (越大越离群) ---
IF <- rownames(M)[arm=="IF"]; AL <- rownames(M)[arm=="AL"]
m2g <- function(s){ g <- if (grepl("IF",s)) IF else AL; mean(M[s, setdiff(g,s)]) }
out <- data.frame(sample=rownames(M), arm=arm, mean_within_group_BC=round(sapply(rownames(M), m2g),3))
out <- out[order(-out$mean_within_group_BC),]
write.table(out, file.path(DIV,"sample_outlier_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# --- 敏感性: 去离群后 between/within 距离比 ---
ratio <- function(drop){
  keep <- setdiff(rownames(M), drop)
  IFk <- keep[grepl("IF",keep)]; ALk <- keep[grepl("AL",keep)]
  wi  <- c(M[IFk,IFk][upper.tri(M[IFk,IFk])], M[ALk,ALk][upper.tri(M[ALk,ALk])])
  bt  <- as.vector(M[IFk,ALk])
  c(within=mean(wi), between=mean(bt), between_within_ratio=mean(bt)/mean(wi))
}
sens <- rbind(`all_10`=ratio(character(0)),
              `drop_IF_4_03_11`=ratio("HFD_IF_4_03_11"),
              `drop_AL_4_02_25`=ratio("HFD_AL_4_02_25"),
              `drop_both`=ratio(c("HFD_IF_4_03_11","HFD_AL_4_02_25")))
write.table(round(sens,3), file.path(DIV,"sample_outlier_sensitivity.tsv"), sep="\t", quote=FALSE, col.names=NA)
cat("DONE: fig9_sample_distance_heatmap + sample_outlier_summary.tsv + sample_outlier_sensitivity.tsv\n")
print(out); cat("\n"); print(round(sens,3))
