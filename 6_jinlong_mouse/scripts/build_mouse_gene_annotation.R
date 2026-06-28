library(biomaRt)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(KEGGREST)
library(dplyr)
library(writexl)

OUT_FILE <- "/home/gao/projects_2026H1/Genes/mouse_Gene_annotation_20260628.xlsx"
cat("=== Building mouse gene annotation ===\n")

# ── 1. biomaRt: gene_id / gene_symbol / gene_type / description / GO / UniProt ──
cat("[1/5] Connecting to Ensembl biomaRt...\n")
mart <- useEnsembl("genes", dataset = "mmusculus_gene_ensembl",
                   mirror = "useast")

cat("[2/5] Querying gene info + GO + UniProt...\n")
bm <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "gene_biotype",
    "description",
    "go_id",
    "uniprot_gn_id"
  ),
  mart = mart
)

cat("  Raw rows from biomaRt:", nrow(bm), "\n")

# ── 2. Collapse multi-value columns per gene ──
cat("[3/5] Collapsing multi-value columns...\n")

collapse_unique <- function(x) {
  vals <- unique(x[!is.na(x) & x != ""])
  if (length(vals) == 0) NA_character_ else paste(vals, collapse = ";")
}

gene_base <- bm %>%
  group_by(ensembl_gene_id, external_gene_name, gene_biotype) %>%
  summarise(
    Description = collapse_unique(sub(" \\[Source:.*\\]", "", description)),
    GO          = collapse_unique(go_id),
    UniProt     = collapse_unique(uniprot_gn_id),
    .groups = "drop"
  ) %>%
  rename(
    gene_id     = ensembl_gene_id,
    gene_symbol = external_gene_name,
    gene_type   = gene_biotype
  )

cat("  Unique genes:", nrow(gene_base), "\n")

# ── 3. KEGG via org.Mm.eg.db + KEGGREST ──
cat("[4/5] Adding KEGG annotations...\n")

# Ensembl → Entrez
ens2entrez <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = gene_base$gene_id,
  columns = c("ENSEMBL", "ENTREZID"),
  keytype = "ENSEMBL"
) %>% distinct(ENSEMBL, .keep_all = TRUE)

# Entrez → KEGG pathway IDs (mmu...)
entrez2kegg <- tryCatch({
  AnnotationDbi::select(
    org.Mm.eg.db,
    keys    = na.omit(ens2entrez$ENTREZID),
    columns = c("ENTREZID", "PATH"),
    keytype = "ENTREZID"
  ) %>%
    filter(!is.na(PATH)) %>%
    mutate(KEGG_ID = paste0("mmu", PATH)) %>%
    group_by(ENTREZID) %>%
    summarise(KEGG_ID = paste(unique(KEGG_ID), collapse = ";"), .groups = "drop")
}, error = function(e) {
  cat("  KEGG PATH lookup failed:", conditionMessage(e), "\n")
  data.frame(ENTREZID = character(), KEGG_ID = character())
})

# KEGG pathway descriptions via KEGGREST (batch, rate-limited)
kegg_ids_all <- unique(unlist(strsplit(entrez2kegg$KEGG_ID, ";")))
cat("  Fetching descriptions for", length(kegg_ids_all), "KEGG pathways...\n")

kegg_desc_map <- list()
batch_size <- 10
for (i in seq(1, length(kegg_ids_all), by = batch_size)) {
  batch <- kegg_ids_all[i:min(i + batch_size - 1, length(kegg_ids_all))]
  tryCatch({
    res <- keggGet(batch)
    for (r in res) {
      kegg_desc_map[[r$ENTRY]] <- r$NAME
    }
  }, error = function(e) NULL)
  Sys.sleep(0.3)   # KEGG API rate limit
}

# Build KEGG description string per gene
entrez2kegg$KEGG_Description <- sapply(entrez2kegg$KEGG_ID, function(ids) {
  if (is.na(ids)) return(NA_character_)
  descs <- sapply(strsplit(ids, ";")[[1]], function(kid) {
    entry <- sub("mmu", "", kid)
    kegg_desc_map[[entry]] %||% NA_character_
  })
  collapse_unique(descs)
})

# KO entries
entrez2ko <- tryCatch({
  AnnotationDbi::select(
    org.Mm.eg.db,
    keys    = na.omit(ens2entrez$ENTREZID),
    columns = c("ENTREZID", "UNIPROT"),
    keytype = "ENTREZID"
  ) %>% distinct()
}, error = function(e) NULL)

# EC numbers via biomaRt
cat("  Fetching EC numbers from biomaRt...\n")
ec_raw <- tryCatch(
  getBM(attributes = c("ensembl_gene_id", "ec"),
        filters    = "biotype",
        values     = "protein_coding",
        mart       = mart),
  error = function(e) data.frame(ensembl_gene_id = character(), ec = character())
)
ec_df <- ec_raw %>%
  filter(!is.na(ec) & ec != "") %>%
  group_by(ensembl_gene_id) %>%
  summarise(EC = paste(unique(ec), collapse = ";"), .groups = "drop") %>%
  rename(gene_id = ensembl_gene_id)

# KO via KEGGREST (gene → KO)
cat("  Fetching KO entries...\n")
entrez_ids <- na.omit(unique(ens2entrez$ENTREZID))
ko_map <- list()
ko_batch <- 10
for (i in seq(1, min(length(entrez_ids), 500), by = ko_batch)) {
  batch <- entrez_ids[i:min(i + ko_batch - 1, length(entrez_ids))]
  tryCatch({
    query_ids <- paste0("mmu:", batch)
    res <- keggGet(query_ids)
    for (r in res) {
      if (!is.null(r$ORTHOLOGY)) {
        eid <- sub("mmu:", "", r$ENTRY)
        ko_map[[eid]] <- paste(names(r$ORTHOLOGY), collapse = ";")
      }
    }
  }, error = function(e) NULL)
  Sys.sleep(0.3)
}

ens2entrez$KO_ENTRY <- sapply(ens2entrez$ENTREZID, function(eid) {
  if (is.na(eid)) return(NA_character_)
  ko_map[[eid]] %||% NA_character_
})

# ── 4. Merge everything ──
cat("[5/5] Merging all columns...\n")

`%||%` <- function(a, b) if (!is.null(a)) a else b

annotation <- gene_base %>%
  left_join(ens2entrez %>% select(ENSEMBL, ENTREZID, KO_ENTRY),
            by = c("gene_id" = "ENSEMBL")) %>%
  left_join(entrez2kegg %>% select(ENTREZID, KEGG_ID, KEGG_Description),
            by = "ENTREZID") %>%
  left_join(ec_df, by = "gene_id") %>%
  select(gene_id, gene_symbol, gene_type, Description,
         GO, UniProt, KEGG_ID, KEGG_Description, KO_ENTRY, EC) %>%
  arrange(gene_id)

cat("Final annotation rows:", nrow(annotation), "\n")
cat("gene_type breakdown:\n")
print(sort(table(annotation$gene_type), decreasing = TRUE)[1:15])

# ── 5. Save ──
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
write_xlsx(annotation, OUT_FILE)
cat("\n✅ Saved:", OUT_FILE, "\n")
cat("   Rows:", nrow(annotation), "  Columns:", ncol(annotation), "\n")
