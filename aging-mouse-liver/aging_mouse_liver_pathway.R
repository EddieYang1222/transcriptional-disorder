# Pathway Enrichment of Transcriptional Dyscoordination (Aging Mouse Liver)

# Online links
# https://www.gsea-msigdb.org/gsea/msigdb
# https://maayanlab.cloud/Enrichr/

# Set up working directory
# setwd("path/to/working/directory")

library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)
library(fgsea)
library(msigdbr)
library(enrichR)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# MSigDB mouse pathway sets: C2:CP + KEGG (legacy + medicus) + C5:GO
msig_cp <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")
msig_kegg_legacy  <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY")
msig_kegg_medicus <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_MEDICUS")
msig_go <- msigdbr(species = "Mus musculus", category = "C5")
msig_pathways <- rbind(msig_cp[, c("gs_name", "gene_symbol")], msig_kegg_legacy[, c("gs_name", "gene_symbol")],
                       msig_kegg_medicus[, c("gs_name", "gene_symbol")], msig_go[, c("gs_name", "gene_symbol")])
msig_pathways <- split(msig_pathways$gene_symbol, msig_pathways$gs_name)

######################################################
# 1. Genes gaining gene-level dyscoordination with age, per zone (Fig 7E)
# Per-zone gene-level dyscoordination LFC from the zonation-stratified run.
liver_hep_only_gene_level_zonation <- read.csv("liver_data_hepatocyte_only_estimated_dispersion_zonation_SAVER.csv")
liver_hep_only_gene_level_zonation_lfc <- liver_hep_only_gene_level_zonation %>%
  select(Gene, Age, cell_type, Gene_level_deviation) %>%
  pivot_wider(names_from = Age, values_from = Gene_level_deviation) %>%
  mutate(logFC = log((old) / (young))) %>%
  drop_na() %>% filter(is.finite(logFC))

fgsea_by_zone <- function(zone, lfc_table, pathways) {
  ranks <- lfc_table %>%
    filter(cell_type == zone) %>%
    arrange(desc(logFC)) %>%
    with(setNames(logFC, Gene))
  fgsea(pathways = pathways, stats = ranks)
}

zones <- unique(liver_hep_only_gene_level_zonation_lfc$cell_type)
liver_aging_hep_only_zonation_fgsea_results <- lapply(zones, function(z) {
  res <- fgsea_by_zone(z, liver_hep_only_gene_level_zonation_lfc, msig_pathways)
  res$zone <- z
  res
})
liver_aging_hep_only_zonation_fgsea_results <- na.omit(bind_rows(liver_aging_hep_only_zonation_fgsea_results))
liver_aging_hep_only_zonation_fgsea_results_sig <- liver_aging_hep_only_zonation_fgsea_results[liver_aging_hep_only_zonation_fgsea_results$pval < 0.05, ]
liver_aging_hep_only_zonation_fgsea_results_sig <- liver_aging_hep_only_zonation_fgsea_results_sig %>%
  mutate(is_prefixed = grepl("^[A-Z]+_", pathway),
         Source = ifelse(is_prefixed, str_extract(pathway, "^[^_]+"), "KEGG"),
         Term = ifelse(is_prefixed, str_remove(pathway, "^[^_]+_"),
                       str_remove(pathway, " - Mus musculus \\(house mouse\\)")),
         Term = str_replace_all(Term, "_", " "),
         Term = str_to_sentence(Term))
liver_aging_hep_only_zonation_fgsea_results_sig <- liver_aging_hep_only_zonation_fgsea_results_sig %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ", ")))
liver_aging_hep_only_zonation_fgsea_results_sig <- liver_aging_hep_only_zonation_fgsea_results_sig[
  liver_aging_hep_only_zonation_fgsea_results_sig$NES > 0 &
  liver_aging_hep_only_zonation_fgsea_results_sig$Source != "HP", ]
write.csv(liver_aging_hep_only_zonation_fgsea_results_sig, "liver_data_hepatocyte_only_pathway_by_zone_gene_level.csv", row.names = FALSE)

# enrichR on zone-specific genes that gain dyscoordination (optional)
liver_hep_only_gene_level_zonation_lfc_pos <- liver_hep_only_gene_level_zonation_lfc[liver_hep_only_gene_level_zonation_lfc$logFC > 0, ]
liver_hep_only_gene_level_zonation_lfc_pos$cell_type <- factor(liver_hep_only_gene_level_zonation_lfc_pos$cell_type, levels = c("Periportal", "Midlobular", "Pericentral"))
enrich_results <- list()
for (zone in levels(liver_hep_only_gene_level_zonation_lfc_pos$cell_type)) {
  sig_genes <- liver_hep_only_gene_level_zonation_lfc_pos[liver_hep_only_gene_level_zonation_lfc_pos$cell_type == zone, ]$Gene
  if (length(sig_genes) > 0) {
    enriched <- enrichr(sig_genes, c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "KEGG_2021_Mouse"))
    for (db in names(enriched)) {
      df <- enriched[[db]]
      if (nrow(df) == 0) next
      df$Zone <- zone
      df$Database <- db
      enrich_results[[length(enrich_results) + 1]] <- df
    }
  }
}
enrich_results <- bind_rows(enrich_results)
enrich_results <- enrich_results[enrich_results$Adjusted.P.value < 0.05, ]
write.csv(enrich_results, "liver_data_hepatocyte_only_pathway_by_zone_gene_level_v2.csv", row.names = FALSE)

######################################################
# 2. Genes activated in high-dyscoordination cells, per zone (Fig 7F)
# Rebuild the zone-annotated hepatocyte Seurat object with cell-level
# dyscoordination attached (bins from the zonation script).
load("liver_data_hepatocyte_only_expr.RData")
liver_aging_hep_only_cell_level_filtered_meta <- read.csv("Peter_liver_aging_zonation_score_bins.csv", stringsAsFactors = FALSE)

liver_aging_hep_only_cell_level_filtered_meta_cells <- liver_aging_hep_only_cell_level_filtered_meta[
  !is.na(liver_aging_hep_only_cell_level_filtered_meta$zonation_region), ]$Cell_barcode
liver_hep_only_expr_filtered <- liver_hep_only_expr[, liver_aging_hep_only_cell_level_filtered_meta_cells]
liver_hep_only_filtered <- CreateSeuratObject(counts = liver_hep_only_expr_filtered,
                                              meta.data = as.data.frame(liver_aging_hep_only_cell_level_filtered_meta[
                                                !is.na(liver_aging_hep_only_cell_level_filtered_meta$zonation_region), ]))
DefaultAssay(liver_hep_only_filtered) <- "RNA"
liver_hep_only_filtered <- NormalizeData(liver_hep_only_filtered)
liver_hep_only_filtered <- FindVariableFeatures(liver_hep_only_filtered)
liver_hep_only_filtered <- ScaleData(liver_hep_only_filtered, features = VariableFeatures(liver_hep_only_filtered))
Idents(liver_hep_only_filtered) <- factor(liver_hep_only_filtered$zonation_region, levels = c("Periportal", "Midlobular", "Pericentral"))
liver_hep_only_filtered$Cell_level_deviation <- liver_aging_hep_only_cell_level_filtered_meta$Cell_level_deviation[
  match(colnames(liver_hep_only_filtered), liver_aging_hep_only_cell_level_filtered_meta$Cell_barcode)]

# High vs low cell-level dyscoordination within each zone (top/bottom quartile)
liver_hep_only_filtered$Dyscoordination_group <- NA
for (zone in levels(Idents(liver_hep_only_filtered))) {
  zone_cells <- WhichCells(liver_hep_only_filtered, idents = zone)
  dev_vals <- liver_hep_only_filtered$Cell_level_deviation[zone_cells]
  q25 <- quantile(dev_vals, 0.25, na.rm = TRUE)
  q75 <- quantile(dev_vals, 0.75, na.rm = TRUE)
  liver_hep_only_filtered$Dyscoordination_group[zone_cells] <- ifelse(
    dev_vals <= q25, "Low", ifelse(dev_vals >= q75, "High", NA))
}

de_results <- list()
fgsea_results <- list()
for (zone in levels(Idents(liver_hep_only_filtered))) {
  zone_cells <- WhichCells(liver_hep_only_filtered, idents = zone)
  valid_cells <- zone_cells[!is.na(liver_hep_only_filtered$Dyscoordination_group[zone_cells])]
  sub_obj <- subset(liver_hep_only_filtered, cells = valid_cells)
  Idents(sub_obj) <- factor(sub_obj$Dyscoordination_group, levels = c("Low", "High"))
  markers <- FindMarkers(sub_obj, ident.1 = "High", ident.2 = "Low", test.use = "wilcox", min.pct = 0, logfc.threshold = 0)
  markers$Zone <- zone
  markers$gene <- rownames(markers)
  de_results[[zone]] <- markers
  ranks <- markers$avg_log2FC
  names(ranks) <- markers$gene
  ranks <- sort(ranks, decreasing = TRUE)
  fgseaRes <- fgsea(pathways = msig_pathways, stats = ranks)
  fgseaRes$Zone <- zone
  fgsea_results[[zone]] <- fgseaRes
}

de_results <- bind_rows(de_results)
fgsea_results <- bind_rows(fgsea_results)
fgsea_results_sig <- fgsea_results[fgsea_results$padj < 0.05, ]
fgsea_results_sig <- fgsea_results_sig %>%
  mutate(is_prefixed = grepl("^[A-Z]+_", pathway),
         Source = ifelse(is_prefixed, str_extract(pathway, "^[^_]+"), "KEGG"),
         Term = ifelse(is_prefixed, str_remove(pathway, "^[^_]+_"),
                       str_remove(pathway, " - Mus musculus \\(house mouse\\)")),
         Term = str_replace_all(Term, "_", " "),
         Term = str_to_sentence(Term))
fgsea_results_sig <- fgsea_results_sig %>% mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ", ")))
fgsea_results_sig <- fgsea_results_sig[fgsea_results_sig$NES > 0 & fgsea_results_sig$Source != "HP", ]
write.csv(fgsea_results_sig, "liver_data_hepatocyte_only_pathway_by_zone.csv", row.names = FALSE)

# enrichR on High-vs-Low DE genes per zone (optional)
de_results <- list()
enrich_results <- list()
for (zone in levels(Idents(liver_hep_only_filtered))) {
  zone_cells <- WhichCells(liver_hep_only_filtered, idents = zone)
  valid_cells <- zone_cells[!is.na(liver_hep_only_filtered$Dyscoordination_group[zone_cells])]
  sub_obj <- subset(liver_hep_only_filtered, cells = valid_cells)
  Idents(sub_obj) <- factor(sub_obj$Dyscoordination_group, levels = c("Low", "High"))
  markers <- FindMarkers(sub_obj, ident.1 = "High", ident.2 = "Low", test.use = "wilcox", min.pct = 0.2, logfc.threshold = 0.25)
  markers$Zone <- zone
  markers$gene <- rownames(markers)
  de_results[[zone]] <- markers
  sig_genes <- rownames(markers[markers$p_val_adj < 0.05, ])
  if (length(sig_genes) > 0) {
    enriched <- enrichr(sig_genes, c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "KEGG_2021_Mouse"))
    for (db in names(enriched)) {
      df <- enriched[[db]]
      if (nrow(df) == 0) next
      df$Zone <- zone
      df$Database <- db
      enrich_results[[length(enrich_results) + 1]] <- df
    }
  }
}
de_results <- bind_rows(de_results)
de_results <- de_results[de_results$p_val_adj < 0.05, ]
enrich_results <- bind_rows(enrich_results)
enrich_results <- enrich_results[enrich_results$Adjusted.P.value < 0.05, ]
write.csv(de_results, "liver_data_hepatocyte_only_de_gene_by_zone.csv", row.names = FALSE)
write.csv(enrich_results, "liver_data_hepatocyte_only_pathway_by_zone_v2.csv", row.names = FALSE)

######################################################
