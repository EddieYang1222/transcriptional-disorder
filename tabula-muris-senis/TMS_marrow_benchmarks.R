# Aging-score benchmarks for TMS bone marrow (SenMayo, SCENT, CytoTRACE)
# The remaining benchmarks in the paper -- mean Euclidean distance to the cell
# type / tissue average (decibel) and Scallop membership -- are computed in the
# repo-root decibel_scallop_analysis.ipynb, not here.

# Online links
# https://www.nature.com/articles/s41467-022-32552-1
# https://bioconductor.org/packages/release/bioc/html/escape.html
# https://github.com/aet21/SCENT
# https://github.com/gunetti/CytoTRACE

# Set up working directory
# setwd("path/to/working/directory")

library(escape)
library(Seurat)
library(Matrix)
library(dplyr)
library(parallel)
library(biomaRt)
library(CytoTRACE)
library(ggplot2)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Load TMS marrow counts and the fitted manifold
# TMS_marrow.RData / TMS_marrow_manifold.RData are written by TMS_marrow_deviation.R
load("TMS_marrow.RData")            # dataset.counts, dataset.age(.levels), dataset.celltype(.levels)
load("TMS_marrow_manifold.RData")   # dataset.saver (SAVER object, mu.out = manifold)

# Match the manifold: drop genes with zero counts
dataset.counts <- dataset.counts[rowSums(dataset.counts) != 0, ]
cell_names <- colnames(dataset.counts)
age_vec    <- dataset.age.levels[dataset.age]

######################################################
# 1. SenMayo senescence score (escape GSEA on the manifold)
# Enrichment is scored on the SAVER manifold (mu.out), not the raw counts.
dataset.metadata <- data.frame(row.names = colnames(dataset.counts), ages = dataset.age)
dataset.seurat <- CreateSeuratObject(counts = dataset.saver$mu.out, meta.data = dataset.metadata)
Idents(object = dataset.seurat) <- "ages"

# SenMayo signature (Saul et al. 2022, mouse symbols)
gene.sets <- list(SenMayo_signature = c(
  "Acvr1b", "Ang", "Angpt1", "Angptl4", "Areg", "Axl", "Bex3", "Bmp2", "Bmp6",
  "C3", "Ccl1", "Ccl13", "Ccl16", "Ccl2", "Ccl20", "Ccl24", "Ccl26", "Ccl3",
  "Ccl3l1", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Cd55", "Cd9", "Csf1", "Csf2",
  "Csf2rb", "Cst4", "Ctnnb1", "Ctsb", "Cxcl1", "Cxcl10", "Cxcl12", "Cxcl16",
  "Cxcl2", "Cxcl3", "Cxcl8", "Cxcr2", "Dkk1", "Edn1", "Egf", "Egfr", "Ereg",
  "Esm1", "Ets2", "Fas", "Fgf1", "Fgf2", "Fgf7", "Gdf15", "Gem", "Gmfg", "Hgf",
  "Hmgb1", "Icam1", "Icam3", "Igf1", "Igfbp1", "Igfbp2", "Igfbp3", "Igfbp4",
  "Igfbp5", "Igfbp6", "Igfbp7", "Il10", "Il13", "Il15", "Il18", "Il1a", "Il1b",
  "Il2", "Il32", "Il6", "Il6st", "Il7", "Inha", "Iqgap2", "Itga2", "Itpka",
  "Jun", "Kitlg", "Lcp1", "Mif", "Mmp1", "Mmp10", "Mmp12", "Mmp13", "Mmp14",
  "Mmp2", "Mmp3", "Mmp9", "Nap1l4", "Nrg1", "Pappa", "Pecam1", "Pgf", "Pigf",
  "Plat", "Plau", "Plaur", "Ptbp1", "Ptger2", "Ptges", "Rps6ka5", "Scamp4",
  "Selplg", "Sema3f", "Serpinb4", "Serpine1", "Serpine2", "Spp1", "Spx", "Timp2",
  "Tnf", "Tnfrsf10c", "Tnfrsf11b", "Tnfrsf1a", "Tnfrsf1b", "Tubgcp2", "Vegfa",
  "Vegfc", "Vgf", "Wnt16", "Wnt2"))

ES.seurat <- enrichIt(obj = dataset.seurat, gene.sets = gene.sets,
                      groups = 1000, cores = 4, min.size = 5)

senmayo_df <- data.frame(Cell_barcode = rownames(ES.seurat),
                         SenMayo = ES.seurat[["SenMayo_signature"]])
write.csv(senmayo_df, "TMS_marrow_senmayo.csv", row.names = FALSE)

######################################################
# 2. SCENT signaling entropy rate (SR)
# "Signaling entropy" is SCENT's own term for this method; it is not our measure.
# SCENT functions and PPI network are sourced directly from the SCENT repo.
scent_url <- "https://raw.githubusercontent.com/aet21/SCENT/master/R/"
scent_fns <- c("DoIntegPPI.R", "CompSRana.R")
tmp <- tempdir()
for (fn in scent_fns) {
  dest <- file.path(tmp, fn)
  download.file(paste0(scent_url, fn), destfile = dest, quiet = TRUE)
  source(dest)
}
# PPI network: net17Jan16.rda loads as net17Jan16.m (Entrez IDs)
ppi_dest <- file.path(tmp, "net17Jan16.rda")
download.file("https://raw.githubusercontent.com/aet21/SCENT/master/data/net17Jan16.rda",
              destfile = ppi_dest, quiet = TRUE, mode = "wb")
load(ppi_dest)

# Library-normalize and log2-transform (pseudocount 1.1 as SCENT expects)
size.factor <- colSums(dataset.counts) / median(colSums(dataset.counts))
norm_mat    <- sweep(as.matrix(dataset.counts), 2, size.factor, "/")
norm_mat    <- log2(norm_mat + 1.1)

# Mouse -> human one-to-one orthologs, then human Entrez IDs, via biomaRt
mouse_genes_ls <- rownames(norm_mat)
mouse_mart <- tryCatch(
  useEnsembl("genes", dataset = "mmusculus_gene_ensembl", mirror = "useast"),
  error = function(e) useEnsembl("genes", dataset = "mmusculus_gene_ensembl", mirror = "uswest"))
mouse2human <- getBM(
  attributes = c("external_gene_name", "hsapiens_homolog_ensembl_gene",
                 "hsapiens_homolog_orthology_type"),
  filters    = "external_gene_name",
  values     = mouse_genes_ls,
  mart       = mouse_mart)
mouse2human <- mouse2human[
  mouse2human$hsapiens_homolog_ensembl_gene != "" &
  mouse2human$hsapiens_homolog_orthology_type == "ortholog_one2one", ]

human_mart <- tryCatch(
  useEnsembl("genes", dataset = "hsapiens_gene_ensembl", mirror = "useast"),
  error = function(e) useEnsembl("genes", dataset = "hsapiens_gene_ensembl", mirror = "uswest"))
human_anno <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  filters    = "ensembl_gene_id",
  values     = unique(mouse2human$hsapiens_homolog_ensembl_gene),
  mart       = human_mart)

mapping_ls <- merge(mouse2human, human_anno,
                    by.x = "hsapiens_homolog_ensembl_gene", by.y = "ensembl_gene_id")
mapping_ls <- mapping_ls[mapping_ls$hgnc_symbol != "" & !is.na(mapping_ls$entrezgene_id), ]
mapping_ls <- mapping_ls[!duplicated(mapping_ls$external_gene_name) &
                         !duplicated(mapping_ls$hgnc_symbol) &
                         !duplicated(mapping_ls$entrezgene_id), ]
write.csv(mapping_ls, "TMS_marrow_mouse_to_human_gene_mapping.csv", row.names = FALSE)

# Build Entrez-ID expression matrix for SCENT
common_mouse <- intersect(rownames(norm_mat), mapping_ls$external_gene_name)
norm_mat_sub <- norm_mat[common_mouse, , drop = FALSE]
mapping_ordered        <- mapping_ls[match(rownames(norm_mat_sub), mapping_ls$external_gene_name), ]
rownames(norm_mat_sub) <- as.character(mapping_ordered$entrezgene_id)
rm(norm_mat)
gc()

# DoIntegPPI + chunked CompSRana
integ <- DoIntegPPI(exp.m = norm_mat_sub, ppiA.m = net17Jan16.m)
rm(norm_mat_sub, net17Jan16.m)
gc()

# PPI integration keeps all cells (only genes are subset), so SR is in cell order
chunk_size  <- 2000
n_cells     <- ncol(integ$expMC)
cell_chunks <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / chunk_size))
sr_vals <- unlist(lapply(seq_along(cell_chunks), function(i) {
  idx               <- cell_chunks[[i]]
  integ_chunk       <- integ
  integ_chunk$expMC <- integ$expMC[, idx, drop = FALSE]
  res <- CompSRana(integ_chunk)
  if ("SR" %in% names(res)) res$SR else res[[1]]
}))
rm(integ)
gc()

scent_df <- data.frame(Cell_barcode = cell_names, SR = sr_vals)
write.csv(scent_df, "TMS_marrow_scent_SR.csv", row.names = FALSE)

######################################################
# 3. CytoTRACE developmental potency
# Runs on genes x cells raw counts; species-agnostic.
# install.packages("CytoTRACE_0.3.3.tar.gz", repos = NULL, type = "source")
ct_res <- CytoTRACE(mat = as.matrix(dataset.counts), enableFast = FALSE)
scores <- ct_res$CytoTRACE
cytotrace_df <- data.frame(
  Cell_barcode = names(scores),
  CytoTRACE    = as.numeric(scores),
  Age          = age_vec[match(names(scores), cell_names)])
write.csv(cytotrace_df, "TMS_marrow_cytotrace_scores.csv", row.names = FALSE)

######################################################
# 4. Benchmark comparison against transcriptional dyscoordination
# Join each score to the cell-level dyscoordination on Cell_barcode and rank-correlate.
cell_level_dyscoordination <- read.csv("TMS_marrow_cellular_dispersion_SAVER.csv")
cell_level_dyscoordination$Cell_level_deviation <- remove_outliers(cell_level_dyscoordination$Cell_level_deviation)
cell_level_dyscoordination <- na.omit(cell_level_dyscoordination)

bench <- cell_level_dyscoordination %>%
  left_join(senmayo_df, by = "Cell_barcode") %>%
  left_join(scent_df, by = "Cell_barcode") %>%
  left_join(cytotrace_df[, c("Cell_barcode", "CytoTRACE")], by = "Cell_barcode")

benchmark_cor <- data.frame(
  Benchmark = c("SenMayo", "SCENT_SR", "CytoTRACE"),
  Spearman_rho = c(
    cor(bench$Cell_level_deviation, bench$SenMayo,   method = "spearman", use = "complete.obs"),
    cor(bench$Cell_level_deviation, bench$SR,        method = "spearman", use = "complete.obs"),
    cor(bench$Cell_level_deviation, bench$CytoTRACE, method = "spearman", use = "complete.obs")))
write.csv(benchmark_cor, "TMS_marrow_benchmark_correlations.csv", row.names = FALSE)
cat("Benchmark Spearman correlations vs cell-level dyscoordination:\n")
print(benchmark_cor)

######################################################
