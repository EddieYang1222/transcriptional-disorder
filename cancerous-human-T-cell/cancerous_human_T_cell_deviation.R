# Transcriptional dyscoordination analysis for cancerous human T cells (CD8 exhaustion)

# Online links
# https://www.biorxiv.org/content/10.1101/2025.09.01.673503v1

# Set up working directory
# setwd("S:/Penn Dropbox/Eddie Yang/Aging/Scripts/T_cell_exhaustion")

library(Seurat)
library(dplyr)
library(Matrix)
library(pbapply)
library(SAVER)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Load and prepare data
data_dir <- 'path/to/T_cell_data'
dataset.counts.sparse <- ReadMtx(
  mtx = file.path(data_dir, 'CD830dim_sparse/matrix.mtx'),
  cells = file.path(data_dir, 'CD830dim_sparse/barcodes.tsv'),
  features = file.path(data_dir, 'CD830dim_sparse/features.tsv'),
  feature.column = 1
)
meta <- read.delim(file.path(data_dir, 'CD8_cell_meta.tsv'), sep = '\t', row.names = 1)

# Filter for cells to keep
common.cells <- intersect(colnames(dataset.counts.sparse), rownames(meta))
dataset.counts.sparse <- dataset.counts.sparse[, common.cells]
meta <- meta[common.cells, , drop = FALSE]
exclude_types <- c('Prolif T', 'ISG T')
keep.cells <- !meta$cell_type %in% exclude_types
dataset.counts.sparse <- dataset.counts.sparse[, keep.cells]
meta <- meta[keep.cells, , drop = FALSE]

# Apply QC filters
keep.cells.qc <- Matrix::colSums(dataset.counts.sparse > 0) >= 200 &
                 Matrix::colSums(dataset.counts.sparse) >= 500
keep.genes.qc <- Matrix::rowSums(dataset.counts.sparse > 0) >= 5
dataset.counts.sparse <- dataset.counts.sparse[keep.genes.qc, keep.cells.qc]
meta <- meta[keep.cells.qc, , drop = FALSE]

# Build condition string (Tumor.Grade x Treatment), keeping Normal grade-only
dataset.celltype <- meta$cell_type
dataset.tumor.grade <- meta$Tumor.Grade
dataset.treatment <- meta$Treatment
dataset.condition_group <- dplyr::case_when(
  dataset.tumor.grade == 'Normal' ~ 'Normal',
  dataset.treatment == 'Treated' ~ paste0(dataset.tumor.grade, ' (treated)'),
  TRUE ~ as.character(dataset.tumor.grade)
)

dataset.counts <- as.matrix(dataset.counts.sparse)
size.factor <- colSums(dataset.counts) / mean(colSums(dataset.counts))
save(dataset.counts, dataset.celltype, dataset.tumor.grade, dataset.treatment,
     dataset.condition_group, size.factor,
     file = 'cancerous_human_T_cell_data.RData')
load('cancerous_human_T_cell_data.RData')

# Run SAVER for manifold fitting
# This step is typically computationally heavy and takes at least several hours for >= 10,000 cells
# We recommend running with at least 128 GB of RAM and saving the manifold separately
dataset.saver <- saver(dataset.counts, ncores = 4)

# Load and save pre-computed manifold
# save(dataset.saver, file = "cancerous_human_T_cell_manifold.RData")
# load("cancerous_human_T_cell_manifold.RData")

rownames(dataset.saver$estimate) <- rownames(dataset.counts)
colnames(dataset.saver$estimate) <- colnames(dataset.counts)
rownames(dataset.saver$mu.out) <- rownames(dataset.counts)
colnames(dataset.saver$mu.out) <- colnames(dataset.counts)

######################################################
# Find dispersion model with highest likelihood for each gene
dataset.saver.mu <- dataset.saver$mu.out
dataset.saver.var.models <- get_var_model(dataset.counts, dataset.saver.mu, size.factor)
# save(dataset.saver.var.models, file = "cancerous_human_T_cell_variance_models_SAVER.RData")
# load("cancerous_human_T_cell_variance_models_SAVER.RData")

######################################################
# Compute gene-level and cell-level deviation
dataset.celltype.levels <- levels(factor(dataset.celltype))
dataset.celltype <- as.integer(factor(dataset.celltype, levels = dataset.celltype.levels))
dataset.condition.levels <- c('Normal', 'II', 'III', 'IV', 'IV (treated)')
dataset.condition <- as.integer(factor(dataset.condition_group, levels = dataset.condition.levels))

# Parallel maps from condition index to grade/treatment label for output
condition.tumor.grade <- c('Normal', 'II', 'III', 'IV', 'IV')
condition.treatment <- c('Untreated', 'Untreated', 'Untreated', 'Untreated', 'Treated')

all_temp_gene_level <- data.frame()
all_temp_cell_level <- data.frame()

for (i in min(dataset.celltype):max(dataset.celltype)) {
  for (j in min(dataset.condition, na.rm = TRUE):max(dataset.condition, na.rm = TRUE)) {
    index <- rep(0, ncol(dataset.counts))
    for (k in 1:ncol(dataset.counts)) {
      # Create an index for each subset in the conditions
      index[k] <- dataset.celltype[k] == i & !is.na(dataset.condition[k]) & dataset.condition[k] == j
    }
    if (sum(index) > 9) { # Only analyze strata with >= 10 cells
      # Subset counts and normalize to match original manifold fitting
      dataset.celltype.counts.condition.norm <- sweep(dataset.counts[, index == 1], 2, size.factor[index == 1], '/')
      # Subset the manifold
      dataset.celltype.mu.condition <- dataset.saver$mu.out[rownames(dataset.saver$mu.out) %in% rownames(dataset.counts), index == 1]
      # Run SAVER on the cell type and condition specific stratum
      dataset.saver.celltype.condition <- saver(dataset.celltype.counts.condition.norm, mu = dataset.celltype.mu.condition, ncores = 4)
      assign(paste0("data.size",i,j), sum(index))
      assign(paste0("cell.barcode",i,j), colnames(dataset.celltype.counts.condition.norm))
      assign(paste0("data.disp.a",i,j), dataset.saver.celltype.condition$a)
      assign(paste0("data.disp.b",i,j), dataset.saver.celltype.condition$b)
      assign(paste0("data.disp.k",i,j), dataset.saver.celltype.condition$k)
      assign(paste0("estimate",i,j), dataset.saver.celltype.condition$estimate)
      assign(paste0("mu.out",i,j), dataset.saver.celltype.condition$mu.out)
    } else {
      assign(paste0("data.size",i,j), 0)
      assign(paste0("cell.barcode",i,j), vector())
      assign(paste0("data.disp.a",i,j), vector())
      assign(paste0("data.disp.b",i,j), vector())
      assign(paste0("data.disp.k",i,j), vector())
    }
  }
}

# Gene-level deviation
if (max(dataset.condition, na.rm = TRUE) > 1) {
  for (i in min(dataset.celltype):max(dataset.celltype)) {
    for (j in min(dataset.condition, na.rm = TRUE):max(dataset.condition, na.rm = TRUE)) {
      if (get(paste0("data.size",i,j)) != 0) {
        # Initialize variables
        Gene <- rownames(dataset.counts)
        Condition <- rep(dataset.condition.levels[j], nrow(dataset.counts))
        Tumor.Grade <- rep(condition.tumor.grade[j], nrow(dataset.counts))
        Treatment <- rep(condition.treatment[j], nrow(dataset.counts))
        Var_model <- names(dataset.saver.var.models)
        Dispersion <- dataset.saver.var.models
        Dispersion_cCV <- 1/get(paste0("data.disp.a",i,j))
        Dispersion_cFF <- 1/get(paste0("data.disp.b",i,j))
        Dispersion_cVar <- get(paste0("data.disp.k",i,j))
        Var_avg <- vector()
        Gene_level_deviation <- vector()
        # Calculate gene-level deviation
        estimate <- get(paste0("estimate",i,j))
        mu.out <- get(paste0("mu.out",i,j))

        for (g in 1:nrow(estimate)) {
          delta <- numeric(ncol(estimate))
          nu <- numeric(ncol(estimate))
          model <- names(dataset.saver.var.models)[g]
          for (c in 1:ncol(estimate)) {
            if (model == 'cCV') {  # cCV model
              nu_c <- mu.out[g, c]^2 / dataset.saver.var.models[g]
            } else if (model == 'cFF') {  # cFF model
              nu_c <- mu.out[g, c] / dataset.saver.var.models[g]
            } else {  # cVar model
              nu_c <- dataset.saver.var.models[g]
            }
            # Compute delta for gene g and cell c
            nu[c] <- nu_c
            delta[c] <- (estimate[g, c] - mu.out[g, c])^2 / nu_c
          }
          # Average delta across all cells to get gene-level deviation
          Var_avg <- append(Var_avg, mean(nu))
          Gene_level_deviation <- append(Gene_level_deviation, mean(delta))
        }
        # Combine data from all cell types
        temp <- data.frame(Gene, Condition, Tumor.Grade, Treatment, Dispersion_cCV, Dispersion_cFF, Dispersion_cVar, Var_model, Dispersion, Var_avg, Gene_level_deviation)
        temp$cell_type <- dataset.celltype.levels[i]
        all_temp_gene_level <- rbind(all_temp_gene_level, temp)
      }
    }
  }

  # Cell-level deviation
  for (i in min(dataset.celltype):max(dataset.celltype)) {
    for (j in min(dataset.condition, na.rm = TRUE):max(dataset.condition, na.rm = TRUE)) {
      if (get(paste0("data.size",i,j)) != 0) {
        # Initialize variables
        Cell_barcode <- get(paste0("cell.barcode",i,j))
        Condition <- rep(dataset.condition.levels[j], length(get(paste0("cell.barcode",i,j))))
        Tumor.Grade <- rep(condition.tumor.grade[j], length(get(paste0("cell.barcode",i,j))))
        Treatment <- rep(condition.treatment[j], length(get(paste0("cell.barcode",i,j))))
        Var_avg <- vector()
        Cell_level_deviation <- vector()
        # Calculate cell-level deviation
        estimate <- get(paste0("estimate",i,j))
        mu.out <- get(paste0("mu.out",i,j))

        for (c in 1:ncol(estimate)) {
          delta <- numeric(nrow(estimate))
          nu <- numeric(nrow(estimate))
          for (g in 1:nrow(estimate)) {
            model <- names(dataset.saver.var.models)[g]
            if (model == 'cCV') {  # cCV model
              nu_g <- mu.out[g, c]^2 / dataset.saver.var.models[g]
            } else if (model == 'cFF') {  # cFF model
              nu_g <- mu.out[g, c] / dataset.saver.var.models[g]
            } else {  # cVar model
              nu_g <- dataset.saver.var.models[g]
            }
            # Compute delta for gene g and cell c
            nu[g] <- nu_g
            delta[g] <- (estimate[g, c] - mu.out[g, c])^2 / nu_g
          }
          # Average delta across all cells to get gene-level deviation
          Var_avg <- append(Var_avg, mean(nu, na.rm = TRUE))
          Cell_level_deviation <- append(Cell_level_deviation, mean(delta, na.rm = TRUE))
        }
        # Combine data from all cell types
        temp <- data.frame(Cell_barcode, Condition, Tumor.Grade, Treatment, Cell_level_deviation)
        temp$cell_type <- dataset.celltype.levels[i]
        all_temp_cell_level <- rbind(all_temp_cell_level, temp)
      }
    }
  }
}

write.csv(all_temp_gene_level, "cancerous_human_T_cell_estimated_dispersion_SAVER.csv", row.names = FALSE)
write.csv(all_temp_cell_level, "cancerous_human_T_cell_cellular_dispersion_SAVER.csv", row.names = FALSE)

######################################################
