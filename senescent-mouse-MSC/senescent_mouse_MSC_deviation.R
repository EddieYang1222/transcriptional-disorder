# Transcriptional dyscoordination analysis for senescent mouse MSC (SC and FAP)

# Online links
# https://www.cell.com/iscience/fulltext/S2589-0042(22)00118-3
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169531

library(Seurat)
library(dplyr)
library(Matrix)
library(pbapply)
library(SAVER)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Cell types: SC (stromal cells),  FAP (fibro-adipogenic progenitors)
# Conditions: PBS (vehicle control), DOXO (doxorubicin-induced senescence), DABT (dasatinib + ABT-263 senolytic clearance)

# Load and prepare data
for (celltype_name in c('SC', 'FAP')) {
  obj <- readRDS(paste0('mouse_MSC_', celltype_name, '_processed.rds'))
  DefaultAssay(obj) <- 'RNA'
  counts.sparse <- GetAssayData(obj, assay = 'RNA', layer = 'counts')
  condition <- as.character(obj@meta.data$condition)
  
  # Apply QC filters
  keep.cells.qc <- colSums(counts.sparse > 0) >= 250 & colSums(counts.sparse) >= 2500 
  keep.genes.qc <- rowSums(counts.sparse > 0) >= 5
  counts.sparse <- counts.sparse[keep.genes.qc, keep.cells.qc]
  condition <- condition[keep.cells.qc]
  
  counts <- as.matrix(counts.sparse)
  rm(counts.sparse)
  size.factor <- colSums(counts) / mean(colSums(counts))

  save(counts, condition, size.factor,
       file = paste0('mouse_MSC_', celltype_name, '_data.RData'))
  }

# Run SAVER for manifold fitting
# As suggested in the original experiment, SC and FAP are processed separately
# This step is typically computationally heavy and takes at least several hours for >= 10,000 cells
# We recommend running with at least 64 GB of RAM and saving the manifold separately
for (celltype_name in c('SC', 'FAP')) {
  load(paste0('mouse_MSC_', celltype_name, '_data.RData'))
  dataset.saver <- saver(counts, size.factor = size.factor, ncores = 4)
  save(dataset.saver, file = paste0('mouse_MSC_', celltype_name, '_SAVER.RData'))
}

######################################################
# Find dispersion model with highest likelihood for each gene
for (celltype_name in c('SC', 'FAP')) {
  load(paste0('mouse_MSC_', celltype_name, '_data.RData'))
  load(paste0('mouse_MSC_', celltype_name, '_SAVER.RData'))
  dataset.saver.var.models <- get_var_model(counts, dataset.saver$mu.out, size.factor)
  save(dataset.saver.var.models, file = paste0('mouse_MSC_', celltype_name, '_variance_models_SAVER.RData'))
}

######################################################
# Compute gene-level and cell-level deviation
# Process each cell type (SC, FAP) independently
condition_levels <- c('PBS', 'DOXO', 'DABT')

for (celltype_name in c('SC', 'FAP')) {
  load(paste0('mouse_MSC_', celltype_name, '_data.RData'))
  load(paste0('mouse_MSC_', celltype_name, '_SAVER.RData'))
  load(paste0('mouse_MSC_', celltype_name, '_variance_models_SAVER.RData'))

  rownames(dataset.saver$estimate) <- rownames(counts)
  colnames(dataset.saver$estimate) <- colnames(counts)
  rownames(dataset.saver$mu.out) <- rownames(counts)
  colnames(dataset.saver$mu.out) <- colnames(counts)

  dataset.counts <- counts
  dataset.condition <- as.integer(factor(condition, levels = condition_levels))

  all_temp_gene_level <- data.frame()
  all_temp_cell_level <- data.frame()

  for (j in min(dataset.condition):max(dataset.condition)) {
    index <- rep(0, ncol(dataset.counts))
    for (k in 1:ncol(dataset.counts)) {
      # Create an index for each subset in the conditions
      index[k] <- dataset.condition[k] == j
    }
    if (sum(index) > 9) { # Only analyze strata with >= 10 cells
      # Subset counts and normalize to match original manifold fitting
      dataset.celltype.counts.age.norm <- sweep(dataset.counts[, index == 1], 2, size.factor[index == 1], '/')
      # Subset the manifold
      dataset.celltype.mu.age <- dataset.saver$mu.out[rownames(dataset.saver$mu.out) %in% rownames(dataset.counts), index == 1]
      # Run SAVER on the cell type and condition specific stratum
      dataset.saver.celltype.age <- saver(dataset.celltype.counts.age.norm, mu = dataset.celltype.mu.age, ncores = 4)
      assign(paste0("data.size",j), sum(index))
      assign(paste0("cell.barcode",j), colnames(dataset.celltype.counts.age.norm))
      assign(paste0("data.disp.a",j), dataset.saver.celltype.age$a)
      assign(paste0("data.disp.b",j), dataset.saver.celltype.age$b)
      assign(paste0("data.disp.k",j), dataset.saver.celltype.age$k)
      assign(paste0("estimate", j), dataset.saver.celltype.age$estimate)
      assign(paste0("mu.out", j), dataset.saver.celltype.age$mu.out)
    } else {
      assign(paste0("data.size",j), 0)
      assign(paste0("cell.barcode",j), vector())
      assign(paste0("data.disp.a",j), vector())
      assign(paste0("data.disp.b",j), vector())
      assign(paste0("data.disp.k",j), vector())
    }
  }

  # Gene-level deviation
  if (max(dataset.condition) > 1) {
    for (j in min(dataset.condition):max(dataset.condition)) {
      if (get(paste0("data.size",j)) != 0) {
        # Initialize variables
        Gene <- rownames(dataset.counts)
        Condition <- rep(condition_levels[j], nrow(dataset.counts))
        Var_model <- names(dataset.saver.var.models)
        Dispersion <- dataset.saver.var.models
        Dispersion_cCV <- 1/get(paste0("data.disp.a",j))
        Dispersion_cFF <- 1/get(paste0("data.disp.b",j))
        Dispersion_cVar <- get(paste0("data.disp.k",j))
        Var_avg <- vector()
        Gene_level_deviation <- vector()
        # Calculate gene-level deviation
        estimate <- get(paste0("estimate", j))
        mu.out <- get(paste0("mu.out", j))

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
        temp <- data.frame(Gene, Condition, Dispersion_cCV, Dispersion_cFF, Dispersion_cVar, Var_model, Dispersion, Var_avg, Gene_level_deviation)
        temp$cell_type <- celltype_name
        all_temp_gene_level <- rbind(all_temp_gene_level, temp)
      }
    }

    # Cell-level deviation
    for (j in min(dataset.condition):max(dataset.condition)) {
      if (get(paste0("data.size",j)) != 0) {
        # Initialize variables
        Cell_barcode <- get(paste0("cell.barcode",j))
        Condition <- rep(condition_levels[j], length(get(paste0("cell.barcode",j))))
        Var_avg <- vector()
        Cell_level_deviation <- vector()
        # Calculate cell-level deviation
        estimate <- get(paste0("estimate", j))
        mu.out <- get(paste0("mu.out", j))

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
        temp <- data.frame(Cell_barcode, Condition, Cell_level_deviation)
        temp$cell_type <- celltype_name
        all_temp_cell_level <- rbind(all_temp_cell_level, temp)
      }
    }
  }

  write.csv(all_temp_gene_level, paste0('senescent_mouse_MSC_', celltype_name, '_estimated_dispersion_SAVER.csv'), row.names = FALSE)
  write.csv(all_temp_cell_level, paste0('senescent_mouse_MSC_', celltype_name, '_cellular_dispersion_SAVER.csv'), row.names = FALSE)
}

######################################################