# Transcriptional dyscoordination analysis for senescent mouse aortic cells

# Online links
# https://www.nature.com/articles/s43587-025-00889-z
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239591

# Set up working directory
# setwd("S:/Penn Dropbox/Eddie Yang/Aging/Scripts/Mouse_VSMC")

library(Seurat)
library(dplyr)
library(Matrix)
library(pbapply)
library(SAVER)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Load and prepare data
# Conditions: ND (normal diet), HFD (high-fat diet), ABT (HFD + ABT-737 senolytic)
aorta <- readRDS('mouse_VSMC_aorta_processed.rds')
DefaultAssay(aorta) <- 'RNA'

dataset.counts.sparse <- GetAssayData(aorta, assay = 'RNA', layer = 'counts')
dataset.celltype.raw <- as.character(aorta@meta.data$celltype)
dataset.condition.raw <- as.character(aorta@meta.data$condition)

# Apply QC filters
keep.cells.qc <- Matrix::colSums(dataset.counts.sparse > 0) >= 250 &
                 Matrix::colSums(dataset.counts.sparse) >= 2500
keep.genes.qc <- Matrix::rowSums(dataset.counts.sparse > 0) >= 5
dataset.counts.sparse <- dataset.counts.sparse[keep.genes.qc, keep.cells.qc]
dataset.celltype.raw <- dataset.celltype.raw[keep.cells.qc]
dataset.condition.raw <- dataset.condition.raw[keep.cells.qc]

dataset.counts <- as.matrix(dataset.counts.sparse)
rm(dataset.counts.sparse, aorta)

dataset.celltype.levels <- sort(unique(dataset.celltype.raw))
dataset.celltype <- as.integer(factor(dataset.celltype.raw, levels = dataset.celltype.levels))
condition_levels <- c('ND', 'HFD', 'ABT')
dataset.condition <- as.integer(factor(dataset.condition.raw, levels = condition_levels))

# Save and load input data
save(dataset.counts, dataset.celltype, dataset.celltype.levels,
     dataset.condition, condition_levels, file = 'mouse_aortic_cell_data.RData')
load('mouse_aortic_cell_data.RData')

# Remove genes with zero counts
message("The initial matrix size is ", nrow(dataset.counts), " genes and ", ncol(dataset.counts), " cells.")
dataset.counts <- dataset.counts[rowSums(dataset.counts) != 0, ]
message("After removing genes with zero counts, the new matrix size is ",
        nrow(dataset.counts), " genes and ", ncol(dataset.counts), " cells.")

# Run SAVER for manifold fitting
# This step is typically computationally heavy and takes at least several hours for >= 10,000 cells
# We recommend running with at least 64 GB of RAM and saving the manifold separately
dataset.saver <- saver(dataset.counts, ncores = 4)

# Load and save pre-computed manifold
# save(dataset.saver, file = 'mouse_aortic_cell_manifold.RData')
# load('mouse_aortic_cell_manifold.RData')

######################################################
# Additional pre-processing
size.factor <- colSums(dataset.counts) / mean(colSums(dataset.counts))

# Find dispersion model with highest likelihood for each gene
dataset.saver.mu <- dataset.saver$mu.out
dataset.saver.var.models <- get_var_model(dataset.counts, dataset.saver.mu, size.factor)
# save(dataset.saver.var.models, file = 'mouse_aortic_cell_variance_models_SAVER.RData')
# load('mouse_aortic_cell_variance_models_SAVER.RData')

######################################################
# Compute gene-level and cell-level deviation
all_temp_gene_level <- data.frame()
all_temp_cell_level <- data.frame()

for (i in min(dataset.celltype):max(dataset.celltype)) {
  for (j in min(dataset.condition):max(dataset.condition)) {
    index <- rep(0, ncol(dataset.counts))
    for (k in 1:ncol(dataset.counts)) {
      # Create an index for each subset in the conditions
      index[k] <- dataset.celltype[k] == i & dataset.condition[k] == j
    }
    if (sum(index) > 9) { # Only analyze strata with >= 10 cells
      # Subset counts and normalize to match original manifold fitting
      dataset.celltype.counts.age.norm <- sweep(dataset.counts[, index == 1], 2, size.factor[index == 1], '/')
      # Subset the manifold
      dataset.celltype.mu.age <- dataset.saver$mu.out[rownames(dataset.saver$mu.out) %in% rownames(dataset.counts), index == 1]
      # Run SAVER on the cell type and condition specific stratum
      dataset.saver.celltype.age <- saver(dataset.celltype.counts.age.norm, mu = dataset.celltype.mu.age, ncores = 4)
      assign(paste0("data.size",i,j), sum(index))
      assign(paste0("cell.barcode",i,j), colnames(dataset.celltype.counts.age.norm))
      assign(paste0("data.disp.a",i,j), dataset.saver.celltype.age$a)
      assign(paste0("data.disp.b",i,j), dataset.saver.celltype.age$b)
      assign(paste0("data.disp.k",i,j), dataset.saver.celltype.age$k)
      assign(paste0("estimate",i,j), dataset.saver.celltype.age$estimate)
      assign(paste0("mu.out",i,j), dataset.saver.celltype.age$mu.out)
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
if (max(dataset.condition) > 1) {
  for (i in min(dataset.celltype):max(dataset.celltype)) {
    for (j in min(dataset.condition):max(dataset.condition)) {
      if (get(paste0("data.size",i,j)) != 0) {
        # Initialize variables
        Gene <- rownames(dataset.counts)
        Condition <- rep(condition_levels[j], nrow(dataset.counts))
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
        temp <- data.frame(Gene, Condition, Dispersion_cCV, Dispersion_cFF, Dispersion_cVar, Var_model, Dispersion, Var_avg, Gene_level_deviation)
        temp$cell_type <- dataset.celltype.levels[i]
        all_temp_gene_level <- rbind(all_temp_gene_level, temp)
      }
    }
  }

  # Cell-level deviation
  for (i in min(dataset.celltype):max(dataset.celltype)) {
    for (j in min(dataset.condition):max(dataset.condition)) {
      if (get(paste0("data.size",i,j)) != 0) {
        # Initialize variables
        Cell_barcode <- get(paste0("cell.barcode",i,j))
        Condition <- rep(condition_levels[j], length(get(paste0("cell.barcode",i,j))))
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
        temp <- data.frame(Cell_barcode, Condition, Cell_level_deviation)
        temp$cell_type <- dataset.celltype.levels[i]
        all_temp_cell_level <- rbind(all_temp_cell_level, temp)
      }
    }
  }
}

write.csv(all_temp_gene_level, "senescent_mouse_aortic_cell_estimated_dispersion_SAVER.csv", row.names = FALSE)
write.csv(all_temp_cell_level, "senescent_mouse_aortic_cell_cellular_dispersion_SAVER.csv", row.names = FALSE)
