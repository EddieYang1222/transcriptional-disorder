# Transcriptional dyscoordination analysis for aging rat kidney (PT cells only)

library(Seurat)
library(dplyr)
library(Matrix)
library(pbapply)
library(SAVER)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Load and prepare data
load("Parker_kidney_aging_PT_only.RData")

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
# save(dataset.saver, file = "aging_rat_kidney_PT_cell_manifold.RData")
# load("aging_rat_kidney_PT_cell_manifold.RData")

######################################################
# Find dispersion model with highest likelihood for each gene
size.factor <- colSums(dataset.counts) / mean(colSums(dataset.counts))
dataset.saver.mu <- dataset.saver$mu.out
dataset.saver.var.models <- get_var_model(dataset.counts, dataset.saver.mu, size.factor)
# save(dataset.saver.var.models, file = "aging_rat_kidney_PT_cell_variance_models_SAVER.RData")
# load("aging_rat_kidney_PT_cell_variance_models_SAVER.RData")

######################################################
# Compute gene-level and cell-level deviation
all_temp_gene_level <- data.frame()
all_temp_cell_level <- data.frame()

for (i in min(dataset.celltype):max(dataset.celltype)) {
  for (j in min(dataset.age):max(dataset.age)) {
    index <- rep(0, ncol(dataset.counts))
    for (k in 1:ncol(dataset.counts)) {
      # Create an index for each subset in the older ages
      index[k] <- dataset.celltype[k] == i && dataset.age[k] == j
    }
    if (sum(index) > 9) { # Only analyze strata with >= 10 cells
      # Subset counts and normalize to match original manifold fitting
      dataset.celltype.counts.age.norm <- sweep(dataset.counts[, index == 1], 2, size.factor[index == 1], "/")
      # Subset the manifold
      dataset.celltype.mu.age <- dataset.saver$mu.out[rownames(dataset.saver$mu.out) %in% rownames(dataset.counts), index == 1]
      # Run SAVER on the cell type and age specific stratum
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
  if (max(dataset.age) > 1) {
    for (j in min(dataset.age):max(dataset.age)) {
      if (get(paste0("data.size",j)) != 0) {
        # Initialize variables
        Gene <- rownames(dataset.counts)
        Age <- rep(dataset.age.levels[j], nrow(dataset.counts))
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
            if (model == "cCV") {  # cCV model
              nu_c <- mu.out[g, c]^2 / dataset.saver.var.models[g]
            } else if (model == "cFF") {  # cFF model
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
        temp <- data.frame(Gene, Age, Dispersion_cCV, Dispersion_cFF, Dispersion_cVar, Var_model, Dispersion, Var_avg, Gene_level_deviation)
        temp$cell_type <- dataset.celltype.levels[i]
        all_temp_gene_level <- rbind(all_temp_gene_level, temp)
      }
    }

    # Cell-level deviation
    for (j in min(dataset.age):max(dataset.age)) {
      if (get(paste0("data.size",j)) != 0) {
        # Initialize variables
        Cell_barcode <- get(paste0("cell.barcode",j))
        Age <- rep(dataset.age.levels[j], length(get(paste0("cell.barcode",j))))
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
            if (model == "cCV") {  # cCV model
              nu_g <- mu.out[g, c]^2 / dataset.saver.var.models[g]
            } else if (model == "cFF") {  # cFF model
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
        temp <- data.frame(Cell_barcode, Age, Cell_level_deviation)
        temp$cell_type <- dataset.celltype.levels[i]
        all_temp_cell_level <- rbind(all_temp_cell_level, temp)
      }
    }
  }
}

write.csv(all_temp_gene_level, "aging_rat_kidney_PT_cell_estimated_dispersion_SAVER.csv", row.names = FALSE)
write.csv(all_temp_cell_level, "aging_rat_kidney_PT_cell_cellular_dispersion_SAVER.csv",  row.names = FALSE)

######################################################
