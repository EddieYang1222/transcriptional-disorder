# Transcriptional dyscoordination analysis for aging human kidney (tubular cells)

# Online links
# https://www.science.org/doi/10.1126/sciadv.adg8287
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211786

library(Seurat)
library(dplyr)
library(Matrix)
library(pbapply)
library(SAVER)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Load and prepare data
data_dir <- 'S:/Penn Dropbox/Eddie Yang/Aging/Data/Human_Kidney'
obj <- readRDS(file.path(data_dir, 'multi_rna.rds'))

# Join clinical metadata
meta <- read.csv(file.path(data_dir, 'clinical_meta.csv'))
colnames(meta)[1] <- 'library_id'
meta_idx <- match(sub("_multi$", "", obj@meta.data$sample), meta$library_id)
obj@meta.data$tissue_category <- meta$tissue_category[meta_idx]
obj@meta.data$age_group <- meta$age_group[meta_idx]
obj@meta.data$age <- meta$age[meta_idx]

# Subset to Control tubular cells
ctrl <- subset(obj, subset = tissue_category == 'Control')
tubular_types <- c('PST', 'PCT', 'PT_VCAM1', 'PT_PROM1', 'TAL', 'DCT1', 'DCT2_PC', 'ICA', 'ICB')
target <- subset(ctrl, subset = celltype %in% tubular_types)

# Create age bins
age_group_col <- target@meta.data$age_group
age_binned <- ifelse(age_group_col %in% c('20-29', '30-39'), '20-39',
              ifelse(age_group_col == '40-49', '40-49',
              ifelse(age_group_col == '50-59', '50-59',
              ifelse(age_group_col %in% c('60-69', '70-79', '80-89', '90-99'), '60+', NA))))
age_bin_levels <- c('20-39', '40-49', '50-59', '60+')
target@meta.data$age_bin <- factor(age_binned, levels = age_bin_levels)
target <- target[, !is.na(target@meta.data$age_bin)]

# Extract counts and metadata from filtered object
DefaultAssay(target) <- 'RNA'
dataset.counts.sparse <- GetAssayData(target, assay = 'RNA', layer = 'counts')
dataset.celltype <- target@meta.data$celltype
dataset.celltype.levels <- tubular_types
dataset.age <- target@meta.data$age_bin
dataset.age.levels <- age_bin_levels

# Remove unannotated cells
keep.cells <- !is.na(dataset.celltype) & !is.na(dataset.age)
dataset.counts.sparse <- dataset.counts.sparse[, keep.cells]
dataset.celltype <- dataset.celltype[keep.cells]
dataset.age <- dataset.age[keep.cells]

# Apply QC filters
keep.cells.qc <- Matrix::colSums(dataset.counts.sparse > 0) >= 250 &
                 Matrix::colSums(dataset.counts.sparse) >= 2500
keep.genes.qc <- Matrix::rowSums(dataset.counts.sparse > 0) >= 10
dataset.counts.sparse <- dataset.counts.sparse[keep.genes.qc, keep.cells.qc]
dataset.celltype <- dataset.celltype[keep.cells.qc]
dataset.age <- dataset.age[keep.cells.qc]
dataset.counts <- as.matrix(dataset.counts.sparse)
size.factor <- colSums(dataset.counts) / mean(colSums(dataset.counts))

rm(dataset.counts.sparse, obj, ctrl, target, meta)
save(dataset.counts, dataset.celltype, dataset.celltype.levels,
     dataset.age, dataset.age.levels, size.factor,
     file = 'aging_human_kidney_data.RData')
load('aging_human_kidney_data.RData')

# Run SAVER for manifold fitting
# This step is typically computationally heavy and takes at least several hours for >= 10,000 cells
# We recommend running with at least 128 GB of RAM and saving the manifold separately
dataset.saver <- saver(dataset.counts, ncores = 8)

# Load and save pre-computed manifold
# save(dataset.saver, file = "aging_human_kidney_manifold.RData")
# load("aging_human_kidney_manifold.RData")

rownames(dataset.saver$estimate) <- rownames(dataset.counts)
colnames(dataset.saver$estimate) <- colnames(dataset.counts)
rownames(dataset.saver$mu.out) <- rownames(dataset.counts)
colnames(dataset.saver$mu.out) <- colnames(dataset.counts)

######################################################
# Find dispersion model with highest likelihood for each gene
dataset.saver.mu <- dataset.saver$mu.out
dataset.saver.var.models <- get_var_model(dataset.counts, dataset.saver.mu, size.factor)
# save(dataset.saver.var.models, file = "aging_human_kidney_variance_models_SAVER.RData")
# load("aging_human_kidney_variance_models_SAVER.RData")

######################################################
# Compute gene-level and cell-level deviation
dataset.celltype <- as.integer(factor(dataset.celltype, levels = dataset.celltype.levels))
dataset.age <- as.integer(dataset.age)

all_temp_gene_level <- data.frame()
all_temp_cell_level <- data.frame()

for (i in min(dataset.celltype):max(dataset.celltype)) {
  for (j in min(dataset.age, na.rm = TRUE):max(dataset.age, na.rm = TRUE)) {
    index <- rep(0, ncol(dataset.counts))
    for (k in 1:ncol(dataset.counts)) {
      # Create an index for each subset in the conditions
      index[k] <- dataset.celltype[k] == i & !is.na(dataset.age[k]) & dataset.age[k] == j
    }
    if (sum(index) > 9) { # Only analyze strata with >= 10 cells
      # Subset counts and normalize to match original manifold fitting
      dataset.celltype.counts.age.norm <- sweep(dataset.counts[, index == 1], 2, size.factor[index == 1], '/')
      # Subset the manifold
      dataset.celltype.mu.age <- dataset.saver$mu.out[rownames(dataset.saver$mu.out) %in% rownames(dataset.counts), index == 1]
      # Run SAVER on the cell type and age specific stratum
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
if (max(dataset.age, na.rm = TRUE) > 1) {
  for (i in min(dataset.celltype):max(dataset.celltype)) {
    for (j in min(dataset.age, na.rm = TRUE):max(dataset.age, na.rm = TRUE)) {
      if (get(paste0("data.size",i,j)) != 0) {
        # Initialize variables
        Gene <- rownames(dataset.counts)
        Age <- rep(dataset.age.levels[j], nrow(dataset.counts))
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
        temp <- data.frame(Gene, Age, Dispersion_cCV, Dispersion_cFF, Dispersion_cVar, Var_model, Dispersion, Var_avg, Gene_level_deviation)
        temp$cell_type <- dataset.celltype.levels[i]
        all_temp_gene_level <- rbind(all_temp_gene_level, temp)
      }
    }
  }

  # Cell-level deviation
  for (i in min(dataset.celltype):max(dataset.celltype)) {
    for (j in min(dataset.age, na.rm = TRUE):max(dataset.age, na.rm = TRUE)) {
      if (get(paste0("data.size",i,j)) != 0) {
        # Initialize variables
        Cell_barcode <- get(paste0("cell.barcode",i,j))
        Age <- rep(dataset.age.levels[j], length(get(paste0("cell.barcode",i,j))))
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
        temp <- data.frame(Cell_barcode, Age, Cell_level_deviation)
        temp$cell_type <- dataset.celltype.levels[i]
        all_temp_cell_level <- rbind(all_temp_cell_level, temp)
      }
    }
  }
}

write.csv(all_temp_gene_level, "aging_human_kidney_estimated_dispersion_SAVER.csv", row.names = FALSE)
write.csv(all_temp_cell_level, "aging_human_kidney_cellular_dispersion_SAVER.csv",  row.names = FALSE)

######################################################