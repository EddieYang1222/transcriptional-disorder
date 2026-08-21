# Transcriptional dyscoordination in healthy-aging PBMC T cells (Wang aging)

# Online links
# https://www.synapse.org/#!Synapse:syn61609846
# Wang et al. 2025, single-cell lifespan PBMC atlas (raw FASTQ at GSA-Human HRA009014)

# Set up working directory
# setwd("path/to/working/directory")

library(SAVER)
library(Matrix)
library(pbapply)
library(Seurat)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Load and prepare data
# Subset: all 28 adult samples (30-40y / 50-60y / 70-80y / 90y) restricted to 12
# T cell subtypes (5 CD4 + 6 CD8 + 1 gdT). Strata = secondary_type x eight_group.
T_SUBTYPES <- c('CD4_Naive_CCR7','CD4_TCM_AQP3','CD4_TEM_ANXA1','CD4_TEM_GNLY','CD4_Treg_FOXP3',
                'CD8_Naive_LEF1','CD8_TCM_HAVCR2','CD8_TEM_CMC1','CD8_TEM_GNLY',
                'CD8_TEM_ZNF683','CD8_MAIT_SLC4A10','gdT')
ADULT_BINS <- c('30-40y','50-60y','70-80y','90y')

data_dir <- 'path/to/Wang_aging'
meta <- read.csv(file.path(data_dir, 'wang_aging_metadata.csv'), stringsAsFactors = FALSE)
rownames(meta) <- meta$cellid

keep_cell <- meta$secondary_type %in% T_SUBTYPES &
             meta$eight_group    %in% ADULT_BINS
meta_sub <- meta[keep_cell, , drop = FALSE]

# Extract counts for the picked T cells from the labeled Seurat object
obj <- readRDS(file.path(data_dir, 'scRNA-seqProcessedLabelledObject.rds'))
if (inherits(obj, 'Seurat')) {
  counts <- GetAssayData(obj, assay = 'RNA', layer = 'counts')
  if (is.null(counts) || ncol(counts) == 0)
    counts <- GetAssayData(obj, assay = 'RNA', slot = 'counts')
} else if (inherits(obj, c('Assay5', 'Assay'))) {
  counts <- tryCatch(LayerData(obj, layer = 'counts'), error = function(e) obj@layers$counts)
} else {
  counts <- obj
}
rm(obj)
gc()

common <- intersect(colnames(counts), rownames(meta_sub))
counts <- counts[, common, drop = FALSE]
meta_sub <- meta_sub[common, , drop = FALSE]

# Apply QC filters (per-cell: features >= 200, counts >= 500; per-gene: cells >= 10)
cell_pass <- Matrix::colSums(counts > 0) >= 200 & Matrix::colSums(counts) >= 500
counts <- counts[, cell_pass, drop = FALSE]
meta_sub <- meta_sub[cell_pass, , drop = FALSE]
counts <- counts[Matrix::rowSums(counts > 0) >= 10, ]
size.factor <- colSums(counts) / mean(colSums(counts))

######################################################
# Run SAVER for manifold fitting
# This step is typically computationally heavy and takes at least several hours for >= 10,000 cells
# We recommend running with at least 128 GB of RAM and saving the manifold separately
dataset.counts <- as.matrix(counts)
dataset.saver <- saver(dataset.counts, ncores = 8)

# Load and save pre-computed manifold
# save(dataset.saver, file = "Wang_aging_T_cell_manifold.RData")
# load("Wang_aging_T_cell_manifold.RData")

rownames(dataset.saver$estimate) <- rownames(counts)
colnames(dataset.saver$estimate) <- colnames(counts)
rownames(dataset.saver$mu.out)   <- rownames(counts)
colnames(dataset.saver$mu.out)   <- colnames(counts)

######################################################
# Find dispersion model with highest likelihood for each gene
dataset.saver.var.models <- get_var_model(dataset.counts, dataset.saver$mu.out, size.factor)
# save(dataset.saver.var.models, file = "Wang_aging_T_cell_variance_models_SAVER.RData")
# load("Wang_aging_T_cell_variance_models_SAVER.RData")

######################################################
# Compute gene-level and cell-level deviation per stratum (secondary_type x eight_group)
cell_types <- intersect(T_SUBTYPES, sort(unique(meta_sub$secondary_type)))
bins       <- intersect(ADULT_BINS, sort(unique(meta_sub$eight_group)))

all_temp_gene_level <- data.frame()
all_temp_cell_level <- data.frame()

for (ct in cell_types) {
  for (ag in bins) {
    idx <- which(meta_sub$secondary_type == ct & meta_sub$eight_group == ag)
    if (length(idx) < 10) next  # Only analyze strata with >= 10 cells

    # Subset counts, normalize to match the manifold, and refit SAVER using the global mu as prior
    cnt_norm <- sweep(counts[, idx, drop = FALSE], 2, size.factor[idx], '/')
    mu_sub   <- dataset.saver$mu.out[rownames(counts), idx, drop = FALSE]
    saver_sub <- saver(as.matrix(cnt_norm), mu = mu_sub, ncores = 8)
    estimate <- saver_sub$estimate
    mu.out   <- saver_sub$mu.out

    # Gene-level deviation
    Gene_level_deviation <- numeric(nrow(counts))
    Var_avg              <- numeric(nrow(counts))
    for (g in seq_len(nrow(estimate))) {
      model <- names(dataset.saver.var.models)[g]
      delta <- numeric(ncol(estimate))
      nu <- numeric(ncol(estimate))
      for (c in seq_len(ncol(estimate))) {
        nu_c <- if (model == 'cCV') mu.out[g, c]^2 / dataset.saver.var.models[g]
                else if (model == 'cFF') mu.out[g, c] / dataset.saver.var.models[g]
                else dataset.saver.var.models[g]
        nu[c] <- nu_c
        delta[c] <- (estimate[g, c] - mu.out[g, c])^2 / nu_c
      }
      Var_avg[g] <- mean(nu)
      Gene_level_deviation[g] <- mean(delta)
    }
    gene_df <- data.frame(
      Gene = rownames(counts),
      Condition = paste0(ct, '_', ag),
      cell_type = ct,
      Age_group = ag,
      Dispersion_cCV = 1 / saver_sub$a,
      Dispersion_cFF = 1 / saver_sub$b,
      Dispersion_cVar = saver_sub$k,
      Var_model = names(dataset.saver.var.models),
      Dispersion = dataset.saver.var.models,
      Var_avg,
      Gene_level_deviation
    )
    all_temp_gene_level <- rbind(all_temp_gene_level, gene_df)

    # Cell-level deviation
    Cell_level_deviation <- numeric(ncol(estimate))
    for (c in seq_len(ncol(estimate))) {
      delta <- numeric(nrow(estimate))
      nu <- numeric(nrow(estimate))
      for (g in seq_len(nrow(estimate))) {
        model <- names(dataset.saver.var.models)[g]
        nu_g <- if (model == 'cCV') mu.out[g, c]^2 / dataset.saver.var.models[g]
                else if (model == 'cFF') mu.out[g, c] / dataset.saver.var.models[g]
                else dataset.saver.var.models[g]
        nu[g] <- nu_g
        delta[g] <- (estimate[g, c] - mu.out[g, c])^2 / nu_g
      }
      Cell_level_deviation[c] <- mean(delta, na.rm = TRUE)
    }
    cell_df <- data.frame(
      Cell_barcode = colnames(cnt_norm),
      nCount_RNA = Matrix::colSums(counts[, idx, drop = FALSE]),
      cell_type = ct,
      Age_group = ag,
      age = meta_sub$age[idx],
      sampleName = meta_sub$sampleName[idx],
      Cell_level_deviation
    )
    all_temp_cell_level <- rbind(all_temp_cell_level, cell_df)

    rm(saver_sub, estimate, mu.out, cnt_norm, mu_sub)
    gc()
  }
}

write.csv(all_temp_gene_level, "Wang_aging_T_cell_estimated_dispersion_SAVER.csv", row.names = FALSE)
write.csv(all_temp_cell_level, "Wang_aging_T_cell_cellular_dispersion_SAVER.csv",  row.names = FALSE)

######################################################
