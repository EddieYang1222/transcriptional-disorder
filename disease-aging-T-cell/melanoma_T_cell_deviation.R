# Transcriptional dyscoordination in ICB-treated melanoma CD8 T cells (Melanoma)

# Online links
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272993
# Wang et al. 2024, longitudinal CD8 T cell response to combination immune-checkpoint blockade

# Set up working directory
# setwd("path/to/working/directory")

library(Seurat)
library(Matrix)
library(pbapply)
library(SAVER)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Load and prepare data
# Cohort: combination therapy (ipilimumab + nivolumab) arm, 8 of 9 patients
# (the 29-year-old outlier 17-2864 is dropped for age consistency). Strata =
# metaclusters_exh x timepoint over the four-timepoint core schedule (Weeks 0/3/6/9).
TIMEPOINTS <- c('Baseline', 'Follow Up 1', 'Follow Up 2', 'Follow Up 3')
KEEP_PATIENTS <- c('18-3313', '17-2796', '17-2631', '18-3253', '14-1319',
                   '17-2865', '15-1471', '17-3041')

data_dir <- 'path/to/GSE272993'
meta <- read.csv(file.path(data_dir, 'wang_melanoma_metadata.csv'), stringsAsFactors = FALSE)
rownames(meta) <- meta$cell_id

keep <- meta$treatment    == 'Combination' &
        meta$patient_alias %in% KEEP_PATIENTS &
        meta$timepoint     %in% TIMEPOINTS    &
        !is.na(meta$metaclusters_exh)
meta_sub <- meta[keep, , drop = FALSE]

# Extract counts for the picked cells from the labeled Seurat object
obj <- readRDS(file.path(data_dir, 'GSE272993_cd8_nn_labeled_FINAL.RDS'))
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
dataset.saver <- saver(dataset.counts, ncores = 4)

# Load and save pre-computed manifold
# save(dataset.saver, file = "melanoma_T_cell_manifold.RData")
# load("melanoma_T_cell_manifold.RData")

rownames(dataset.saver$estimate) <- rownames(counts)
colnames(dataset.saver$estimate) <- colnames(counts)
rownames(dataset.saver$mu.out)   <- rownames(counts)
colnames(dataset.saver$mu.out)   <- colnames(counts)

######################################################
# Find dispersion model with highest likelihood for each gene
dataset.saver.var.models <- get_var_model(dataset.counts, dataset.saver$mu.out, size.factor)
# save(dataset.saver.var.models, file = "melanoma_T_cell_variance_models_SAVER.RData")
# load("melanoma_T_cell_variance_models_SAVER.RData")

######################################################
# Compute gene-level and cell-level deviation per stratum (metaclusters_exh x timepoint)
cell_types <- sort(unique(meta_sub$metaclusters_exh))
tps        <- TIMEPOINTS[TIMEPOINTS %in% unique(meta_sub$timepoint)]

all_temp_gene_level <- data.frame()
all_temp_cell_level <- data.frame()

for (ct in cell_types) {
  for (tp in tps) {
    idx <- which(meta_sub$metaclusters_exh == ct & meta_sub$timepoint == tp)
    if (length(idx) < 10) next  # Only analyze strata with >= 10 cells

    # Subset counts, normalize to match the manifold, and refit SAVER using the global mu as prior
    cnt_norm <- sweep(counts[, idx, drop = FALSE], 2, size.factor[idx], '/')
    mu_sub   <- dataset.saver$mu.out[rownames(counts), idx, drop = FALSE]
    saver_sub <- saver(as.matrix(cnt_norm), mu = mu_sub, ncores = 4)
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
      Condition = paste0(ct, '_', gsub(' ', '_', tp)),
      cell_type = ct,
      Timepoint = tp,
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
      Timepoint = tp,
      patient_alias = meta_sub$patient_alias[idx],
      age = meta_sub$age[idx],
      gender = meta_sub$gender[idx],
      response = meta_sub$response[idx],
      Cell_level_deviation
    )
    all_temp_cell_level <- rbind(all_temp_cell_level, cell_df)

    rm(saver_sub, estimate, mu.out, cnt_norm, mu_sub)
    gc()
  }
}

write.csv(all_temp_gene_level, "melanoma_T_cell_estimated_dispersion_SAVER.csv", row.names = FALSE)
write.csv(all_temp_cell_level, "melanoma_T_cell_cellular_dispersion_SAVER.csv",  row.names = FALSE)

######################################################
