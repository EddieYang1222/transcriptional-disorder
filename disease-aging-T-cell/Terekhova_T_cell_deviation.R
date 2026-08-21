# Transcriptional dyscoordination in healthy-aging blood T cells (Terekhova)

# Online links
# https://www.synapse.org/#!Synapse:syn49637038
# Terekhova et al. 2023, single-cell atlas of healthy human blood across the lifespan

# Set up working directory
# setwd("path/to/working/directory")

library(SAVER)
library(Matrix)
library(pbapply)
library(zellkonverter)
library(SingleCellExperiment)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Load and prepare data
# Subset: 35 T cell types = 16 CD4 fine + 12 CD8 fine + 5 gd fine + 1 MAIT
# (collapsed) + 1 DN_T, over 20 tubes (4 donors per Age_group A-E). Strata =
# cell_type x Age_group.
data_dir <- 'path/to/Terekhova_atlas'
picks <- read.csv(file.path(data_dir, 'picked_samples_v8.csv'), stringsAsFactors = FALSE)
picked_tubes <- picks$Tube_id

# Load per-lineage fine metadata, filter to picks, prefix labels by lineage
load_lineage_meta <- function(path, lineage) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  colnames(df)[1] <- 'cell_id'
  df <- df[df$Tube_id %in% picked_tubes, , drop = FALSE]
  df$Lineage <- lineage
  df
}
cd4 <- load_lineage_meta(file.path(data_dir, 'cd4_metadata.csv'), 'CD4')
cd8 <- load_lineage_meta(file.path(data_dir, 'conventional_cd8_metadata.csv'), 'CD8')
gd  <- load_lineage_meta(file.path(data_dir, 'gd_t_cells_metadata.csv'), 'gd')

# Combined fine-label dataframe with cell_type = Lineage_<Cluster_names>
keep_cols <- c('cell_id','Tube_id','Donor_id','Age_group','Age','Sex','Batch',
               'nCount_RNA','nFeature_RNA','percent.mt','percent.ribo',
               'Lineage','Cluster_names')
fine_meta <- rbind(cd4[, keep_cols], cd8[, keep_cols], gd[, keep_cols])
fine_meta$cell_type <- paste0(fine_meta$Lineage, '_', fine_meta$Cluster_names)

# Add MAIT and DN_T cells from the top-level all_pbmcs metadata
ap <- read.csv(file.path(data_dir, 'all_pbmcs_metadata.csv'), stringsAsFactors = FALSE)
colnames(ap)[1] <- 'cell_id'
ap <- ap[ap$Tube_id %in% picked_tubes, , drop = FALSE]
mait_dn <- ap[ap$Cluster_names %in% c('MAIT cells','DN T cells'), , drop = FALSE]
mait_dn$Lineage   <- ifelse(mait_dn$Cluster_names == 'MAIT cells', 'MAIT', 'DN_T')
mait_dn$cell_type <- mait_dn$Lineage
keep_in_ap <- intersect(keep_cols, colnames(mait_dn))
mait_dn <- mait_dn[, c(keep_in_ap, 'cell_type')]
for (col in setdiff(colnames(fine_meta), colnames(mait_dn)))
  mait_dn[[col]] <- NA
mait_dn <- mait_dn[, colnames(fine_meta), drop = FALSE]

meta <- rbind(fine_meta, mait_dn)
meta <- meta[!duplicated(meta$cell_id), , drop = FALSE]
rownames(meta) <- meta$cell_id

# Load the raw-count h5ad on-disk and subset to the picked cells
sce <- readH5AD(file.path(data_dir, 'all_pbmcs_rna.h5ad'), use_hdf5 = TRUE, raw = FALSE,
                X_name = 'counts')
common <- intersect(colnames(sce), rownames(meta))
sce <- sce[, common]
meta <- meta[common, , drop = FALSE]
counts <- as(assay(sce, 1), 'CsparseMatrix')
rm(sce)
gc()

# Apply QC filters (per-cell: features >= 200, counts >= 500). At ~187k cells the
# usual >=10-cell gene floor is effectively no filter, so keep a gene only if it is
# detected in >= 1% of QC-passing cells: this removes rarely-detected genes whose
# manifold is poorly estimated and whose variance term nu explodes (the source of
# unstable, outlier-driven gene-level dyscoordination). The manifold is refit on the
# filtered gene set, not post-hoc subset.
DETECTION_FRAC <- 0.01
cell_pass <- Matrix::colSums(counts > 0) >= 200 & Matrix::colSums(counts) >= 500
counts <- counts[, cell_pass, drop = FALSE]
meta   <- meta[cell_pass, , drop = FALSE]
gene_min_cells <- max(10, ceiling(DETECTION_FRAC * ncol(counts)))
counts <- counts[Matrix::rowSums(counts > 0) >= gene_min_cells, ]
size.factor <- colSums(counts) / mean(colSums(counts))

######################################################
# Run SAVER for manifold fitting
# This step is typically computationally heavy and takes at least several hours for >= 10,000 cells
# We recommend running with at least 128 GB of RAM and saving the manifold separately
dataset.counts <- as.matrix(counts)
dataset.saver <- saver(dataset.counts, ncores = 8)

# Load and save pre-computed manifold
# save(dataset.saver, file = "Terekhova_T_cell_manifold.RData")
# load("Terekhova_T_cell_manifold.RData")

rownames(dataset.saver$estimate) <- rownames(counts)
colnames(dataset.saver$estimate) <- colnames(counts)
rownames(dataset.saver$mu.out)   <- rownames(counts)
colnames(dataset.saver$mu.out)   <- colnames(counts)

######################################################
# Find dispersion model with highest likelihood for each gene
dataset.saver.var.models <- get_var_model(dataset.counts, dataset.saver$mu.out, size.factor)
# save(dataset.saver.var.models, file = "Terekhova_T_cell_variance_models_SAVER.RData")
# load("Terekhova_T_cell_variance_models_SAVER.RData")

######################################################
# Compute gene-level and cell-level deviation per stratum (cell_type x Age_group)
cell_types <- sort(unique(meta$cell_type))
age_levels <- c('A','B','C','D','E')
ages       <- age_levels[age_levels %in% unique(meta$Age_group)]

all_temp_gene_level <- data.frame()
all_temp_cell_level <- data.frame()

for (ct in cell_types) {
  for (ag in ages) {
    idx <- which(meta$cell_type == ct & meta$Age_group == ag)
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
      Age = meta$Age[idx],
      Sex = meta$Sex[idx],
      Donor_id = meta$Donor_id[idx],
      Tube_id = meta$Tube_id[idx],
      Batch = meta$Batch[idx],
      Cell_level_deviation
    )
    all_temp_cell_level <- rbind(all_temp_cell_level, cell_df)

    rm(saver_sub, estimate, mu.out, cnt_norm, mu_sub)
    gc()
  }
}

write.csv(all_temp_gene_level, "Terekhova_T_cell_estimated_dispersion_SAVER.csv", row.names = FALSE)
write.csv(all_temp_cell_level, "Terekhova_T_cell_cellular_dispersion_SAVER.csv",  row.names = FALSE)

######################################################
