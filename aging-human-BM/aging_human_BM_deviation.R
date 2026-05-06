# Transcriptional dyscoordination analysis for aging human bone marrow

# Online links
# https://insight.jci.org/articles/view/124928
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120446

library(Seurat)
library(dplyr)
library(Matrix)
library(pbapply)
library(SAVER)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Load and prepare data
sample_info <- data.frame(
  sample = c("A","B","C1","Ck","C2","E","F","G","H","J","K","L","M","N","O","P","Q","R","S1","Sk1","S2","Sk2","T","U","W"),
  gsm = c("GSM3396161","GSM3396162","GSM3396163","GSM3396164","GSM3396165",
          "GSM3396166","GSM3396167","GSM3396168","GSM3396169","GSM3396170",
          "GSM3396171","GSM3396172","GSM3396173","GSM3396174","GSM3396175",
          "GSM3396176","GSM3396177","GSM3396178","GSM3396179","GSM3396180",
          "GSM3396181","GSM3396182","GSM3396183","GSM3396184","GSM3396185"),
  donor = c("A","B","C","C","C","E","F","G","H","J","K","L","M","N","O","P","Q","R","S","Sk","S","Sk","T","U","W"),
  age = c(59,47,60,59,60,30,41,58,50,43,84,57,60,67,50,58,66,31,56,55,56,55,24,46,28),
  sex = c("F","M","F","F","F","M","F","M","F","F","M","M","M","M","M","F","M","M","F","F","F","F","F","F","F"),
  stringsAsFactors = FALSE
)

data_dir <- 'path/to/GSE120446_raw_data'
seurat_list <- list()
for (i in 1:nrow(sample_info)) {
  s <- sample_info$sample[i]; gsm <- sample_info$gsm[i]
  mat <- ReadMtx(mtx = file.path(data_dir, paste0(gsm, '_matrix_', s, '.mtx.gz')),
                 cells = file.path(data_dir, paste0(gsm, '_barcodes_', s, '.tsv.gz')),
                 features = file.path(data_dir, paste0(gsm, '_genes_', s, '.tsv.gz')))
  obj <- CreateSeuratObject(counts = mat, project = s, min.cells = 3, min.features = 200)
  obj@meta.data$sample <- s
  obj@meta.data$donor <- sample_info$donor[i]
  obj@meta.data$age <- sample_info$age[i]
  obj@meta.data$sex <- sample_info$sex[i]
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = 2000)
  seurat_list[[s]] <- obj
}
anchors <- FindIntegrationAnchors(seurat_list, dims = 1:30, reduction = 'cca')
bm <- IntegrateData(anchors, dims = 1:30)

# Cluster and annotate cell types
# Cell type labels assigned manually based on canonical marker gene expression
DefaultAssay(bm) <- 'integrated'
bm <- ScaleData(bm) %>% 
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  RunTSNE(dims = 1:20)

saveRDS(bm, 'human_BM_Oetjen_annotated.rds')

# Extract counts and metadata from annotated object
bm <- readRDS('human_BM_Oetjen_annotated.rds')
DefaultAssay(bm) <- 'RNA'

dataset.counts.sparse <- GetAssayData(bm, assay = 'RNA', layer = 'counts')
dataset.celltype <- bm@meta.data$celltype
dataset.age <- bm@meta.data$age

# Remove unannotated cells
keep.cells <- !is.na(dataset.celltype)
dataset.counts.sparse <- dataset.counts.sparse[, keep.cells]
dataset.celltype <- dataset.celltype[keep.cells]
dataset.age <- dataset.age[keep.cells]

# Apply QC filters
keep.cells.qc <- Matrix::colSums(dataset.counts.sparse > 0) >= 250 &
                 Matrix::colSums(dataset.counts.sparse) >= 2500
keep.genes.qc <- Matrix::rowSums(dataset.counts.sparse > 0) >= 5
dataset.counts.sparse <- dataset.counts.sparse[keep.genes.qc, keep.cells.qc]
dataset.celltype <- dataset.celltype[keep.cells.qc]
dataset.age <- dataset.age[keep.cells.qc]

# Create age bins for older samples
dataset.age_group <- dplyr::case_when(
  dataset.age >= 46 & dataset.age <= 50 ~ "46-50",
  dataset.age >= 55 & dataset.age <= 57 ~ "55-57",
  dataset.age >= 58 & dataset.age <= 60 ~ "58-60",
  dataset.age >= 66 & dataset.age <= 67 ~ "66-67",
  dataset.age == 84 ~ "84",
  TRUE ~ NA_character_
)

dataset.counts <- as.matrix(dataset.counts.sparse)
size.factor <- colSums(dataset.counts) / mean(colSums(dataset.counts))
save(dataset.counts, dataset.celltype, dataset.age, dataset.age_group,
     size.factor, file = 'human_BM_data.RData')
load('human_BM_data.RData')

# Run SAVER for manifold fitting
# This step is typically computationally heavy and takes at least several hours for >= 10,000 cells
# We recommend running with at least 128 GB of RAM and saving the manifold separately
dataset.saver <- saver(dataset.counts, ncores = 4)

# Load and save pre-computed manifold
# save(dataset.saver, file = "human_BM_manifold.RData")
# load("human_BM_manifold.RData")

rownames(dataset.saver$estimate) <- rownames(dataset.counts)
colnames(dataset.saver$estimate) <- colnames(dataset.counts)
rownames(dataset.saver$mu.out) <- rownames(dataset.counts)
colnames(dataset.saver$mu.out) <- colnames(dataset.counts)

######################################################
# Find dispersion model with highest likelihood for each gene
dataset.saver.mu <- dataset.saver$mu.out
dataset.saver.var.models <- get_var_model(dataset.counts, dataset.saver.mu, size.factor)
# save(dataset.saver.var.models, file = "human_BM_variance_models_SAVER.RData")
# load("human_BM_variance_models_SAVER.RData")

######################################################
# Compute gene-level and cell-level deviation
dataset.celltype.levels <- levels(factor(dataset.celltype))
dataset.celltype <- as.integer(factor(dataset.celltype, levels = dataset.celltype.levels))
dataset.age.levels <- c('46-50', '55-57', '58-60', '66-67', '84')
dataset.age <- as.integer(factor(dataset.age_group, levels = dataset.age.levels))

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

write.csv(all_temp_gene_level, "aging_human_BM_estimated_dispersion_SAVER.csv", row.names = FALSE)
write.csv(all_temp_cell_level, "aging_human_BM_cellular_dispersion_SAVER.csv", row.names = FALSE)

######################################################