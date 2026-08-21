# Transcriptional dyscoordination analysis for TMS bone marrow
# This analysis can be directly applied to all other TMS tissues (lung, limb muscle, adipose tissue)

# Online links
# https://www.nature.com/articles/s41586-020-2496-1
# https://figshare.com/articles/dataset/Processed_files_to_use_with_scanpy_/8273102/2

# Set up working directory
# setwd("path/to/working/directory")

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(Matrix)
library(stringr)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

# Depth covariates default to loess (almost always non-linear at single-cell scale).
loess_covs <- c("nCount_RNA", "nFeature_RNA")

######################################################
# 1. Transcriptional dyscoordination analysis
# Load in and clean results
gene_level_dyscoordination <- read.csv("Dispersion_final/TMS_marrow_estimated_dispersion_SAVER.csv")
cell_level_dyscoordination <- read.csv("Dispersion_final/TMS_marrow_cellular_dispersion_SAVER.csv")

gene_level_dyscoordination$Gene_level_deviation <- remove_outliers(gene_level_dyscoordination$Gene_level_deviation)
cell_level_dyscoordination$Cell_level_deviation <- remove_outliers(cell_level_dyscoordination$Cell_level_deviation)
gene_level_dyscoordination_cleaned <- na.omit(gene_level_dyscoordination)
cell_level_dyscoordination_cleaned <- na.omit(cell_level_dyscoordination)

# Capitalize cell type names
cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
gene_level_dyscoordination_cleaned$cell_type <- cap_first(gene_level_dyscoordination_cleaned$cell_type)
cell_level_dyscoordination_cleaned$cell_type <- cap_first(cell_level_dyscoordination_cleaned$cell_type)

age_levels <- c("01m", "03m", "18m", "21m", "24m", "30m")
old_ages <- age_levels[-1]
cell_types <- sort(unique(gene_level_dyscoordination_cleaned$cell_type))

# Compute log-fold changes (LFCs) relative to the youngest age group per cell type
lfc_data <- data.frame()
for (ct in cell_types) {
  df_ct <- gene_level_dyscoordination_cleaned %>% filter(cell_type == ct)
  # Skip if less than 2 age groups present
  if (length(unique(df_ct$Age)) < 2) next
  # Add LFC values and cell type labels
  ct_baseline_age <- sort(unique(as.character(df_ct$Age)))[1]
  baseline_df <- df_ct %>% 
    filter(Age == ct_baseline_age) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    select(Gene, Baseline_deviation = Gene_level_deviation)
  df_ct_old <- df_ct %>% 
    filter(Age != ct_baseline_age) %>%
    left_join(baseline_df, by = "Gene") %>%
    mutate(LFC = log(Gene_level_deviation / Baseline_deviation),
           Baseline_Age = ct_baseline_age,
           Cell_type = paste0(ct, " (", ct_baseline_age, ")")) %>%
    filter(is.finite(LFC))
  lfc_data <- rbind(lfc_data, df_ct_old)
}

# Violin plots for transcriptional dyscoordination
signif_comparisons <- lapply(seq_along(old_ages[-length(old_ages)]),
                             function(i) c(old_ages[i], old_ages[i + 1]))

p_gene_level_dyscoordination <- lfc_data %>%
  group_by(Cell_type) %>%
  filter(LFC >= quantile(LFC, 0.025, na.rm = TRUE),
         LFC <= quantile(LFC, 0.975, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = Age, y = LFC, fill = Age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Blues") +
  facet_wrap(~ Cell_type, ncol = 3) +
  theme_minimal() +
  labs(x = "Age", y = "Log-fold change of gene-level dyscoordination") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "less"), step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

ggsave("TMS_marrow_gene_level_dyscoordination.png", p_gene_level_dyscoordination, width = 9.5, height = 11)

signif_comparisons <- lapply(seq_along(age_levels[-length(age_levels)]),
                             function(i) c(age_levels[i], age_levels[i + 1]))

p_cell_level_dyscoordination <- cell_level_dyscoordination_cleaned %>%
  group_by(cell_type) %>%
  filter(Cell_level_deviation >= quantile(Cell_level_deviation, 0.025, na.rm = TRUE),
         Cell_level_deviation <= quantile(Cell_level_deviation, 0.975, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = Age, y = log(Cell_level_deviation), fill = Age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Reds") +
  facet_wrap(~ cell_type, ncol = 3) +
  theme_minimal() +
  labs(x = "Age", y = "Cell-level dyscoordination (log-transformed)") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "less"), step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("TMS_marrow_cell_level_dyscoordination.png", p_cell_level_dyscoordination, width = 9.5, height = 11)

######################################################
# 2. Technical covariate analysis
cov_list <- c(nCount_RNA = "Sequencing depth",
              nFeature_RNA = "Genes detected",
              percent.ribo = "Ribosomal reads (%)",
              S.Score = "Cell cycle S score",
              G2M.Score = "Cell cycle G2/M score")

# Load in data
load('./TMS_marrow.RData')
  
ages <- data.frame(index = dataset.age)
age_map <- data.frame(index = c(1:length(dataset.age.levels)),
                      age = dataset.age.levels)
celltypes_df <- data.frame(index = dataset.celltype)
celltype_map <- data.frame(index = c(1:length(dataset.celltype.levels)),
                           celltype = cap_first(dataset.celltype.levels))
ages <- left_join(ages, age_map, by = "index")
celltypes_df <- left_join(celltypes_df, celltype_map, by = "index")
  
# Create and preprocess Seurat object
TMS_marrow <- CreateSeuratObject(counts = dataset.counts)
TMS_marrow$age <- factor(ages$age, levels = dataset.age.levels)
TMS_marrow$celltype <- factor(celltypes_df$celltype, levels = cap_first(dataset.celltype.levels))
Idents(TMS_marrow) <- TMS_marrow$celltype
TMS_marrow <- NormalizeData(TMS_marrow)
TMS_marrow <- FindVariableFeatures(TMS_marrow, selection.method = "vst", nfeatures = 2000)
TMS_marrow <- ScaleData(TMS_marrow)
TMS_marrow <- RunPCA(TMS_marrow, features = VariableFeatures(object = TMS_marrow))
TMS_marrow <- FindNeighbors(TMS_marrow, dims = 1:30)
TMS_marrow <- RunUMAP(TMS_marrow, dims = 1:30)
  
# Compute QC metrics
TMS_marrow[["percent.ribo"]] <- PercentageFeatureSet(TMS_marrow, pattern = "^Rp[sl]")
s_genes <- intersect(str_to_title(cc.genes.updated.2019$s.genes), rownames(TMS_marrow))
g2m_genes <- intersect(str_to_title(cc.genes.updated.2019$g2m.genes), rownames(TMS_marrow))
TMS_marrow <- CellCycleScoring(TMS_marrow, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

# UMAP of cell types included in downstream analyses
p_umap <- DimPlot(TMS_marrow, reduction = "umap", group.by = "celltype",
                  label = TRUE, repel = TRUE, label.size = 4) +
  theme_minimal() +
  labs(title = "TMS marrow cell types (UMAP)")
ggsave("TMS_marrow_UMAP_celltype.png", p_umap, width = 8, height = 8)

# Cell-count table — Celltype rows × Age columns, values = n_cells, with margins
cell_counts <- TMS_marrow@meta.data %>%
  as.data.frame() %>%
  filter(age %in% age_levels) %>%
  count(celltype, age) %>%
  rename(Celltype = celltype, Age = age) %>%
  mutate(Age = factor(Age, levels = age_levels)) %>%
  tidyr::pivot_wider(names_from = Age, values_from = n, values_fill = 0) %>%
  arrange(Celltype)
cell_counts <- cell_counts %>% mutate(Total = rowSums(dplyr::select(., -Celltype)))
cell_counts <- dplyr::bind_rows(
  cell_counts,
  cell_counts %>% summarise(Celltype = "Total", dplyr::across(-Celltype, sum)))
write.csv(cell_counts, "TMS_marrow_cell_counts.csv", row.names = FALSE)
cat("Cell counts (Celltype x Age) with margins:\n")
print(cell_counts)

TMS_marrow_metadata <- TMS_marrow@meta.data %>%
    as.data.frame() %>%
    select(nCount_RNA, nFeature_RNA, percent.ribo, S.Score, G2M.Score) %>%
    mutate(Cell_barcode = rownames(TMS_marrow@meta.data))

rm(TMS_marrow, dataset.counts)
  
cell_level_cov <- cell_level_dyscoordination_cleaned %>%
    filter(Cell_level_deviation > 0) %>%
    inner_join(TMS_marrow_metadata, by = "Cell_barcode") %>%
    filter(!is.na(Cell_level_deviation)) %>%
    mutate(log_deviation = log(Cell_level_deviation),
           Age = factor(Age, levels = age_levels))
  
# Create scatter plots for covariates vs cell-level transcriptional dyscoordination
plot_cov_scatter <- function(df, x_var, x_label) {
  x_lo <- quantile(df[[x_var]], 0.01, na.rm = TRUE)
  x_hi <- quantile(df[[x_var]], 0.99, na.rm = TRUE)
  ggplot(df %>% filter(.data[[x_var]] >= x_lo & .data[[x_var]] <= x_hi),
         aes_string(x = x_var, y = "log_deviation")) +
    geom_point(aes(fill = Age), shape = 21, color = "black",
               stroke = 0.2, alpha = 0.5, size = 1.2) +
    geom_smooth(method = "lm", linetype = "dashed", color = "black", linewidth = 0.5) +
    stat_cor(method = "spearman", size = 2.5, cor.coef.name = "rho") +
    facet_wrap(~ cell_type, ncol = 3) +
    scale_fill_brewer(palette = "Blues") +
    theme_minimal() +
    labs(x = x_label, y = "Cell-level dyscoordination (log-transformed)")
  }
  
for (v in names(cov_list)) {
  ggsave(paste0("TMS_marrow_covariate_scatter_", gsub("\\.", "_", v), ".png"),
         plot_cov_scatter(cell_level_cov, v, cov_list[[v]]),
         width = 8, height = 10)
}

# Compute corrected cell-level transcriptional dyscoordination with the following covariates
fit <- loess(log_deviation ~ nCount_RNA, data = cell_level_cov, span = 1)
cell_level_cov$log_deviation_corrected <- cell_level_cov$log_deviation - predict(fit, newdata = cell_level_cov)
na_idx <- is.na(cell_level_cov$log_deviation_corrected)
cell_level_cov$log_deviation_corrected[na_idx] <- cell_level_cov$log_deviation[na_idx] - mean(cell_level_cov$log_deviation, na.rm = TRUE) 

# Violin plot for corrected cell-level transcriptional dyscoordination
p_cell_level_corrected <- cell_level_cov %>%
  filter(!is.na(log_deviation_corrected),
         log_deviation_corrected >= quantile(log_deviation_corrected, 0.025, na.rm = TRUE),
         log_deviation_corrected <= quantile(log_deviation_corrected, 0.975, na.rm = TRUE)) %>%
  ggplot(aes(x = Age, y = log_deviation_corrected, fill = Age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Reds") +
  facet_wrap(~ cell_type, ncol = 3) +
  theme_minimal() +
  labs(x = "Age", y = "Cell-level dyscoordination (log-transformed, corrected)") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "less"),
              step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("TMS_marrow_cell_level_dyscoordination_corrected.png", p_cell_level_corrected, width = 9.5, height = 11)