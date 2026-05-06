# Transcriptional dyscoordination analysis for cancerous human T cells

# Online links
# https://www.biorxiv.org/content/10.1101/2025.09.01.673503v1

# Set up working directory
# setwd("S:/Penn Dropbox/Eddie Yang/Aging/Scripts/T_cell_exhaustion")

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

######################################################
# 1. Transcriptional dyscoordination analysis
# Load in and clean results
gene_level_dyscoordination <- read.csv("cancerous_human_T_cell_estimated_dispersion_SAVER.csv")
cell_level_dyscoordination <- read.csv("cancerous_human_T_cell_cellular_dispersion_SAVER.csv")

gene_level_dyscoordination$Gene_level_deviation <- remove_outliers(gene_level_dyscoordination$Gene_level_deviation)
cell_level_dyscoordination$Cell_level_deviation <- remove_outliers(cell_level_dyscoordination$Cell_level_deviation)
gene_level_dyscoordination_cleaned <- na.omit(gene_level_dyscoordination)
cell_level_dyscoordination_cleaned <- na.omit(cell_level_dyscoordination)

exhaustion_levels <- c("Tnaive", "IL7R+ Tm", "GZMK+ Tem", "GZMK+ Trm",
                        "NME1+ T", "NR4A2+ Teff", "IL7R+ Temra",
                        "Activated Trm", "Temra", "NK-like", "Term Teff")

grade_levels <- c("Normal", "II", "III", "IV")
treatment_levels <- c("Untreated", "Treated")
condition_levels <- c(grade_levels, paste0(tail(grade_levels, 1), " (treated)"))
gene_level_dyscoordination_cleaned <- gene_level_dyscoordination_cleaned %>%
  filter(cell_type %in% exhaustion_levels) %>%
  mutate(cell_type = factor(cell_type, levels = exhaustion_levels),
         Condition = factor(ifelse(Treatment == "Treated", paste0(Tumor.Grade, " (treated)"),
                                   Tumor.Grade), levels = condition_levels))
cell_level_dyscoordination_cleaned <- cell_level_dyscoordination_cleaned %>%
  filter(cell_type %in% exhaustion_levels) %>%
  mutate(cell_type = factor(cell_type, levels = exhaustion_levels),
         Condition = factor(ifelse(Treatment == "Treated", paste0(Tumor.Grade, " (treated)"),
                                   Tumor.Grade), levels = condition_levels))

# Compute log-fold change of gene-level dyscoordination relative to Normal baseline per cell type
cell_types <- sort(unique(gene_level_dyscoordination_cleaned$cell_type))
lfc_data <- data.frame()
for (ct in cell_types) {
  df_ct <- gene_level_dyscoordination_cleaned %>% filter(cell_type == ct)
  if (length(unique(df_ct$Condition)) < 2) next
  baseline_df <- df_ct %>%
    filter(Condition == "Normal") %>%
    distinct(Gene, .keep_all = TRUE) %>%
    select(Gene, Baseline_deviation = Gene_level_deviation)
  df_ct_non_baseline <- df_ct %>%
    filter(Condition != "Normal") %>%
    left_join(baseline_df, by = "Gene") %>%
    mutate(LFC = log(Gene_level_deviation / Baseline_deviation),
           Cell_type = paste0(ct, " (Normal)")) %>%
    filter(is.finite(LFC))
  lfc_data <- rbind(lfc_data, df_ct_non_baseline)
}

# Violin plots for transcriptional dyscoordination
non_baseline_conditions <- condition_levels[-1]
signif_comparisons <- lapply(seq_along(non_baseline_conditions[-length(non_baseline_conditions)]),
                             function(i) c(non_baseline_conditions[i], non_baseline_conditions[i + 1]))

p_gene_level_dyscoordination <- lfc_data %>%
  filter(LFC >= quantile(LFC, 0.01, na.rm = TRUE),
         LFC <= quantile(LFC, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Condition, y = LFC, fill = Condition)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Blues") +
  facet_wrap(~ Cell_type, ncol = 3) +
  theme_minimal() +
  labs(x = "Condition", y = "Log-fold change of gene-level dyscoordination") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "less"), step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

ggsave("cancerous_human_T_cell_gene_level_dyscoordination.png",
       p_gene_level_dyscoordination, width = 8, height = 8)

signif_comparisons <- lapply(seq_along(condition_levels[-length(condition_levels)]),
                             function(i) c(condition_levels[i], condition_levels[i + 1]))

p_cell_level_dyscoordination <- cell_level_dyscoordination_cleaned %>%
  filter(Cell_level_deviation >= quantile(Cell_level_deviation, 0.01, na.rm = TRUE),
         Cell_level_deviation <= quantile(Cell_level_deviation, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Condition, y = log(Cell_level_deviation), fill = Condition)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Reds") +
  facet_wrap(~ cell_type, ncol = 3) +
  theme_minimal() +
  labs(x = "Condition", y = "Cell-level dyscoordination (log-transformed)") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "less"), step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("cancerous_human_T_cell_cell_level_dyscoordination.png",
       p_cell_level_dyscoordination, width = 8, height = 8)

######################################################
# 2. Technical covariate analysis
cov_list <- c(nCount_RNA = "Sequencing depth",
              nFeature_RNA = "Genes detected",
              percent.mt = "Mitochondrial reads (%)",
              percent.ribo = "Ribosomal reads (%)",
              S.Score = "Cell cycle S score",
              G2M.Score = "Cell cycle G2/M score")

# Build Seurat object
# data_dir <- 'S:/Penn Dropbox/Eddie Yang/Aging/Data/T_cell'
sparse_counts <- Matrix::readMM(file.path(data_dir, "CD830dim_sparse/matrix.mtx"))
barcodes_all <- readLines(file.path(data_dir, "CD830dim_sparse/barcodes.tsv"))
features_all <- read.delim(file.path(data_dir, "CD830dim_sparse/features.tsv"), header = FALSE)$V1
rownames(sparse_counts) <- features_all
colnames(sparse_counts) <- barcodes_all

common_cells <- intersect(colnames(sparse_counts), cell_level_dyscoordination_cleaned$Cell_barcode)
sparse_counts <- sparse_counts[, common_cells]

seu <- CreateSeuratObject(counts = sparse_counts, min.cells = 0, min.features = 0)
rm(sparse_counts)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
seu <- NormalizeData(seu)
s_genes <- intersect(cc.genes.updated.2019$s.genes, rownames(seu))
g2m_genes <- intersect(cc.genes.updated.2019$g2m.genes, rownames(seu))
seu <- CellCycleScoring(seu, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

t_cell_metadata <- seu@meta.data %>%
  as.data.frame() %>%
  select(nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, S.Score, G2M.Score) %>%
  mutate(Cell_barcode = rownames(seu@meta.data))

rm(seu)

cell_level_cov <- cell_level_dyscoordination_cleaned %>%
  inner_join(t_cell_metadata, by = "Cell_barcode") %>%
  filter(Cell_level_deviation > quantile(Cell_level_deviation, 0.01, na.rm = TRUE),
         Cell_level_deviation < quantile(Cell_level_deviation, 0.99, na.rm = TRUE)) %>%
  mutate(log_deviation = log(Cell_level_deviation))

# Create scatter plots for covariates vs cell-level transcriptional dyscoordination
plot_cov_scatter <- function(df, x_var, x_label) {
  x_lo <- quantile(df[[x_var]], 0.01, na.rm = TRUE)
  x_hi <- quantile(df[[x_var]], 0.99, na.rm = TRUE)
  ggplot(df %>% filter(.data[[x_var]] >= x_lo & .data[[x_var]] <= x_hi),
         aes_string(x = x_var, y = "log_deviation")) +
    geom_point(aes(fill = Condition), shape = 21, color = "black",
               stroke = 0.2, alpha = 0.5, size = 1.2) +
    geom_smooth(method = "lm", linetype = "dashed", color = "black", linewidth = 0.5) +
    stat_cor(method = "spearman", size = 2.5, cor.coef.name = "rho") +
    facet_wrap(~ cell_type, ncol = 3) +
    scale_fill_brewer(palette = "Blues") +
    theme_minimal() +
    labs(x = x_label, y = "Cell-level dyscoordination (log-transformed)")
}

for (v in names(cov_list)) {
  ggsave(paste0("cancerous_human_T_cell_covariate_scatter_", gsub("_RNA", "", v), ".png"),
         plot_cov_scatter(cell_level_cov, v, cov_list[[v]]),
         width = 8, height = 8)
}

# Compute corrected cell-level transcriptional dyscoordination with the following covariates
correct_covariates <- c("nCount_RNA", "nFeature_RNA")
fmla <- as.formula(paste("log_deviation ~", paste(correct_covariates, collapse = "+")))

cell_level_cov$log_deviation_corrected <- unsplit(
  lapply(split(cell_level_cov, cell_level_cov$cell_type), function(g) {
    if (nrow(g) < max(5, length(correct_covariates) + 2))
      return(rep(NA_real_, nrow(g)))
    residuals(lm(fmla, data = g))
  }),
  cell_level_cov$cell_type
)

p_cell_level_corrected <- cell_level_cov %>%
  filter(!is.na(log_deviation_corrected),
         log_deviation_corrected >= quantile(log_deviation_corrected, 0.01, na.rm = TRUE),
         log_deviation_corrected <= quantile(log_deviation_corrected, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Condition, y = log_deviation_corrected, fill = Condition)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Reds") +
  facet_wrap(~ cell_type, ncol = 3) +
  theme_minimal() +
  labs(x = "Condition", y = "Cell-level dyscoordination (log-transformed, corrected)") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "less"),
              step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("cancerous_human_T_cell_cell_level_dyscoordination_corrected.png",
       p_cell_level_corrected, width = 8, height = 8)

######################################################