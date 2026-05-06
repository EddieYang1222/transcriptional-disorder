# Transcriptional dyscoordination analysis for senescent mouse aortic cells

# Online links
# https://www.nature.com/articles/s43587-025-00889-z
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239591

# Set up working directory
# setwd("S:/Penn Dropbox/Eddie Yang/Aging/Scripts/Mouse_VSMC")

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
gene_level_dyscoordination <- read.csv("senescent_mouse_aortic_cell_estimated_dispersion_SAVER.csv")
cell_level_dyscoordination <- read.csv("senescent_mouse_aortic_cell_cellular_dispersion_SAVER.csv")

gene_level_dyscoordination$Gene_level_deviation <- remove_outliers(gene_level_dyscoordination$Gene_level_deviation)
cell_level_dyscoordination$Cell_level_deviation <- remove_outliers(cell_level_dyscoordination$Cell_level_deviation)
gene_level_dyscoordination_cleaned <- na.omit(gene_level_dyscoordination)
cell_level_dyscoordination_cleaned <- na.omit(cell_level_dyscoordination)

# Capitalize cell type names
cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

cond_levels <- c("ND", "HFD", "ABT")
baseline <- "ND"
non_base <- setdiff(cond_levels, baseline)
cond_labels <- c("ND" = "ND", "HFD" = "HFD", "ABT" = "HFD + ABT-737")
display_levels <- unname(cond_labels[cond_levels])

gene_level_dyscoordination_cleaned <- gene_level_dyscoordination_cleaned %>%
  mutate(Celltype = cap_first(Celltype),
         Condition = recode(Condition, "ABT" = "HFD + ABT-737"),
         Condition = factor(Condition, levels = display_levels))
cell_level_dyscoordination_cleaned <- cell_level_dyscoordination_cleaned %>%
  mutate(Celltype = cap_first(Celltype),
         Condition = recode(Condition, "ABT" = "HFD + ABT-737"),
         Condition = factor(Condition, levels = display_levels))

cell_types <- sort(unique(gene_level_dyscoordination_cleaned$Celltype))

non_baseline_labels <- cond_labels[non_base]

# Compute log-fold change of gene-level dyscoordination relative to ND baseline per cell type
lfc_data <- data.frame()
for (ct in cell_types) {
  df_ct <- gene_level_dyscoordination_cleaned %>% filter(Celltype == ct)
  # Skip if baseline condition not present
  if (!baseline %in% as.character(df_ct$Condition)) next
  # Add LFC values and cell type labels
  baseline_df <- df_ct %>%
    filter(Condition == baseline) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    select(Gene, Baseline_deviation = Gene_level_deviation)
  df_ct_old <- df_ct %>%
    filter(Condition != baseline) %>%
    left_join(baseline_df, by = "Gene") %>%
    mutate(LFC = log(Gene_level_deviation / Baseline_deviation),
           Cell_type = paste0(ct, " (", baseline, ")"),
           Condition = factor(as.character(Condition), levels = display_levels)) %>%
    filter(is.finite(LFC))
  lfc_data <- rbind(lfc_data, df_ct_old)
}

# Violin plots for transcriptional dyscoordination
signif_comparisons <- list(c("HFD", "HFD + ABT-737"))

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
              test = function(x, y) wilcox.test(x, y, alternative = "greater"), step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("senescent_mouse_aortic_cell_gene_level_dyscoordination.png", p_gene_level_dyscoordination, width = 8, height = 6)

p_cell_level_dyscoordination <- cell_level_dyscoordination_cleaned %>%
  filter(Cell_level_deviation >= quantile(Cell_level_deviation, 0.01, na.rm = TRUE),
         Cell_level_deviation <= quantile(Cell_level_deviation, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Condition, y = log(Cell_level_deviation), fill = Condition)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Reds") +
  facet_wrap(~ Celltype, ncol = 3) +
  theme_minimal() +
  labs(x = "Condition", y = "Cell-level dyscoordination (log-transformed)") +
  # ND -> HFD: disease induction, dyscoordination increases
  geom_signif(comparisons = list(c("ND", "HFD")),
              test = function(x, y) wilcox.test(x, y, alternative = "less"),
              map_signif_level = TRUE, step_increase = 0.1) +
  # HFD -> HFD+ABT: senolytic reversal, dyscoordination decreases
  geom_signif(comparisons = list(c("HFD", "HFD + ABT-737")),
              test = function(x, y) wilcox.test(x, y, alternative = "greater"),
              map_signif_level = TRUE, step_increase = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("senescent_mouse_aortic_cell_cell_level_dyscoordination.png", p_cell_level_dyscoordination, width = 8, height = 6)

######################################################
# 2. Technical covariate analysis
cov_list <- c(nCount_RNA = "Sequencing depth",
              nFeature_RNA = "Genes detected",
              percent.ribo = "Ribosomal reads (%)",
              S.Score = "Cell cycle S score",
              G2M.Score = "Cell cycle G2/M score")

# Load in data
data_dir <- "S:/Penn Dropbox/Eddie Yang/Aging/Data/Mouse_VSMC"
seu <- readRDS(file.path(data_dir, "mouse_VSMC_aorta_processed.rds"))

seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^Rp[sl]")
s_genes_mm <- str_to_title(cc.genes.updated.2019$s.genes)
g2m_genes_mm <- str_to_title(cc.genes.updated.2019$g2m.genes)
s_found <- intersect(s_genes_mm, rownames(seu))
g2m_found <- intersect(g2m_genes_mm, rownames(seu))
seu <- CellCycleScoring(seu, s.features = s_found, g2m.features = g2m_found, set.ident = FALSE)

aortic_cell_metadata <- seu@meta.data %>%
  as.data.frame() %>%
  mutate(Cell_barcode = rownames(seu@meta.data),
         Celltype = cap_first(as.character(seu$celltype))) %>%
  select(Cell_barcode, Celltype, nCount_RNA, nFeature_RNA, percent.ribo, S.Score, G2M.Score)

rm(seu)

cell_level_cov <- cell_level_dyscoordination_cleaned %>%
  filter(Cell_level_deviation > 0) %>%
  inner_join(aortic_cell_metadata, by = c("Cell_barcode", "Celltype")) %>%
  filter(!is.na(Cell_level_deviation)) %>%
  mutate(log_deviation = log(Cell_level_deviation),
         Condition = factor(Condition, levels = unname(cond_labels[cond_levels])))

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
    facet_wrap(~ Celltype, ncol = 3) +
    scale_fill_brewer(palette = "Blues") +
    theme_minimal() +
    labs(x = x_label, y = "Cell-level dyscoordination (log-transformed)") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

for (v in names(cov_list)) {
  ggsave(paste0("senescent_mouse_aortic_cell_covariate_scatter_", gsub("\\.", "_", v), ".png"),
         plot_cov_scatter(cell_level_cov, v, cov_list[[v]]),
         width = 8, height = 6)
}

# Compute corrected cell-level transcriptional dyscoordination with the following covariates
correct_covariates <- c("nCount_RNA", "nFeature_RNA")
fmla <- as.formula(paste("log_deviation ~", paste(correct_covariates, collapse = " + ")))

# Regress within each cell type (condition effect is preserved across groups)
cell_level_cov$log_deviation_corrected <- unsplit(
  lapply(split(cell_level_cov, cell_level_cov$Celltype), function(g) {
    if (nrow(g) < max(5, length(correct_covariates) + 2))
      return(rep(NA_real_, nrow(g)))
    residuals(lm(fmla, data = g))
  }),
  cell_level_cov$Celltype
)

# Violin plot for corrected cell-level transcriptional dyscoordination
p_cell_level_corrected <- cell_level_cov %>%
  filter(!is.na(log_deviation_corrected),
         log_deviation_corrected >= quantile(log_deviation_corrected, 0.01, na.rm = TRUE),
         log_deviation_corrected <= quantile(log_deviation_corrected, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Condition, y = log_deviation_corrected, fill = Condition)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Reds") +
  facet_wrap(~ Celltype, ncol = 3) +
  theme_minimal() +
  labs(x = "Condition", y = "Cell-level dyscoordination (log-transformed, corrected)") +
  # ND -> HFD: disease induction, dyscoordination increases
  geom_signif(comparisons = list(c("ND", "HFD")),
              test = function(x, y) wilcox.test(x, y, alternative = "less"),
              map_signif_level = TRUE, step_increase = 0.1) +
  # HFD -> HFD+ABT: senolytic reversal, dyscoordination decreases
  geom_signif(comparisons = list(c("HFD", "HFD + ABT-737")),
              test = function(x, y) wilcox.test(x, y, alternative = "greater"),
              map_signif_level = TRUE, step_increase = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("senescent_mouse_aortic_cell_cell_level_dyscoordination_corrected.png",
       p_cell_level_corrected, width = 8, height = 6)

######################################################