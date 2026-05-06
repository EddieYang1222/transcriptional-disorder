# Transcriptional dyscoordination analysis for aging mouse liver (hepatocyte only)

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
gene_level_dyscoordination <- read.csv("liver_data_hepatocyte_only_estimated_dispersion_SAVER.csv")
cell_level_dyscoordination <- read.csv("liver_data_hepatocyte_only_cellular_dispersion_SAVER.csv")

gene_level_dyscoordination$Gene_level_deviation <- remove_outliers(gene_level_dyscoordination$Gene_level_deviation)
cell_level_dyscoordination$Cell_level_deviation <- remove_outliers(cell_level_dyscoordination$Cell_level_deviation)
gene_level_dyscoordination_cleaned <- na.omit(gene_level_dyscoordination)
cell_level_dyscoordination_cleaned <- na.omit(cell_level_dyscoordination)

age_levels <- c("young", "old")
gene_level_dyscoordination_cleaned$Age <- factor(gene_level_dyscoordination_cleaned$Age, levels = age_levels)
cell_level_dyscoordination_cleaned$Age <- factor(cell_level_dyscoordination_cleaned$Age, levels = age_levels)

signif_comparisons <- list(c("young", "old"))

p_gene_level_dyscoordination <- gene_level_dyscoordination_cleaned %>%
  filter(Gene_level_deviation >= quantile(Gene_level_deviation, 0.01, na.rm = TRUE),
         Gene_level_deviation <= quantile(Gene_level_deviation, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Age, y = log(Gene_level_deviation), fill = Age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Blues") +
  theme_minimal() +
  labs(x = "Age", y = "Gene-level dyscoordination (log-transformed)") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "less"), step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("aging_mouse_liver_gene_level_dyscoordination.png", p_gene_level_dyscoordination, width = 5, height = 4)

p_cell_level_dyscoordination <- cell_level_dyscoordination_cleaned %>%
  filter(Cell_level_deviation >= quantile(Cell_level_deviation, 0.01, na.rm = TRUE),
         Cell_level_deviation <= quantile(Cell_level_deviation, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Age, y = log(Cell_level_deviation), fill = Age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Reds") +
  theme_minimal() +
  labs(x = "Age", y = "Cell-level dyscoordination (log-transformed)") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "less"), step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("aging_mouse_liver_cell_level_dyscoordination.png", p_cell_level_dyscoordination, width = 5, height = 4)

######################################################
# 2. Technical covariate analysis
cov_list <- c(nCount_RNA = "Sequencing depth",
              nFeature_RNA = "Genes detected",
              percent.mt = "Mitochondrial reads (%)",
              percent.ribo = "Ribosomal reads (%)",
              S.Score = "Cell cycle S score",
              G2M.Score = "Cell cycle G2/M score")

# Load in data
load('./liver_data_integrated.RData')

liver_hepatocyte <- subset(liver_data_integrated,
                            subset = doublet == "Singlet" & annotation == "Hepatocyte")
rm(liver_data_integrated)

liver_hepatocyte[["percent.mt"]] <- PercentageFeatureSet(liver_hepatocyte, pattern = "^mt-")
liver_hepatocyte[["percent.ribo"]] <- PercentageFeatureSet(liver_hepatocyte, pattern = "^Rp[sl]")
s_genes <- intersect(str_to_title(cc.genes.updated.2019$s.genes), rownames(liver_hepatocyte))
g2m_genes <- intersect(str_to_title(cc.genes.updated.2019$g2m.genes), rownames(liver_hepatocyte))
liver_hepatocyte <- CellCycleScoring(liver_hepatocyte, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

liver_hepatocyte_metadata <- liver_hepatocyte@meta.data %>%
  as.data.frame() %>%
  select(nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, S.Score, G2M.Score) %>%
  mutate(Cell_barcode = rownames(liver_hepatocyte@meta.data))

rm(liver_hepatocyte)

cell_level_cov <- cell_level_dyscoordination_cleaned %>%
  filter(Cell_level_deviation > 0) %>%
  inner_join(liver_hepatocyte_metadata, by = "Cell_barcode") %>%
  filter(!is.na(Cell_level_deviation)) %>%
  mutate(log_deviation = log(Cell_level_deviation),
         Age = factor(Age, levels = age_levels))

# Create one combined covariate scatter plot
cov_long <- do.call(rbind, lapply(names(cov_list), function(v) {
  x_lo <- quantile(cell_level_cov[[v]], 0.01, na.rm = TRUE)
  x_hi <- quantile(cell_level_cov[[v]], 0.99, na.rm = TRUE)
  cell_level_cov %>%
    filter(.data[[v]] >= x_lo & .data[[v]] <= x_hi) %>%
    select(Cell_barcode, Age, log_deviation, all_of(v)) %>%
    rename(x_value = all_of(v)) %>%
    mutate(Covariate = factor(cov_list[[v]], levels = unname(cov_list)))
}))

p_cov_scatter <- ggplot(cov_long, aes(x = x_value, y = log_deviation)) +
  geom_point(aes(fill = Age), shape = 21, color = "black",
             stroke = 0.2, alpha = 0.5, size = 1.2) +
  geom_smooth(method = "lm", linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", size = 2.5, cor.coef.name = "rho") +
  scale_fill_brewer(palette = "Blues") +
  facet_wrap(~ Covariate, ncol = 3, scales = "free_x") +
  theme_minimal() +
  labs(x = NULL, y = "Cell-level dyscoordination (log-transformed)")

ggsave("aging_mouse_liver_covariate_scatter.png", p_cov_scatter, width = 10, height = 7)

# Compute corrected cell-level transcriptional dyscoordination with the following covariates
correct_covariates <- c("nCount_RNA", "nFeature_RNA")
fmla <- as.formula(paste("log_deviation ~", paste(correct_covariates, collapse = " + ")))

# Regress within each cell type (age effect is preserved across groups)
cell_level_cov$log_deviation_corrected <- residuals(lm(fmla, data = cell_level_cov))

# Violin plot for corrected cell-level transcriptional dyscoordination
p_cell_level_corrected <- cell_level_cov %>%
  filter(!is.na(log_deviation_corrected),
         log_deviation_corrected >= quantile(log_deviation_corrected, 0.01, na.rm = TRUE),
         log_deviation_corrected <= quantile(log_deviation_corrected, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Age, y = log_deviation_corrected, fill = Age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Reds") +
  theme_minimal() +
  labs(x = "Age", y = "Cell-level dyscoordination (log-transformed, corrected)") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "less"),
              step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("aging_mouse_liver_cell_level_dyscoordination_corrected.png", p_cell_level_corrected, width = 5, height = 4.5)

######################################################