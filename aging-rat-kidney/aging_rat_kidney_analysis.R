# Transcriptional dyscoordination analysis for aging rat kidney (PT cells only)

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
gene_level_dyscoordination <- read.csv("kidney_aging_all_samples_PT_only_estimated_dispersion_SAVER.csv")
cell_level_dyscoordination <- read.csv("kidney_aging_all_samples_PT_only_cellular_dispersion_SAVER.csv")

gene_level_dyscoordination$Gene_level_deviation <- remove_outliers(gene_level_dyscoordination$Gene_level_deviation)
cell_level_dyscoordination$Cell_level_deviation <- remove_outliers(cell_level_dyscoordination$Cell_level_deviation)
gene_level_dyscoordination_cleaned <- na.omit(gene_level_dyscoordination)
cell_level_dyscoordination_cleaned <- na.omit(cell_level_dyscoordination)

age_levels <- c("16wk", "30wk", "56wk", "82wk")
old_ages <- age_levels[-1]
cell_types <- c("PT-S1", "PT-S2", "PT-S3")

gene_level_dyscoordination_cleaned <- gene_level_dyscoordination_cleaned %>%
  filter(cell_type %in% cell_types) %>%
  mutate(Age = factor(Age, levels = age_levels))
cell_level_dyscoordination_cleaned <- cell_level_dyscoordination_cleaned %>%
  filter(cell_type %in% cell_types) %>%
  mutate(Age = factor(Age, levels = age_levels))

# Compute log-fold change of gene-level dyscoordination relative to each cell type's youngest age
lfc_data <- data.frame()
for (ct in cell_types) {
  df_ct <- gene_level_dyscoordination_cleaned %>%
    filter(cell_type == ct) %>%
    arrange(Age)
  # Skip if less than 2 age groups present
  if (length(unique(df_ct$Age)) == 1) next
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
  filter(LFC >= quantile(LFC, 0.01, na.rm = TRUE),
         LFC <= quantile(LFC, 0.99, na.rm = TRUE)) %>%
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

ggsave("aging_rat_kidney_gene_level_dyscoordination.png", p_gene_level_dyscoordination, width = 8, height = 4)

signif_comparisons <- lapply(seq_along(age_levels[-length(age_levels)]),
                             function(i) c(age_levels[i], age_levels[i + 1]))

p_cell_level_dyscoordination <- cell_level_dyscoordination_cleaned %>%
  filter(Cell_level_deviation >= quantile(Cell_level_deviation, 0.01, na.rm = TRUE),
         Cell_level_deviation <= quantile(Cell_level_deviation, 0.99, na.rm = TRUE)) %>%
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

ggsave("aging_rat_kidney_cell_level_dyscoordination.png", p_cell_level_dyscoordination, width = 8, height = 4)

######################################################
# 2. Technical covariate analysis
cov_list <- c(nCount_RNA = "Sequencing depth",
              nFeature_RNA = "Genes detected",
              percent.mt = "Mitochondrial reads (%)",
              percent.ribo = "Ribosomal reads (%)",
              S.Score = "Cell cycle S score",
              G2M.Score = "Cell cycle G2/M score")

# Load in data
obj <- readRDS("./step2_anno.rds")
obj <- subset(obj, subset = nFeature_RNA >= 250 & nCount_RNA >= 2500)
pt_cells <- rownames(obj@meta.data)[obj@meta.data$celltype_refined2 %in% cell_types]
obj <- obj[, pt_cells]

# Compute QC metrics
DefaultAssay(obj) <- "RNA"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^[Mm][Tt]-")
obj[["percent.ribo"]] <- obj$percent.rpl + obj$percent.rps
s_genes <- intersect(str_to_title(cc.genes.updated.2019$s.genes), rownames(obj))
g2m_genes <- intersect(str_to_title(cc.genes.updated.2019$g2m.genes), rownames(obj))
obj <- CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

kidney_metadata <- obj@meta.data %>%
  as.data.frame() %>%
  select(nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, S.Score, G2M.Score) %>%
  mutate(Cell_barcode = rownames(obj@meta.data))

rm(obj)

cell_level_cov <- cell_level_dyscoordination_cleaned %>%
  filter(Cell_level_deviation > 0) %>%
  inner_join(kidney_metadata, by = "Cell_barcode") %>%
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
  ggsave(paste0("aging_rat_kidney_covariate_scatter_", gsub("\\.", "_", v), ".png"),
         plot_cov_scatter(cell_level_cov, v, cov_list[[v]]),
         width = 8, height = 4)
}

# Compute corrected cell-level transcriptional dyscoordination with the following covariates
correct_covariates <- c("nCount_RNA", "nFeature_RNA")
fmla <- as.formula(paste("log_deviation ~", paste(correct_covariates, collapse = " + ")))

# Regress within each cell type (age effect is preserved across groups)
cell_level_cov$log_deviation_corrected <- unsplit(
  lapply(split(cell_level_cov, cell_level_cov$cell_type), function(g) {
    if (nrow(g) < max(5, length(correct_covariates) + 2))
      return(rep(NA_real_, nrow(g)))
    residuals(lm(fmla, data = g))
  }),
  cell_level_cov$cell_type 
)

# Violin plot for corrected cell-level transcriptional dyscoordination
p_cell_level_corrected <- cell_level_cov %>%
  filter(!is.na(log_deviation_corrected),
         log_deviation_corrected >= quantile(log_deviation_corrected, 0.01, na.rm = TRUE),
         log_deviation_corrected <= quantile(log_deviation_corrected, 0.99, na.rm = TRUE)) %>%
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

ggsave("aging_rat_kidney_cell_level_dyscoordination_corrected.png", p_cell_level_corrected, width = 8, height = 4)

######################################################