# Transcriptional dyscoordination analysis for aging human bone marrow

# Online links
# https://insight.jci.org/articles/view/124928
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120446

# Set up working directory
# setwd("S:/Penn Dropbox/Eddie Yang/Aging/Scripts/Human_BM")

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
gene_level_dyscoordination <- read.csv("human_BM_estimated_dispersion_SAVER_v2.csv")
cell_level_dyscoordination <- read.csv("human_BM_cellular_dispersion_SAVER_v2.csv")

gene_level_dyscoordination$Gene_level_deviation <- remove_outliers(gene_level_dyscoordination$Gene_level_deviation)
cell_level_dyscoordination$Cell_level_deviation <- remove_outliers(cell_level_dyscoordination$Cell_level_deviation)
gene_level_dyscoordination_cleaned <- na.omit(gene_level_dyscoordination)
cell_level_dyscoordination_cleaned <- na.omit(cell_level_dyscoordination)

# Capitalize cell type names
cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
gene_level_dyscoordination_cleaned$Celltype <- cap_first(gene_level_dyscoordination_cleaned$Celltype)
cell_level_dyscoordination_cleaned$Celltype <- cap_first(cell_level_dyscoordination_cleaned$Celltype)

age_levels <- c("46-50", "55-57", "58-60", "66-67", "84")
old_ages <- age_levels[-1]
cell_types <- sort(unique(gene_level_dyscoordination_cleaned$Celltype))

gene_level_dyscoordination_cleaned <- gene_level_dyscoordination_cleaned %>%
  mutate(Age = factor(Age, levels = age_levels))
cell_level_dyscoordination_cleaned <- cell_level_dyscoordination_cleaned %>%
  mutate(Age = factor(Age, levels = age_levels))

# Compute log-fold change of gene-level dyscoordination relative to each cell type's youngest age
lfc_data <- data.frame()
for (ct in cell_types) {
  df_ct <- gene_level_dyscoordination_cleaned %>% filter(Celltype == ct)
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

ggsave("aging_human_BM_gene_level_dyscoordination.png", p_gene_level_dyscoordination, width = 8, height = 10)

signif_comparisons <- lapply(seq_along(age_levels[-length(age_levels)]),
                             function(i) c(age_levels[i], age_levels[i + 1]))

p_cell_level_dyscoordination <- cell_level_dyscoordination_cleaned %>%
  filter(Cell_level_deviation >= quantile(Cell_level_deviation, 0.01, na.rm = TRUE),
         Cell_level_deviation <= quantile(Cell_level_deviation, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Age, y = log(Cell_level_deviation), fill = Age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Reds") +
  facet_wrap(~ Celltype, ncol = 3) +
  theme_minimal() +
  labs(x = "Age", y = "Cell-level dyscoordination (log-transformed)") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "less"), step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("aging_human_BM_cell_level_dyscoordination.png", p_cell_level_dyscoordination, width = 8, height = 10)

######################################################
# 2. Technical covariate analysis
cov_list <- c(nCount_RNA = "Sequencing depth",
              nFeature_RNA = "Genes detected",
              percent.mt = "Mitochondrial reads (%)",
              percent.ribo = "Ribosomal reads (%)",
              S.Score = "Cell cycle S score",
              G2M.Score = "Cell cycle G2/M score")

# Load in data
bm <- readRDS(file.path(data_dir, "./human_BM_Oetjen_annotated.rds"))
DefaultAssay(bm) <- "RNA"

# Compute QC metrics
bm[["percent.mt"]] <- PercentageFeatureSet(bm, pattern = "^MT-")
bm[["percent.ribo"]] <- PercentageFeatureSet(bm, pattern = "^RP[SL]")
s_genes <- intersect(cc.genes.updated.2019$s.genes, rownames(bm))
g2m_genes <- intersect(cc.genes.updated.2019$g2m.genes, rownames(bm))
bm <- CellCycleScoring(bm, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

bm_metadata <- bm@meta.data %>%
  as.data.frame() %>%
  select(nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, S.Score, G2M.Score) %>%
  mutate(Cell_barcode = rownames(bm@meta.data))

rm(bm)

cell_level_cov <- cell_level_dyscoordination_cleaned %>%
  filter(Cell_level_deviation > 0) %>%
  inner_join(bm_metadata, by = "Cell_barcode") %>%
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
    facet_wrap(~ Celltype, ncol = 3) +
    scale_fill_brewer(palette = "Blues") +
    theme_minimal() +
    labs(x = x_label, y = "Cell-level dyscoordination (log-transformed)")
}

for (v in names(cov_list)) {
  ggsave(paste0("aging_human_BM_covariate_scatter_", gsub("\\.", "_", v), ".png"),
         plot_cov_scatter(cell_level_cov, v, cov_list[[v]]),
         width = 8, height = 10)
}

# Compute corrected cell-level transcriptional dyscoordination with the following covariates
correct_covariates <- c("nCount_RNA", "nFeature_RNA")
fmla <- as.formula(paste("log_deviation ~", paste(correct_covariates, collapse = " + ")))

# Regress within each cell type (age effect is preserved across groups)
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
  ggplot(aes(x = Age, y = log_deviation_corrected, fill = Age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Reds") +
  facet_wrap(~ Celltype, ncol = 3) +
  theme_minimal() +
  labs(x = "Age", y = "Cell-level dyscoordination (log-transformed, corrected)") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "less"),
              step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("aging_human_BM_cell_level_dyscoordination_corrected.png", p_cell_level_corrected, width = 8, height = 10)

######################################################