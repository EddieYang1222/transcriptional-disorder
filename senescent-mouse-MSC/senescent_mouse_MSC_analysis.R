# Transcriptional dyscoordination analysis for senescent mouse MSC (SC and FAP)

# Online links
# https://www.cell.com/iscience/fulltext/S2589-0042(22)00118-3
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169531

# Set up working directory
# setwd("S:/Penn Dropbox/Eddie Yang/Aging/Scripts/Mouse_MSC")

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
gene_level_dyscoordination <- rbind(
  read.csv("mouse_MSC_SC_estimated_dispersion_SAVER.csv")  %>% mutate(Celltype = "SC"),
  read.csv("mouse_MSC_FAP_estimated_dispersion_SAVER.csv") %>% mutate(Celltype = "FAP")
)
cell_level_dyscoordination <- rbind(
  read.csv("mouse_MSC_SC_cellular_dispersion_SAVER.csv")  %>% mutate(Celltype = "SC"),
  read.csv("mouse_MSC_FAP_cellular_dispersion_SAVER.csv") %>% mutate(Celltype = "FAP")
)

gene_level_dyscoordination$Gene_level_deviation <- remove_outliers(gene_level_dyscoordination$Gene_level_deviation)
cell_level_dyscoordination$Cell_level_deviation <- remove_outliers(cell_level_dyscoordination$Cell_level_deviation)
gene_level_dyscoordination_cleaned <- na.omit(gene_level_dyscoordination)
cell_level_dyscoordination_cleaned <- na.omit(cell_level_dyscoordination)

# Capitalize cell type names
cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
gene_level_dyscoordination_cleaned$Celltype <- cap_first(gene_level_dyscoordination_cleaned$Celltype)
cell_level_dyscoordination_cleaned$Celltype <- cap_first(cell_level_dyscoordination_cleaned$Celltype)

cond_levels <- c("PBS", "DOXO", "DABT")
baseline <- "PBS"
non_baseline <- setdiff(cond_levels, baseline)
cell_types <- sort(unique(gene_level_dyscoordination_cleaned$Celltype))

gene_level_dyscoordination_cleaned <- gene_level_dyscoordination_cleaned %>%
  mutate(Condition = factor(Condition, levels = cond_levels))
cell_level_dyscoordination_cleaned <- cell_level_dyscoordination_cleaned %>%
  mutate(Condition = factor(Condition, levels = cond_levels))

# Compute log-fold change of gene-level dyscoordination relative to PBS baseline per cell type
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
           Cell_type = paste0(ct, " (", baseline, ")")) %>%
    filter(is.finite(LFC))
  lfc_data <- rbind(lfc_data, df_ct_old)
}

# Violin plots for transcriptional dyscoordination
signif_comparisons <- lapply(seq_along(non_baseline[-length(non_baseline)]),
                             function(i) c(non_baseline[i], non_baseline[i + 1]))

p_gene_level_dyscoordination <- lfc_data %>%
  filter(LFC >= quantile(LFC, 0.01, na.rm = TRUE),
         LFC <= quantile(LFC, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Condition, y = LFC, fill = Condition)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Blues") +
  facet_wrap(~ Cell_type, ncol = 2) +
  theme_minimal() +
  labs(x = "Condition", y = "Log-fold change of gene-level dyscoordination") +
  geom_signif(comparisons = signif_comparisons, map_signif_level = TRUE,
              test = function(x, y) wilcox.test(x, y, alternative = "greater"), step_increase = 0.1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

ggsave("senescent_mouse_MSC_gene_level_dyscoordination.png", p_gene_level_dyscoordination, width = 6, height = 4)

p_cell_level_dyscoordination <- cell_level_dyscoordination_cleaned %>%
  filter(Cell_level_deviation >= quantile(Cell_level_deviation, 0.01, na.rm = TRUE),
         Cell_level_deviation <= quantile(Cell_level_deviation, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = Condition, y = log(Cell_level_deviation), fill = Condition)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Reds") +
  facet_wrap(~ Celltype, ncol = 2) +
  theme_minimal() +
  labs(x = "Condition", y = "Cell-level dyscoordination (log-transformed)") +
  # PBS -> DOXO: senescence induction, dyscoordination increases
  geom_signif(comparisons = list(c("PBS", "DOXO")),
              test = function(x, y) wilcox.test(x, y, alternative = "less"),
              map_signif_level = TRUE, step_increase = 0.1) +
  # DOXO -> DABT: senolytic reversal, dyscoordination decreases
  geom_signif(comparisons = list(c("DOXO", "DABT")),
              test = function(x, y) wilcox.test(x, y, alternative = "greater"),
              map_signif_level = TRUE, step_increase = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("senescent_mouse_MSC_cell_level_dyscoordination.png", p_cell_level_dyscoordination, width = 6, height = 4)

######################################################
# 2. Technical covariate analysis
cov_list <- c(nCount_RNA = "Sequencing depth",
              nFeature_RNA = "Genes detected",
              percent.mt = "Mitochondrial reads (%)",
              percent.ribo = "Ribosomal reads (%)",
              S.Score = "Cell cycle S score",
              G2M.Score = "Cell cycle G2/M score")

# Load in processed data
# data_dir <- "S:/Penn Dropbox/Eddie Yang/Aging/Data/Mouse_MSC"
sc_seu <- readRDS(file.path(data_dir, "mouse_MSC_SC_processed.rds"))
fap_seu <- readRDS(file.path(data_dir, "mouse_MSC_FAP_processed.rds"))

compute_qc <- function(seu) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^Rp[sl]")
  s_genes_mm <- str_to_title(cc.genes.updated.2019$s.genes)
  g2m_genes_mm <- str_to_title(cc.genes.updated.2019$g2m.genes)
  s_found <- intersect(s_genes_mm,   rownames(seu))
  g2m_found <- intersect(g2m_genes_mm, rownames(seu))
  if (length(s_found) >= 5 && length(g2m_found) >= 5)
    seu <- CellCycleScoring(seu, s.features = s_found, g2m.features = g2m_found, set.ident = FALSE)
  else
    seu$S.Score <- seu$G2M.Score <- NA_real_
  seu
}

sc_seu <- compute_qc(sc_seu)
fap_seu <- compute_qc(fap_seu)

extract_meta <- function(seu, celltype_name) {
  m <- seu@meta.data
  for (col in c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "S.Score", "G2M.Score"))
    if (!col %in% colnames(m)) m[[col]] <- NA_real_
  data.frame(
    Cell_barcode = rownames(m),
    Celltype = celltype_name,
    nCount_RNA = m$nCount_RNA,
    nFeature_RNA = m$nFeature_RNA,
    percent.mt = m$percent.mt,
    percent.ribo = m$percent.ribo,
    S.Score = m$S.Score,
    G2M.Score = m$G2M.Score
  )
}

meta_all <- rbind(extract_meta(sc_seu, "SC"), extract_meta(fap_seu, "FAP"))
rm(sc_seu, fap_seu)

cell_level_cov <- cell_level_dyscoordination_cleaned %>%
  filter(Cell_level_deviation > 0) %>%
  inner_join(meta_all, by = c("Cell_barcode", "Celltype")) %>%
  filter(!is.na(Cell_level_deviation)) %>%
  mutate(log_deviation = log(Cell_level_deviation),
         Condition = factor(Condition, levels = cond_levels))

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
    facet_wrap(~ Celltype, ncol = 2) +
    scale_fill_brewer(palette = "Blues") +
    theme_minimal() +
    labs(x = x_label, y = "Cell-level dyscoordination (log-transformed)")
}

for (v in names(cov_list)) {
  ggsave(paste0("senescent_mouse_MSC_covariate_scatter_", gsub("\\.", "_", v), ".png"),
         plot_cov_scatter(cell_level_cov, v, cov_list[[v]]),
         width = 6, height = 4)
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
  facet_wrap(~ Celltype, ncol = 2) +
  theme_minimal() +
  labs(x = "Condition", y = "Cell-level dyscoordination (log-transformed, corrected)") +
  # PBS -> DOXO: senescence induction, dyscoordination increases
  geom_signif(comparisons = list(c("PBS", "DOXO")),
              test = function(x, y) wilcox.test(x, y, alternative = "less"),
              map_signif_level = TRUE, step_increase = 0.1) +
  # DOXO -> DABT: senolytic reversal, dyscoordination decreases
  geom_signif(comparisons = list(c("DOXO", "DABT")),
              test = function(x, y) wilcox.test(x, y, alternative = "greater"),
              map_signif_level = TRUE, step_increase = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("senescent_mouse_MSC_cell_level_dyscoordination_corrected.png",
       p_cell_level_corrected, width = 6, height = 4)

######################################################