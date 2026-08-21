# Hepatocyte Zonation and Transcriptional Dyscoordination (Aging Mouse Liver)

# Online links
# https://www.gsea-msigdb.org/gsea/msigdb

# Set up working directory
# setwd("path/to/working/directory")

library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(Seurat)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

data_dir <- 'path/to/C57BL6J_liver_hippocampus'

######################################################
# 1. Zonation markers (21 markers) and per-marker-type gene-level dyscoordination (Fig 7C)
zonation_markers <- readLines(file.path(data_dir, "zonation_genes.txt"))
marker_labels <- c(
  "Cyp1a2" = "Pericentral", "Gsta3" = "Pericentral", "Cyp2f2" = "Pericentral",
  "Cyp27a1" = "Pericentral", "Mup17" = "Pericentral", "Nt5e" = "Pericentral",
  "Cyp1a1" = "Pericentral", "Gstm3" = "Pericentral", "Oat" = "Pericentral",
  "Igfbp2" = "Mid", "Hamp" = "Mid", "Hamp2" = "Mid", "Ccnd1" = "Mid",
  "Alb" = "Periportal", "Cyp2e1" = "Periportal", "Asl" = "Periportal",
  "Gls2" = "Periportal", "Cdh1" = "Periportal", "Cps1" = "Periportal",
  "Pck1" = "Periportal", "Sdhd" = "Periportal", "Mup3" = "Periportal"
)
marker_type <- ifelse(zonation_markers %in% names(marker_labels), marker_labels[zonation_markers], "Unknown")
zonation_markers_info <- data.frame(gene = zonation_markers, marker_type = marker_type)
zonation_markers_info$marker_type <- factor(zonation_markers_info$marker_type, levels = c("Periportal", "Mid", "Pericentral"))

liver_aging_hep_only_gene_level <- read.csv("liver_data_hepatocyte_only_estimated_dispersion_SAVER_v2.csv")
liver_aging_hep_only_gene_level_filtered <- merge(liver_aging_hep_only_gene_level, zonation_markers_info, by.x = "Gene", by.y = "gene")
liver_aging_hep_only_gene_level_filtered <- liver_aging_hep_only_gene_level_filtered %>%
  mutate(Age = factor(Age, levels = c("young", "old"))) %>%
  arrange(marker_type, Gene, Age)

bootstrap_mean <- function(v, R = 2000L) {
  v <- na.omit(v)
  n <- length(v)
  if (n <= 1) return(list(mean = mean(v), lo = NA_real_, hi = NA_real_))
  idx <- replicate(R, sample.int(n, n, replace = TRUE))
  means <- colMeans(matrix(v[idx], nrow = n))
  list(mean = mean(v), low = quantile(means, 0.025), high = quantile(means, 0.975))
}

marker_type_ci <- liver_aging_hep_only_gene_level_filtered %>%
  group_by(marker_type, Age) %>%
  summarise(ci = list(bootstrap_mean(log(Gene_level_deviation))), .groups = "drop") %>%
  unnest_wider(ci)

liver_aging_hep_only_gene_level_filtered_labels <- liver_aging_hep_only_gene_level_filtered %>%
  group_by(marker_type, Gene) %>%
  summarise(delta = diff(log(Gene_level_deviation)), .groups = "drop") %>%
  filter(delta > 0) %>%
  left_join(liver_aging_hep_only_gene_level_filtered %>% filter(Age == "old"), by = c("marker_type", "Gene"))

# NOTE: `pal` is used in the faceted plot below; defined here (source relied on a
# commented definition). Colours match the source's intended palette.
pal <- c("Pericentral" = "#2C7BB6", "Mid" = "#7F7F7F", "Periportal" = "#D7191C")

p_zonation_marker_gene_level <- ggplot(liver_aging_hep_only_gene_level_filtered,
       aes(x = Age, y = log(Gene_level_deviation), group = Gene, color = marker_type)) +
  geom_line(alpha = 0.5, linewidth = 0.5) +
  geom_point(alpha = 0.75, size = 1.5) +
  geom_line(data = marker_type_ci, aes(x = Age, y = mean, group = marker_type), inherit.aes = FALSE, linewidth = 1) +
  geom_errorbar(data = marker_type_ci, aes(x = Age, ymin = low, ymax = high), inherit.aes = FALSE, width = 0.1, linewidth = 0.75, alpha = 0.75) +
  geom_point(data = marker_type_ci, aes(x = Age, y = mean), inherit.aes = FALSE, size = 2, alpha = 0.75) +
  ggrepel::geom_text_repel(data = liver_aging_hep_only_gene_level_filtered_labels, aes(label = Gene),
                           nudge_x = 0.25, hjust = 0, direction = "y", force = 1, size = 1.75,
                           box.padding = 0.15, segment.color = NA, show.legend = FALSE) +
  scale_color_manual(values = pal, guide = guide_legend(title = "Zonation marker")) +
  facet_wrap(~ marker_type, scales = "free_y") +
  labs(x = NULL, y = "Gene-level dyscoordination (log-transformed)") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))

ggsave("Peter_liver_aging_zonation_marker_gene_level_dyscoordination_summary.pdf", plot = p_zonation_marker_gene_level, width = 5, height = 3.5)

# Combined overall change with top-5 gainers labelled
marker_ci_overall <- liver_aging_hep_only_gene_level_filtered %>%
  group_by(Age) %>%
  summarise(ci = list(bootstrap_mean(log(Gene_level_deviation))), .groups = "drop") %>%
  unnest_wider(ci)

liver_aging_hep_only_gene_level_filtered_labels <- liver_aging_hep_only_gene_level_filtered %>%
  group_by(Gene) %>%
  summarise(delta = diff(log(Gene_level_deviation)), .groups = "drop") %>%
  arrange(desc(delta)) %>% slice_head(n = 5) %>%
  left_join(liver_aging_hep_only_gene_level_filtered %>% filter(Age == "old"), by = "Gene")

pal <- c("Periportal" = "#1B9E77", "Mid" = "#7570B3", "Pericentral" = "#E6AB02")
p_zonation_marker_gene_level_combined <- ggplot(liver_aging_hep_only_gene_level_filtered,
                                                aes(x = Age, y = log(Gene_level_deviation), group = Gene)) +
  geom_line(alpha = 0.5, linewidth = 0.5, aes(color = marker_type)) +
  geom_point(alpha = 0.75, size = 1.5, aes(color = marker_type)) +
  geom_line(data = marker_ci_overall, aes(x = Age, y = mean, group = 1), inherit.aes = FALSE, linewidth = 1, color = "black") +
  geom_errorbar(data = marker_ci_overall, aes(x = Age, ymin = low, ymax = high), inherit.aes = FALSE, width = 0.1, linewidth = 0.75, alpha = 0.75, color = "black") +
  geom_point(data = marker_ci_overall, aes(x = Age, y = mean), inherit.aes = FALSE, size = 2, alpha = 0.75, color = "black") +
  facet_wrap(~ marker_type, ncol = 3) +
  ggrepel::geom_text_repel(data = liver_aging_hep_only_gene_level_filtered_labels, aes(label = Gene, color = marker_type),
                           nudge_x = 0.25, hjust = 0, direction = "y", force = 1, size = 1.75,
                           box.padding = 0.15, segment.color = NA, show.legend = FALSE) +
  scale_color_manual(values = pal, guide = guide_legend(title = "Zonation marker")) +
  labs(x = NULL, y = "Gene-level dyscoordination (log-transformed)") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave("Peter_liver_aging_zonation_marker_gene_level_dyscoordination_combined.pdf", plot = p_zonation_marker_gene_level_combined, width = 4.5, height = 3.5)

######################################################
# 2. Hepatocyte zonation score per cell
genes_midlobular  <- c("Hamp", "Hamp2", "Ccnd1", "Cyp8b1")
genes_pericentral <- c("Cyp2e1", "Gsta3", "Cyp27a1", "Mup17", "Nt5e", "Axin2", "Cyp1a2", "Oat", "Gstm3", "Axin2", "Lgr5")
genes_periportal  <- c("Alb", "Cyp2f2", "Asl", "Gls2", "Cdh1", "Cps1", "Pck1", "Sdhd")
genes_housekeeping <- c("Gapdh", "Actb", "Rplp0", "Vcl", "L3mbtl2", "Rbck1", "Vamp7", "Wdr55", "Hprt")
genes_immune       <- c("Cd45", "Cd68", "Cd64", "Cd3", "Cd4", "Cd8")
genes_zonation     <- c(genes_midlobular, genes_pericentral, genes_periportal)
genes_marker       <- c(genes_midlobular, genes_pericentral, genes_periportal, genes_housekeeping, genes_immune)

expr_sum <- function(gene_set, expr_mat) {
  gene_set <- intersect(gene_set, rownames(expr_mat))
  colSums(expr_mat[gene_set, ])
}

# Hepatocyte expression matrix (genes x cells), computed upstream
load("liver_data_hepatocyte_only_expr.RData")

expr_midlobular   <- expr_sum(genes_midlobular, liver_hep_only_expr)
expr_hamp2        <- liver_hep_only_expr["Hamp2", ]
expr_pericentral  <- expr_sum(genes_pericentral, liver_hep_only_expr)
expr_periportal   <- expr_sum(genes_periportal, liver_hep_only_expr)
expr_zonation     <- expr_sum(genes_zonation, liver_hep_only_expr)
expr_housekeeping <- expr_sum(genes_housekeeping, liver_hep_only_expr)
expr_immune       <- expr_sum(genes_immune, liver_hep_only_expr)
expr_all_marker   <- expr_sum(genes_marker, liver_hep_only_expr)

liver_aging_hep_only_cell_level <- read.csv("liver_data_hepatocyte_only_cellular_dispersion_SAVER_v2.csv")
liver_aging_hep_only_cell_level <- liver_aging_hep_only_cell_level[match(colnames(liver_hep_only_expr), liver_aging_hep_only_cell_level$Cell_barcode), ]
liver_aging_hep_only_cell_level_filtered <- liver_aging_hep_only_cell_level %>%
  mutate(expr_midlobular = expr_midlobular, expr_hamp2 = expr_hamp2,
         expr_pericentral = expr_pericentral, expr_periportal = expr_periportal,
         expr_zonation = expr_zonation, expr_housekeeping = expr_housekeeping,
         expr_immune = expr_immune, expr_all_marker = expr_all_marker,
         expr_zonation_score = expr_pericentral / (expr_periportal + expr_pericentral))

# Bin cells by zonation score, assign Periportal / Midlobular / Pericentral regions
liver_aging_hep_only_cell_level_filtered <- liver_aging_hep_only_cell_level_filtered[
  liver_aging_hep_only_cell_level_filtered$expr_zonation_score > 0 &
  liver_aging_hep_only_cell_level_filtered$expr_zonation_score < 1, ]
zonation_quantile_bins <- quantile(liver_aging_hep_only_cell_level_filtered$expr_zonation_score,
                                   probs = seq(0, 1, length.out = 10 + 1), na.rm = TRUE)

liver_aging_hep_only_cell_level_filtered_meta <- liver_aging_hep_only_cell_level_filtered %>%
  mutate(zonation_bin = cut(expr_zonation_score, breaks = zonation_quantile_bins,
                            include.lowest = TRUE, labels = paste0(seq(10, 100, by = 10), "%")),
         zonation_region = case_when(zonation_bin %in% c("10%", "20%") ~ "Periportal",
                                     zonation_bin %in% c("50%", "60%") ~ "Midlobular",
                                     zonation_bin %in% c("90%", "100%") ~ "Pericentral",
                                     TRUE ~ NA_character_))
write.csv(liver_aging_hep_only_cell_level_filtered_meta, "Peter_liver_aging_zonation_score_bins.csv", row.names = FALSE)

######################################################
# 3. Binned summary of cell-level dyscoordination vs zonation score (Fig 7C)
high_dyscoordination_cutoff <- quantile(liver_aging_hep_only_cell_level_filtered$Cell_level_deviation, 0.95)
liver_aging_hep_only_cell_level_filtered_binned_summary <- liver_aging_hep_only_cell_level_filtered %>%
  filter(Cell_level_deviation <= 50) %>%
  mutate(zonation_bin = cut(expr_zonation_score, breaks = zonation_quantile_bins,
                            include.lowest = TRUE, labels = paste0(seq(10, 100, by = 10), "%")),
         high_dyscoordination = Cell_level_deviation > high_dyscoordination_cutoff) %>%
  group_by(zonation_bin) %>%
  summarise(mean_dyscoordination = mean(Cell_level_deviation, na.rm = TRUE),
            percent_high_dyscoordination = mean(high_dyscoordination, na.rm = TRUE) * 100,
            n = n(), .groups = "drop")

lo <- min(liver_aging_hep_only_cell_level_filtered_binned_summary$mean_dyscoordination, na.rm = TRUE) - 0.15
hi <- max(liver_aging_hep_only_cell_level_filtered_binned_summary$mean_dyscoordination, na.rm = TRUE) + 0.15
y_base <- lo
p_max <- max(liver_aging_hep_only_cell_level_filtered_binned_summary$percent_high_dyscoordination, na.rm = TRUE)
k <- 0.8 * (hi - lo) / p_max

p_binned_summary <- ggplot(liver_aging_hep_only_cell_level_filtered_binned_summary, aes(x = zonation_bin)) +
  geom_linerange(aes(ymin = y_base, ymax = y_base + k * percent_high_dyscoordination, color = "Mean dyscoordination"),
                 linewidth = 10, alpha = 1) +
  geom_line(aes(y = mean_dyscoordination, color = "% high dyscoordination", group = 1), linewidth = 1) +
  geom_point(aes(y = mean_dyscoordination, color = "% high dyscoordination"), size = 1.8) +
  scale_color_manual(name = NULL, values = c("Mean dyscoordination" = "darkgrey", "% high dyscoordination" = "#009E73")) +
  scale_y_continuous(name = "Mean cell-level dyscoordination", limits = c(lo, hi),
                     sec.axis = sec_axis(~ (. - y_base) / k, name = "% high dyscoordination")) +
  labs(x = "Zonation score bin") +
  theme_minimal()

ggsave("Peter_liver_aging_zonation_score_cell_level_dyscoordination_summary.png", plot = p_binned_summary, width = 7.5, height = 4, dpi = 300)

######################################################
# 4. Per-zone gene-level dyscoordination change with age (Fig 7D)
# Build the zonation stratum vector for the zonation-stratified dispersion run.
liver_aging_hep_only_cell_level_filtered_v2 <- liver_aging_hep_only_cell_level %>%
  mutate(expr_zonation_score = expr_pericentral / (expr_periportal + expr_pericentral)) %>%
  left_join(liver_aging_hep_only_cell_level_filtered_meta)
liver_aging_hep_only_cell_level_filtered_v2$zonation_region2 <- ifelse(
  is.na(liver_aging_hep_only_cell_level_filtered_v2$zonation_region),
  "Unknown", liver_aging_hep_only_cell_level_filtered_v2$zonation_region)
dataset.celltype <- as.numeric(factor(liver_aging_hep_only_cell_level_filtered_v2$zonation_region2,
                                      levels = c("Periportal", "Midlobular", "Pericentral", "Unknown")))
dataset.celltype.levels <- c("Periportal", "Midlobular", "Pericentral", "Unknown")
# save(dataset.celltype, dataset.celltype.levels, file = "liver_data_zonation_celltype.RData")

# Per-zone gene-level dyscoordination CSV comes from the zonation-stratified run.
liver_hep_only_gene_level_zonation <- read.csv("liver_data_hepatocyte_only_estimated_dispersion_zonation_SAVER.csv")
liver_hep_only_gene_level_zonation_lfc <- liver_hep_only_gene_level_zonation %>%
  select(Gene, Age, cell_type, Gene_level_deviation) %>%
  pivot_wider(names_from = Age, values_from = Gene_level_deviation) %>%
  mutate(logFC = log((old) / (young))) %>%
  drop_na() %>% filter(is.finite(logFC))

liver_hep_only_gene_level_zonation_lfc_pos <- liver_hep_only_gene_level_zonation_lfc[liver_hep_only_gene_level_zonation_lfc$logFC > 0, ]
liver_hep_only_gene_level_zonation_lfc_pos$cell_type <- factor(liver_hep_only_gene_level_zonation_lfc_pos$cell_type, levels = c("Periportal", "Midlobular", "Pericentral"))
write.csv(liver_hep_only_gene_level_zonation_lfc_pos, "liver_data_hepatocyte_only_estimated_dispersion_zonation_LFC_positive.csv", row.names = FALSE)

######################################################
