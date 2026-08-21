# EpiTrace Mitotic Age vs Transcriptional Dyscoordination (Aging Mouse Liver)

# Online links
# https://github.com/MagpiePKU/EpiTrace

# Set up working directory
# setwd("path/to/working/directory")

library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(Matrix)
library(EpiTrace)
library(Seurat)
library(SeuratObject)
library(easyLift)
library(plyranges)
library(Signac)
library(rtracklayer)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Load integrated liver object + shared ATAC peaks
data_dir     <- 'path/to/C57BL6J_liver_hippocampus'
epitrace_dir <- 'path/to/EpiTrace'

load(file.path(data_dir, "liver_data_integrated.RData"))
load(file.path(data_dir, "shared_peaks_PetersLab_liver_ATAC.RData"))

metadata_df <- liver_data_integrated@meta.data
metadata_df$cell <- rownames(metadata_df)
merged_obj$cell <- rownames(merged_obj@meta.data)
merged_obj@meta.data <- merge(merged_obj@meta.data, metadata_df[, c(6:11)], by.x = 'cell', by.y = 'cell')
merged_obj@meta.data <- merged_obj@meta.data[match(colnames(merged_obj@assays$ATAC@counts), merged_obj@meta.data$cell), ]
rm(liver_data_integrated)

# Mouse clock lifted from the human clock (mm10)
human_clock <- plyranges::reduce_ranges(c(clock_gr_list[[1]], clock_gr_list[[2]]))
mouse_clock_by_human <- readRDS(file.path(epitrace_dir, 'mouse_clock_lifted_human_to_mouse_mm10.rds'))

######################################################
# Initialise EpiTrace inputs; keep non-empty peaks/cells
init_mm <- Init_Matrix(peakname = paste(seqnames(merged_obj@assays$ATAC@ranges), start(merged_obj@assays$ATAC@ranges), end(merged_obj@assays$ATAC@ranges), sep = "-"),
                       cellname = colnames(merged_obj@assays$ATAC@counts),
                       matrix = merged_obj@assays$ATAC@counts)
init_gr <- Init_Peakset(merged_obj@assays$ATAC@ranges[rowSums(init_mm) > 0])
init_mm <- init_mm[rowSums(init_mm) > 0, ]
init_mm_hep_only <- init_mm[, merged_obj$annotation == "Hepatocyte"]

merged_metadata <- merged_obj@meta.data[colSums(init_mm) > 0, ]
init_mm <- init_mm[, colSums(init_mm) > 0]
init_mm_hep_only <- init_mm[, merged_metadata$annotation == "Hepatocyte"]

######################################################
# EpiTrace convergence on all cells
aging_liver_epitrace_obj_age_conv_human_clock <- EpiTraceAge_Convergence(peakSet = init_gr,
                                                                         matrix = init_mm,
                                                                         ref_genome = 'mm10',
                                                                         clock_gr = mouse_clock_by_human,
                                                                         iterative_time = 5,
                                                                         min.cutoff = 0,
                                                                         non_standard_clock = T,
                                                                         qualnum = 10,
                                                                         ncore_lim = 1,
                                                                         select_minimal_percentage = 0.0075,
                                                                         mean_error_limit = 0.1)
aging_liver_epitrace_obj_age_conv_human_clock@meta.data <- merge(aging_liver_epitrace_obj_age_conv_human_clock@meta.data, metadata_df[, c(6:11)], by.x = "cell", by.y = "cell", all.x = TRUE)
aging_liver_epitrace_obj_age_conv_human_clock$age <- factor(aging_liver_epitrace_obj_age_conv_human_clock$age, levels = c("young", "old"))
# save(aging_liver_epitrace_obj_age_conv_human_clock, file = "Peter_liver_aging_EpiTrace_human_clock.RData")

######################################################
# EpiTrace convergence on hepatocytes only
aging_liver_epitrace_obj_age_conv_human_clock_hep_only <- EpiTraceAge_Convergence(peakSet = init_gr,
                                                                                 matrix = init_mm_hep_only,
                                                                                 ref_genome = 'mm10',
                                                                                 clock_gr = mouse_clock_by_human,
                                                                                 iterative_time = 5,
                                                                                 min.cutoff = 0,
                                                                                 non_standard_clock = T,
                                                                                 qualnum = 10,
                                                                                 ncore_lim = 1,
                                                                                 select_minimal_percentage = 0.01,
                                                                                 mean_error_limit = 0.1)
aging_liver_epitrace_obj_age_conv_human_clock_hep_only@meta.data <- merge(aging_liver_epitrace_obj_age_conv_human_clock_hep_only@meta.data, metadata_df[, c(6:11)], by.x = "cell", by.y = "cell", all.x = TRUE)
aging_liver_epitrace_obj_age_conv_human_clock_hep_only$age <- factor(aging_liver_epitrace_obj_age_conv_human_clock_hep_only$age, levels = c("young", "old"))
# save(aging_liver_epitrace_obj_age_conv_human_clock_hep_only, file = "Peter_liver_aging_EpiTrace_human_clock_hep_only.RData")
# load("Peter_liver_aging_EpiTrace_human_clock_hep_only.RData")

######################################################
# EpiTrace mitotic age distribution for hepatocytes (young vs old)
p_epitrace_iter <- ggplot(aging_liver_epitrace_obj_age_conv_human_clock_hep_only@meta.data, aes(x = age, y = EpiTraceAge_iterative, fill = age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("young" = "#e41a1c", "old" = "#377eb8")) +
  theme_minimal() +
  xlab('Age') +
  ylab('EpiTrace mitotic age score') +
  labs(fill = "Age") +
  geom_signif(comparisons = list(c("young", "old")), map_signif_level = TRUE, step_increase = 0.1,
              test = function(x, y) wilcox.test(x, y, alternative = "less"))

ggsave("Peter_liver_aging_EpiTrace_iterative_hep_only_filtered.jpg", plot = p_epitrace_iter, width = 6, height = 4.5)

######################################################
# EpiTrace age vs cell-level dyscoordination (Fig 7B)
liver_aging_hep_only_cell_level <- read.csv("liver_data_hepatocyte_only_cellular_dispersion_SAVER_v2.csv")
liver_aging_hep_only_cell_level <- merge(liver_aging_hep_only_cell_level,
                                         aging_liver_epitrace_obj_age_conv_human_clock_hep_only@meta.data[, c(1, 7, 10)],
                                         by.x = "Cell_barcode", by.y = "cell")

cor.test(liver_aging_hep_only_cell_level$Cell_level_deviation, liver_aging_hep_only_cell_level$EpiTraceAge_Clock_initial)
cor.test(liver_aging_hep_only_cell_level$Cell_level_deviation, liver_aging_hep_only_cell_level$EpiTraceAge_iterative)

p_epitrace_init_cell_level_hep_only <- ggplot(liver_aging_hep_only_cell_level[liver_aging_hep_only_cell_level$Cell_level_deviation < 50, ], aes(x = EpiTraceAge_Clock_initial, y = log(Cell_level_deviation), color = Age)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "EpiTrace age (initial)", y = "Cell-level dyscoordination (log-transformed)", color = "Age") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_cor(method = "pearson", size = 4, show.legend = FALSE)

p_epitrace_iter_cell_level_hep_only <- ggplot(liver_aging_hep_only_cell_level[liver_aging_hep_only_cell_level$Cell_level_deviation < 50, ], aes(x = EpiTraceAge_iterative, y = log(Cell_level_deviation), color = Age)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "EpiTrace age (iterative)", y = "Cell-level dyscoordination (log-transformed)", color = "Age") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_cor(method = "pearson", size = 4, show.legend = FALSE)

ggsave("Peter_liver_aging_EpiTrace_initial_cell_level_dyscoordination_hep_only.jpg", plot = p_epitrace_init_cell_level_hep_only, width = 8, height = 6)
ggsave("Peter_liver_aging_EpiTrace_iterative_cell_level_dyscoordination_hep_only.jpg", plot = p_epitrace_iter_cell_level_hep_only, width = 8, height = 6)

######################################################
