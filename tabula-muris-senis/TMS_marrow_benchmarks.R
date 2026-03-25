# Additional benchmark analyses using:
# 1. SenMayo senescence score
# 2. Mean Euclidean distance to cell type average (Decibel)
# 3. Mean Euclidean distance to tissue average using invariant genes (Decibel)
# 4. Scallop membership score

# Online links
# https://www.nature.com/articles/s41467-022-32552-1
# https://elifesciences.org/articles/80380
# https://gitlab.com/olgaibanez/decibel
# https://gitlab.com/olgaibanez/scallop

library(escape)
library(Seurat)
library(ggplot2)

######################################################
# 1. SenMayo senescence score
# Set up working directory
# setwd("S:/Penn Dropbox/Eddie Yang/Aging/Scripts/TMS")

# Load in data
load("./TMS_marrow/TMS_marrow.RData") # TMS marrow data
load("./TMS_marrow/TMS_marrow_manifold.RData") # TMS marrow Manifold

# Remove genes with zero counts, to match the manifold.
message("The initial matrix size is ", nrow(dataset.counts), " genes and ", ncol(dataset.counts), " cells.")
dataset.counts <- dataset.counts[rowSums(dataset.counts) != 0, ] 
message("After removing genes with zero counts, the new matrix size is ",
        nrow(dataset.counts), " genes and ", ncol(dataset.counts)," cells.")

dataset.age.levels <- replace(dataset.age.levels, dataset.age.levels == "1m", "01m")
dataset.age.levels <- replace(dataset.age.levels, dataset.age.levels == "3m", "03m")

# Create a Seurat object
dataset.metadata <- data.frame(row.names = colnames(dataset.counts), ages = dataset.age.levels[dataset.age], celltype = dataset.celltype.levels[dataset.celltype])
dataset.seurat <- CreateSeuratObject(counts = dataset.saver$estimate, meta.data = dataset.metadata)
dataset.seurat <- AddMetaData(dataset.seurat, factor(x= dataset.seurat$ages, levels = dataset.age.levels), col.name = 'ages.factor')
Idents(object = dataset.seurat) <- "celltype"

rm(dataset.counts)
rm(dataset.saver)

######################################################
# Load cell-level deviation
TMS_marrow_cell_level <- read.csv("./TMS_marrow/Dispersion_final/TMS_marrow_cellular_dispersion_SAVER.csv")
TMS_marrow_cell_level <- TMS_marrow_cell_level[match(rownames(dataset.seurat@meta.data), TMS_marrow_cell_level$Cell_barcode), ]
dataset.seurat <- AddMetaData(dataset.seurat, TMS_marrow_cell_level$Cell_level_deviation, col.name = 'entropy')
dataset.seurat.hspc = subset(x = dataset.seurat, idents = "hematopoietic precursor cell")
Idents(object = dataset.seurat.hspc) <- "ages.factor"

# Define the SenMayo gene set
gene.sets <- list(SenMayo_signature = c("Acvr1b", "Ang", "Angpt1", "Angptl4", "Areg", "Axl", "Bex3", "Bmp2",
                                        "Bmp6", "C3", "Ccl1", "Ccl13", "Ccl16", "Ccl2", "Ccl20", "Ccl24", 
                                        "Ccl26", "Ccl3", "Ccl3l1", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Cd55", 
                                        "Cd9", "Csf1", "Csf2", "Csf2rb", "Cst4", "Ctnnb1", "Ctsb", "Cxcl1",
                                        "Cxcl10", "Cxcl12", "Cxcl16", "Cxcl2", "Cxcl3", "Cxcl8", "Cxcr2", "Dkk1",
                                        "Edn1", "Egf", "Egfr", "Ereg", "Esm1", "Ets2", "Fas", "Fgf1",
                                        "Fgf2", "Fgf7", "Gdf15", "Gem", "Gmfg", "Hgf", "Hmgb1", "Icam1",
                                        "Icam3", "Igf1", "Igfbp1", "Igfbp2", "Igfbp3", "Igfbp4", "Igfbp5", "Igfbp6",
                                        "Igfbp7", "Il10", "Il13", "Il15", "Il18", "Il1a", "Il1b", "Il2",
                                        "Il32", "Il6", "Il6st", "Il7", "Inha", "Iqgap2", "Itga2", "Itpka",
                                        "Jun", "Kitlg", "Lcp1", "Mif", "Mmp1", "Mmp10", "Mmp12", "Mmp13",
                                        "Mmp14", "Mmp2", "Mmp3", "Mmp9", "Nap1l4", "Nrg1", "Pappa", "Pecam1", 
                                        "Pgf", "Pigf", "Plat", "Plau", "Plaur", "Ptbp1", "Ptger2", "Ptges", 
                                        "Rps6ka5", "Scamp4", "Selplg", "Sema3f", "Serpinb4", "Serpine1", "Serpine2", "Spp1",
                                        "Spx", "Timp2", "Tnf", "Tnfrsf10c", "Tnfrsf11b", "Tnfrsf1a", "Tnfrsf1b", "Tubgcp2",
                                        "Vegfa", "Vegfc", "Vgf", "Wnt16", "Wnt2"))

# Calculate SenMayo senescence score
ES.seurat <- escape.matrix(dataset.seurat.hspc, gene.sets = gene.sets, groups = 1000, min.size = 5)
dataset.seurat.hspc <- AddMetaData(dataset.seurat.hspc, ES.seurat, col.name = 'Enrichment_Score')
Age = dataset.seurat.hspc$ages
Enrichment_score = dataset.seurat.hspc$Enrichment_Score

######################################################
# SenMayo senescence score by age group
temp <- data.frame(Age, Enrichment_score)
temp$Age <- factor(temp$Age, levels = c("01m", "03m", "18m", "21m", "24m", "30m"),
                   labels = c("1mo", "3mo", "18mo", "21mo", "24mo", "30mo"))

p <- ggplot(temp[temp$Age %in% c("3mo", "18mo", "21mo", "24mo", "30mo"),], 
            aes(x = Age, y = Enrichment_score, fill = Age)) + 
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette="Greens") + 
  ylim(NA, max(temp$Enrichment_score[is.finite(temp$Enrichment_score)], na.rm = TRUE) * 1.2) + 
  theme_minimal(base_size = 10) +
  theme(axis.title = element_blank(), legend.position = "none",
        panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2),
        axis.line.x = element_line(color = "black", linewidth = 0.4),
        axis.line.y = element_line(color = "black", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black", linewidth = 0.3),
        axis.ticks.y = element_line(color = "black", linewidth = 0.3)) +
  geom_signif(comparisons = list(c("18mo", "21mo"), c("24mo", "30mo")),
              map_signif_level = TRUE, tip_length = 0.01,
              step_increase = 0.025, textsize = 3, 
              test = function(x, y) wilcox.test(x, y, alternative = "less"))

pdf(file = "./TMS_marrow/TMS_hspcs_Senmayo_violin.pdf", width = 3.25, height = 2.25)
print(p)
dev.off()

######################################################
# SenMayo senescence score vs. cell-level deviation
Entropy = dataset.seurat.hspc$entropy
temp <- data.frame(Age, Enrichment_score, Entropy)
temp <- temp %>%
  filter(Entropy >= quantile(temp$Entropy, 0.01, na.rm = TRUE), 
         Entropy <= quantile(temp$Entropy, 0.99, na.rm = TRUE))

p <- ggplot(temp, aes(x = Enrichment_score, y = log(Entropy), fill = Age)) + 
  geom_point() + 
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2) + 
  scale_fill_brewer(palette = "Blues") + 
  theme_minimal() + 
  ggtitle("TMS marrow HSPCs: Intrinsic noise vs SenMayo enrichment score")

pdf(file = "./TMS_marrow/TMS_hspcs_Senmayo_scatter.pdf")
print(p)
dev.off()

# Correlation between SenMayo score and cell-level deviation
cor.test(temp[is.finite(temp$Entropy),]$Enrichment_score, log(temp[is.finite(temp$Entropy),]$Entropy))

######################################################
# Proportion of low vs. high entropy cells (below/above the best fitted line)
temp_classified <- temp %>%
  mutate(logEntropy = log(Entropy),
         fitted = intercept + slope * Enrichment_score,
         resid_val = logEntropy - fitted,
         Entropy_group = if_else(logEntropy > fitted, "High entropy", "Low entropy"))

prop_data <- temp_classified %>%
  count(Entropy_group, Age) %>%
  group_by(Entropy_group) %>%
  mutate(prop = n / sum(n))

p_bar <- ggplot(prop_data, aes(x = Entropy_group, y = prop, fill = Age)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = NULL, y = "Proportion of cells", fill = "Age") +
  theme_minimal()

ggsave("./TMS_marrow/TMS_hspcs_Senmayo_scatter_region_bargraph.pdf", plot = p_bar, width = 3, height = 5)

######################################################
# 2. Mean Euclidean distance to cell type average
# 3. Mean Euclidean distance to tissue average using invariant genes
# 4. Scallop membership score
# For detailed workflow to compute these metrics, please refer to "decibel_scallop_analysis.ipynb"

# Load in data
load("TMS_marrow.RData")
TMS_marrow_metadata <- read.csv("adata_TMS_marrow_metadata.csv")
dataset.age.levels <- replace(dataset.age.levels, dataset.age.levels == "1m", "01m")
dataset.age.levels <- replace(dataset.age.levels, dataset.age.levels == "3m", "03m")
TMS_marrow_metadata$Age <- dataset.age.levels[dataset.age]
TMS_marrow_metadata$Age <- factor(TMS_marrow_metadata$Age, levels = c("01m", "03m", "18m", "21m", "24m", "30m"))
TMS_marrow_metadata$cell_type <- dataset.celltype.levels[dataset.celltype]

# Plot score distributions across cell types and age groups
p_marrow_euc_dist <- TMS_marrow_metadata %>%
  ggplot(aes(x = Age, y = log(euc_dist), fill = Age)) + 
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  facet_wrap(~ cell_type, ncol = 3) + 
  scale_fill_brewer(palette = "Greens") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  labs(x = "Age", y = "Euclidean distance (log-transformed)") +
  geom_signif(comparisons = list(c("01m", "03m"), c("03m", "18m"), c("18m", "21m"), c("21m", "24m"), c("24m", "30m")),
              test = function(x, y) wilcox.test(x, y, alternative = "less"), map_signif_level = TRUE, step_increase = 0.1)

p_marrow_euc_dist_tissue_invar <- TMS_marrow_metadata %>%
  ggplot(aes(x = Age, y = log(euc_dist_tissue_invar), fill = Age)) + 
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  facet_wrap(~ cell_type, ncol = 3) +
  scale_fill_brewer(palette = "Greens") +
  theme_minimal() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  labs(x = "Age", y = "Euclidean distance using invariant genes (log-transformed)") +
  geom_signif(comparisons = list(c("01m", "03m"), c("03m", "18m"), c("18m", "21m"), c("21m", "24m"), c("24m", "30m")),
              test = function(x, y) wilcox.test(x, y, alternative = "less"), map_signif_level = TRUE, step_increase = 0.1)

p_marrow_scallop_score <- TMS_marrow_metadata %>%
  ggplot(aes(x = Age, y = scallop_score, fill = Age)) + 
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  facet_wrap(~ cell_type, ncol = 3) +
  scale_fill_brewer(palette = "Greens") +
  theme_minimal() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), breaks = seq(0, 1, by = 0.2)) + 
  labs(x = "Age", y = "Scallop membership score") +
  geom_signif(comparisons = list(c("01m", "03m"), c("03m", "18m"), c("18m", "21m"), c("21m", "24m"), c("24m", "30m")),
              test = function(x, y) wilcox.test(x, y, alternative = "greater"), map_signif_level = TRUE, step_increase = 0.1)

ggsave("TMS_marrow_euc_dist.jpg", plot = p_marrow_euc_dist, width = 9, height = 11)
ggsave("TMS_marrow_euc_dist_tissue_invar.jpg", plot = p_marrow_euc_dist_tissue_invar, width = 9, height = 11)
ggsave("TMS_marrow_scallop.jpg", plot = p_marrow_scallop_score, width = 9, height = 11)

######################################################
# Plot score vs. cell-level deviation (optional)
TMS_marrow_cell_level <- read.csv("./TMS_marrow/Dispersion_final/TMS_marrow_cellular_dispersion_SAVER.csv")
TMS_marrow_cell_level_cleaned <- na.omit(distinct(TMS_marrow_cell_level))
TMS_marrow_cell_level_cleaned <- inner_join(TMS_marrow_cell_level_cleaned, 
                                            TMS_marrow_metadata[, c("X", "euc_dist", "euc_dist_tissue_invar", "scallop_score")], 
                                            by = c("Cell_barcode" = "X"))

p_scatter_cell_level_euc_dist <- ggplot(TMS_marrow_cell_level_cleaned, aes(x = log(euc_dist), y = log(Cell_level_deviation))) +
  geom_point(aes(color = Age), alpha = 0.7) +
  geom_smooth(method = "lm", linetype = "dashed", se = TRUE) +
  facet_wrap(~ cell_type, ncol = 6) + 
  labs(
    x = "Log-transformed Euclidean distance to cell-type average",
    y = "Log-transformed cell-level entropy",
    color = "Age"
  ) + theme_minimal(base_size = 10) +
  stat_cor(method = "pearson", size = 4, show.legend = FALSE)

p_scatter_cell_level_euc_dist_tissue_invar <- ggplot(TMS_marrow_cell_level_cleaned, aes(x = log(euc_dist_tissue_invar), y = log(Cell_level_deviation))) +
  geom_point(aes(color = Age), alpha = 0.7) +
  geom_smooth(method = "lm", linetype = "dashed", se = TRUE) +
  facet_wrap(~ cell_type, ncol = 6) + 
  labs(
    x = "Log-transformed Euclidean distance to tissue average using invariant genes",
    y = "Log-transformed cell-level entropy",
    color = "Age"
  ) + theme_minimal(base_size = 10) +
  stat_cor(method = "pearson", size = 4, show.legend = FALSE)

p_scatter_cell_level_scallop_score <- ggplot(TMS_marrow_cell_level_cleaned, aes(x = scallop_score, y = log(Cell_level_deviation))) +
  geom_point(aes(color = Age), alpha = 0.7) +
  geom_smooth(method = "lm", linetype = "dashed", se = TRUE) +
  facet_wrap(~ cell_type, ncol = 6) +
  labs(
    x = "Scallop membership score",
    y = "Log-transformed cell-level entropy",
    color = "Age"
  ) + theme_minimal(base_size = 10) +
  stat_cor(method = "pearson", size = 4, show.legend = FALSE)

######################################################