# Transcriptional dyscoordination analysis for ionizing radiated ESCC cells

# Online links
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6831193/
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81812

# Set up working directory
# setwd("S:/Penn Dropbox/Eddie Yang/Aging/Scripts/Code/radiation/ESCC_cell_line")

library(SAVER)
library(pbapply)
library(ggplot2)
library(ggsignif)

######################################################
# Load in data
dataset <- read.csv("./GSE81812_Normalized_counts.txt", header = TRUE, row.names = 1)
dataset <- dataset[, 0:133]
dataset.rad0 <- dataset[,0:41]  # 0 Gy
dataset.rad1 <- dataset[,42:133] # 12 Gy

# Remove cells and genes with zero counts
message("The initial matrix size is ", nrow(dataset), " genes and ", ncol(dataset), " cells.")
dataset <- dataset[rowSums(dataset) != 0, colSums(dataset) != 0] 
message("After removing cells and genes with zero counts, the new matrix size is ",
        nrow(dataset), " genes and ", ncol(dataset), " cells.")

######################################################
# Run SAVER
dataset.saver <- saver(dataset, ncores = 4)
# save(dataset.saver, file="ESCC_radiation_manifold.RData")

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

# Find dispersion model with highest likelihood for each gene
size.factor <- colSums(dataset) / mean(colSums(dataset))
dataset.saver.mu <- dataset.saver$mu.out
dataset.saver.var.models <- get_var_model(as.matrix(dataset), dataset.saver.mu, size.factor)
# save(dataset.saver.var.models, file="ESCC_radiation_variance_models_SAVER.RData")

# Load pre-computed data
# load("ESCC_radiation_manifold.RData")
# load("ESCC_radiation_variance_models_SAVER.RData")

######################################################
# Normalize each cell group
dataset.rad0 <- dataset[,0:41] / dataset.saver$info$size.factor[0:41] # 0 Gy
dataset.rad1 <- dataset[,42:133] / dataset.saver$info$size.factor[42:133] # 12 Gy

# Subset the manifold
dataset.rad0.mu <- dataset.saver$mu.out[,0:41] # 0 Gy
dataset.rad1.mu <- dataset.saver$mu.out[,42:133] # 12 Gy

# Calculate dispersion
dataset.saver.rad0 <- saver(dataset.rad0, mu = dataset.rad0.mu)
dataset.saver.rad1 <- saver(dataset.rad1, mu = dataset.rad1.mu)

# Prepare inputs
dataset.age.levels <- c("0 Gy", "12 Gy")
dataset.saver.rad.a <- list(dataset.saver.rad0$a, dataset.saver.rad1$a)
dataset.saver.rad.b <- list(dataset.saver.rad0$b, dataset.saver.rad1$b)
dataset.saver.rad.k <- list(dataset.saver.rad0$k, dataset.saver.rad1$k)
dataset.saver.rad.est <- list(dataset.saver.rad0$estimate, dataset.saver.rad1$estimate)
dataset.saver.rad.mu <- list(dataset.saver.rad0$mu.out, dataset.saver.rad1$mu.out)

######################################################
# Gene-level deviation
all_temp_gene_level <- data.frame()

for (j in c(1:2)) {
  # Initialize variables
  Gene <- rownames(dataset)
  Age <- rep(dataset.age.levels[j], nrow(dataset))
  Var_model <- names(dataset.saver.var.models)
  Dispersion <- dataset.saver.var.models
  Dispersion_cCV <- ifelse(dataset.saver.rad.a[[j]] == 0, NA_real_, 1 / dataset.saver.rad.a[[j]])
  Dispersion_cFF <- ifelse(dataset.saver.rad.b[[j]] == 0, NA_real_, 1 / dataset.saver.rad.b[[j]])
  Dispersion_cVar <- ifelse(dataset.saver.rad.k[[j]] == 0, NA_real_, dataset.saver.rad.k[[j]])
  Var_avg <- vector()
  Gene_level_deviation <- vector()
  
  # Calculate gene-level deviation
  estimate <- dataset.saver.rad.est[[j]]
  mu.out <- dataset.saver.rad.mu[[j]]
  
  for (g in 1:nrow(estimate)) {
    delta <- numeric(ncol(estimate))
    nu <- numeric(ncol(estimate))
    model <- names(dataset.saver.var.models)[g]
    for (c in 1:ncol(estimate)) {
      if (model == "cCV") {  # cCV model
        nu_c <- mu.out[g, c]^2 / dataset.saver.var.models[g]
      } else if (model == "cFF") {  # cFF model
        nu_c <- mu.out[g, c] / dataset.saver.var.models[g]
      } else if (model == "cVar") {  # cVar model
        nu_c <- dataset.saver.var.models[g]
      } else {
        nu_c = NA  
      }
      # Compute delta for gene g and cell c
      if (!is.na(nu_c)) {
        nu[c] <- nu_c
        delta[c] <- (estimate[g, c] - mu.out[g, c])^2 / nu_c
      } else {
        nu[c] <- NA
        delta[c] <- NA
      }
    }
    # Average delta across all cells to get gene-level deviation
    Var_avg <- append(Var_avg, mean(nu, na.rm = TRUE))
    Gene_level_deviation <- append(Gene_level_deviation, mean(delta, na.rm = TRUE))
  }
  # Combine data from all cell types
  temp <- data.frame(Gene, Age, Dispersion_cCV, Dispersion_cFF, Dispersion_cVar, Var_model, Dispersion, Var_avg, Gene_level_deviation)
  all_temp_gene_level <- rbind(all_temp_gene_level, temp)
}

all_temp_gene_level$Gene_level_deviation <- remove_outliers(all_temp_gene_level$Gene_level_deviation)

######################################################
# Radiation dosage vs. gene-level deviation
p <- ggplot(all_temp_gene_level[!is.nan(all_temp_gene_level$Gene_level_deviation) & !is.na(all_temp_gene_level$Gene_level_deviation) & 
                                  all_temp_gene_level$Gene_level_deviation <= quantile(all_temp_gene_level$Gene_level_deviation, 0.99, na.rm = TRUE) &
                                  all_temp_gene_level$Gene_level_deviation >= quantile(all_temp_gene_level$Gene_level_deviation, 0.01, na.rm = TRUE),], 
            aes(x = Age, y = log(Gene_level_deviation), fill = Age)) + 
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  labs(x = "Radiation dosage", y = "Gene-level entropy (log-transformed)") + 
  scale_fill_brewer(palette="Blues") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2),
        axis.line.x = element_line(color = "black", linewidth = 0.4),
        axis.line.y = element_line(color = "black", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black", linewidth = 0.3),
        axis.ticks.y = element_line(color = "black", linewidth = 0.3)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_signif(comparisons = list(c("0 Gy", "12 Gy")), map_signif_level = TRUE, step_increase = 0.1, 
              test = function(x, y) wilcox.test(x, y, alternative = "less"))

ggsave("disp_by_dosage_gene_level.pdf", plot = p)

######################################################
# Cell-level deviation
all_temp_cell_level <- data.frame()

for (j in c(1:2)) {
  # Initialize variables
  Cell_barcode <- colnames(dataset.saver.rad.est[[j]])
  Age <- rep(dataset.age.levels[j], length(Cell_barcode))
  Var_avg <- vector()
  Cell_level_deviation <- vector()
  
  # Calculate cell-level deviation
  estimate <- dataset.saver.rad.est[[j]]
  mu.out <- dataset.saver.rad.mu[[j]]
  
  for (c in 1:ncol(estimate)) {
    delta <- numeric(nrow(estimate))
    nu <- numeric(nrow(estimate))
    for (g in 1:nrow(estimate)) {
      model <- names(dataset.saver.var.models)[g]
      if (model == "cCV") {  # cCV model
        nu_g <- mu.out[g, c]^2 / dataset.saver.var.models[g]
      } else if (model == "cFF") {  # cFF model
        nu_g <- mu.out[g, c] / dataset.saver.var.models[g]
      } else if (model == "cVar") {  # cVar model
        nu_g <- dataset.saver.var.models[g]
      } else {
        nu_g <- NA  
      }
      # Compute delta for gene g and cell c
      if (!is.na(nu_g)) {
        nu[g] <- nu_g
        delta[g] <- (estimate[g, c] - mu.out[g, c])^2 / nu_g
      } else {
        nu[g] <- NA
        delta[g] <- NA
      }
    }
    # Average delta across all cells to get gene-level deviation
    Var_avg <- append(Var_avg, mean(nu, na.rm=TRUE))
    Cell_level_deviation <- append(Cell_level_deviation, mean(delta, na.rm = TRUE))
  }
  # Combine data from all cell types
  temp <- data.frame(Cell_barcode, Age, Var_avg, Cell_level_deviation)
  all_temp_cell_level <- rbind(all_temp_cell_level, temp)
}

######################################################
# Radiation dosage vs. cell-level deviation
p <- ggplot(all_temp_cell_level[all_temp_cell_level$Cell_level_deviation <= quantile(all_temp_cell_level$Cell_level_deviation, 0.99, na.rm = TRUE) & 
                                  all_temp_cell_level$Cell_level_deviation >= quantile(all_temp_cell_level$Cell_level_deviation, 0.01, na.rm = TRUE),], 
            aes(x = Age, y = log(Cell_level_deviation), fill = Age)) + 
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  labs(x = "Radiation dosage", y = "Cell-level entropy (log-transformed)") + 
  scale_fill_brewer(palette="Reds") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2),
        axis.line.x = element_line(color = "black", linewidth = 0.4),
        axis.line.y = element_line(color = "black", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black", linewidth = 0.3),
        axis.ticks.y = element_line(color = "black", linewidth = 0.3)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_signif(comparisons = list(c("0 Gy", "12 Gy")), map_signif_level = TRUE, step_increase = 0.1, 
              test = function(x, y) wilcox.test(x, y, alternative = "less"))

ggsave("disp_by_dosage_cell_level.pdf", plot = p)

######################################################
# Radiation dosage vs. gene-level dispersion (under 3 models)
temp <- data.frame(dosage = c(rep('0 Gy',nrow(dataset.rad0)), rep('12 Gy',nrow(dataset.rad1))),
                   value_cCV = c(log(1/dataset.saver.rad0$a), log(1/dataset.saver.rad1$a)),
                   value_cFF = c(log(1/dataset.saver.rad0$b), log(1/dataset.saver.rad1$b)),
                   value_cVar = c(log(dataset.saver.rad0$k), log(dataset.saver.rad1$k)))

p <- ggplot(temp, aes(x = dosage, y = value_cCV, fill = dosage)) + 
  geom_violin(trim=FALSE) + stat_summary(fun= mean, geom= "point", shape= 23, size = 2, fill = "white") +
  labs(x = "Radiation dosage", y = "Log-transformed dispersion (cCV)") + 
  ggtitle("Gene-level dispersion (cCV) vs radiation dosage") + 
  scale_fill_brewer(palette="Blues") +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black"))

ggsave("disp_by_dosage_cCV.pdf", plot = p)

p <- ggplot(temp, aes(x = dosage, y = value_cFF, fill = dosage)) + 
  geom_violin(trim=FALSE) + stat_summary(fun= mean, geom= "point", shape= 23, size = 2, fill = "white") +
  labs(x = "Radiation dosage", y = "Log-transformed dispersion (cFF)") + 
  ggtitle("Gene-level dispersion (cFF) vs radiation dosage") + 
  scale_fill_brewer(palette="Blues") +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black"))

ggsave("disp_by_dosage_cFF.pdf", plot = p)

p <- ggplot(temp, aes(x = dosage, y = value_cVar, fill = dosage)) + 
  geom_violin(trim=FALSE) + stat_summary(fun= mean, geom= "point", shape= 23, size = 2, fill = "white") +
  labs(x = "Radiation dosage", y = "Log-transformed dispersion (cVar)") + 
  ggtitle("Gene-level dispersion (cVar) vs radiation dosage") + 
  scale_fill_brewer(palette="Blues") +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(color = "black"))

ggsave("disp_by_dosage_cVar.pdf", plot = p)

######################################################