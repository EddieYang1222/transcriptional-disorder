# Allele Specific Expression Analysis (Ochiai)

# Online links
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7299619/
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132589

# Set up working directory
# setwd("S:/Penn Dropbox/Eddie Yang/Aging/Scripts/Code/allele_specific_analysis/ochiai")

######################################################
# Load in data
# 117,762 genes and 419 cells
counts_129 <- read.csv("/GSE132589_ASEcount_G1_129.txt", header = TRUE, row.names = 1, sep = "")
counts_129 <- as.matrix(counts_129)

counts_cast <- read.csv("/GSE132589_ASEcount_G1_CAST.txt", header = TRUE, row.names = 1, sep = "")
counts_cast <- as.matrix(counts_cast)

# Remove genes with zero counts
# 97,175 genes and 419 cells
index_genes <- rowSums(counts_129 + counts_cast) != 0
counts_129 <- counts_129[index_genes,]
counts_cast <- counts_cast[index_genes,]

######################################################
# Filter out genes with low expression
Aallele <- counts_129
Ballele <- counts_cast

MIN.TOT.COUNT = 4000
index_genes <- (rowSums(Aallele) > MIN.TOT.COUNT)
Aallele <- Aallele[index_genes,]
Ballele <- Ballele[index_genes,]
index_genes <- (rowSums(Ballele) > MIN.TOT.COUNT)
Aallele <- Aallele[index_genes,]
Ballele <- Ballele[index_genes,]

# Percent of genes remaining
nrow(Aallele)/nrow(counts_129) # 0.05971701

######################################################
# Normalize the data to diploid counts
celltot_diploid = 0.5*(colSums(Aallele + Ballele))
Aallele = (10000)*(Aallele/matrix(nrow=nrow(Aallele), ncol=ncol(Aallele), data=celltot_diploid, byrow=TRUE))
Ballele = (10000)*(Ballele/matrix(nrow=nrow(Ballele), ncol=ncol(Ballele), data=celltot_diploid, byrow=TRUE))

ABallele <- (0.5)*(Aallele + Ballele)

# Clean up intermediates
rm(MIN.TOT.COUNT)
rm(index_genes)
rm(celltot_diploid)

######################################################
# Run SAVER
library(SAVER)
library(pbapply)

ABallele.saver <- saver(ABallele, ncores = 10, size.factor = 1)
# save(ABallele.saver, file = "ABallele_Ochiai_SAVER.RData")
# load("ABallele_Ochiai_SAVER.RData")

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

# Find dispersion model with highest likelihood for each gene
size.factor <- colSums(ABallele) / mean(colSums(ABallele))
ABallele.saver.mu <- ABallele.saver$mu.out
ABallele.saver.var.models <- get_var_model(as.matrix(ABallele), ABallele.saver.mu, size.factor)
# save(ABallele.saver.var.models, file = "ABallele_variance_models_SAVER.RData")

######################################################
all_temp_gene_level <- data.frame()

# Initialize variables
Gene <- rownames(ABallele)
Var_model <- names(ABallele.saver.var.models)
Dispersion <- ABallele.saver.var.models
Dispersion_cCV <- ifelse(ABallele.saver$a == 0, NA_real_, 1 / ABallele.saver$a)
Dispersion_cFF <- ifelse(ABallele.saver$b == 0, NA_real_, 1 / ABallele.saver$b)
Dispersion_cVar <- ifelse(ABallele.saver$k == 0, NA_real_, ABallele.saver$k)
Var_avg <- vector()
Gene_level_deviation <- vector()

# Calculate gene-level deviation
estimate <- ABallele.saver$estimate
mu.out <- ABallele.saver$mu.out

for (g in 1:nrow(estimate)) {
  delta <- numeric(ncol(estimate))
  nu <- numeric(ncol(estimate))
  model <- names(ABallele.saver.var.models)[g]
  for (c in 1:ncol(estimate)) {
    if (model == "cCV") {  # cCV model
      nu_c <- mu.out[g, c]^2 / ABallele.saver.var.models[g]
    } else if (model == "cFF") {  # cFF model
      nu_c <- mu.out[g, c] / ABallele.saver.var.models[g]
    } else if (model == "cVar") {  # cVar model
      nu_c <- ABallele.saver.var.models[g]
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
temp <- data.frame(Gene, Dispersion_cCV, Dispersion_cFF, Dispersion_cVar, Var_model, Dispersion, Var_avg, Gene_level_deviation)
all_temp_gene_level <- rbind(all_temp_gene_level, temp)
all_temp_gene_level$Gene_level_deviation <- remove_outliers(all_temp_gene_level$Gene_level_deviation)

######################################################
all_temp_cell_level <- data.frame()

# Initialize variables
Cell_barcode <- colnames(estimate)
Var_avg <- vector()
Cell_level_deviation <- vector()

# Calculate cell-level deviation
estimate <- ABallele.saver$estimate
mu.out <- ABallele.saver$mu.out

for (c in 1:ncol(estimate)) {
  delta <- numeric(nrow(estimate))
  nu <- numeric(nrow(estimate))
  for (g in 1:nrow(estimate)) {
    model <- names(ABallele.saver.var.models)[g]
    if (model == "cCV") {  # cCV model
      nu_g <- mu.out[g, c]^2 / ABallele.saver.var.models[g]
    } else if (model == "cFF") {  # cFF model
      nu_g <- mu.out[g, c] / ABallele.saver.var.models[g]
    } else if (model == "cVar") {  # cVar model
      nu_g <- ABallele.saver.var.models[g]
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
temp <- data.frame(Cell_barcode, Cell_level_deviation)
all_temp_cell_level <- rbind(all_temp_cell_level, temp)

######################################################
# Allelic discordance

# Initialize variables
G <- nrow(ABallele)
C <- ncol(ABallele)
allelic_discordance <- numeric(G)

# Calculate allelic discordance for each cell
for (g in 1:G) {
  discordance_sum <- 0
  for (j in 1:C) {
    Y_A <- Aallele[g, j]
    Y_B <- Ballele[g, j]
    alpha_c <- (Y_A + Y_B)
    
    # Avoid division by zero
    if (alpha_c == 0) alpha_c <- 1
    
    numerator <- (Y_A / alpha_c - Y_B / alpha_c)^2
    denominator <- (Y_A / (2 * alpha_c) + Y_B / (2 * alpha_c))
    
    # Avoid division by zero in denominator
    if (denominator == 0) denominator <- 1
    
    discordance_sum <- discordance_sum + 2 * (numerator / denominator)
  }
  allelic_discordance[g] <- discordance_sum / C
}

######################################################
# a = constant CV (quadratic)
# b = constant FF (linear)
# k = constant variance (independent)

# Constant CV (quadratic)
Aallele.disp <- matrix(nrow=nrow(Aallele), ncol=ncol(Aallele))
for (i in 1:nrow(Aallele)) {
  Aallele.disp[i, ] <- ((Aallele[i, ] - Ballele[i, ])^2) / (2*((ABallele.saver$mu.out[i, ])^2))
}
Aallele.disp.final <- (1/ncol(Aallele.disp))*(rowSums(Aallele.disp))
Aallele.disp.final.cv <- Aallele.disp.final

# Constant FF (linear)
Aallele.disp <- matrix(nrow=nrow(Aallele), ncol=ncol(Aallele))
for (i in 1:nrow(Aallele)) {
  Aallele.disp[i, ] <- ((Aallele[i, ] - Ballele[i, ])^2)/ (2*(ABallele.saver$mu.out[i, ]))
}
Aallele.disp.final <- (1/ncol(Aallele.disp))*(rowSums(Aallele.disp))
Aallele.disp.final.ff <- Aallele.disp.final

# Constant variance (independent)
Aallele.disp <- matrix(nrow=nrow(Aallele), ncol=ncol(Aallele))
for (i in 1:nrow(Aallele)) {
  Aallele.disp[i, ] <- (Aallele[i, ] - Ballele[i, ])^2 / (2)
}
Aallele.disp.final <- (1/ncol(Aallele.disp))*(rowSums(Aallele.disp))
Aallele.disp.final.var <- Aallele.disp.final

######################################################

# Gene-level deviation vs. allelic discordance
all_temp_gene_level$allelic_discordance <- allelic_discordance

p <- ggplot(all_temp_gene_level, aes(x = log(Gene_level_deviation), y = log(allelic_discordance))) + 
  geom_point(shape = 21, size = 2, fill = "forestgreen", color = "black") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(x = "Log-transformed gene-level entropy", y = "Log-transformed allelic discordance") + 
  ggtitle("Gene-level entropy vs allelic discordance") +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        plot.background = element_blank(), axis.line = element_line(color = "black")) + 
  stat_cor() +
  ylim(-8, 2)

ggsave("entropy_vs_discordance.pdf", plot = p)

# Gene-level dispersion (under 3 models) vs. allelic discordance
plot_df <- data.frame(CV = log(Aallele.disp.final.cv), FF = log(Aallele.disp.final.ff),
                      Var = log(Aallele.disp.final.var), allelic_discordance = log(allelic_discordance))

p <- ggplot(plot_df, aes(x = CV, y = allelic_discordance)) + 
  geom_point(shape = 21, size = 2, fill = "forestgreen", color = "black") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(x = "Log-transformed dispersion (cCV)", y = "Log-transformed allelic discordance") + 
  ggtitle("Gene-level dispersion (cCV) vs allelic discordance") +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        plot.background = element_blank(), axis.line = element_line(color = "black"))

ggsave("disp_vs_discordance_cCV.pdf", plot = p)

p <- ggplot(plot_df, aes(x = FF, y = allelic_discordance)) + geom_point(shape = 21, size = 2, fill = "forestgreen", color = "black") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(x = "Log-transformed dispersion (cFF)", y = "Log-transformed allelic discordance") + 
  ggtitle("Gene-level dispersion (cFF) vs allelic discordance") +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        plot.background = element_blank(), axis.line = element_line(color = "black"))

ggsave("disp_vs_discordance_cFF.pdf", plot = p)

p <- ggplot(plot_df, aes(x = Var, y = allelic_discordance)) + geom_point(shape = 21, size = 2, fill = "forestgreen", color = "black") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(x = "Log-transformed dispersion (cVar)", y = "Log-transformed allelic discordance") + 
  ggtitle("Gene-level dispersion (cVar) vs allelic discordance") +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        plot.background = element_blank(), axis.line = element_line(color = "black"))

ggsave("disp_vs_discordance_cVar.pdf", plot = p)

######################################################