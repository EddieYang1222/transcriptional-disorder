# Allele Specific Expression Analysis (Weber)

# Online links
# https://www.biorxiv.org/content/10.1101/2026.04.02.716195
# https://data.igvf.org/

# Set up working directory
# setwd("path/to/working/directory")

library(Matrix)
library(SAVER)
library(pbapply)
library(ggplot2)
library(ggpubr)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Load in data
# Liver snRNA-seq of B6 crossed with seven strains. Each cross is analyzed
# separately (its allele-specific cis-bias is cross-specific), so allelic
# discordance is centered within cross. Analysis is restricted to hepatocytes.
data_dir <- 'path/to/IGVF_liver_allele_matrices'
cell_labels <- read.delim(file.path(data_dir, "liver_cell_labels.tsv"), stringsAsFactors = FALSE)
cell_labels$barcode <- sub("^.*\\|", "", cell_labels$barcode)
hepatocyte_types <- c("Hepatocyte", "Hepatocyte-periportal")
crosses <- setdiff(unique(cell_labels$cross), c("unassigned", "cross", NA, ""))

######################################################
# Compute dyscoordination and allelic discordance per cross
# Crosses with fewer than 50 hepatocytes are skipped.
MIN.READS = 5
FRAC.CELLS = 0.20

all_temp_gene_level <- data.frame()

for (cr in crosses) {
  cross_dir <- file.path(data_dir, "allele_matrices", cr)
  Aallele <- as.matrix(Matrix::readMM(list.files(cross_dir, "^Aallele.*\\.mtx$", full.names = TRUE)[1]))
  Ballele <- as.matrix(Matrix::readMM(list.files(cross_dir, "^Ballele.*\\.mtx$", full.names = TRUE)[1]))
  genes    <- read.delim(file.path(cross_dir, "genes.tsv"), header = FALSE)
  barcodes <- readLines(file.path(cross_dir, "barcodes.tsv"))
  rownames(Aallele) <- rownames(Ballele) <- make.unique(genes[[ncol(genes)]])
  colnames(Aallele) <- colnames(Ballele) <- barcodes

  cells <- intersect(cell_labels$barcode[cell_labels$cross == cr & cell_labels$cell_type %in% hepatocyte_types], barcodes)
  if (length(cells) < 50) next
  Aallele <- Aallele[, cells, drop = FALSE]
  Ballele <- Ballele[, cells, drop = FALSE]

  # Gene filter: covered in at least 20% of the cross's hepatocytes
  Nallele <- Aallele + Ballele
  index_genes <- rowSums(Nallele >= MIN.READS) >= FRAC.CELLS * ncol(Nallele)
  Aallele <- Aallele[index_genes, ]
  Ballele <- Ballele[index_genes, ]

  ######################################################
  # Normalize the data to diploid counts
  celltot_diploid = 0.5 * (colSums(Aallele + Ballele))
  Aallele = (10000) * (Aallele / matrix(nrow = nrow(Aallele), ncol = ncol(Aallele), data = celltot_diploid, byrow = TRUE))
  Ballele = (10000) * (Ballele / matrix(nrow = nrow(Ballele), ncol = ncol(Ballele), data = celltot_diploid, byrow = TRUE))
  ABallele <- (0.5) * (Aallele + Ballele)

  ######################################################
  # Run SAVER on the allele-sum expression
  ABallele.saver <- saver(ABallele, ncores = 10, size.factor = 1)
  # save(ABallele.saver, file = paste0("ABallele_Weber_", cr, "_SAVER.RData"))

  # Find dispersion model with highest likelihood for each gene
  size.factor <- colSums(ABallele) / mean(colSums(ABallele))
  ABallele.saver.mu <- ABallele.saver$mu.out
  ABallele.saver.var.models <- get_var_model(as.matrix(ABallele), ABallele.saver.mu, size.factor)

  ######################################################
  # Calculate gene-level transcriptional dyscoordination
  Gene <- rownames(ABallele)
  Var_model <- names(ABallele.saver.var.models)
  Gene_level_deviation <- vector()
  estimate <- ABallele.saver$estimate
  mu.out <- ABallele.saver$mu.out

  for (g in 1:nrow(estimate)) {
    delta <- numeric(ncol(estimate))
    model <- names(ABallele.saver.var.models)[g]
    for (c in 1:ncol(estimate)) {
      if (model == "cCV") {  # cCV model
        nu_c <- mu.out[g, c]^2 / ABallele.saver.var.models[g]
      } else if (model == "cFF") {  # cFF model
        nu_c <- mu.out[g, c] / ABallele.saver.var.models[g]
      } else if (model == "cVar") {  # cVar model
        nu_c <- ABallele.saver.var.models[g]
      } else {
        nu_c <- NA
      }
      delta[c] <- if (!is.na(nu_c)) (estimate[g, c] - mu.out[g, c])^2 / nu_c else NA
    }
    Gene_level_deviation <- append(Gene_level_deviation, mean(delta, na.rm = TRUE))
  }

  ######################################################
  # Allelic discordance: allele fraction centered on the gene's baseline ratio
  # p_hat and standardized by the binomial variance (expectation 1 under sampling)
  G <- nrow(ABallele)
  allelic_discordance <- numeric(G)
  for (g in 1:G) {
    Y_A <- Aallele[g, ]
    Y_B <- Ballele[g, ]
    alpha_c <- Y_A + Y_B
    alpha_c[alpha_c == 0] <- 1
    p_hat <- sum(Y_A) / sum(Y_A + Y_B)
    allelic_discordance[g] <- mean((Y_A - alpha_c * p_hat)^2 / (alpha_c * p_hat * (1 - p_hat)), na.rm = TRUE)
  }

  # Log mean expression, used as the control in the comparison below
  log_mean <- log(rowMeans(mu.out))

  temp <- data.frame(Gene, cross = cr, Var_model, Gene_level_deviation,
                     allelic_discordance, log_mean)
  all_temp_gene_level <- rbind(all_temp_gene_level, temp)
}

all_temp_gene_level$Gene_level_deviation <- remove_outliers(all_temp_gene_level$Gene_level_deviation)

######################################################
# Gene-level dyscoordination vs. allelic discordance, per cross
# Expression-controlled partial Spearman: rank-regress both on log mean and
# correlate the residuals within each cross.
partial_rho_by_cross <- sapply(split(all_temp_gene_level, all_temp_gene_level$cross), function(df) {
  ok <- with(df, is.finite(Gene_level_deviation) & is.finite(allelic_discordance) & is.finite(log_mean))
  df <- df[ok, ]
  if (nrow(df) < 10) return(NA_real_)
  res_dev  <- resid(lm(rank(df$Gene_level_deviation) ~ rank(df$log_mean)))
  res_disc <- resid(lm(rank(df$allelic_discordance) ~ rank(df$log_mean)))
  cor(res_dev, res_disc)
})

write.csv(data.frame(cross = names(partial_rho_by_cross), partial_rho = partial_rho_by_cross),
          "dyscoordination_vs_discordance_Weber_by_cross.csv", row.names = FALSE)

######################################################
