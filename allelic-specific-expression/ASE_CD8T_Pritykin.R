# Allele Specific Expression Analysis (Pritykin)

# Online links
# https://pubmed.ncbi.nlm.nih.gov/33862018/
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164978

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
# Allele-specific UMI counts (B6 = Aallele, SPRET/EiJ = Ballele) aggregated to
# genes, with the CD8 T cell type annotation. The uninfected F1 sample is used.
data_dir <- 'path/to/GSE164978_allele_matrices/uninfected_F1'
Aallele <- as.matrix(Matrix::readMM(list.files(data_dir, "^Aallele.*\\.mtx$", full.names = TRUE)[1]))
Ballele <- as.matrix(Matrix::readMM(list.files(data_dir, "^Ballele.*\\.mtx$", full.names = TRUE)[1]))
genes    <- read.delim(file.path(data_dir, "genes.tsv"), header = FALSE)
barcodes <- readLines(file.path(data_dir, "barcodes.tsv"))
rownames(Aallele) <- rownames(Ballele) <- make.unique(genes[[ncol(genes)]])
colnames(Aallele) <- colnames(Ballele) <- barcodes

cell_labels <- read.delim(file.path(data_dir, "pritykin_cell_labels.tsv"), stringsAsFactors = FALSE)
cell_labels <- cell_labels[cell_labels$sample == "uninfected_F1", ]
cell_labels$barcode <- sub("^.*\\|", "", cell_labels$barcode)   # strip the sample| prefix
cell_type <- setNames(cell_labels$cell_type, cell_labels$barcode)

######################################################
# Compute dyscoordination and allelic discordance per cell type
# Cell types with fewer than 150 cells are skipped.
MIN.READS = 5
FRAC.CELLS = 0.20

all_temp_gene_level <- data.frame()

for (ct in sort(unique(cell_type))) {
  cells <- intersect(names(cell_type)[cell_type == ct], colnames(Aallele))
  if (length(cells) < 150) next
  Aallele_ct <- Aallele[, cells, drop = FALSE]
  Ballele_ct <- Ballele[, cells, drop = FALSE]

  # Gene filter: covered in at least 20% of the cell type's cells
  Nallele <- Aallele_ct + Ballele_ct
  index_genes <- rowSums(Nallele >= MIN.READS) >= FRAC.CELLS * ncol(Nallele)
  Aallele_ct <- Aallele_ct[index_genes, ]
  Ballele_ct <- Ballele_ct[index_genes, ]

  ######################################################
  # Normalize the data to diploid counts
  celltot_diploid = 0.5 * (colSums(Aallele_ct + Ballele_ct))
  Aallele_ct = (10000) * (Aallele_ct / matrix(nrow = nrow(Aallele_ct), ncol = ncol(Aallele_ct), data = celltot_diploid, byrow = TRUE))
  Ballele_ct = (10000) * (Ballele_ct / matrix(nrow = nrow(Ballele_ct), ncol = ncol(Ballele_ct), data = celltot_diploid, byrow = TRUE))
  ABallele <- (0.5) * (Aallele_ct + Ballele_ct)

  ######################################################
  # Run SAVER on the allele-sum expression
  ABallele.saver <- saver(ABallele, ncores = 10, size.factor = 1)
  # save(ABallele.saver, file = paste0("ABallele_Pritykin_", gsub("[/ ]", "-", ct), "_SAVER.RData"))

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
    Y_A <- Aallele_ct[g, ]
    Y_B <- Ballele_ct[g, ]
    alpha_c <- Y_A + Y_B
    alpha_c[alpha_c == 0] <- 1
    p_hat <- sum(Y_A) / sum(Y_A + Y_B)
    allelic_discordance[g] <- mean((Y_A - alpha_c * p_hat)^2 / (alpha_c * p_hat * (1 - p_hat)), na.rm = TRUE)
  }

  # Log mean expression, used as the control in the comparison below
  log_mean <- log(rowMeans(mu.out))

  temp <- data.frame(Gene, cell_type = ct, Var_model, Gene_level_deviation,
                     allelic_discordance, log_mean)
  all_temp_gene_level <- rbind(all_temp_gene_level, temp)
}

all_temp_gene_level$Gene_level_deviation <- remove_outliers(all_temp_gene_level$Gene_level_deviation)

######################################################
# Gene-level dyscoordination vs. allelic discordance
# Expression-controlled partial Spearman: rank-regress both on log mean and
# correlate the residuals.
ok <- with(all_temp_gene_level, is.finite(Gene_level_deviation) & is.finite(allelic_discordance) & is.finite(log_mean))
df <- all_temp_gene_level[ok, ]
res_dev  <- resid(lm(rank(df$Gene_level_deviation) ~ rank(df$log_mean)))
res_disc <- resid(lm(rank(df$allelic_discordance) ~ rank(df$log_mean)))
partial_rho <- cor(res_dev, res_disc)

p <- ggplot(data.frame(x = res_dev, y = res_disc), aes(x = x, y = y)) +
  geom_point(shape = 21, size = 2, fill = "forestgreen", color = "black") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(x = "Partial residual of gene-level dyscoordination",
       y = "Partial residual of allelic discordance") +
  ggtitle(sprintf("Pritykin: partial rho = %.2f", partial_rho)) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        plot.background = element_blank(), axis.line = element_line(color = "black"))

ggsave("dyscoordination_vs_discordance_Pritykin.pdf", plot = p)

######################################################
