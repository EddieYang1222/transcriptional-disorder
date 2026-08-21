# Pathway Enrichment of Transcriptional Dyscoordination (Human Kidney)

# Online links
# https://www.gsea-msigdb.org/gsea/msigdb

# Set up working directory
# setwd("path/to/working/directory")

library(dplyr)
library(fgsea)
library(msigdbr)
library(data.table)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

set.seed(42)

######################################################
# MSigDB C2:CP* + C5:GO* gene sets (robust to msigdbr column-name changes)
msig <- as.data.frame(msigdbr(species = "Homo sapiens"))
cn <- colnames(msig)
catcol <- if ("gs_cat" %in% cn) "gs_cat" else "gs_collection"
subcol <- if ("gs_subcat" %in% cn) "gs_subcat" else "gs_subcollection"
symcol <- if ("gene_symbol" %in% cn) "gene_symbol" else "gene_symbol"
keep <- (msig[[catcol]] == "C2" & grepl("^CP", msig[[subcol]])) |
        (msig[[catcol]] == "C5" & grepl("^GO", msig[[subcol]]))
msig <- msig[keep, ]
cat("gene-set rows:", nrow(msig), " unique sets:", length(unique(msig$gs_name)), "\n")
pathways <- split(msig[[symcol]], msig$gs_name)
src_map <- msig %>% distinct(gs_name, .keep_all = TRUE) %>% select(gs_name, !!subcol) %>%
  mutate(Source = sub("^CP:", "", sub("^GO:", "GO_", .data[[subcol]])))

######################################################
# H1: genes gaining gene-level dyscoordination with age
# Rank each gene by Spearman(ordinal age, Gene_level_deviation) pooled over the 3
# regrouped types PT(=PST+PCT)/DCT(=DCT1+DCT2_PC)/TAL; fgsea over the ranking.
age_ord <- c("30-39" = 1, "40-49" = 2, "50-59" = 3, "60+" = 4)

g <- read.csv("human_kidney_aging_PT_DCT_TAL_regrouped_estimated_dispersion_SAVER_geneFilt0p1.csv")
g$Gene_level_deviation <- remove_outliers(g$Gene_level_deviation)
g <- g %>% mutate(Age = ifelse(as.character(Age) == "20-39", "30-39", as.character(Age)),
                  aord = age_ord[Age]) %>%
  filter(!is.na(Gene_level_deviation), !is.na(aord))

# per-gene Spearman(age, dev) pooling the 3 cell types (up to 12 points/gene)
rho_tab <- g %>% group_by(Gene) %>%
  summarise(rho = if (length(unique(aord)) >= 3 && sd(Gene_level_deviation) > 0)
                    suppressWarnings(cor(aord, Gene_level_deviation, method = "spearman")) else NA_real_,
            npts = dplyr::n(), .groups = "drop") %>%
  filter(!is.na(rho))
cat("genes ranked:", nrow(rho_tab), "\n")
write.csv(rho_tab %>% arrange(desc(rho)) %>% select(Gene, rho),
          "hk_H1_regrouped_gene_rho.csv", row.names = FALSE)

ranks <- setNames(rho_tab$rho, rho_tab$Gene)
ranks <- ranks[!duplicated(names(ranks))]
ranks <- sort(ranks, decreasing = TRUE)

resH1 <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500, eps = 0)
resH1 <- resH1 %>% as.data.frame() %>%
  left_join(src_map %>% select(gs_name, Source), by = c("pathway" = "gs_name")) %>%
  arrange(padj)
resH1$leadingEdge <- vapply(resH1$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
write.csv(resH1, "gseaH1_regrouped_PT_DCT_TAL.csv", row.names = FALSE)
cat("H1 significant (padj<0.05, NES>0):", sum(resH1$NES > 0 & resH1$padj < 0.05, na.rm = TRUE), "\n")

######################################################
# H2: genes activated in high-dyscoordination cells
# Input hk_H2_regrouped_gene_r.csv (Gene, r, n_nz) is the per-gene Pearson
# correlation between gene expression and cell-level dyscoordination, computed
# upstream. Keep only genes FDR<0.05 (rat-consistent), fgsea over their r.
r2 <- read.csv("hk_H2_regrouped_gene_r.csv")   # Gene, r, n_nz
n <- r2$n_nz
rr <- r2$r
t <- rr * sqrt((n - 2) / (1 - rr^2))
r2$p <- 2 * pt(-abs(t), df = n - 2)
r2$fdr <- p.adjust(r2$p, "fdr")
sig <- r2 %>% filter(fdr < 0.05)               # enrichment input
cat("H2 FDR<0.05 genes (enrichment input):", nrow(sig), "\n")

ranks2 <- setNames(sig$r, sig$Gene)
ranks2 <- ranks2[!duplicated(names(ranks2))]
ranks2 <- sort(ranks2, decreasing = TRUE)
resH2 <- fgsea(pathways = pathways, stats = ranks2, minSize = 5, maxSize = 500, eps = 0) %>%
  as.data.frame() %>% arrange(pval)
resH2$leadingEdge <- vapply(resH2$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
resH2 <- resH2[!is.na(resH2$pval), ]
write.csv(resH2, "gseaH2_regrouped_FDRsub_PT_DCT_TAL.csv", row.names = FALSE)
cat("H2 pathways pval<0.05:", sum(resH2$pval < 0.05, na.rm = TRUE), "\n")

######################################################
