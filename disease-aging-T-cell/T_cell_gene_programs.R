# Cell-level dyscoordination vs T-cell gene-program scores (Figure 4C)

# Set up working directory
# setwd("path/to/working/directory")

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')
source('T_cell_common.R')

DATASETS <- c("glioma", "melanoma", "Terekhova", "Wang_aging")
MIN_CELLS_DONOR <- 100     # donors need this many cells to contribute a per-donor rho

######################################################
# Gene programs. Score per cell: Seurat AddModuleScore (ctrl=100 bins) for the Seurat
# objects; manual mean-of-per-gene-z for the glioma log-normalized matrix.
PROGRAMS <- list(
  GZMK_inflammaging = c("GZMK","EOMES","CXCR3","CCL5","CST7","GZMA"),
  cytotoxicity      = c("GZMB","PRF1","NKG7","GNLY","KLRG1"),
  exhaustion        = c("PDCD1","TOX","LAG3","HAVCR2","TIGIT","ENTPD1"),
  progenitor_memory = c("TCF7","IL7R","CCR7","LEF1","SELL"),
  proliferation     = c("MKI67","TOP2A","STMN1"),
  interferon        = c("ISG15","IFIT1","IFIT3","MX1","OAS1","STAT1"),
  nfkb_inflammatory = c("NFKB1","NFKBIA","REL","TNF","IL1B","CXCL8"),
  mito_oxidative    = c("SOD1","TXN","NQO1","NDUFA1","GPX1","MT-CO1","MT-ND1"),
  generic_activation= c("CD69","FOS","JUN","JUNB","NR4A1","EGR1","DUSP1"))

score_seurat <- function(obj) {
  genes <- rownames(obj)
  feats <- lapply(PROGRAMS, function(g) intersect(g, genes))
  feats <- feats[sapply(feats, length) >= 2]
  obj <- AddModuleScore(obj, features = feats, name = "PROG_", nbin = 24, ctrl = 100, seed = 42)
  sc <- obj@meta.data[, paste0("PROG_", seq_along(feats)), drop = FALSE]
  colnames(sc) <- names(feats)
  data.frame(Cell_barcode = rownames(obj@meta.data),
             nCount_RNA = obj@meta.data$nCount_RNA, nFeature_RNA = obj@meta.data$nFeature_RNA,
             sc, check.names = FALSE)
}
score_matrix <- function(mat) {  # genes x cells (log-normalized)
manual mean of per-gene z
  genes <- rownames(mat)
  out <- sapply(PROGRAMS, function(g) {
    gg <- intersect(g, genes)
    if (length(gg) < 2) return(rep(NA_real_, ncol(mat)))
    z <- t(scale(t(as.matrix(mat[gg, , drop = FALSE]))))   # z-score each gene across cells
    colMeans(z, na.rm = TRUE)
  })
  data.frame(Cell_barcode = colnames(mat), nCount_RNA = Matrix::colSums(mat),
             nFeature_RNA = Matrix::colSums(mat > 0), out, check.names = FALSE)
}

######################################################
# Score once per dataset and cache to module_scores_<ds>.csv.gz.
# Expression objects (Seurat rds / GBM log-normalized matrix) come from data_dir.
data_dir <- 'path/to/T_cell_expression'
SEURAT_RDS <- c(melanoma = "melanoma_T_cell_processed_seurat.rds",
                Terekhova = "Terekhova_T_cell_processed_seurat.rds",
                Wang_aging = "Wang_aging_T_cell_processed_seurat.rds")

for (ds in DATASETS) {
  out_f <- sprintf("module_scores_%s.csv.gz", ds)
  if (file.exists(out_f)) { cat(ds, "scores exist, skip\n")
  next }
  cat("scoring", ds, "...\n")
  if (ds == "glioma") {
    mat <- readRDS(file.path(data_dir, "CD8_exprs_norm.rds"))
    sc  <- score_matrix(mat)
    rm(mat)
    gc()
  } else {
    obj <- readRDS(file.path(data_dir, SEURAT_RDS[[ds]]))
    sc  <- score_seurat(obj)
    rm(obj)
    gc()
  }
  write_csv(sc, out_f)
  cat("  wrote", out_f, "\n")
}

######################################################
# Join scores to dyscoordination + donor, and add the two dyscoordination modes.
aux_dir <- 'path/to/T_cell_aux'
build_master <- function(ds) {
  sc <- read_csv(sprintf("module_scores_%s.csv.gz", ds), show_col_types = FALSE)
  if (ds == "glioma") {
    ent <- read_dyscoordination("glioma_T_cell_cellular_dispersion_SAVER.csv") %>%
      select(Cell_barcode, cell_type, dyscoordination, log_dyscoordination)
    obs <- read.csv(file.path(aux_dir, "GBM_CD8_obs.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    names(obs)[1] <- "Cell_barcode"
    ent <- ent %>% inner_join(obs %>% transmute(Cell_barcode, donor = Patient_x), by = "Cell_barcode")
    # prefer obs depth for glioma (matches the covariate table)
    m <- ent %>% inner_join(sc %>% select(-nCount_RNA, -nFeature_RNA) %>%
                              left_join(obs %>% transmute(Cell_barcode, nCount_RNA = nCount_RNA,
                                                          nFeature_RNA = nFeature_RNA), by = "Cell_barcode"),
                            by = "Cell_barcode")
  } else {
    csv       <- sprintf("%s_T_cell_cellular_dispersion_SAVER.csv", ds)
    donor_col <- c(melanoma = "patient_alias", Terekhova = "Donor_id", Wang_aging = "sampleName")[[ds]]
    ent <- read_dyscoordination(csv) %>% rename(donor = !!donor_col) %>%
      select(Cell_barcode, cell_type, donor, dyscoordination, log_dyscoordination)
    m <- ent %>% inner_join(sc, by = "Cell_barcode")
  }
  m %>% filter(nCount_RNA > 0)
}

######################################################
# Donor-level inference: per-donor Spearman(dysc_used, program score) -> median across
# donors, number positive, signed-rank p. Pooled rho reported for reference. Both modes.
rows <- list()
donor_rows <- list()
for (ds in DATASETS) {
  m0    <- build_master(ds)
  progs <- intersect(PROG, names(m0))
  for (mode in MODES) {
    m <- add_dysc_used(m0, mode) %>% filter(is.finite(dysc_used))
    for (p in progs) {
      pooled <- suppressWarnings(cor(m$dysc_used, m[[p]], method = "spearman", use = "complete.obs"))
      dd <- m %>% filter(is.finite(.data[[p]])) %>% group_by(donor) %>%
        filter(dplyr::n() >= MIN_CELLS_DONOR) %>%
        summarise(rho = suppressWarnings(cor(dysc_used, .data[[p]], method = "spearman")),
                  n = dplyr::n(), .groups = "drop") %>% filter(is.finite(rho))
      wt <- if (nrow(dd) >= 3) suppressWarnings(wilcox.test(dd$rho, mu = 0)$p.value) else NA_real_
      rows[[length(rows) + 1]] <- tibble(dataset = ds, mode = mode, program = p,
                                         pooled_rho = pooled, donor_med_rho = median(dd$rho),
                                         donor_pos = sum(dd$rho > 0), donor_n = nrow(dd), donor_p = wt)
      donor_rows[[length(donor_rows) + 1]] <- dd %>% mutate(dataset = ds, mode = mode, program = p)
    }
  }
}
A  <- bind_rows(rows)
DR <- bind_rows(donor_rows)
write_csv(A,  "gene_program_dyscoordination_correlations.csv")
write_csv(DR, "gene_program_per_donor_rho.csv")

cat("\n### donor-level program correlations (ncount_corrected) ###\n")
print(as.data.frame(A %>% filter(mode == "ncount_corrected") %>%
  transmute(dataset, program, pooled = round(pooled_rho, 3), donor_med = round(donor_med_rho, 3),
            donors = sprintf("%d/%d", donor_pos, donor_n), p = signif(donor_p, 2))), right = FALSE)

######################################################
# Figure 4C: six-program forest (drop mito/oxidative, interferon, NF-kB/inflammatory).
# Grey points = per-donor rho; coloured point = donor median (dark if signed-rank p<0.05).
DROP <- c("mito_oxidative", "interferon", "nfkb_inflammatory")
PROG_LAB <- c(cytotoxicity = "Cytotoxicity", GZMK_inflammaging = "GZMK-associated inflammaging",
              exhaustion = "Exhaustion", proliferation = "Proliferation",
              generic_activation = "Generic activation", progenitor_memory = "Progenitor/memory")
DS_LAB <- c(glioma = "Glioma (Wang et al. 2025)", melanoma = "Melanoma (Wang et al. 2024)",
            Terekhova = "Aging (Terekhova et al. 2023)", Wang_aging = "Aging (Wang et al. 2025)")
ds_fac <- function(x) factor(unname(DS_LAB[x]), levels = unname(DS_LAB[DATASETS]))

prog_panel <- function(mode) {
  Am <- A  %>% filter(mode == !!mode, !program %in% DROP)
  Dm <- DR %>% filter(mode == !!mode, !program %in% DROP)
  ord <- Am %>% group_by(program) %>% summarise(m = mean(donor_med_rho, na.rm = TRUE)) %>%
    arrange(m) %>% pull(program)
  fac <- function(x) x %>% mutate(program = factor(unname(PROG_LAB[program]), levels = unname(PROG_LAB[ord])),
                                  dataset = ds_fac(dataset))
  ggplot(fac(Dm), aes(program, rho)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey45", linewidth = 0.3) +
    geom_jitter(width = 0.13, height = 0, size = 0.7, alpha = 0.45, colour = "grey45", stroke = 0) +
    geom_point(data = fac(Am), aes(program, donor_med_rho, colour = donor_p < 0.05), size = 1.6) +
    scale_colour_manual(values = c(`TRUE` = "#08306B", `FALSE` = "grey55"),
                        labels = c(`TRUE` = "donor median, signed-rank p < 0.05", `FALSE` = "donor median, n.s."),
                        breaks = c("TRUE", "FALSE"), name = NULL) +
    coord_flip() + facet_wrap(~ dataset, nrow = 1) +
    labs(x = NULL, y = "Spearman correlation between cell-level dyscoordination and T cell program score") +
    theme(legend.position = "bottom")
}

p4C <- prog_panel("ncount_corrected")
save_plot(p4C, "Fig4C_gene_programs.png", 6.7, 2.6)
ggsave("Fig4C_gene_programs.pdf", p4C, width = 6.7, height = 2.6)   # base pdf(): links cleanly into Illustrator

cat("\nT_cell_gene_programs done.\n")
