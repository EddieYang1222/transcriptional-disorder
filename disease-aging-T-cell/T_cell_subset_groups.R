# Transcriptional dyscoordination across coarse T-cell subset groups (Figure 4B)

# Set up working directory
# setwd("path/to/working/directory")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')
source('T_cell_common.R')

DATASETS  <- c("glioma", "melanoma", "Terekhova", "Wang_aging")
MIN_CELLS  <- 10      # cells per donor x cell_type stratum
MIN_DONORS <- 5       # donors a subset must appear in to be tested/plotted

# Dataset labels used in Figure 4
DS_LAB <- c(glioma = "Glioma (Wang et al. 2025)", melanoma = "Melanoma (Wang et al. 2024)",
            Terekhova = "Aging (Terekhova et al. 2023)", Wang_aging = "Aging (Wang et al. 2025)")

######################################################
# Coarse, investigator-defined subset grouping (a SUMMARY device, NOT a developmental or
# antigen-exposure trajectory). It orders subsets 0->3 by increasing cytotoxic/effector
# character: A naive | B central/stem memory | C effector-memory/helper/activated |
# D cytotoxic/NK-like/innate-like. Tregs, IFN, proliferative, DN and NME1+ subsets are
# ranked but left ungrouped (NA).
STAGE <- list(
  glioma = c(
    "Tnaive" = 0, "IL7R+ Tm" = 1,
    "GZMK+ Tem" = 2, "GZMK+ Trm" = 2, "Activated Trm" = 2,
    "NR4A2+ Teff" = 3, "Term Teff" = 3, "Temra" = 3, "IL7R+ Temra" = 3, "NK-like" = 3,
    "NME1+ T" = NA),
  melanoma = c(
    "SCM" = 1, "CM" = 1,
    "Early Activated" = 2, "Activated" = 2, "Early Effector" = 2,
    "Effector" = 3, "Exhausted" = 3, "NK-like" = 3,
    "IFN" = NA),
  Terekhova = c(
    "CD4_Naive" = 0, "CD8_Naive" = 0, "CD4_Naive-IFN" = 0, "CD8_Naive-IFN" = 0, "gd_gd naive" = 0,
    "CD8_Tcm CCR4+" = 1, "CD8_Tcm CCR4-" = 1,
    "CD4_Tfh" = 2, "CD4_Th1" = 2, "CD4_Th2" = 2, "CD4_Th17" = 2, "CD4_Th22" = 2,
    "CD4_Th1/Th17" = 2, "CD4_HLA-DR+ memory" = 2, "CD8_HLA-DR+" = 2,
    "CD4_Exhausted-like memory" = 2, "CD8_Trm" = 2, "CD8_Tmem KLRC2+" = 2,
    "CD8_Tem GZMK+" = 3, "CD8_Tem GZMB+" = 3, "CD8_Temra" = 3, "CD4_Temra" = 3,
    "CD4_Terminal effector" = 3, "CD8_NKT-like" = 3, "MAIT" = 3,
    "gd_Vd2 GZMK+" = 3, "gd_Vd2 GZMB+" = 3, "gd_Vd1 GZMB+" = 3, "gd_Vd1 GZMK+" = 3,
    "CD4_Treg naive" = NA, "CD4_Treg memory" = NA, "CD4_Treg KLRB1+RORC+" = NA,
    "CD4_Treg cytotoxic" = NA, "DN_T" = NA, "CD8_Proliferative" = NA),
  Wang_aging = c(
    "CD4_Naive_CCR7" = 0, "CD8_Naive_LEF1" = 0,
    "CD4_TCM_AQP3" = 1, "CD8_TCM_HAVCR2" = 1,
    "CD4_TEM_ANXA1" = 2, "CD8_TEM_CMC1" = 2, "CD8_TEM_ZNF683" = 2,
    "CD8_TEM_GNLY" = 3, "CD4_TEM_GNLY" = 3, "CD8_MAIT_SLC4A10" = 3, "gdT" = 3,
    "CD4_Treg_FOXP3" = NA)
)
STAGE_LAB <- c("A naive", "B central/stem mem", "C eff-mem/helper/act", "D cytotox/NK-like/innate")

######################################################
# Per-cell tables: Cell_barcode, cell_type, donor, dyscoordination (+ log), nCount_RNA.
# nCount_RNA is written into the cellular CSV by the *_deviation.R scripts. Glioma donor
# and depth come from GBM_CD8_obs.csv.
data_dir <- 'path/to/T_cell_aux'
build_cells <- function(ds) {
  if (ds == "glioma") {
    ent <- read_dyscoordination("glioma_T_cell_cellular_dispersion_SAVER.csv") %>%
      select(Cell_barcode, cell_type, dyscoordination, log_dyscoordination)
    obs <- read.csv(file.path(data_dir, "GBM_CD8_obs.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    names(obs)[1] <- "Cell_barcode"
    ent %>% inner_join(obs %>% transmute(Cell_barcode, donor = Patient_x, nCount_RNA = nCount_RNA),
                       by = "Cell_barcode")
  } else {
    csv   <- sprintf("%s_T_cell_cellular_dispersion_SAVER.csv", ds)
    donor_col <- c(melanoma = "patient_alias", Terekhova = "Donor_id", Wang_aging = "sampleName")[[ds]]
    read_dyscoordination(csv) %>% rename(donor = !!donor_col) %>%
      select(Cell_barcode, cell_type, donor, nCount_RNA, dyscoordination, log_dyscoordination)
  }
}

# donor x cell_type medians of dysc_used, within-donor centred, with stage label
donor_ct_table <- function(m, ds) {
  st <- STAGE[[ds]]
  m %>% filter(is.finite(dysc_used)) %>%
    group_by(donor, cell_type) %>%
    summarise(n_cells = dplyr::n(), med = median(dysc_used), .groups = "drop") %>%
    filter(n_cells >= MIN_CELLS) %>%
    group_by(donor) %>% mutate(med_c = med - mean(med)) %>% ungroup() %>%
    mutate(stage = unname(st[cell_type]),
           stage_lab = ifelse(is.na(stage), "not grouped", STAGE_LAB[stage + 1]))
}

wilcox_1s <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(c(n = length(x), med = NA, p = NA, pos = NA))
  c(n = length(x), med = median(x), p = suppressWarnings(wilcox.test(x, mu = 0)$p.value), pos = sum(x > 0))
}
per_donor_stage_rho <- function(tab) {
  tab %>% filter(!is.na(stage)) %>% group_by(donor) %>%
    filter(dplyr::n_distinct(stage) >= 3, dplyr::n() >= 4) %>%
    summarise(rho = suppressWarnings(cor(stage, med, method = "spearman")), .groups = "drop") %>%
    filter(is.finite(rho))
}
paired_contrast <- function(tab, hi = 3, lo = c(0, 1)) {
  tab %>% filter(!is.na(stage)) %>% group_by(donor) %>%
    summarise(hi_m = mean(med[stage %in% hi]), lo_m = mean(med[stage %in% lo]),
              n_hi = sum(stage %in% hi), n_lo = sum(stage %in% lo), .groups = "drop") %>%
    filter(n_hi >= 1, n_lo >= 1) %>% mutate(delta = hi_m - lo_m)
}
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 2e-16, "<2e-16", sprintf("%.1e", p)))

######################################################
# Compute (both modes) and collect the coarse-group tables
all_tabs <- list()
summ <- list()
for (ds in DATASETS) {
  m0 <- build_cells(ds) %>% filter(nCount_RNA > 0)
  for (mode in MODES) {
    m   <- add_dysc_used(m0, mode)
    tab <- donor_ct_table(m, ds) %>% mutate(dataset = ds, mode = mode)
    all_tabs[[paste(ds, mode)]] <- tab

    keep <- tab %>% count(cell_type) %>% filter(n >= MIN_DONORS) %>% pull(cell_type)
    tk   <- tab %>% filter(cell_type %in% keep)

    # ranking of subsets by within-donor-centred median
    rk <- tk %>% group_by(cell_type) %>%
      summarise(n_donors = dplyr::n(), n_cells = sum(n_cells), mean_c = mean(med_c),
                p_vs0 = unname(wilcox_1s(med_c)["p"]),
                stage = dplyr::first(stage), stage_lab = dplyr::first(stage_lab), .groups = "drop") %>%
      mutate(fdr = p.adjust(p_vs0, "BH")) %>% arrange(desc(mean_c)) %>%
      mutate(rank = row_number(), dataset = ds, mode = mode)
    write_csv(rk, sprintf("subset_ranking_%s_%s.csv", ds, mode))

    # T1 cell-type effect: Friedman on the complete donor x subset matrix
    wide <- tk %>% select(donor, cell_type, med) %>% pivot_wider(names_from = cell_type, values_from = med)
    cc <- wide %>% select(-donor)
    keep_ct <- names(cc)[colSums(!is.na(cc)) >= 0.8 * nrow(cc)]
    mat <- wide %>% select(donor, all_of(keep_ct)) %>% tidyr::drop_na()
    fried <- if (ncol(mat) >= 3 && nrow(mat) >= 5) friedman.test(as.matrix(mat[, -1])) else NULL

    # T2 ordered trend across coarse groups
    on_ax  <- tk %>% filter(!is.na(stage))
    s_pool <- spear(on_ax$stage, on_ax$med_c)
    pdr    <- per_donor_stage_rho(tk)
    w_pdr <- wilcox_1s(pdr$rho)

    # T3 paired contrast: group D vs groups A+B
    pc <- paired_contrast(tk, hi = 3, lo = c(0, 1))
    w3 <- wilcox_1s(pc$delta)

    cat(sprintf("\n== %s [%s] : %d donors, %d subsets ==\n", ds, mode, n_distinct(tab$donor), length(keep)))
    if (!is.null(fried)) cat(sprintf("  T1 Friedman: chi2=%.1f, p=%s (%d donors x %d subsets)\n",
                                     unname(fried$statistic), fmt_p(fried$p.value), nrow(mat), ncol(mat) - 1))
    cat(sprintf("  T2 stage vs centred median: rho=%.3f, p=%s; per-donor %d/%d donors >0, signed-rank p=%s\n",
                s_pool["rho"], fmt_p(s_pool["p"]), w_pdr["pos"], w_pdr["n"], fmt_p(w_pdr["p"])))
    cat(sprintf("  T3 delta [D] - [A+B]: median=%+.3f, %d/%d donors >0, p=%s\n",
                w3["med"], w3["pos"], w3["n"], fmt_p(w3["p"])))

    summ[[paste(ds, mode)]] <- tibble(
      dataset = ds, mode = mode, n_donors = n_distinct(tab$donor), n_subsets = length(keep),
      friedman_chi2 = if (is.null(fried)) NA else unname(fried$statistic),
      friedman_p = if (is.null(fried)) NA else fried$p.value,
      stage_rho = s_pool["rho"], stage_p = s_pool["p"],
      pdonor_pos = w_pdr["pos"], pdonor_n = w_pdr["n"], pdonor_p = w_pdr["p"],
      d31_med = w3["med"], d31_pos = w3["pos"], d31_n = w3["n"], d31_p = w3["p"])
  }
}
S  <- bind_rows(summ)
TT <- bind_rows(all_tabs)
write_csv(S, "subset_groups_summary_both_modes.csv")

######################################################
# Figure 4B: cross-dataset coarse-group panel (nCount-corrected, within-donor centred)
# Blue ramp, dark (A) -> light (D). ASCII text; base pdf() (no cairo).
RANK_LAB <- c("Naive", "Central/stem memory",
              "Effector memory/helper/activated", "Cytotoxic/NK-like/innate-like")
BLUE <- c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF")
RANK_COLS <- setNames(BLUE, RANK_LAB)
grp_fac <- function(x) factor(RANK_LAB[match(x, STAGE_LAB)], levels = RANK_LAB)
ds_fac  <- function(x) factor(unname(DS_LAB[x]), levels = unname(DS_LAB[DATASETS]))

group_panel <- function(mode, centred = TRUE) {
  d <- TT %>% filter(mode == !!mode, !is.na(stage)) %>%
    mutate(grp = grp_fac(stage_lab), dataset = ds_fac(dataset),
           v = if (centred) med_c else med)
  p <- ggplot(d, aes(grp, v))
  if (centred) p <- p + geom_hline(yintercept = 0, colour = "grey45", linetype = 2, linewidth = 0.3)
  p +
    geom_boxplot(aes(fill = grp), outlier.shape = NA, alpha = 0.9, linewidth = 0.3,
                 colour = "grey20", width = 0.68) +
    geom_jitter(width = 0.17, height = 0, size = 0.7, alpha = 0.5, colour = "grey15", stroke = 0) +
    scale_fill_manual(values = RANK_COLS, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    facet_wrap(~ dataset, nrow = 1) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(x = NULL, y = "Donor- and subset-level median\ndyscoordination (centered, corrected)", fill = NULL) +
    theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

p4B <- group_panel("ncount_corrected", centred = TRUE)
save_plot(p4B, "Fig4B_subset_group_centred.png", 7.3, 3.0)
ggsave("Fig4B_subset_group_centred.pdf", p4B, width = 7.3, height = 3.0)   # base pdf(): links cleanly into Illustrator

cat("\nT_cell_subset_groups done.\n")
