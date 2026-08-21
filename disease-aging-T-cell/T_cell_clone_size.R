# Clone size vs clone-level transcriptional dyscoordination across four T-cell datasets (Figure 4A)

# Set up working directory
# setwd("path/to/working/directory")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')
source('T_cell_common.R')

######################################################
# Inputs per dataset
# Each dataset contributes a per-cell table: Cell_barcode, cell_type, donor,
# dyscoordination (+ log), nCount_RNA (written into the cellular CSV by the
# *_deviation.R scripts), clone_id, has_clone. The clonotype is joined from an
# auxiliary table extracted once from the source objects or VDJ downloads:
#   *_tcr.tsv : Cell_barcode, TRA_cdr3, TRB_cdr3 (paired alpha/beta CDR3)
# The glioma per-cell obs table (GBM_CD8_obs.csv) carries donor, cdr3 and depth.
# clone_id = donor | clonotype, so clones are within-donor by construction.
data_dir <- 'path/to/T_cell_aux'
DATASETS <- c("glioma", "melanoma", "Terekhova", "Wang_aging")
bad_cdr3 <- c("", "NA", "None", "NoneNone")

build_cells <- function(ds) {
  if (ds == "glioma") {
    ent <- read_dyscoordination("glioma_T_cell_cellular_dispersion_SAVER.csv") %>%
      select(Cell_barcode, cell_type, dyscoordination, log_dyscoordination)
    obs <- read.csv(file.path(data_dir, "GBM_CD8_obs.csv"), check.names = FALSE, stringsAsFactors = FALSE)
    names(obs)[1] <- "Cell_barcode"
    obs <- obs %>% transmute(Cell_barcode, donor = Patient_x, nCount_RNA = nCount_RNA,
                             cdr3 = as.character(cdr3))
    ent %>% inner_join(obs, by = "Cell_barcode") %>%
      mutate(has_clone = !is.na(cdr3) & !(cdr3 %in% bad_cdr3),
             clone_id = ifelse(has_clone, paste0(donor, "|", cdr3), NA_character_))

  } else if (ds == "melanoma") {
    ent <- read_dyscoordination("melanoma_T_cell_cellular_dispersion_SAVER.csv") %>%
      rename(donor = patient_alias) %>%
      select(Cell_barcode, cell_type, donor, nCount_RNA, dyscoordination, log_dyscoordination)
    tcr <- read_tsv(file.path(data_dir, "melanoma_T_cell_tcr.tsv"), show_col_types = FALSE)
    ent %>%
      left_join(tcr %>% select(Cell_barcode, TRB_cdr3), by = "Cell_barcode") %>%
      mutate(has_clone = !is.na(TRB_cdr3) & !(TRB_cdr3 %in% bad_cdr3),
             clone_id = ifelse(has_clone, paste0(donor, "|", TRB_cdr3), NA_character_))

  } else if (ds == "Terekhova") {
    ent <- read_dyscoordination("Terekhova_T_cell_cellular_dispersion_SAVER.csv") %>%
      rename(donor = Donor_id) %>%
      select(Cell_barcode, cell_type, donor, nCount_RNA, dyscoordination, log_dyscoordination)
    tcr <- read_tsv(file.path(data_dir, "Terekhova_T_cell_tcr.tsv"), show_col_types = FALSE)
    ent %>%
      left_join(tcr %>% select(Cell_barcode, TRA_cdr3, TRB_cdr3), by = "Cell_barcode") %>%
      mutate(has_clone = !is.na(TRB_cdr3) & !(TRB_cdr3 %in% bad_cdr3),
             clone_id = ifelse(has_clone, paste0(donor, "|", TRB_cdr3), NA_character_))

  } else if (ds == "Wang_aging") {
    ent <- read_dyscoordination("Wang_aging_T_cell_cellular_dispersion_SAVER.csv") %>%
      rename(donor = sampleName) %>%
      select(Cell_barcode, cell_type, donor, nCount_RNA, dyscoordination, log_dyscoordination)
    tcr <- read_tsv(file.path(data_dir, "Wang_aging_T_cell_tcr.tsv"), show_col_types = FALSE)
    ent %>%
      left_join(tcr %>% select(Cell_barcode, TRB_cdr3), by = "Cell_barcode") %>%
      mutate(has_clone = !is.na(TRB_cdr3) & !(TRB_cdr3 %in% bad_cdr3),
             clone_id = ifelse(has_clone, paste0(donor, "|", TRB_cdr3), NA_character_))
  }
}

######################################################
# Per-donor Spearman(clone size, clone-level dyscoordination), both modes
# Clone-level dyscoordination = mean of dysc_used over the clone's cells; clone size =
# number of cells. Pooled Spearman mixes a between-donor confound (age, grade constant
# within donor) with the within-donor component; within_group_rho() group-mean-centres
# per donor to isolate the deconfounded within-donor relationship.
combo <- list()
for (ds in DATASETS) {
  m0 <- build_cells(ds) %>% filter(nCount_RNA > 0)
  for (mode in MODES) {
    m   <- add_dysc_used(m0, mode)
    agg <- clone_agg(m) %>% filter(size >= 2)

    pooled <- spear(agg$log_size, agg$mean_dysc)
    wp     <- within_group_rho(agg, "log_size", "mean_dysc", gcol = "donor")
    summ <- bind_rows(
      tibble(stratum = "ALL pooled (size>=2)", rho = pooled["rho"], p = pooled["p"],
             n = pooled["n"], n_donors = NA_real_, frac_pos = NA_real_, wilcox_p = NA_real_),
      tibble(stratum = "within-donor partial", rho = wp$partial_rho, p = wp$p, n = wp$n,
             n_donors = wp$n_groups, frac_pos = wp$frac_pos, wilcox_p = wp$wilcox_p))
    write_csv(summ, sprintf("clone_size_dyscoordination_%s_%s.csv", ds, mode))
    cat(sprintf("%-11s [%-16s] pooled rho=%+.3f (p=%.1e) | within-donor partial=%+.3f (p=%.1e), %.0f%% donors +\n",
                ds, mode, pooled["rho"], pooled["p"], wp$partial_rho, wp$p, 100 * wp$frac_pos))

    combo[[paste(ds, mode)]] <- agg %>% transmute(Dataset = ds, mode = mode,
                                                  size, mean_dysc, log_size, size_bin)
  }
}

######################################################
# Figure 4A: clone-level dyscoordination by expansion stratum, all datasets (nCount-corrected)
for (mode in MODES) {
  ag <- bind_rows(combo[grepl(paste0(" ", mode, "$"), names(combo))])
  ag$Dataset <- factor(ag$Dataset, levels = DATASETS)
  save_plot(violin_lin(ag, "mean_dysc",
                       sprintf("Clone-level dyscoordination by expansion stratum - all datasets [%s]", mode),
                       facet_var = "Dataset", ncol = 2, point_col = NULL, ylab = ylab_of(mode)),
            sprintf("Fig4A_clone_size_dyscoordination_%s.png", mode), 12, 9)
}

cat("\nT_cell_clone_size done.\n")
