# Glioma patient-level transcriptional dyscoordination vs tumor grade (Figure 4D)

# Set up working directory
# setwd("path/to/working/directory")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')
source('T_cell_common.R')

# Inference unit = patient (donor). Grade and Age are patient-level. Question: does grade
# raise (depth-corrected) dyscoordination after adjusting for clone size and age, or is it
# absorbed by differentiation state (effector programs rise with grade, and effector cells
# are intrinsically high-dyscoordination)? nCount-corrected mode is the headline; raw kept.
GRADES <- c("Normal", "II", "III", "IV")
EFFECTOR_CT <- c("GZMK+ Tem","Temra","Term Teff","NR4A2+ Teff","NK-like","Activated Trm","GZMK+ Trm","IL7R+ Temra")
MEMNAIVE_CT <- c("IL7R+ Tm","Tnaive")

######################################################
# Build the glioma per-cell master: dyscoordination + depth + clone + program scores +
# patient covariates. Depth/clone/covariates come from GBM_CD8_obs.csv; program scores
# from module_scores_glioma.csv.gz (written by T_cell_gene_programs.R).
data_dir <- 'path/to/T_cell_aux'
bad_cdr3 <- c("", "NA", "NoneNone", "None")

ent <- read_dyscoordination("glioma_T_cell_cellular_dispersion_SAVER.csv") %>%
  select(Cell_barcode, cell_type, dyscoordination, log_dyscoordination)
obs <- read.csv(file.path(data_dir, "GBM_CD8_obs.csv"), check.names = FALSE, stringsAsFactors = FALSE)
names(obs)[1] <- "Cell_barcode"
obs <- obs %>% transmute(Cell_barcode, donor = Patient_x, nCount_RNA = nCount_RNA,
                         Grade = Tumor.Grade, Age = suppressWarnings(as.numeric(Age)), Sex,
                         cdr3 = as.character(cdr3))
sc  <- read_csv("module_scores_glioma.csv.gz", show_col_types = FALSE) %>%
  select(Cell_barcode, cytotoxicity, progenitor_memory)

m0 <- ent %>% inner_join(obs, by = "Cell_barcode") %>% left_join(sc, by = "Cell_barcode") %>%
  filter(!is.na(Grade), nCount_RNA > 0) %>%
  mutate(has_clone = !is.na(cdr3) & !(cdr3 %in% bad_cdr3),
         clone_id = ifelse(has_clone, paste0(donor, "|", cdr3), NA_character_),
         Grade = factor(Grade, levels = GRADES), grade_num = as.integer(Grade))
m0 <- m0 %>% group_by(clone_id) %>%
  mutate(clone_size = ifelse(is.na(clone_id), NA_integer_, dplyr::n())) %>% ungroup()

######################################################
# Grade (and sex) analysis, both modes
for (mode in MODES) {
  m <- add_dysc_used(m0, mode)
  cat(sprintf("\n== glioma [%s] grade + sex ==\n", mode))

  # ---- per-patient summary (donor-level inference unit) ----
  don <- m %>% group_by(donor) %>%
    summarise(Grade = first(Grade), grade_num = first(grade_num), Age = first(Age), Sex = first(Sex),
              n = dplyr::n(), mean_dysc = mean(dysc_used),
              mean_log_clonesize = mean(log(clone_size[has_clone]), na.rm = TRUE),
              mean_cyto = mean(cytotoxicity, na.rm = TRUE),
              mean_prog = mean(progenitor_memory, na.rm = TRUE),
              eff_frac  = mean(cell_type %in% EFFECTOR_CT), .groups = "drop")
  write_csv(don, sprintf("grade_patient_summary_%s.csv", mode))

  # ---- grade -> dyscoordination, crude vs sequential adjustment (donor-level) ----
  sp_crude  <- spear(don$grade_num, don$mean_dysc)
  kruskal_p <- suppressWarnings(kruskal.test(mean_dysc ~ Grade, don)$p.value)
  b_crude <- coef(lm(mean_dysc ~ grade_num, don))["grade_num"]
  b_age   <- coef(lm(mean_dysc ~ grade_num + Age, don))["grade_num"]
  b_full  <- coef(lm(mean_dysc ~ grade_num + Age + mean_log_clonesize, don))["grade_num"]
  b_diff  <- coef(lm(mean_dysc ~ grade_num + Age + mean_log_clonesize + mean_cyto + mean_prog, don))["grade_num"]
  p_full  <- summary(lm(mean_dysc ~ grade_num + Age + mean_log_clonesize, don))$coef["grade_num","Pr(>|t|)"]
  p_diff  <- summary(lm(mean_dysc ~ grade_num + Age + mean_log_clonesize + mean_cyto + mean_prog, don))$coef["grade_num","Pr(>|t|)"]
  sp_eff  <- spear(don$grade_num, don$eff_frac)
  sp_cyto <- spear(don$grade_num, don$mean_cyto)

  donor_stats <- tibble(
    quantity = c("Spearman grade~mean_dysc (crude)","Kruskal grade (p)",
                 "lm grade coef: crude","lm grade coef: +Age","lm grade coef: +Age+cloneSize",
                 "lm grade coef: +Age+cloneSize+differentiation",
                 "p(grade | Age+cloneSize)","p(grade | Age+cloneSize+diff)",
                 "Spearman grade~effector_fraction","Spearman grade~mean_cytotoxicity"),
    value = c(sp_crude["rho"], kruskal_p, b_crude, b_age, b_full, b_diff, p_full, p_diff,
              sp_eff["rho"], sp_cyto["rho"]))
  write_csv(donor_stats, sprintf("grade_donor_level_stats_%s.csv", mode))
  cat(sprintf("  grade~dyscoordination Spearman rho=%.2f (Kruskal p=%.3f, n=%d patients); grade slope crude=%.3f -> +diff=%.3f\n",
              sp_crude["rho"], kruskal_p, nrow(don), b_crude, b_diff))

  # ---- per-cell-type grade trend, adjusted for age + clone size ----
  dct <- m %>% group_by(donor, cell_type) %>% filter(dplyr::n() >= 20) %>%
    summarise(Grade = first(Grade), grade_num = first(grade_num), Age = first(Age), Sex = first(Sex),
              med = median(dysc_used),
              mlcs = mean(log(clone_size[has_clone]), na.rm = TRUE), .groups = "drop")
  cts <- sort(unique(dct$cell_type))
  per_ct <- lapply(cts, function(ct) {
    d <- dct %>% filter(cell_type == ct)
    if (dplyr::n_distinct(d$donor) < 8 || dplyr::n_distinct(d$grade_num) < 2) return(NULL)
    a  <- spear(d$grade_num, d$med)
    b0 <- tryCatch(coef(lm(med ~ grade_num, d))["grade_num"], error = function(e) NA)
    d2 <- d %>% filter(is.finite(mlcs))
    b1 <- tryCatch(coef(lm(med ~ grade_num + Age + mlcs, d2))["grade_num"], error = function(e) NA)
    p1 <- tryCatch(summary(lm(med ~ grade_num + Age + mlcs, d2))$coef["grade_num","Pr(>|t|)"], error = function(e) NA)
    data.frame(cell_type = ct, n_donors = nrow(d), rho_grade_crude = a["rho"], p_grade = a["p"],
               coef_crude = b0, coef_adj_age_clone = b1, p_adj = p1)
  }) %>% bind_rows()
  per_ct$p_grade_fdr <- p.adjust(per_ct$p_grade, "BH")
  write_csv(per_ct, sprintf("grade_by_celltype_%s.csv", mode))

  # ---- plots ----
  yl <- as.numeric(quantile(don$mean_dysc, c(.02, .98), na.rm = TRUE))
  # (Figure 4D) patient dyscoordination by tumor grade
  p1 <- ggplot(don, aes(Grade, mean_dysc)) +
    geom_violin(fill = "grey80", colour = "grey30", scale = "width", alpha = .5, draw_quantiles = .5) +
    geom_jitter(aes(colour = Age), width = .15, size = 2.4) + scale_colour_viridis_c() +
    coord_cartesian(ylim = yl) +
    labs(title = sprintf("Glioma [%s]: patient mean dyscoordination by tumor grade", mode),
         subtitle = sprintf("Spearman grade~dyscoordination rho=%.2f; Kruskal p=%.3f (n=%d patients)",
                            sp_crude["rho"], kruskal_p, nrow(don)),
         x = "tumor grade", y = ylab_of(mode))
  save_plot(p1, sprintf("Fig4D_dyscoordination_by_grade_%s.png", mode), 6, 5)
  ggsave(sprintf("Fig4D_dyscoordination_by_grade_%s.pdf", mode), p1, width = 6, height = 5)

  # grade slope, sequential adjustment (shrinkage at last step = absorbed by differentiation)
  att <- tibble(model = factor(c("crude","+Age","+Age +cloneSize","+Age +cloneSize +differentiation"),
                               levels = c("crude","+Age","+Age +cloneSize","+Age +cloneSize +differentiation")),
                grade_coef = c(b_crude, b_age, b_full, b_diff))
  p2 <- ggplot(att, aes(model, grade_coef, fill = model)) + geom_col(colour = "grey25", show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = 2) + geom_text(aes(label = sprintf("%.3f", grade_coef)), vjust = -0.4, size = 3) +
    labs(title = sprintf("Glioma [%s]: grade slope on dyscoordination, sequential adjustment", mode),
         x = NULL, y = "grade slope (per grade step)") + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  save_plot(p2, sprintf("grade_coef_attenuation_%s.png", mode), 7, 5)

  # per-cell-type grade trend (donor x cell-type medians)
  dct$cell_type <- factor(dct$cell_type, levels = cts)
  p3 <- ggplot(dct, aes(grade_num, med)) + geom_jitter(aes(colour = Age), width = .12, size = 1.3, alpha = .8) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = .6) + scale_colour_viridis_c() +
    scale_x_continuous(breaks = 1:4, labels = GRADES) + facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
    labs(title = sprintf("Glioma [%s]: dyscoordination vs grade by cell type", mode),
         x = "tumor grade", y = ylab_of(mode)) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p3, sprintf("grade_by_celltype_facet_%s.png", mode), 11, 8)

  # differentiation tracks grade: effector fraction + cytotoxicity by grade
  comp <- don %>% select(donor, Grade, eff_frac, mean_cyto) %>%
    pivot_longer(c(eff_frac, mean_cyto), names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric, eff_frac = "effector-like cell fraction", mean_cyto = "mean cytotoxicity score"))
  p4 <- ggplot(comp, aes(Grade, value)) + geom_violin(fill = "grey80", scale = "width", alpha = .5, draw_quantiles = .5) +
    geom_jitter(width = .15, size = 1.8, colour = "firebrick") + facet_wrap(~ metric, scales = "free_y") +
    labs(title = sprintf("Glioma [%s]: differentiation state by grade", mode),
         subtitle = sprintf("Spearman grade~effector frac=%.2f; grade~cytotoxicity=%.2f", sp_eff["rho"], sp_cyto["rho"]),
         x = "tumor grade", y = NULL)
  save_plot(p4, sprintf("differentiation_by_grade_%s.png", mode), 9, 5)

  # ---- sex stratification (secondary) ----
  sx_overall <- suppressWarnings(wilcox.test(mean_dysc ~ Sex, don)$p.value)
  sx_ct <- lapply(cts, function(ct) {
    d <- dct %>% filter(cell_type == ct)
    if (any(table(d$Sex) < 3) || dplyr::n_distinct(d$Sex) < 2) return(NULL)
    p <- suppressWarnings(wilcox.test(med ~ Sex, d)$p.value)
    medF <- median(d$med[d$Sex == "F"], na.rm = TRUE)
    medM <- median(d$med[d$Sex == "M"], na.rm = TRUE)
    data.frame(cell_type = ct, n_F = sum(d$Sex == "F"), n_M = sum(d$Sex == "M"),
               med_F = medF, med_M = medM, delta_F_minus_M = medF - medM, wilcox_p = p)
  }) %>% bind_rows()
  if (nrow(sx_ct)) sx_ct$p_fdr <- p.adjust(sx_ct$wilcox_p, "BH")
  write_csv(sx_ct, sprintf("sex_by_celltype_%s.csv", mode))
  write_csv(tibble(quantity = "donor mean_dysc ~ Sex (Wilcoxon p)", value = sx_overall,
                   n_F = sum(don$Sex == "F"), n_M = sum(don$Sex == "M")), sprintf("sex_overall_%s.csv", mode))
  p5 <- ggplot(don, aes(Sex, mean_dysc, fill = Sex)) +
    geom_violin(scale = "width", alpha = .5, draw_quantiles = .5, colour = "grey30") +
    geom_jitter(width = .12, size = 2.2, alpha = .8) +
    scale_fill_manual(values = c(F = "#d7191c", M = "#2c7bb6")) + coord_cartesian(ylim = yl) +
    labs(title = sprintf("Glioma [%s]: patient mean dyscoordination by sex", mode),
         subtitle = sprintf("Wilcoxon p=%.2f (%dF / %dM)", sx_overall, sum(don$Sex == "F"), sum(don$Sex == "M")),
         x = NULL, y = ylab_of(mode)) + theme(legend.position = "none")
  save_plot(p5, sprintf("dyscoordination_by_sex_%s.png", mode), 5.5, 5)
}

cat("\nT_cell_tumor_grade done.\n")
