# Longitudinal within-clone dyscoordination dynamics under ICB (Figure 4E)

# Set up working directory
# setwd("path/to/working/directory")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')
source('T_cell_common.R')

# Clones are tracked across timepoints by TCR (clone_id = patient | TRB_cdr3). Question:
# is the post-ICB dyscoordination change WITHIN-clone (scenario a, reversible) or
# COMPOSITIONAL (scenario b, depletion/selection)? Analyses A-D + a symmetric
# Marshall-Edgeworth shift-share decomposition adjudicate, for every consecutive
# transition. Raw dyscoordination on a log axis with rank-based tests throughout.
TP_LEVELS <- c("Baseline", "Follow Up 1", "Follow Up 2", "Follow Up 3")
TP_TAG    <- c("Baseline" = "BL", "Follow Up 1" = "FU1", "Follow Up 2" = "FU2", "Follow Up 3" = "FU3")
TRANSITIONS <- list(c("Baseline", "Follow Up 1"), c("Follow Up 1", "Follow Up 2"),
                    c("Follow Up 2", "Follow Up 3"))
MIN_CELLS <- 3L
set.seed(42)
tag_of <- function(t0, t1) paste(TP_TAG[[t0]], TP_TAG[[t1]], sep = "_")

######################################################
# Violin-on-points over an ordered x with a log y axis (the approved clone-expansion
# styling; used for the fate panels B2/C3). Heavy-tailed dyscoordination -> log y,
# display clipped to the central bulk (no data dropped).
violin_bin <- function(df, ycol, title, xvar, xcont, xlab, point_col = "dominant_cell_type") {
  yv  <- df[[ycol]]
  yv <- yv[is.finite(yv) & yv > 0]
  qlo <- as.numeric(quantile(yv, 0.002))
  qhi <- as.numeric(quantile(yv, 0.99))
  cdf <- df %>% group_by(.data[[xvar]]) %>% summarise(n = dplyr::n(), .groups = "drop") %>% mutate(.y = qlo)
  f <- function(x) { x <- x[is.finite(x[[xcont]]) & x[[ycol]] > 0, ]
    if (nrow(x) < 5) return(sprintf("n=%d", nrow(x)))
    s <- suppressWarnings(cor.test(x[[xcont]], x[[ycol]], method = "spearman"))
    sprintf("rho=%.2f, p=%.0e\nn=%d", unname(s$estimate), s$p.value, nrow(x)) }
  set.seed(42)
  ggplot(df, aes(.data[[xvar]], .data[[ycol]])) +
    geom_jitter(aes(colour = .data[[point_col]]), width = 0.32, height = 0, size = 0.5, alpha = 0.22) +
    geom_violin(fill = "grey75", colour = "grey15", linewidth = 0.55, scale = "width",
                draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.35) +
    stat_summary(aes(group = 1), fun = median, geom = "line", colour = "black", linewidth = 0.7) +
    stat_summary(fun = median, geom = "point", colour = "black", size = 1.3) +
    geom_text(data = cdf, aes(.data[[xvar]], .y, label = n), inherit.aes = FALSE, vjust = 1.4, size = 2.4, colour = "grey30") +
    geom_text(data = tibble(label = f(df)), aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.05, vjust = 1.2, size = 3.0, fontface = "bold") +
    scale_y_log10() + coord_cartesian(ylim = c(qlo, qhi)) +
    labs(x = xlab, y = paste0(ycol, " (log10, clipped)"), colour = "dominant\ncell type", title = title) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
}

######################################################
# Step 1: per-cell table (raw dyscoordination + TCR)
data_dir <- 'path/to/T_cell_aux'
ent <- read_dyscoordination("melanoma_T_cell_cellular_dispersion_SAVER.csv") %>%
  rename(patient = patient_alias) %>% select(Cell_barcode, cell_type, Timepoint, patient, dyscoordination)
tcr <- read_tsv(file.path(data_dir, "melanoma_T_cell_tcr.tsv"), show_col_types = FALSE) %>%
  select(Cell_barcode, TRB_cdr3)
bad <- c("", "NA", "None", "NoneNone")
cells <- ent %>% left_join(tcr, by = "Cell_barcode") %>%
  filter(!is.na(TRB_cdr3), !(TRB_cdr3 %in% bad)) %>%
  mutate(clone_id = paste0(patient, "|", TRB_cdr3),
         Timepoint = factor(Timepoint, levels = TP_LEVELS, ordered = TRUE))

pop <- cells %>% group_by(Timepoint) %>%
  summarise(n = dplyr::n(), mean_dyscoordination = mean(dyscoordination),
            median_dyscoordination = median(dyscoordination),
            mean_log_dyscoordination = mean(log(dyscoordination)), .groups = "drop")
write_csv(pop, "longitudinal_population_mean_by_timepoint.csv")
cat("Population dyscoordination by timepoint:\n")
print(as.data.frame(pop), digits = 4)

######################################################
# Step 2: clone x timepoint table
ctt <- cells %>% group_by(clone_id, Timepoint) %>%
  summarise(patient = dplyr::first(patient), size = dplyr::n(),
            mean_dyscoordination = mean(dyscoordination),
            mean_log_dyscoordination = mean(log(dyscoordination)),
            median_dyscoordination = median(dyscoordination),
            dominant_cell_type = { tb <- sort(table(cell_type), decreasing = TRUE)
            names(tb)[1] },
            .groups = "drop") %>%
  mutate(log_size = log10(size))
write_csv(ctt, "longitudinal_clone_timepoint_table.csv")

pair_clones <- function(tab, t0, t1) {
  a <- tab %>% filter(Timepoint == t0) %>%
    select(clone_id, size, mean_dyscoordination, mean_log_dyscoordination, median_dyscoordination, dominant_cell_type)
  b <- tab %>% filter(Timepoint == t1) %>%
    select(clone_id, size, mean_dyscoordination, mean_log_dyscoordination, median_dyscoordination)
  full_join(a, b, by = "clone_id", suffix = c("_t0", "_t1")) %>%
    mutate(patient = sub("\\|.*", "", clone_id),
           size_t0 = tidyr::replace_na(size_t0, 0L), size_t1 = tidyr::replace_na(size_t1, 0L),
           presence = dplyr::case_when(size_t0 > 0 & size_t1 > 0 ~ "both",
                                       size_t0 > 0 ~ "lost", TRUE ~ "new"))
}
wp_spearman <- function(d, xcol, ycol, min_clones = 8L) {
  per <- d %>% group_by(patient) %>% filter(dplyr::n() >= min_clones) %>%
    summarise(rho = suppressWarnings(cor(.data[[xcol]], .data[[ycol]], method = "spearman")),
              n = dplyr::n(), .groups = "drop") %>% filter(is.finite(rho))
  list(n_patients = nrow(per), median_rho = if (nrow(per)) median(per$rho) else NA_real_,
       frac_neg = if (nrow(per)) mean(per$rho < 0) else NA_real_,
       wilcox_p = if (nrow(per) >= 3) suppressWarnings(wilcox.test(per$rho, mu = 0)$p.value) else NA_real_)
}

# tracking power (all transitions)
power <- lapply(TRANSITIONS, function(tt) {
  pw <- pair_clones(ctt, tt[1], tt[2])
  data.frame(transition = paste(tt, collapse = " -> "),
             clones_t0_ge3 = sum(pw$size_t0 >= MIN_CELLS),
             persist_ge1 = sum(pw$size_t0 >= MIN_CELLS & pw$size_t1 >= 1),
             persist_ge3 = sum(pw$size_t0 >= MIN_CELLS & pw$size_t1 >= MIN_CELLS),
             disappear = sum(pw$size_t0 >= MIN_CELLS & pw$size_t1 == 0))
}) %>% bind_rows()
write_csv(power, "longitudinal_tracking_power.csv")

######################################################
# Symmetric Marshall-Edgeworth shift-share decomposition of the population-mean change:
# WITHIN (scenario a) vs BETWEEN / NEW / LOST (scenario b). No interaction term.
decomp <- function(size0, size1, m0, m1) {
  m0[size0 == 0] <- 0
  m1[size1 == 0] <- 0
  N0 <- sum(size0)
  N1 <- sum(size1)
  w0 <- size0 / N0
  w1 <- size1 / N1
  M0 <- sum(w0 * m0)
  M1 <- sum(w1 * m1)
  both <- size0 > 0 & size1 > 0
  lost <- size0 > 0 & size1 == 0
  new <- size0 == 0 & size1 > 0
  NEW <- sum((w1 * m1)[new])
  LOST <- sum((w0 * m0)[lost])
  wbar <- (w0 + w1) / 2
  mbar <- (m0 + m1) / 2
  dw <- w1 - w0
  dm <- m1 - m0
  WITHIN <- sum((wbar * dm)[both])
  BETWEEN <- sum((mbar * dw)[both])
  c(M0 = M0, M1 = M1, total = M1 - M0, WITHIN = WITHIN, BETWEEN = BETWEEN,
    NEW = NEW, LOST = LOST, resid = (M1 - M0) - (WITHIN + BETWEEN + NEW - LOST))
}
# patient-stratified clone bootstrap for the decomposition components
boot_decomp <- function(d, c0, c1, B = 2000) {
  idx <- split(seq_len(nrow(d)), sub("\\|.*", "", d$clone_id))
  comp <- matrix(NA_real_, B, 5, dimnames = list(NULL, c("WITHIN", "BETWEEN", "NEW", "LOST", "total")))
  for (b in seq_len(B)) {
    s <- unlist(lapply(idx, function(ix) sample(ix, length(ix), replace = TRUE)))
    dd <- d[s, ]
    comp[b, ] <- decomp(dd$size_t0, dd$size_t1, ifelse(is.na(dd[[c0]]), 0, dd[[c0]]),
                        ifelse(is.na(dd[[c1]]), 0, dd[[c1]]))[c("WITHIN", "BETWEEN", "NEW", "LOST", "total")]
  }
  comp
}
m_cols <- list(raw = c("mean_dyscoordination_t0", "mean_dyscoordination_t1"),
               median = c("median_dyscoordination_t0", "median_dyscoordination_t1"),
               log = c("mean_log_dyscoordination_t0", "mean_log_dyscoordination_t1"))

######################################################
# Per-transition analysis (A within-clone, B expansion/fate, C depletion, D decomposition)
run_transition <- function(t0, t1) {
  tag <- tag_of(t0, t1)
  lab <- paste(t0, "->", t1)
  wide <- pair_clones(ctt, t0, t1)

  ## ---- A. within-clone change (both >=3 cells) ----
  pwA <- wide %>% filter(presence == "both", size_t0 >= MIN_CELLS, size_t1 >= MIN_CELLS) %>%
    mutate(delta = mean_dyscoordination_t1 - mean_dyscoordination_t0,
           delta_log = mean_log_dyscoordination_t1 - mean_log_dyscoordination_t0,
           dir = ifelse(delta < 0, "decrease", "increase"))
  A <- data.frame(transition = lab, n_clones = nrow(pwA),
                  frac_decreasing = NA, median_ratio = NA, wilcox_p_raw = NA, wilcox_p_log = NA,
                  median_delta_log = NA, n_patients = NA, patients_neg = NA)
  if (nrow(pwA) >= 5) {
    wr <- suppressWarnings(wilcox.test(pwA$mean_dyscoordination_t1, pwA$mean_dyscoordination_t0, paired = TRUE))
    wl <- suppressWarnings(wilcox.test(pwA$mean_log_dyscoordination_t1, pwA$mean_log_dyscoordination_t0, paired = TRUE))
    pp <- pwA %>% group_by(patient) %>% summarise(md = median(delta_log), n = dplyr::n(), .groups = "drop") %>% filter(n >= 5)
    A <- data.frame(transition = lab, n_clones = nrow(pwA), frac_decreasing = mean(pwA$delta < 0),
                    median_ratio = median(pwA$mean_dyscoordination_t1 / pwA$mean_dyscoordination_t0),
                    wilcox_p_raw = wr$p.value, wilcox_p_log = wl$p.value,
                    median_delta_log = median(pwA$delta_log), n_patients = nrow(pp), patients_neg = sum(pp$md < 0))
    longA <- bind_rows(
      pwA %>% transmute(clone_id, patient, Timepoint = t0, dyscoordination = mean_dyscoordination_t0, dir),
      pwA %>% transmute(clone_id, patient, Timepoint = t1, dyscoordination = mean_dyscoordination_t1, dir)) %>%
      mutate(Timepoint = factor(Timepoint, levels = TP_LEVELS))
    pA <- ggplot(longA, aes(Timepoint, dyscoordination, group = clone_id)) +
      geom_line(aes(colour = dir), alpha = 0.25, linewidth = 0.3) +
      stat_summary(aes(group = 1), fun = median, geom = "line", colour = "black", linewidth = 1.1) +
      stat_summary(aes(group = 1), fun = median, geom = "point", colour = "black", size = 2) +
      scale_y_log10() + scale_colour_manual(values = c(decrease = "#2c7bb6", increase = "#d7191c")) +
      labs(title = sprintf("A. Within-clone dyscoordination, %s (>=3 cells both)", lab),
           subtitle = sprintf("Wilcoxon p(raw)=%.1e; %.0f%% of %d clones decrease; median ratio=%.2f",
                              wr$p.value, 100 * A$frac_decreasing, A$n_clones, A$median_ratio),
           y = "clone mean dyscoordination (log10)", colour = "direction")
    save_plot(pA, sprintf("Fig4E_A1_spaghetti_%s.png", tag), 6.5, 6)
    ggsave(sprintf("Fig4E_A1_spaghetti_%s.pdf", tag), pA, width = 6.5, height = 6)
    pA2 <- ggplot(pwA, aes(delta_log)) + geom_histogram(bins = 40, fill = "grey70", colour = "grey30") +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_vline(xintercept = median(pwA$delta_log), colour = "#2c7bb6", linewidth = 1) +
      labs(title = sprintf("A2. Within-clone delta log-dyscoordination, %s", lab),
           x = "delta log-dyscoordination (t1 - t0)", y = "clones")
    save_plot(pA2, sprintf("A2_delta_log_%s.png", tag), 6.5, 4.5)
  }

  ## ---- B. expansion vs prior dyscoordination + fate (t0 >=3) ----
  pwB <- wide %>% filter(size_t0 >= MIN_CELLS) %>%
    mutate(expansion = log2((size_t1 + 1) / (size_t0 + 1)),
           fate = factor(dplyr::case_when(size_t1 == 0 ~ "Disappeared", expansion < -1 ~ "Contracted",
                                          expansion > 1 ~ "Expanded", TRUE ~ "Stable"),
                         levels = c("Expanded", "Stable", "Contracted", "Disappeared")),
           fate_rank = as.integer(fate), mean_dyscoordination = mean_dyscoordination_t0)
  sp <- if (nrow(pwB) >= 5) suppressWarnings(cor.test(pwB$mean_dyscoordination_t0, pwB$expansion, method = "spearman")) else NULL
  wp <- wp_spearman(pwB, "mean_dyscoordination_t0", "expansion")
  kw <- if (nlevels(droplevels(pwB$fate)) >= 2) suppressWarnings(kruskal.test(mean_dyscoordination_t0 ~ fate, pwB)) else NULL
  B <- data.frame(transition = lab, n = nrow(pwB),
                  spearman = if (!is.null(sp)) unname(sp$estimate) else NA, spearman_p = if (!is.null(sp)) sp$p.value else NA,
                  wp_median_rho = wp$median_rho, kruskal_p = if (!is.null(kw)) kw$p.value else NA,
                  med_Expanded = median(pwB$mean_dyscoordination_t0[pwB$fate == "Expanded"]),
                  med_Disappeared = median(pwB$mean_dyscoordination_t0[pwB$fate == "Disappeared"]))
  if (nrow(pwB) >= 20) {
    pB2 <- violin_bin(pwB, "mean_dyscoordination", sprintf("B2. Dyscoordination at earlier tp by clone fate, %s", lab),
                      xvar = "fate", xcont = "fate_rank", xlab = sprintf("clone fate %s", lab))
    save_plot(pB2, sprintf("B2_fate_%s.png", tag), 7.5, 6)
  }

  ## ---- C. depletion: disappear vs persist baseline dyscoordination (t0 >=3) ----
  pwC <- pwB %>% mutate(fate = factor(ifelse(size_t1 == 0, "Disappeared", "Persisted"),
                                      levels = c("Persisted", "Disappeared")), fate_rank = as.integer(fate))
  wC <- if (all(table(pwC$fate) > 0)) suppressWarnings(wilcox.test(mean_dyscoordination ~ fate, pwC)) else NULL
  Cc <- data.frame(transition = lab, n_persist = sum(pwC$fate == "Persisted"),
                   n_disappear = sum(pwC$fate == "Disappeared"),
                   med_persist = median(pwC$mean_dyscoordination[pwC$fate == "Persisted"]),
                   med_disappear = median(pwC$mean_dyscoordination[pwC$fate == "Disappeared"]),
                   dyscoordination_wilcox_p = if (!is.null(wC)) wC$p.value else NA)
  if (!is.null(wC) && nrow(pwC) >= 20) {
    pC3 <- violin_bin(pwC, "mean_dyscoordination",
                      sprintf("C3. Earlier-tp dyscoordination: persist vs disappear, %s (p=%.1e)", lab, wC$p.value),
                      xvar = "fate", xcont = "fate_rank", xlab = sprintf("fate %s", lab))
    save_plot(pC3, sprintf("C3_persist_vs_disappear_%s.png", tag), 6.5, 6)
  }

  ## ---- D. decomposition (raw/median/log x global/shared>=3) + bootstrap ----
  Dg <- lapply(names(m_cols), function(s) { cc <- m_cols[[s]]
    data.frame(transition = lab, scale = s, variant = "global_all",
               t(decomp(wide$size_t0, wide$size_t1, ifelse(is.na(wide[[cc[1]]]), 0, wide[[cc[1]]]),
                        ifelse(is.na(wide[[cc[2]]]), 0, wide[[cc[2]]])))) }) %>% bind_rows()
  wrb <- wide %>% filter(size_t0 >= MIN_CELLS, size_t1 >= MIN_CELLS)
  Drb <- lapply(names(m_cols), function(s) { cc <- m_cols[[s]]
    data.frame(transition = lab, scale = s, variant = "shared_ge3",
               t(decomp(wrb$size_t0, wrb$size_t1, wrb[[cc[1]]], wrb[[cc[2]]]))) }) %>% bind_rows()
  set.seed(42)
  bc <- boot_decomp(wide, "mean_dyscoordination_t0", "mean_dyscoordination_t1", B = 2000)
  ci <- apply(bc, 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
  Dboot <- data.frame(transition = lab, component = colnames(bc), lo = ci[1, ], med = ci[2, ], hi = ci[3, ])
  Draw <- Dg[Dg$scale == "raw", ]
  comp_df <- data.frame(component = c("WITHIN", "BETWEEN", "NEW", "LOST"),
                        value = c(Draw$WITHIN, Draw$BETWEEN, Draw$NEW, -Draw$LOST)) %>%
    left_join(Dboot %>% transmute(component, lo = ifelse(component == "LOST", -hi, lo),
                                  hi = ifelse(component == "LOST", -lo, hi)), by = "component") %>%
    mutate(component = factor(component, levels = c("WITHIN", "BETWEEN", "NEW", "LOST")),
           sign = ifelse(value < 0, "decrease", "increase"))
  pD1 <- ggplot(comp_df, aes(component, value, fill = sign)) +
    geom_col(colour = "grey20", width = 0.7) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = c(decrease = "#2c7bb6", increase = "#d7191c"), guide = "none") +
    labs(title = sprintf("D1. Decomposition of the %s mean-dyscoordination change (raw)", lab),
         subtitle = sprintf("M_t0=%.2f -> M_t1=%.2f (total %.2f). WITHIN=scenario a; BETWEEN/NEW/LOST=scenario b.",
                            Draw$M0, Draw$M1, Draw$total),
         y = "contribution to change in population mean dyscoordination", x = NULL)
  save_plot(pD1, sprintf("Fig4E_D1_decomposition_%s.png", tag), 8.5, 5)
  ggsave(sprintf("Fig4E_D1_decomposition_%s.pdf", tag), pD1, width = 8.5, height = 5)

  within_share <- Draw$WITHIN / Draw$total
  list(A = A, B = B, C = Cc, Dg = Dg, Drb = Drb, Dboot = Dboot,
       headline = data.frame(transition = lab, total = Draw$total, WITHIN = Draw$WITHIN,
                             within_share = within_share, BETWEEN = Draw$BETWEEN, NEW = Draw$NEW, LOST = Draw$LOST,
                             within_lo = Dboot$lo[Dboot$component == "WITHIN"],
                             within_hi = Dboot$hi[Dboot$component == "WITHIN"]))
}

res <- lapply(TRANSITIONS, function(tt) run_transition(tt[1], tt[2]))

A_tab <- bind_rows(lapply(res, `[[`, "A"))
write_csv(A_tab, "longitudinal_analysisA_within_clone_paired.csv")
B_tab <- bind_rows(lapply(res, `[[`, "B"))
write_csv(B_tab, "longitudinal_analysisB_expansion.csv")
C_tab <- bind_rows(lapply(res, `[[`, "C"))
write_csv(C_tab, "longitudinal_analysisC_depletion.csv")
D_tab <- bind_rows(lapply(res, function(r) bind_rows(r$Dg, r$Drb)))
write_csv(D_tab, "longitudinal_analysisD_decomposition.csv")
Dboot_tab <- bind_rows(lapply(res, `[[`, "Dboot"))
write_csv(Dboot_tab, "longitudinal_analysisD_bootstrap_ci.csv")
head_tab <- bind_rows(lapply(res, `[[`, "headline"))
write_csv(head_tab, "longitudinal_decomposition_by_transition.csv")
cat("\nAnalysis A (within-clone paired):\n")
print(A_tab, digits = 3)
cat("\nDecomposition headline by transition:\n")
print(head_tab, digits = 3)

# Top-20 Baseline clones: size trajectories coloured by Baseline dyscoordination
bl <- ctt %>% filter(Timepoint == "Baseline")
top20 <- bl %>% slice_max(size, n = 20) %>% pull(clone_id)
traj <- ctt %>% filter(clone_id %in% top20) %>%
  left_join(bl %>% transmute(clone_id,
              grp = ifelse(mean_dyscoordination >= median(bl$mean_dyscoordination), "high", "low")), by = "clone_id")
pC2 <- ggplot(traj, aes(Timepoint, size, group = clone_id, colour = grp)) +
  geom_line(alpha = 0.8) + geom_point(size = 1) + scale_y_log10() +
  scale_colour_manual(values = c(high = "#d7191c", low = "#2c7bb6")) +
  labs(title = "Top-20 Baseline clones: size trajectory (colour = Baseline dyscoordination)",
       y = "clone size (log10)", colour = "Baseline\ndyscoordination")
save_plot(pC2, "C2_top_baseline_clone_trajectories.png", 7.5, 5)

cat("\nT_cell_longitudinal_ICB done.\n")
