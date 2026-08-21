# Shared helpers for the T-cell transcriptional dyscoordination analyses (Figure 4)

# Consolidated from the per-cell dyscoordination measures. Every downstream Figure-4
# script sources this file (and Transcriptional_dyscoordination_functions.R) and reads
# the *_cellular_dispersion_SAVER.csv tables written by the *_deviation.R scripts.
#
# Two conventions used throughout:
#   * inference unit = donor / patient (cells are pseudoreplicates); rank-based stats
#     (Spearman / Wilcoxon) on the heavy-tailed, strictly-positive Cell_level_deviation.
#   * depth is the dominant technical confound (nCount_RNA). Every result is produced in
#     two modes:
#       raw               dysc_used = log(Cell_level_deviation)
#       ncount_corrected  dysc_used = residual of a GLOBAL loess(log dyscoordination ~ nCount_RNA)
#     Both are real-valued and analysed identically.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

theme_set(theme_bw(base_size = 11))

MODES     <- c("raw", "ncount_corrected")
sz_breaks <- c(0, 1, 2, 5, 10, 20, 50, 100, Inf)          # clone-size (expansion) strata
sz_labs   <- c("1", "2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100")

# T-cell gene programs used for the program-correlation panel (Figure 4C)
PROG   <- c("GZMK_inflammaging", "cytotoxicity", "exhaustion", "progenitor_memory",
            "proliferation", "interferon", "nfkb_inflammatory", "mito_oxidative",
            "generic_activation")
INFLAM <- c("GZMK_inflammaging", "cytotoxicity", "exhaustion", "interferon", "nfkb_inflammatory")

######################################################
# Read a cellular-dispersion table and standardize the per-cell measure column.
# Keeps the on-disk column Cell_level_deviation but adds a working `dyscoordination`
# (> 0) and its log for the two analysis modes.
read_dyscoordination <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    rename(dyscoordination = Cell_level_deviation) %>%
    filter(dyscoordination > 0, is.finite(dyscoordination)) %>%
    mutate(log_dyscoordination = log(dyscoordination))
}

######################################################
# Add dysc_used for a mode. ncount_corrected removes library size only (global loess
# across all cells), preserving cell-type / state biology; requires an nCount_RNA column.
add_dysc_used <- function(m, mode) {
  if (mode == "raw") {
    m$dysc_used <- m$log_dyscoordination
  } else {
    fit <- loess(log_dyscoordination ~ nCount_RNA, data = m, span = 1)
    m$dysc_used <- m$log_dyscoordination - predict(fit, newdata = m)
    na <- is.na(m$dysc_used)
    m$dysc_used[na] <- m$log_dyscoordination[na] - mean(m$log_dyscoordination, na.rm = TRUE)
  }
  m
}

ylab_of <- function(mode)
  if (mode == "raw") "dyscoordination (log scale)" else "dyscoordination (nCount-corrected residual)"

######################################################
# Rank-based correlation summary (returns rho, p, n)
spear <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(c(rho = NA, p = NA, n = sum(ok)))
  s <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  c(rho = unname(s$estimate), p = s$p.value, n = sum(ok))
}

######################################################
# Within-group (group-mean-centred) partial Spearman: rank both axes WITHIN each group
# (donor / patient), centre, pool. This removes ALL group-level confounders (age, grade,
# sex, batch, depth-per-donor) exactly, since they are constant within a group. Also
# returns the distribution of per-group Spearman rho and a signed-rank test on them.
within_group_rho <- function(df, xcol, ycol, gcol = "donor", min_g = 8L) {
  d <- df %>%
    filter(is.finite(.data[[xcol]]), is.finite(.data[[ycol]])) %>%
    group_by(.data[[gcol]]) %>% filter(dplyr::n() >= min_g) %>%
    mutate(rx = rank(.data[[xcol]]) - mean(rank(.data[[xcol]])),
           ry = rank(.data[[ycol]]) - mean(rank(.data[[ycol]]))) %>% ungroup()
  if (nrow(d) < 5)
    return(tibble(n = 0L, n_groups = 0L, partial_rho = NA_real_, p = NA_real_,
                  median_group_rho = NA_real_, frac_pos = NA_real_, wilcox_p = NA_real_))
  pe  <- suppressWarnings(cor.test(d$rx, d$ry))                     # pooled centred-rank corr
  per <- d %>% group_by(.data[[gcol]]) %>%
    summarise(rho = suppressWarnings(cor(.data[[xcol]], .data[[ycol]], method = "spearman")),
              .groups = "drop") %>% filter(is.finite(rho))
  wt <- if (nrow(per) >= 3) suppressWarnings(wilcox.test(per$rho, mu = 0)$p.value) else NA_real_
  tibble(n = nrow(d), n_groups = nrow(per), partial_rho = unname(pe$estimate), p = pe$p.value,
         median_group_rho = median(per$rho), frac_pos = mean(per$rho > 0), wilcox_p = wt)
}

######################################################
# Clone-level aggregate on dysc_used. `cells` needs clone_id, has_clone, donor,
# cell_type, dysc_used. One row per clone, pooled across cell types (clones straddle
# several subsets); dominant_cell_type recorded. Clone size = number of cells with a
# dyscoordination value that carry the clonotype.
clone_agg <- function(cells) {
  cells %>%
    filter(has_clone, is.finite(dysc_used)) %>%
    group_by(clone_id) %>%
    summarise(donor = dplyr::first(donor), size = dplyr::n(), mean_dysc = mean(dysc_used),
              dominant_cell_type = { tb <- sort(table(cell_type), decreasing = TRUE)
              names(tb)[1] },
              .groups = "drop") %>%
    mutate(log_size = log10(size), size_bin = cut(size, sz_breaks, sz_labs))
}

######################################################
# Violin-on-points on a LINEAR clipped axis (handles residuals). x is an ordered factor
# `xvar`; y is `ycol`; point colour `point_col` (NULL = constant). Prints per-bin n and a
# Spearman(continuous size, y) label. This is the approved clone-expansion styling.
violin_lin <- function(df, ycol, title, xvar = "size_bin", xcont = "log_size", facet_var = NULL,
                       point_col = "dominant_cell_type", ncol = NULL,
                       xlab = "clone size (expansion stratum)", ylab = "dyscoordination") {
  d <- df %>% filter(is.finite(.data[[ycol]]))
  ylim <- as.numeric(quantile(d[[ycol]], c(0.01, 0.99), na.rm = TRUE))
  cdf <- d %>% group_by(across(all_of(c(facet_var, xvar)))) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>% mutate(.y = ylim[1])
  labf <- function(x) {
    x <- x[is.finite(x[[xcont]]) & is.finite(x[[ycol]]), ]
    if (nrow(x) < 5) return(sprintf("n=%d", nrow(x)))
    s <- suppressWarnings(cor.test(x[[xcont]], x[[ycol]], method = "spearman"))
    sprintf("rho=%.2f, p=%.0e\nn=%d", unname(s$estimate), s$p.value, nrow(x)) }
  lab <- if (is.null(facet_var)) tibble(label = labf(d)) else
    d %>% group_by(across(all_of(facet_var))) %>% group_modify(~ tibble(label = labf(.x))) %>% ungroup()
  set.seed(42)
  jit <- if (is.null(point_col))
    geom_jitter(width = 0.32, height = 0, size = 0.5, alpha = 0.25, colour = "steelblue4") else
    geom_jitter(aes(colour = .data[[point_col]]), width = 0.32, height = 0, size = 0.5, alpha = 0.25)
  p <- ggplot(d, aes(.data[[xvar]], .data[[ycol]])) + jit +
    geom_violin(fill = "grey75", colour = "grey15", linewidth = 0.55, scale = "width",
                alpha = 0.35, draw_quantiles = 0.5) +
    stat_summary(aes(group = 1), fun = median, geom = "line", colour = "black", linewidth = 0.7) +
    stat_summary(fun = median, geom = "point", colour = "black", size = 1.3) +
    geom_text(data = cdf, aes(.data[[xvar]], .y, label = n), inherit.aes = FALSE,
              vjust = 1.3, size = 2.3, colour = "grey30") +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.05, vjust = 1.2, size = 2.9, fontface = "bold") +
    coord_cartesian(ylim = ylim) +
    labs(x = xlab, y = ylab, colour = "dominant\ncell type", title = title) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = if (is.null(point_col)) "none" else "bottom")
  if (!is.null(facet_var)) p <- p + facet_wrap(stats::as.formula(paste("~", facet_var)),
                                               scales = "free_y", ncol = ncol)
  p
}

######################################################
# Save a plot to a bare relative path (PNG). Base pdf() is used elsewhere for PDFs
# per the plotting convention; here dpi is fine for PNG previews.
save_plot <- function(p, path, w = 8, h = 6) {
  ggsave(path, p, width = w, height = h, dpi = 150)
  cat("  wrote", path, "\n")
}
