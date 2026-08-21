# In-silico spike-in simulation figures (Supplementary Note 1)
# DE-gene count vs dyscoordination ratio across three substrates, plus the
# prevalence and gene-count learnability sweeps. Reads the replicate metrics
# written by spike_in_simulation.R.

# Online links
# https://www.biorxiv.org/content/10.64898/2026.01.24.701460v1

# Set up working directory
# setwd("path/to/working/directory")

library(ggplot2)
library(dplyr)
library(patchwork)

rd <- function(f) read.csv(file.path("results", f))
# Mean with a t-based CI across replicates
ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- mean(x)
  s <- sd(x) / sqrt(max(1, n))
  tc <- qt(0.975, max(1, n - 1))
  c(m = m, lo = m - tc * s, hi = m + tc * s)
}
pal  <- c(Coordinated = "#2c7fb8", Dyscoordinated = "#c0392b", Both = "#6a51a3")
YLAB <- "Replicate-level dyscoordination ratio"

######################################################
# 1. DE-gene count vs dyscoordination ratio, per arm, across the three substrates
ds <- c(`Kidney PT cell (rat)`      = "replicate_metrics_B_primary.csv",
        `Liver hepatocyte (mouse)`  = "replicate_metrics_BLIVER_primary.csv",
        `Bone marrow HPC (mouse)`   = "replicate_metrics_BTMS_primary.csv")
A <- bind_rows(lapply(names(ds), function(nm) { d <- rd(ds[[nm]])
d$dataset <- nm
d }))
A$dataset <- factor(A$dataset, levels = names(ds))
A$arm <- factor(A$arm, levels = c("coord", "dyscoord", "both"),
                labels = c("Coordinated", "Dyscoordinated", "Both"))

set.seed(1)
pA <- ggplot(A, aes(n_de_spike, dys_ratio, color = arm)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey75") +
  geom_point(position = position_jitter(width = 3, height = 0), size = 1.5, alpha = 0.75) +
  facet_wrap(~ dataset, ncol = 1) +
  scale_color_manual(values = pal) +
  scale_x_continuous(breaks = c(0, 50, 100, 150), limits = c(-6, 162)) +
  labs(x = "Spike-in genes called DE at FDR 5% (of 150)", y = YLAB, color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("spike_in_de_vs_dyscoord.pdf", pA, width = 4.3, height = 8)

######################################################
# 2. Learnability sweeps (rat kidney): prevalence and spike-in gene count
# CI drawn as vertical lines; coordinated and dyscoordinated arms only.

# Prevalence sweep: percentage of cells carrying the coordinated program
prev <- c(`5`  = "replicate_metrics_B_prev05.csv", `10` = "replicate_metrics_B_prev10.csv",
          `20` = "replicate_metrics_B_prev20.csv", `30` = "replicate_metrics_B_prev30.csv",
          `50` = "replicate_metrics_B_prev50.csv")
B <- bind_rows(lapply(names(prev), function(k) { d <- rd(prev[[k]])
d$x <- as.integer(k)
d })) %>%
  filter(arm %in% c("coord", "dyscoord")) %>%
  group_by(x, arm) %>%
  summarise(m = ci(dys_ratio)["m"], lo = ci(dys_ratio)["lo"], hi = ci(dys_ratio)["hi"], .groups = "drop")

# Gene-count sweep: number of spike-in genes (prev20 = 150-gene midpoint)
gs <- c(`30`  = "replicate_metrics_B_ng20u10d.csv", `75`  = "replicate_metrics_B_ng50u25d.csv",
        `150` = "replicate_metrics_B_prev20.csv",    `300` = "replicate_metrics_B_ng200u100d.csv",
        `600` = "replicate_metrics_B_ng400u200d.csv")
C <- bind_rows(lapply(names(gs), function(k) { d <- rd(gs[[k]])
d$x <- as.integer(k)
d })) %>%
  filter(arm %in% c("coord", "dyscoord")) %>%
  group_by(x, arm) %>%
  summarise(m = ci(dys_ratio)["m"], lo = ci(dys_ratio)["lo"], hi = ci(dys_ratio)["hi"], .groups = "drop")

lab_arm <- function(d) {
  d$arm <- factor(d$arm, levels = c("coord", "dyscoord"),
                  labels = c("Coordinated", "Dyscoordinated"))
  d
}
B <- lab_arm(B)
C <- lab_arm(C)

mk <- function(d, xlab, xbreaks, logx = FALSE) {
  p <- ggplot(d, aes(x, m, color = arm)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
    geom_line(linewidth = 0.7) +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.7) +
    geom_point(size = 1.9) +
    scale_color_manual(values = pal[c("Coordinated", "Dyscoordinated")]) +
    labs(x = xlab, y = YLAB, color = NULL) +
    theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
  if (logx) p <- p + scale_x_log10(breaks = xbreaks) else p <- p + scale_x_continuous(breaks = xbreaks)
  p
}
pB <- mk(B, "Percentage of cells carrying program", c(5, 10, 20, 30, 50))
pC <- mk(C, "Number of spike-in genes", c(30, 75, 150, 300, 600), logx = TRUE)

combined <- (pB / pC) + plot_layout(guides = "collect") & theme(legend.position = "top")
ggsave("spike_in_sweeps.pdf", combined, width = 5.4, height = 6.4)

######################################################
