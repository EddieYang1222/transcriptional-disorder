# Allele Specific Expression Summary (Figure 2C and Supplementary Figure 8)

# Set up working directory
# setwd("path/to/working/directory")

library(ggplot2)
library(grid)

# Point estimates and gene-bootstrap 95% confidence intervals are assembled from
# the per-dataset ASE analyses (all_rho / h_rho = partial Spearman over all
# covered genes / the top-500 highly variable genes). Base pdf() only, ASCII text.

######################################################
# Figure 2C: six-dataset forest of partial Spearman (dyscoordination vs. discordance)
# Single-cross datasets first (Larsson, Ochiai, van der Veeken, Pritykin), then the
# multi-cross datasets (Weber, Medina-Cano), ordered by sequence divergence.
d <- read.csv(text = "
dataset,disp,ord,all_rho,all_lo,all_hi,h_rho,h_lo,h_hi
Larsson et al. 2019,C57BL/6J x CAST/EiJ,1,0.404,0.320,0.482,0.407,0.325,0.482
Ochiai et al. 2020,129 x CAST/EiJ,2,0.230,0.201,0.258,0.167,0.060,0.258
van der Veeken et al. 2019,SPRET/EiJ,3,0.079,0.044,0.113,0.355,0.234,0.473
Pritykin et al. 2021,SPRET/EiJ (pooled),4,0.134,0.089,0.176,0.350,0.199,0.483
Weber et al. 2026,129S1/SvImJ,6,0.284,0.236,0.325,0.256,0.121,0.387
Weber et al. 2026,A/J,7,0.267,0.223,0.312,0.256,0.123,0.397
Weber et al. 2026,NOD/ShiLtJ,8,0.276,0.229,0.318,0.220,0.078,0.353
Weber et al. 2026,NZO/HlLtJ,9,0.287,0.245,0.328,0.319,0.180,0.451
Weber et al. 2026,WSB/EiJ,10,0.312,0.273,0.349,0.301,0.196,0.396
Weber et al. 2026,CAST/EiJ,11,0.288,0.258,0.317,0.443,0.356,0.523
Weber et al. 2026,PWK/PhJ,12,0.308,0.283,0.331,0.507,0.427,0.583
Medina-Cano et al. 2025,MOLF/EiJ,13,0.198,0.144,0.249,0.302,0.073,0.504
Medina-Cano et al. 2025,CAST/EiJ,14,0.116,0.062,0.172,0.222,-0.014,0.447
Medina-Cano et al. 2025,PWK/PhJ,15,0.223,0.169,0.274,0.335,0.085,0.536
Medina-Cano et al. 2025,SPRET/EiJ,16,0.185,0.147,0.221,0.458,0.298,0.586
", stringsAsFactors = FALSE, strip.white = TRUE)

ds_levels <- c("Larsson et al. 2019", "Ochiai et al. 2020", "van der Veeken et al. 2019",
               "Pritykin et al. 2021", "Weber et al. 2026", "Medina-Cano et al. 2025")
d$dataset <- factor(d$dataset, levels = ds_levels)
d$key <- paste(sprintf("%02d", d$ord), d$disp)
d$key <- factor(d$key, levels = d$key[order(-d$ord)])

long <- rbind(
  data.frame(dataset = d$dataset, key = d$key, series = "All genes",
             rho = d$all_rho, lo = d$all_lo, hi = d$all_hi),
  data.frame(dataset = d$dataset, key = d$key, series = "Top 500 HVGs",
             rho = d$h_rho, lo = d$h_lo, hi = d$h_hi))
long$series <- factor(long$series, levels = c("All genes", "Top 500 HVGs"))

labmap <- setNames(sub("^[0-9]+ ", "", levels(d$key)), levels(d$key))
sepdf  <- data.frame(dataset = factor(head(ds_levels, -1), levels = ds_levels))
pd     <- position_dodge(width = 0.62)

p_forest <- ggplot(long, aes(rho, key, color = series)) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.45) +
  geom_hline(data = sepdf, aes(yintercept = 0.5), color = "grey80", linewidth = 0.3) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = pd, linewidth = 0.6) +
  geom_point(position = pd, size = 2.1) +
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = labmap, expand = expansion(add = 0.6)) +
  scale_x_continuous(breaks = seq(0, 0.6, 0.2)) +
  coord_cartesian(xlim = c(-0.02, 0.6)) +
  scale_color_manual(values = c("All genes" = "grey55", "Top 500 HVGs" = "#C0392B"), name = NULL) +
  labs(x = "Partial Spearman rho (dyscoordination vs. discordance)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey65", linewidth = 0.35),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.spacing.y = unit(0, "pt"), strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 9.5),
        axis.ticks.y = element_blank(), axis.text.y = element_text(size = 9.5, color = "grey25"),
        legend.position = "top", legend.margin = margin(0, 0, 2, 0),
        plot.margin = margin(5, 10, 5, 4))

pdf("ASE_summary_forest.pdf", width = 10.2, height = 5.8)
grid.draw(ggplotGrob(p_forest))
dev.off()

######################################################
# Supplementary Figure 8: top-500-HVG partial Spearman per cross vs the parental strain's
# genome-wide sequence divergence from C57BL/6J (approx. Mouse Genomes Project
# SNVs, in millions), for the two multi-cross datasets.
d2 <- read.csv(text = "
dataset,strain,snv,rho,lo,hi
Weber et al. 2026,129S1/SvImJ,4.6,0.256,0.121,0.387
Weber et al. 2026,A/J,4.7,0.256,0.123,0.397
Weber et al. 2026,NOD/ShiLtJ,5.0,0.220,0.078,0.353
Weber et al. 2026,NZO/HlLtJ,5.3,0.319,0.180,0.451
Weber et al. 2026,WSB/EiJ,6.2,0.301,0.196,0.396
Weber et al. 2026,CAST/EiJ,17.6,0.443,0.356,0.523
Weber et al. 2026,PWK/PhJ,17.1,0.507,0.427,0.583
Medina-Cano et al. 2025,MOLF/EiJ,16.5,0.302,0.073,0.504
Medina-Cano et al. 2025,PWK/PhJ,17.1,0.335,0.085,0.536
Medina-Cano et al. 2025,CAST/EiJ,17.6,0.222,-0.014,0.447
Medina-Cano et al. 2025,SPRET/EiJ,35.2,0.458,0.298,0.586
", stringsAsFactors = FALSE, strip.white = TRUE)
d2$dataset <- factor(d2$dataset, levels = c("Weber et al. 2026", "Medina-Cano et al. 2025"))

ct  <- cor.test(log10(d2$snv), d2$rho, method = "pearson")
ttl <- sprintf("Pearson R = %.2f, p = %.3f", unname(ct$estimate), ct$p.value)

p_divergence <- ggplot(d2, aes(snv, rho)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey45", linetype = "dashed", linewidth = 0.6) +
  geom_linerange(aes(ymin = lo, ymax = hi, color = dataset), linewidth = 0.5, alpha = 0.9) +
  geom_point(aes(color = dataset), size = 2.6) +
  ggtitle(ttl) +
  scale_x_log10(breaks = c(5, 10, 20, 40)) +
  scale_color_manual(values = c("Weber et al. 2026" = "#2C7FB8", "Medina-Cano et al. 2025" = "#D95F02"), name = NULL) +
  labs(x = "Sequence divergence from C57BL/6J\n(approx. SNVs, millions; log scale)",
       y = "Partial Spearman rho\n(Top 500 HVGs)") +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(), axis.line = element_line(color = "black", linewidth = 0.4),
        axis.title = element_text(size = 10), plot.title = element_text(size = 11, hjust = 0.5))

ggsave("ASE_summary_divergence.pdf", p_divergence, width = 5.2, height = 4.2)

######################################################
