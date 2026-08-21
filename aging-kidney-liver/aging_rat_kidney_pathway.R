# Pathway and Marker Analysis of Transcriptional Dyscoordination (Aging Rat Kidney)

# Online links
# https://www.gsea-msigdb.org/gsea/msigdb

# Set up working directory
# setwd("path/to/working/directory")

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(ggpubr)
library(openxlsx)
library(fgsea)
library(KEGGREST)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# 1. Segment markers (Fig 5D) - S1/S2/S3 PT segment genes gaining dyscoordination
kidney_PT_segment_markers <- c("Gpx3", "Miox", "Slc34a1", "Cyp4b1", "Odc1", "Ass1",
                               "Timp3", "Pah", "Akr1a1", "Gpx1", "Pck1", "Igfbp4",
                               "Dnase1", "Tcn2", "Inmt", "Slc13a3", "Cndp2", "Cyp2e1",
                               "Alpl", "Errfi1", "Ndrg1", "Slc5a12", "Gatm", # SpotGLM
                               "Aqp1", "Pdzk1", "Lrp2", "Cubn", "Cyp2e1", "Clcn5",
                               "Gcm1" # GPT
)
segment_assignment <- c(
  "Gpx3" = "S1/S2", "Miox" = "S3", "Slc34a1" = "S1/S2", "Cyp4b1" = "S3", "Odc1" = "Unclear",
  "Ass1" = "S1/S2", "Timp3" = "Unclear", "Pah" = "Unclear", "Akr1a1" = "S1/S2", "Gpx1" = "S1/S2",
  "Pck1" = "S1/S2", "Igfbp4" = "Unclear", "Dnase1" = "S1/S2", "Tcn2" = "Unclear", "Inmt" = "S3",
  "Slc13a3" = "S3", "Cndp2" = "S1/S2", "Cyp2e1" = "S3", "Alpl" = "S1/S2", "Errfi1" = "Unclear",
  "Ndrg1" = "Unclear", "Slc5a12" = "S3", "Gatm" = "S1/S2", "Aqp1" = "S1/S2", "Pdzk1" = "S1/S2",
  "Lrp2" = "S1/S2", "Cubn" = "S1/S2", "Clcn5" = "S1/S2", "Gcm1" = "S3"
)
kidney_PT_segment_markers <- data.frame(
  Gene = kidney_PT_segment_markers,
  Segment = segment_assignment[kidney_PT_segment_markers],
  stringsAsFactors = FALSE
)

kidney_aging_PT_only_gene_level <- read.csv("aging_rat_kidney_all_samples_PT_only_estimated_dispersion_SAVER.csv")
kidney_aging_PT_only_gene_level_filtered <- kidney_aging_PT_only_gene_level[kidney_aging_PT_only_gene_level$Gene %in% kidney_PT_segment_markers$Gene, ]
kidney_aging_PT_only_gene_level_filtered <- kidney_aging_PT_only_gene_level_filtered %>%
  filter(cell_type %in% c("PT-S1", "PT-S2", "PT-S3", "PT-injured1", "PT_injured2")) %>%
  mutate(Age = factor(Age, levels = c("16wk", "30wk", "56wk", "82wk"))) %>%
  arrange(Gene, Age, cell_type)
kidney_aging_PT_only_gene_level_filtered <- na.omit(kidney_aging_PT_only_gene_level_filtered)

# Violin of gene-level dyscoordination by PT subtype
p_segment_marker_dyscoordination_by_cell_type <- kidney_aging_PT_only_gene_level_filtered[kidney_aging_PT_only_gene_level_filtered$cell_type %in% c("PT-S1", "PT-S2", "PT-S3"), ] %>%
  ggplot(aes(x = Age, y = log(Gene_level_deviation), fill = Age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Blues") +
  facet_wrap(~ cell_type, ncol = 3) +
  theme_minimal(base_size = 10) +
  theme(axis.title = element_blank(), legend.position = "none", panel.grid = element_blank(),
        plot.margin = margin(2, 2, 2, 2),
        axis.line.x = element_line(color = "black", linewidth = 0.4),
        axis.line.y = element_line(color = "black", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black", linewidth = 0.3),
        axis.ticks.y = element_line(color = "black", linewidth = 0.3)) +
  ylim(-8, max(log(kidney_aging_PT_only_gene_level_filtered$Gene_level_deviation[kidney_aging_PT_only_gene_level_filtered$cell_type %in% c("PT-S1", "PT-S2", "PT-S3")])) + 2) +
  labs(x = "Age", y = "Gene-level dyscoordination (log-transformed)") +
  geom_signif(data = subset(kidney_aging_PT_only_gene_level_filtered, cell_type == "PT-S1"),
              comparisons = list(c("30wk", "56wk")),
              test = function(x, y) wilcox.test(x, y, alternative = "less"), map_signif_level = TRUE, step_increase = 0.1) +
  geom_signif(data = subset(kidney_aging_PT_only_gene_level_filtered, cell_type == "PT-S2"),
              comparisons = list(c("16wk", "30wk"), c("30wk", "56wk")),
              test = function(x, y) wilcox.test(x, y, alternative = "less"), map_signif_level = TRUE, step_increase = 0.1) +
  geom_signif(data = subset(kidney_aging_PT_only_gene_level_filtered, cell_type == "PT-S3"),
              comparisons = list(c("16wk", "30wk"), c("30wk", "56wk")),
              test = function(x, y) wilcox.test(x, y, alternative = "less"), map_signif_level = TRUE, step_increase = 0.1)

ggsave("Parker_kidney_aging_segment_marker_gene_level_dyscoordination_by_cell_type.pdf",
       plot = p_segment_marker_dyscoordination_by_cell_type, width = 6.25, height = 2.75)

# Line plot: per-gene mean dyscoordination trajectory, highlighting select genes
kidney_aging_PT_only_gene_level_filtered_mean <- kidney_aging_PT_only_gene_level_filtered %>%
  filter(cell_type %in% c("PT-S1", "PT-S2", "PT-S3")) %>%
  mutate(Age_num = as.numeric(str_extract(Age, "\\d+")), log_dev = log(Gene_level_deviation)) %>%
  group_by(Gene, Age, Age_num) %>%
  summarize(mean_log_dev = mean(log_dev, na.rm = TRUE), .groups = "drop")

highlight_genes <- c("Alpl", "Ndrg1", "Igfbp4", "Inmt", "Pah")
mean_labels <- kidney_aging_PT_only_gene_level_filtered_mean %>%
  filter(Gene %in% highlight_genes) %>%
  group_by(Gene) %>% slice_max(order_by = Age_num, n = 1, with_ties = FALSE) %>% ungroup()

p_segment_marker_dyscoordination_lines <- ggplot() +
  geom_line(data = kidney_aging_PT_only_gene_level_filtered_mean %>% filter(!Gene %in% highlight_genes),
            aes(x = Age, y = mean_log_dev, group = Gene), color = "darkgrey", alpha = 0.75, linewidth = 0.5) +
  geom_point(data = kidney_aging_PT_only_gene_level_filtered_mean %>% filter(!Gene %in% highlight_genes),
             aes(x = Age, y = mean_log_dev), size = 0.25, color = "darkgrey") +
  geom_line(data = kidney_aging_PT_only_gene_level_filtered_mean %>% filter(Gene %in% highlight_genes),
            aes(x = Age, y = mean_log_dev, group = Gene), linewidth = 0.5) +
  geom_point(data = kidney_aging_PT_only_gene_level_filtered_mean %>% filter(Gene %in% highlight_genes),
             aes(x = Age, y = mean_log_dev), size = 0.5) +
  geom_text_repel(data = mean_labels %>% filter(Gene %in% highlight_genes),
                  aes(x = Age, y = mean_log_dev, label = Gene),
                  segment.size = 0.3, segment.linetype = "dashed", size = 1.75, nudge_x = 0.25, hjust = 0) +
  theme_minimal(base_size = 8) +
  labs(x = "Age", y = "Mean gene-level dyscoordination (log-transformed)") +
  theme(legend.position = "none", panel.grid.minor = element_blank())

ggsave("Parker_kidney_aging_segment_marker_gene_level_dyscoordination_lineplot.pdf",
       plot = p_segment_marker_dyscoordination_lines, width = 2.5, height = 3)

######################################################
# 2. De-differentiation and inflammation markers (Fig 5E; incl. Cd44, Ccl2)
kidney_marker_genes <- data.frame(
  Gene = c(
    # Dedifferentiation markers
    "Sox9", "Ccng1", "Cdh6", "Vcam1", "Cd24", "Pax2", "Cd44", "Havcr1", "Foxm1",
    # Inflammation markers
    "Il6", "Il1a", "Il8", "Tnfrsf1a", "Tnfrsf1b", "Cxcl1", "Cxcl2", "Ccl2", "Ccl5", "Mif", "Cd74", "Hla-dqb1", "Hla-drb1",
    "Serpine1", "Tgfb1", "Col1a1", "Acta2", "Mmp7", "Stat3", "Rela", "Vegfa",
    # Stress-related markers
    "Qprt", "Krt8", "Krt18", "Vim"
  ),
  Category = c(rep("Dedifferentiation", 9), rep("Inflammation", 21), rep("Stress", 4)))

kidney_aging_PT_only_gene_level <- read.csv("kidney_aging_PT_only_estimated_dispersion_SAVER.csv")
kidney_aging_PT_only_gene_level_filtered <- merge(kidney_aging_PT_only_gene_level, kidney_marker_genes, by.x = "Gene", by.y = "Gene")
kidney_aging_PT_only_gene_level_filtered <- kidney_aging_PT_only_gene_level_filtered %>%
  filter(cell_type %in% c("PT-S1", "PT-S2", "PT-S3", "PT-injured1", "PT_injured2")) %>%
  mutate(Age = factor(Age, levels = c("16wk", "30wk", "56wk", "82wk"))) %>%
  arrange(Category, Gene, Age, cell_type)
kidney_aging_PT_only_gene_level_filtered <- na.omit(kidney_aging_PT_only_gene_level_filtered)

# Gene-level dyscoordination violin by PT subtype
p_dediff_marker_dyscoordination_by_cell_type <- kidney_aging_PT_only_gene_level_filtered[kidney_aging_PT_only_gene_level_filtered$cell_type %in% c("PT-S1", "PT-S2", "PT-S3"), ] %>%
  ggplot(aes(x = Age, y = log(Gene_level_deviation), fill = Age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_brewer(palette = "Blues") +
  facet_wrap(~ cell_type, ncol = 3) +
  theme_minimal(base_size = 10) +
  theme(axis.title = element_blank(), legend.position = "none", panel.grid = element_blank(),
        plot.margin = margin(2, 2, 2, 2),
        axis.line.x = element_line(color = "black", linewidth = 0.4),
        axis.line.y = element_line(color = "black", linewidth = 0.4),
        axis.ticks.x = element_line(color = "black", linewidth = 0.3),
        axis.ticks.y = element_line(color = "black", linewidth = 0.3)) +
  ylim(min(log(kidney_aging_PT_only_gene_level_filtered$Gene_level_deviation[kidney_aging_PT_only_gene_level_filtered$cell_type %in% c("PT-S1", "PT-S2", "PT-S3")])),
       max(log(kidney_aging_PT_only_gene_level_filtered$Gene_level_deviation[kidney_aging_PT_only_gene_level_filtered$cell_type %in% c("PT-S1", "PT-S2", "PT-S3")])) + 2) +
  labs(x = "Age", y = "Gene-level dyscoordination (log-transformed)") +
  geom_signif(data = subset(kidney_aging_PT_only_gene_level_filtered, cell_type == "PT-S2"),
              comparisons = list(c("30wk", "56wk")),
              test = function(x, y) wilcox.test(x, y, alternative = "less"), map_signif_level = TRUE, step_increase = 0.1) +
  geom_signif(data = subset(kidney_aging_PT_only_gene_level_filtered, cell_type == "PT-S3"),
              comparisons = list(c("16wk", "30wk"), c("30wk", "56wk")),
              test = function(x, y) wilcox.test(x, y, alternative = "less"), map_signif_level = TRUE, step_increase = 0.1)

ggsave("Parker_kidney_aging_dedifferentiation_marker_gene_level_dyscoordination_by_cell_type.pdf",
       plot = p_dediff_marker_dyscoordination_by_cell_type, width = 6.25, height = 2.75)

# Volcano of per-gene correlation with cell-level dyscoordination (highlights
# markers such as Cd44 / Ccl2 activated in high-dyscoordination cells).
cor_results <- read.csv("Parker_kidney_aging_PT_only_cell_level_expr_corr_adjusted.csv")
cor_results <- cor_results %>%
  mutate(Sig = case_when(FDR < 0.05 & Correlation > 0 ~ "Up",
                         FDR < 0.05 & Correlation < 0 ~ "Down",
                         TRUE ~ "Not significant"))
selected_genes <- cor_results %>%
  filter(!grepl("^LOC", Gene), -log10(FDR) > 4) %>%
  arrange(Correlation) %>% slice_head(n = 5) %>%
  bind_rows(cor_results %>% filter(!grepl("^LOC", Gene), -log10(FDR) > 4) %>%
              arrange(desc(Correlation)) %>% slice_head(n = 5))

p_cell_level_corr_volcano <- ggplot(cor_results, aes(x = Correlation, y = -log10(FDR), color = Sig)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("Up" = "#009E73", "Down" = "#D55E00", "Not significant" = "gray70")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(x = "Correlation with cell-level dyscoordination", y = "-log10(FDR)", color = NULL) +
  theme_minimal(base_size = 10) +
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 0.5))) +
  geom_point(data = selected_genes, aes(x = Correlation, y = -log10(FDR)), size = 1) +
  geom_text_repel(data = selected_genes, aes(x = Correlation, y = -log10(FDR), label = Gene),
                  size = 3, color = "black", max.overlaps = Inf)

ggsave("kidney_aging_PT_only_cell_level_corr_volcano.png",
       p_cell_level_corr_volcano, width = 6.5, height = 5, dpi = 300)

######################################################
# 3. Gene-set enrichment (Fig 5F)
# MSigDB mouse CP + GO + KEGG gene sets.
msigdb_dir <- 'path/to/MSigDB_mouse'
# kegg_mouse computed once from KEGGREST, then cached:
# kegg_ids  <- keggList("pathway", "mmu")
# kegg_mouse <- lapply(kegg_ids, function(pid) keggGet(pid)[[1]]$GENE)
# names(kegg_mouse) <- kegg_ids
# kegg_mouse <- kegg_mouse[!sapply(kegg_mouse, is.null)]
# saveRDS(kegg_mouse, "kegg.v2025.1.Mm.symbols.RDS")
cp_mouse   <- gmtPathways(file.path(msigdb_dir, "m2.cp.v2025.1.Mm.symbols.gmt"))
go_mouse   <- gmtPathways(file.path(msigdb_dir, "m5.go.v2025.1.Mm.symbols.gmt"))
kegg_mouse <- readRDS("kegg.v2025.1.Mm.symbols.RDS")
combined_pathways <- c(cp_mouse, go_mouse, kegg_mouse)

# 3a. Genes activated in high-dyscoordination cells (cell-level).
# cor_results = per-gene Pearson corr of expression with cell-level dyscoordination
# (depth-adjusted), computed upstream. Rank FDR<0.05 genes by correlation.
cor_results <- read.csv("Parker_kidney_aging_PT_only_cell_level_expr_corr_adjusted.csv")
ranked_genes <- cor_results[cor_results$FDR < 0.05, ]$Correlation
names(ranked_genes) <- cor_results[cor_results$FDR < 0.05, ]$Gene
ranked_genes <- sort(na.omit(ranked_genes), decreasing = TRUE)

kidney_aging_PT_only_fgsea_results <- fgsea(pathways = combined_pathways, stats = ranked_genes)
kidney_aging_PT_only_fgsea_results_sig <- kidney_aging_PT_only_fgsea_results[kidney_aging_PT_only_fgsea_results$pval < 0.05, ]
kidney_aging_PT_only_fgsea_results_sig <- kidney_aging_PT_only_fgsea_results_sig %>%
  mutate(is_prefixed = grepl("^[A-Z]+_", pathway),
         Source = ifelse(is_prefixed, str_extract(pathway, "^[^_]+"), "KEGG"),
         Term = ifelse(is_prefixed, str_remove(pathway, "^[^_]+_"),
                       str_remove(pathway, " - Mus musculus \\(house mouse\\)")),
         Term = str_replace_all(Term, "_", " "),
         Term = str_to_sentence(Term))
kidney_aging_PT_only_fgsea_results_sig <- as.data.frame(kidney_aging_PT_only_fgsea_results_sig)
list_cols <- sapply(kidney_aging_PT_only_fgsea_results_sig, is.list)
kidney_aging_PT_only_fgsea_results_sig[, list_cols] <- lapply(
  kidney_aging_PT_only_fgsea_results_sig[, list_cols, drop = FALSE],
  function(col) sapply(col, function(x) paste(x, collapse = ";")))
write.csv(kidney_aging_PT_only_fgsea_results_sig,
          "Parker_kidney_aging_PT_only_corr_sig_gene_pathway_results_fgsea.csv", row.names = FALSE)

# 3b. Genes gaining gene-level dyscoordination with age.
# Rank each gene by Spearman(ordinal age, Gene_level_deviation); FDR<=0.05; fgsea.
kidney_aging_PT_only_gene_level <- read.csv("kidney_aging_PT_only_estimated_dispersion_SAVER.csv")
kidney_aging_PT_only_gene_level$Gene_level_deviation <- remove_outliers(kidney_aging_PT_only_gene_level$Gene_level_deviation)
kidney_aging_PT_only_gene_level_cleaned <- na.omit(distinct(kidney_aging_PT_only_gene_level))
kidney_aging_PT_only_gene_level_cleaned_summary <- kidney_aging_PT_only_gene_level_cleaned %>%
  group_by(Gene) %>%
  summarise(test = list(cor.test(as.numeric(factor(Age, ordered = TRUE)), Gene_level_deviation,
                                 method = "spearman", exact = FALSE)), .groups = "drop") %>%
  mutate(rho = purrr::map_dbl(test, "estimate"), p = purrr::map_dbl(test, "p.value")) %>%
  select(Gene, rho, p) %>% drop_na()

kidney_aging_PT_only_gene_level_cleaned_summary_filtered <- kidney_aging_PT_only_gene_level_cleaned_summary %>%
  mutate(FDR = p.adjust(p, method = "BH")) %>%
  filter(FDR <= 0.05) %>% arrange(desc(rho))

ranked_genes <- kidney_aging_PT_only_gene_level_cleaned_summary_filtered$rho
names(ranked_genes) <- kidney_aging_PT_only_gene_level_cleaned_summary_filtered$Gene
ranked_genes <- sort(na.omit(ranked_genes), decreasing = TRUE)

kidney_aging_PT_only_gene_level_fgsea_results <- fgsea(pathways = combined_pathways, stats = ranked_genes)
kidney_aging_PT_only_gene_level_fgsea_results_sig <- kidney_aging_PT_only_gene_level_fgsea_results[kidney_aging_PT_only_gene_level_fgsea_results$pval < 0.05, ]
kidney_aging_PT_only_gene_level_fgsea_results_sig <- kidney_aging_PT_only_gene_level_fgsea_results_sig %>%
  mutate(is_prefixed = grepl("^[A-Z]+_", pathway),
         Source = ifelse(is_prefixed, str_extract(pathway, "^[^_]+"), "KEGG"),
         Term = ifelse(is_prefixed, str_remove(pathway, "^[^_]+_"),
                       str_remove(pathway, " - Mus musculus \\(house mouse\\)")),
         Term = str_replace_all(Term, "_", " "),
         Term = str_to_sentence(Term))
kidney_aging_PT_only_gene_level_fgsea_results_sig <- as.data.frame(kidney_aging_PT_only_gene_level_fgsea_results_sig)
list_cols <- sapply(kidney_aging_PT_only_gene_level_fgsea_results_sig, is.list)
kidney_aging_PT_only_gene_level_fgsea_results_sig[, list_cols] <- lapply(
  kidney_aging_PT_only_gene_level_fgsea_results_sig[, list_cols, drop = FALSE],
  function(col) sapply(col, function(x) paste(x, collapse = ";")))
write.csv(kidney_aging_PT_only_gene_level_fgsea_results_sig,
          "Parker_kidney_aging_PT_only_gene_level_corr_sig_gene_pathway_results_fgsea.csv", row.names = FALSE)

######################################################
# 3c. Selected-pathway dot plot (gene-level vs cell-level hypotheses)
kidney_aging_PT_only_fgsea_results_sig <- read.csv("Parker_kidney_aging_PT_only_corr_sig_gene_pathway_results_fgsea.csv")
kidney_aging_PT_only_fgsea_results_sig <- kidney_aging_PT_only_fgsea_results_sig[kidney_aging_PT_only_fgsea_results_sig$ES > 0, ]
kidney_aging_PT_only_gene_level_fgsea_results_sig <- read.csv("Parker_kidney_aging_PT_only_gene_level_corr_sig_gene_pathway_results_fgsea.csv")
kidney_aging_PT_only_gene_level_fgsea_results_sig <- kidney_aging_PT_only_gene_level_fgsea_results_sig[kidney_aging_PT_only_gene_level_fgsea_results_sig$ES > 0, ]

gene_level_pathways <- c(
  "P27 pathway", "P53 signaling pathway", "Caspase pathway",
  "Senescence associated secretory phenotype sasp",
  "Mitochondrial calcium ion homeostasis",
  "Transcriptional activation of mitochondrial biogenesis",
  "Citric acid cycle tca cycle", "Nucleotide excision repair", "Homologous recombination")
cell_level_pathways <- c(
  "Antigen processing and presentation", "Cytokine activity",
  "Regulation of leukocyte chemotaxis", "Extracellular matrix organization")

kidney_aging_PT_only_gene_level_fgsea_results_sig_selected <- kidney_aging_PT_only_gene_level_fgsea_results_sig %>%
  filter(Term %in% gene_level_pathways & ES >= 0.3) %>% mutate(Level = "Gene-level")
kidney_aging_PT_only_fgsea_results_sig_selected <- kidney_aging_PT_only_fgsea_results_sig %>%
  filter(Term %in% cell_level_pathways) %>% mutate(Level = "Cell-level")

kidney_aging_PT_only_fgsea_results_plot <- bind_rows(kidney_aging_PT_only_gene_level_fgsea_results_sig_selected, kidney_aging_PT_only_fgsea_results_sig_selected) %>%
  mutate(minuslog10pval = -log10(pval),
         Level = factor(Level, levels = c("Gene-level", "Cell-level"))) %>%
  group_by(Level) %>% arrange(desc(NES), .by_group = TRUE) %>% ungroup() %>%
  mutate(Term = factor(Term, levels = unique(Term)))

p_selected_pathways <- ggplot(kidney_aging_PT_only_fgsea_results_plot,
       aes(x = NES, y = Term, size = minuslog10pval, color = Level)) +
  geom_point() +
  scale_size_continuous(name = "-log10(p-value)") +
  scale_color_manual(values = c("Gene-level" = "steelblue", "Cell-level" = "#E64B35")) +
  labs(x = "NES", y = NULL) +
  facet_wrap(~ Level, scales = "free_y", ncol = 2) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9), strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("kidney_aging_PT_only_selected_pathways.pdf", plot = p_selected_pathways, width = 12, height = 2.5)

######################################################
