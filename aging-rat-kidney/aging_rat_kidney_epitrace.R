# EpiTrace Mitotic Age vs Transcriptional Dyscoordination (Aging Rat Kidney)

# Online links
# https://github.com/MagpiePKU/EpiTrace

# Set up working directory
# setwd("path/to/working/directory")

library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(ggplot2)
library(Matrix)
library(EpiTrace)
library(Seurat)
library(SeuratObject)
library(ggsignif)
library(ggpubr)
library(ggExtra)
library(patchwork)
library(Signac)
library(rtracklayer)
library(easyLift)

# Load helper functions
source('Transcriptional_dyscoordination_functions.R')

######################################################
# Clock liftover: human -> mouse (mm10) -> rat (rn7)
epitrace_dir <- 'path/to/EpiTrace'
mouse_clock_by_human <- readRDS(file.path(epitrace_dir, 'mouse_clock_lifted_human_to_mouse_mm10.rds'))
rat_clock_by_human   <- easyLiftOver(mouse_clock_by_human,
                                     map = file.path(epitrace_dir, 'mm10ToRn7.over.chain'))

######################################################
# Sequence name remapping: NCBI accession -> chr notation (rat GRCr8)
sequence_report_GRCr8 <- read.delim("sequence_report.tsv", header = TRUE, sep = "\t",
                                    stringsAsFactors = FALSE)
sequence_report_GRCr8 <- sequence_report_GRCr8[, c("RefSeq.seq.accession", "Sequence.name")]
colnames(sequence_report_GRCr8) <- c("ncbi", "chr")
sequence_report_GRCr8 <- sequence_report_GRCr8[!is.na(sequence_report_GRCr8$chr), ]

update_seqnames <- function(peaks_df, mapping_df) {
  peaks_df$chr <- ifelse(peaks_df$chr %in% mapping_df$ncbi,
                         mapping_df$chr[match(peaks_df$chr, mapping_df$ncbi)],
                         peaks_df$chr)
  return(peaks_df)
}

######################################################
# Sample manifest - all 22 ATAC libraries matching step2_anno RNA libraries
#   sample_key (list name) is unique per ATAC library and prefixes cell_remapped.
#   Multiome-ATAC anchors (C03/L12/Q17/U21) carry RNA labels via barcode
#   conversion; all DEFND-ATAC libraries get labels via NN inference.
data_dir <- 'path/to/RAGE24_cellranger'

sample_manifest <- list(
  # Multiome-ATAC anchor samples (RNA label transfer via barcode conversion)
  C03   = list(folder = "RAGE24-C03-KYC-LN-01", atac_id = "RAGE24-C03-KYC-LN-01-Multiome-ATAC", rat_id = "C03", age_wks = "16wk", is_anchor = TRUE),
  L12   = list(folder = "RAGE24-L12-KYC-LN-01", atac_id = "RAGE24-L12-KYC-LN-01-Multiome-ATAC", rat_id = "L12", age_wks = "30wk", is_anchor = TRUE),
  Q17   = list(folder = "RAGE24-Q17-KYC-LN-01", atac_id = "RAGE24-Q17-KYC-LN-01-Multiome-ATAC", rat_id = "Q17", age_wks = "56wk", is_anchor = TRUE),
  U21   = list(folder = "RAGE24-U21-KYC-LN-01", atac_id = "RAGE24-U21-KYC-LN-01-Multiome-ATAC", rat_id = "U21", age_wks = "82wk", is_anchor = TRUE),
  # DEFND-ATAC for anchor rats (separate aliquot, NN label inference)
  C03_D = list(folder = "RAGE24-C03-KYC-LN-01", atac_id = "RAGE24-C03-KYC-LN-01-DEFND-ATAC", rat_id = "C03", age_wks = "16wk", is_anchor = FALSE),
  L12_D = list(folder = "RAGE24-L12-KYC-LN-01", atac_id = "RAGE24-L12-KYC-LN-01-DEFND-ATAC", rat_id = "L12", age_wks = "30wk", is_anchor = FALSE),
  Q17_D = list(folder = "RAGE24-Q17-KYC-LN-01", atac_id = "RAGE24-Q17-KYC-LN-01-DEFND-ATAC", rat_id = "Q17", age_wks = "56wk", is_anchor = FALSE),
  U21_D = list(folder = "RAGE24-U21-KYC-LN-01", atac_id = "RAGE24-U21-KYC-LN-01-DEFND-ATAC", rat_id = "U21", age_wks = "82wk", is_anchor = FALSE),
  # DEFND-ATAC only rats (16 wk)
  B02   = list(folder = "RAGE24-B02-KYC-LN-01", atac_id = "RAGE24-B02-KYC-LN-01-DEFND-ATAC", rat_id = "B02", age_wks = "16wk", is_anchor = FALSE),
  D04   = list(folder = "RAGE24-D04-KYC-LN-01", atac_id = "RAGE24-D04-KYC-LN-01-DEFND-ATAC", rat_id = "D04", age_wks = "16wk", is_anchor = FALSE),
  E05   = list(folder = "RAGE24-E05-KYC-LN-01", atac_id = "RAGE24-E05-KYC-LN-01-DEFND-ATAC", rat_id = "E05", age_wks = "16wk", is_anchor = FALSE),
  F06   = list(folder = "RAGE24-F06-KYC-LN-01", atac_id = "RAGE24-F06-KYC-LN-01-DEFND-ATAC", rat_id = "F06", age_wks = "16wk", is_anchor = FALSE),
  # DEFND-ATAC only rats (30 wk)
  G07   = list(folder = "RAGE24-G07-KYC-LN-01", atac_id = "RAGE24-G07-KYC-LN-01-DEFND-ATAC", rat_id = "G07", age_wks = "30wk", is_anchor = FALSE),
  H08   = list(folder = "RAGE24-H08-KYC-LN-01", atac_id = "RAGE24-H08-KYC-LN-01-DEFND-ATAC", rat_id = "H08", age_wks = "30wk", is_anchor = FALSE),
  I09   = list(folder = "RAGE24-I09-KYC-LN-01", atac_id = "RAGE24-I09-KYC-LN-01-DEFND-ATAC", rat_id = "I09", age_wks = "30wk", is_anchor = FALSE),
  # J10: no ATAC data in cellranger folder - excluded
  # DEFND-ATAC only rats (56 wk)
  N14   = list(folder = "RAGE24-N14-KYC-LN-01", atac_id = "RAGE24-N14-KYC-LN-01-DEFND-ATAC", rat_id = "N14", age_wks = "56wk", is_anchor = FALSE),
  O15   = list(folder = "RAGE24-O15-KYC-LN-02", atac_id = "RAGE24-O15-KYC-LN-02-DEFND-ATAC", rat_id = "O15", age_wks = "56wk", is_anchor = FALSE),
  P16   = list(folder = "RAGE24-P16-KYC-LN-01", atac_id = "RAGE24-P16-KYC-LN-01-DEFND-ATAC", rat_id = "P16", age_wks = "56wk", is_anchor = FALSE),
  # DEFND-ATAC only rats (82 wk)
  T20   = list(folder = "RAGE24-T20-KYC-LN-02", atac_id = "RAGE24-T20-KYC-LN-02-DEFND-ATAC", rat_id = "T20", age_wks = "82wk", is_anchor = FALSE),
  V22   = list(folder = "RAGE24-V22-KYC-LN-01", atac_id = "RAGE24-V22-KYC-LN-01-DEFND-ATAC", rat_id = "V22", age_wks = "82wk", is_anchor = FALSE),
  W23   = list(folder = "RAGE24-W23-KYC-LN-01", atac_id = "RAGE24-W23-KYC-LN-01-DEFND-ATAC", rat_id = "W23", age_wks = "82wk", is_anchor = FALSE)
)

######################################################
# Load peaks / barcodes / count matrices for all samples
read_atac_sample <- function(s) {
  base  <- file.path(data_dir, s$folder, s$atac_id, "outs", "filtered_peak_bc_matrix")
  peaks <- read_tsv(file.path(base, "peaks.bed"),    col_names = c("chr", "start", "end"), show_col_types = FALSE)
  cells <- read_tsv(file.path(base, "barcodes.tsv"), col_names = c("cell"), show_col_types = FALSE)
  mm    <- readMM(file.path(base, "matrix.mtx"))
  list(peaks = peaks, cells = cells, mm = mm)
}

message("[1/10] Loading ATAC data for ", length(sample_manifest), " libraries...")
atac_data <- lapply(names(sample_manifest), function(nm) {
  message("  Reading: ", nm)
  read_atac_sample(sample_manifest[[nm]])
})
names(atac_data) <- names(sample_manifest)

# Remap NCBI seqnames -> chr notation
message("[2/10] Remapping NCBI sequence names to chr notation...")
for (nm in names(atac_data)) {
  atac_data[[nm]]$peaks <- update_seqnames(atac_data[[nm]]$peaks, sequence_report_GRCr8)
}

######################################################
# Build merged peak set across all samples
message("[3/10] Building merged peak set...")
gr_list <- lapply(atac_data, function(d) {
  GRanges(seqnames = d$peaks$chr,
          ranges   = IRanges(start = d$peaks$start, end = d$peaks$end))
})
merged_peaks <- suppressWarnings(reduce(do.call(c, unname(gr_list))))

######################################################
# Align each sample's count matrix to the merged peak set
align_matrix_to_peaks <- function(mat, orig_peaks, merged_peaks) {
  overlap_idx <- findOverlaps(orig_peaks, merged_peaks, select = "first")
  aligned_mat <- Matrix(0, nrow = length(merged_peaks), ncol = ncol(mat), sparse = TRUE)
  rownames(aligned_mat) <- paste(seqnames(merged_peaks), start(merged_peaks), end(merged_peaks), sep = "-")
  colnames(aligned_mat) <- colnames(mat)
  valid_idx <- which(!is.na(overlap_idx))
  aligned_mat[overlap_idx[valid_idx], ] <- mat[valid_idx, ]
  return(aligned_mat)
}

message("[4/10] Aligning count matrices to merged peak set...")
realigned_list <- lapply(names(atac_data), function(nm) {
  d   <- atac_data[[nm]]
  mat <- d$mm
  colnames(mat) <- d$cells$cell
  align_matrix_to_peaks(mat, gr_list[[nm]], merged_peaks)
})
names(realigned_list) <- names(atac_data)

merged_matrix_realigned <- do.call(cbind, realigned_list)

# Flat cell table preserving cbind order (sample_key unique per ATAC library)
cells_flat <- do.call(rbind, lapply(names(sample_manifest), function(nm) {
  data.frame(cell       = atac_data[[nm]]$cells$cell,
             sample_key = nm,
             rat_id     = sample_manifest[[nm]]$rat_id,
             age_wks    = sample_manifest[[nm]]$age_wks,
             is_anchor  = sample_manifest[[nm]]$is_anchor,
             stringsAsFactors = FALSE)
}))
cells_flat$cell_remapped <- paste(cells_flat$sample_key, cells_flat$cell, sep = ":")
colnames(merged_matrix_realigned) <- cells_flat$cell_remapped

merged_metadata <- data.frame(
  cell       = cells_flat$cell,
  sample_key = cells_flat$sample_key,
  rat_id     = cells_flat$rat_id,
  age_wks    = cells_flat$age_wks,
  is_anchor  = cells_flat$is_anchor,
  stringsAsFactors = FALSE
)

######################################################
# Cell type label transfer via barcode conversion (anchor samples only)
message("[5/10] Loading barcode conversion table and RNA anchor metadata...")
barcode_conversion_df <- read.csv("barcode_conversion_10XMultiome.csv", stringsAsFactors = FALSE)

# RNA anchor object (step2_anno) - the 4 Multiome anchor rats
kidney_anchor <- readRDS("step2_anno.rds")
kidney_anchor <- subset(kidney_anchor, subset = rat_id %in% c("C03", "L12", "Q17", "U21"))

ordered_celltypes   <- c("PT-S1","PT-S2","PT-S3","PT-MT","PT-injured1","PT-injured2",
                         "PEC","M-TAL","C-TAL","DCT","CNT","CNT-PC","PC","IC-B",
                         "EC-PTC","EC-GC","EC-AEA_DVR","EC-AVR","EC-LYM","FIB","VSMC/P",
                         "Lym","Mono","POD")
ordered_celltypes_2 <- c("PT","PT","PT","PT","PT","PT",
                         "PEC","M-TAL","C-TAL","DCT","CNT","CNT-PC","PC","IC-B",
                         "EC","EC","EC","EC","EC","FIB","VSMC/P",
                         "Lym","Mono","POD")
full_names <- c(
  "proximal tubule epithelial cell (segment 1)",
  "proximal tubule epithelial cell (segment 2)",
  "proximal tubule epithelial cell (segment 3)",
  "proximal tubule epithelial cell (mixed type)",
  "proximal tubule epithelial cell (injured 1)",
  "proximal tubule epithelial cell (injured 2)",
  "parietal epithelial cell","medullary thick ascending limb cell",
  "cortical thick ascending limb cell","distal convoluted tubule cell",
  "connecting tubule cell","connecting tubule principal cell",
  "principal cell","intercalated cell type B",
  "peritubular capillary endothelial cell",
  "glomerular capillary endothelial cell",
  "afferent/efferent arteriole & descending vasa recta endothelial cell",
  "ascending vasa recta endothelial cell","lymphatic endothelial cell",
  "fibroblast","vascular smooth muscle cell/pericyte",
  "lymphocyte","monocyte","podocyte"
)
full_names_2 <- c(
  rep("proximal tubule epithelial cell", 6),
  "parietal epithelial cell","medullary thick ascending limb cell",
  "cortical thick ascending limb cell","distal convoluted tubule cell",
  "connecting tubule cell","connecting tubule principal cell",
  "principal cell","intercalated cell type B",
  rep("endothelial cell", 5),
  "fibroblast","vascular smooth muscle cell/pericyte",
  "lymphocyte","monocyte","podocyte"
)

kidney_anchor$age_wks <- paste0(kidney_anchor$age_wks, "wk")
kidney_anchor$celltype_refined2_regrouped <-
  ordered_celltypes_2[match(as.character(kidney_anchor$celltype_refined2), ordered_celltypes)]
kidney_anchor$celltype_refined2_full <-
  full_names[match(as.character(kidney_anchor$celltype_refined2), ordered_celltypes)]
kidney_anchor$celltype_refined2_regrouped_full <-
  full_names_2[match(as.character(kidney_anchor$celltype_refined2), ordered_celltypes)]

rna_anchor_metadata <- kidney_anchor@meta.data[,
  c("rat_id","age_wks","celltype_refined2","celltype_refined2_full",
    "celltype_refined2_regrouped","celltype_refined2_regrouped_full")]
rna_anchor_metadata$cell_name <- sapply(strsplit(rownames(rna_anchor_metadata), "-"), `[`, 1)

barcode_conversion_df <- merge(barcode_conversion_df, rna_anchor_metadata,
                               by.x = "GEX_bc", by.y = "cell_name")
barcode_conversion_df$ATAC_bc <- paste0(barcode_conversion_df$ATAC_bc, "-1")

# Merge ATAC metadata with barcode conversion (DEFND cells get NA)
merged_metadata <- merge(merged_metadata,
                         barcode_conversion_df[, -1],   # drop GEX_bc
                         by.x  = "cell",
                         by.y  = "ATAC_bc",
                         all.x = TRUE)

# Only assign RNA labels for Multiome anchor cells whose rat_id agrees
merged_metadata <- merged_metadata %>%
  mutate(
    valid_label = is_anchor & !is.na(rat_id.y) & (rat_id.x == rat_id.y),
    celltype_refined2                = if_else(valid_label, celltype_refined2,                NA_character_),
    celltype_refined2_full           = if_else(valid_label, celltype_refined2_full,           NA_character_),
    celltype_refined2_regrouped      = if_else(valid_label, celltype_refined2_regrouped,      NA_character_),
    celltype_refined2_regrouped_full = if_else(valid_label, celltype_refined2_regrouped_full, NA_character_)
  ) %>%
  select(-valid_label)

merged_metadata$rat_id  <- merged_metadata$rat_id.x
merged_metadata$age_wks <- merged_metadata$age_wks.x
merged_metadata$rat_id.x    <- NULL
merged_metadata$age_wks.x   <- NULL
merged_metadata$rat_id.y    <- NULL
merged_metadata$age_wks.y   <- NULL
merged_metadata$cell_remapped <- paste(merged_metadata$sample_key, merged_metadata$cell, sep = ":")

######################################################
# Filter to standard chromosomes (rat rn7: chr1-20, X, Y) + deduplicate metadata
message("[6/10] Filtering to standard chromosomes and deduplicating metadata...")
standard_chr    <- paste0("chr", c(1:20, "X", "Y"))
merged_peaks    <- keepSeqlevels(merged_peaks, standard_chr, pruning.mode = "coarse")
valid_peaknames <- paste(seqnames(merged_peaks), start(merged_peaks), end(merged_peaks), sep = "-")

merged_matrix_realigned <- merged_matrix_realigned[rownames(merged_matrix_realigned) %in% valid_peaknames, ]
merged_matrix_realigned <- merged_matrix_realigned[, colnames(merged_matrix_realigned) %in% merged_metadata$cell_remapped]

filtered_metadata <- data.frame()
for (cell in unique(merged_metadata$cell_remapped)) {
  cell_rows <- merged_metadata[merged_metadata$cell_remapped == cell, ]
  if (nrow(cell_rows) == 1) { filtered_metadata <- rbind(filtered_metadata, cell_rows)
  next }
  if (all(is.na(cell_rows$celltype_refined2_regrouped))) { filtered_metadata <- rbind(filtered_metadata, cell_rows[1, ])
  next }
  unique_ct  <- unique(na.omit(cell_rows$celltype_refined2_regrouped))
  unique_ctf <- unique(na.omit(cell_rows$celltype_refined2_regrouped_full))
  if (length(unique_ct) == 1 && length(unique_ctf) == 1) { filtered_metadata <- rbind(filtered_metadata, cell_rows[1, ])
  next }
  # Conflicting labels -> NA out
  cell_rows$celltype_refined2                <- NA
  cell_rows$celltype_refined2_full           <- NA
  cell_rows$celltype_refined2_regrouped      <- NA
  cell_rows$celltype_refined2_regrouped_full <- NA
  filtered_metadata <- rbind(filtered_metadata, cell_rows[1, ])
}
filtered_metadata <- filtered_metadata[order(match(filtered_metadata$cell_remapped, colnames(merged_matrix_realigned))), ]
rownames(filtered_metadata) <- filtered_metadata$cell_remapped

# save(merged_matrix_realigned, filtered_metadata, merged_peaks,
#      file = "Parker_kidney_aging_all_samples_ATAC_for_EpiTrace.RData")
# load("Parker_kidney_aging_all_samples_ATAC_for_EpiTrace.RData")

######################################################
# Seurat ATAC object -> LSI -> NN label transfer
message("[7/10] Creating Seurat ATAC object and running LSI + NN graph...")
cell_type_annotation <- read.csv("cell_type_annotation.csv", stringsAsFactors = FALSE)

kidney_all_atac_obj <- CreateSeuratObject(
  counts    = merged_matrix_realigned,
  assay     = "ATAC",
  meta.data = filtered_metadata
)
kidney_all_atac_obj <- RunTFIDF(kidney_all_atac_obj)
kidney_all_atac_obj <- FindTopFeatures(kidney_all_atac_obj, min.cutoff = 10)
kidney_all_atac_obj <- RunSVD(kidney_all_atac_obj)
kidney_all_atac_obj <- FindNeighbors(kidney_all_atac_obj, k.param = 31,
                                     reduction = "lsi", dims = 1:30, return.neighbor = TRUE)

kidney_all_atac_obj$annotation_type <- ifelse(
  !is.na(kidney_all_atac_obj$celltype_refined2_regrouped), "original", "inferred")

labelled_cells   <- rownames(kidney_all_atac_obj@meta.data)[!is.na(kidney_all_atac_obj$celltype_refined2_regrouped)]
unlabelled_cells <- rownames(kidney_all_atac_obj@meta.data)[is.na(kidney_all_atac_obj$celltype_refined2_regrouped)]

message("[8/10] Running vectorized NN label inference...")
# Vectorized majority-vote NN label inference
nn_idx <- kidney_all_atac_obj@neighbors$ATAC.nn@nn.idx
all_cell_names  <- rownames(kidney_all_atac_obj@meta.data)
unlabelled_rows <- match(unlabelled_cells, all_cell_names)
all_labels <- setNames(kidney_all_atac_obj$celltype_refined2, seq_along(all_cell_names))

inferred_labels <- apply(nn_idx[unlabelled_rows, -1, drop = FALSE], 1, function(nbr_rows) {
  lbls <- all_labels[as.character(nbr_rows)]
  lbls <- lbls[!is.na(lbls)]
  if (length(lbls) == 0) return(NA_character_)
  names(sort(table(lbls), decreasing = TRUE))[1]
})
names(inferred_labels) <- unlabelled_cells

kidney_all_atac_obj@meta.data[unlabelled_cells, "celltype_refined2"] <- inferred_labels

valid_mask  <- !is.na(inferred_labels)
valid_cells <- unlabelled_cells[valid_mask]
if (length(valid_cells) > 0) {
  ann_matched <- cell_type_annotation[match(inferred_labels[valid_mask], cell_type_annotation$celltype), ]
  kidney_all_atac_obj@meta.data[valid_cells, "celltype_refined2_full"]           <- ann_matched$full_name
  kidney_all_atac_obj@meta.data[valid_cells, "celltype_refined2_regrouped"]      <- ann_matched$celltype_regrouped
  kidney_all_atac_obj@meta.data[valid_cells, "celltype_refined2_regrouped_full"] <- ann_matched$full_name_regrouped
}

# save(kidney_all_atac_obj, file = "Parker_kidney_aging_all_samples_ATAC_for_EpiTrace_obj.RData")

######################################################
# EpiTrace on PT cells only (iterative_time = 1)
message("[9/10] Running EpiTrace on PT cells...")
pt_mask <- !is.na(kidney_all_atac_obj$celltype_refined2_regrouped) &
           kidney_all_atac_obj$celltype_refined2_regrouped == "PT"
filtered_matrix_PT_only <- merged_matrix_realigned[, pt_mask]

init_gr <- Init_Peakset(merged_peaks)
init_mm <- Init_Matrix(
  peakname = paste(seqnames(merged_peaks), start(merged_peaks), end(merged_peaks), sep = "-"),
  cellname = colnames(filtered_matrix_PT_only),
  matrix   = filtered_matrix_PT_only
)

epitrace_all_samples_PT <- EpiTraceAge_Convergence(
  peakSet            = init_gr,
  matrix             = init_mm,
  ref_genome         = 'rn7',
  clock_gr           = rat_clock_by_human,
  iterative_time     = 1,
  min.cutoff         = 0,
  non_standard_clock = FALSE,
  qualnum            = 10,
  ncore_lim          = 1,
  mean_error_limit   = 0.1
)

# Attach full metadata
epitrace_all_samples_PT@meta.data <- merge(
  epitrace_all_samples_PT@meta.data,
  kidney_all_atac_obj@meta.data,
  by.x  = "cell",
  by.y  = "cell_remapped",
  all.x = TRUE
)
epitrace_all_samples_PT@meta.data$orig.iden    <- epitrace_all_samples_PT@meta.data$orig.iden.x
epitrace_all_samples_PT@meta.data$cell.y       <- NULL
epitrace_all_samples_PT@meta.data$orig.ident.x <- NULL
epitrace_all_samples_PT@meta.data$orig.iden.y  <- NULL

# save(epitrace_all_samples_PT, file = "Parker_kidney_aging_all_samples_EpiTrace_human_clock_PT_only.RData")

######################################################
# EpiTrace age distribution across PT subtypes
message("[10/10] Generating plots...")
age_levels <- c("16wk", "30wk", "56wk", "82wk")

plot_df <- epitrace_all_samples_PT@meta.data %>%
  filter(!is.na(celltype_refined2), celltype_refined2 %in% c("PT-S1", "PT-S2", "PT-S3")) %>%
  mutate(Age = factor(age_wks, levels = age_levels))

p_all_epitrace_PT_by_subtype <- plot_df %>%
  ggplot(aes(x = Age, y = EpiTraceAge_iterative, fill = Age)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", fill = "white", linewidth = 0.4) +
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  facet_wrap(~ celltype_refined2, ncol = 3) +
  theme_minimal() +
  labs(x = "Age", y = "EpiTrace mitotic age score") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_signif(
    comparisons      = list(c("16wk", "30wk"), c("30wk", "56wk"), c("56wk", "82wk")),
    test             = function(x, y) wilcox.test(x, y, alternative = "less"),
    map_signif_level = TRUE, step_increase = 0.1) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("Parker_kidney_aging_all_samples_EpiTrace_PT_by_subtype.pdf",
       p_all_epitrace_PT_by_subtype, width = 9, height = 4)

######################################################
# EpiTrace age vs cell-level dyscoordination (Fig 5C)
# Original 4 Multiome samples: PT-only dispersion CSV -> barcode conversion
# (GEX -> ATAC) -> EpiTrace metadata.
load("Parker_kidney_aging_EpiTrace_human_clock_PT_only.RData")

kidney_aging_PT_only_cell_level <- read.csv("kidney_aging_PT_only_cellular_dispersion_SAVER.csv",
                                            stringsAsFactors = FALSE)
kidney_aging_PT_only_cell_level$Cell_barcode <-
  sapply(strsplit(kidney_aging_PT_only_cell_level$Cell_barcode, "-"), `[`, 1)

kidney_aging_PT_only_cell_level <- merge(
  kidney_aging_PT_only_cell_level, barcode_conversion_df,
  by.x = c("Cell_barcode", "Age", "cell_type"),
  by.y = c("GEX_bc",       "age_wks", "celltype_refined2"),
  all.x = TRUE)

kidney_aging_PT_only_cell_level_cleaned <- na.omit(distinct(kidney_aging_PT_only_cell_level))

# Remove ambiguous duplicates
kidney_aging_PT_only_cell_level_cleaned <- kidney_aging_PT_only_cell_level_cleaned[
  !(kidney_aging_PT_only_cell_level_cleaned$Cell_barcode %in%
      kidney_aging_PT_only_cell_level_cleaned$Cell_barcode[
        duplicated(kidney_aging_PT_only_cell_level_cleaned$Cell_barcode)]), ]

kidney_aging_PT_only_cell_level_cleaned$ATAC_bc <-
  paste(kidney_aging_PT_only_cell_level_cleaned$rat_id,
        kidney_aging_PT_only_cell_level_cleaned$ATAC_bc, sep = ":")

kidney_aging_PT_only_cell_level_cleaned <- merge(
  kidney_aging_PT_only_cell_level_cleaned,
  epitrace_obj_age_conv_human_clock_PT_only@meta.data[,
    c("cell", "EpiTraceAge_Clock_initial", "EpiTraceAge_iterative", "celltype_refined2")],
  by.x = "ATAC_bc", by.y = "cell")

kidney_aging_PT_only_cell_level_cleaned$Age <-
  factor(kidney_aging_PT_only_cell_level_cleaned$Age, levels = age_levels)

# Scatter: iterative EpiTrace vs cell-level dyscoordination, faceted by PT subtype
p_epitrace_iter_by_cell_type_PT_only <- ggplot(
  kidney_aging_PT_only_cell_level_cleaned[
    kidney_aging_PT_only_cell_level_cleaned$Cell_level_deviation < 100 &
    kidney_aging_PT_only_cell_level_cleaned$cell_type %in% c("PT-S1", "PT-S2", "PT-S3"), ],
  aes(x = EpiTraceAge_iterative, y = log(Cell_level_deviation), color = Age)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE, linewidth = 0.75) +
  facet_wrap(~ cell_type, ncol = 3) +
  scale_color_viridis_d(option = "plasma", direction = -1) +
  labs(x = "EpiTrace age", y = "Cell-level dyscoordination (log-transformed)", color = "Age") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    plot.margin  = margin(2, 2, 2, 2),
    axis.line.x  = element_line(color = "black", linewidth = 0.4),
    axis.line.y  = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black", linewidth = 0.3),
    axis.ticks.y = element_line(color = "black", linewidth = 0.3)
  )

# Scatter: iterative EpiTrace vs cell-level dyscoordination, all PT, marginal density
p_epitrace_iter_PT_only <- ggplot(
  kidney_aging_PT_only_cell_level_cleaned[
    kidney_aging_PT_only_cell_level_cleaned$Cell_level_deviation < 100, ],
  aes(x = EpiTraceAge_iterative, y = log(Cell_level_deviation), color = Age)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE, linewidth = 0.75) +
  scale_color_viridis_d(option = "plasma", direction = -1) +
  labs(x = "EpiTrace age", y = "Cell-level dyscoordination (log-transformed)", color = "Age") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    plot.margin  = margin(2, 2, 2, 2),
    axis.line.x  = element_line(color = "black", linewidth = 0.4),
    axis.line.y  = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black", linewidth = 0.3),
    axis.ticks.y = element_line(color = "black", linewidth = 0.3)
  )
p_epitrace_iter_PT_only <- ggMarginal(p_epitrace_iter_PT_only, type = "density",
                                      groupColour = TRUE, groupFill = TRUE)

ggsave("Parker_kidney_aging_EpiTrace_iterative_cell_level_dyscoordination_by_cell_type_PT_only.pdf",
       plot = p_epitrace_iter_by_cell_type_PT_only, width = 5, height = 3)
ggsave("Parker_kidney_aging_EpiTrace_iterative_cell_level_dyscoordination_PT_only.pdf",
       plot = p_epitrace_iter_PT_only, width = 3.5, height = 2.5)

######################################################
