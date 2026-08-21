# Transcriptional Dyscoordination

This is the code repository for the [paper](https://www.biorxiv.org/content/10.64898/2026.01.24.701460v1):

**A Measure of Transcriptional Dyscoordination for Quantifying Aging in Single Cells**

## Standard Workflow

The standard workflow to compute transcriptional dyscoordination includes (1) loading and preparing data (the raw count matrix and desired cell type and age annotation vectors), (2) fitting the manifold of gene expression using SAVER, (3) selecting the best-fitting variance model for each gene, and (4) computing the gene-level and cell-level transcriptional dyscoordination measures across cell types and age groups.

Each analysis folder corresponds to one section of the paper. Datasets with a two-step pipeline use a `*_deviation.R` script (heavy compute: builds the SAVER manifold and writes the gene-level and cell-level dyscoordination CSVs) followed by a `*_analysis.R` script (reads those CSVs and produces the figures). One-off datasets are handled by a single self-contained script. External data are referenced through a `data_dir <- 'path/to/...'` placeholder that the user sets to their own local path; intermediate results (dispersion CSVs, `.RData`) are read by their bare relative filenames.

## Repository Structure

### `Transcriptional_dyscoordination_functions.R`

Helper functions for computing transcriptional dyscoordination measures (variance-model selection, manifold fitting, and dispersion-table cleaning).

### `allelic-specific-expression/` Figure 2A-2C

Allele-specific expression analyses that validate transcriptional dyscoordination against classical intrinsic transcriptional noise. In each F1 hybrid dataset, allelic discordance (from the difference between the two alleles) is compared with gene-level dyscoordination (from their sum) by an expression-controlled partial Spearman correlation.

- `ASE_mESC_Larsson.R` — B6 x CAST mouse embryonic stem cells ([Larsson et al. 2019](https://pubmed.ncbi.nlm.nih.gov/30602787/))
- `ASE_mESC_Ochiai.R` — 129 x CAST mouse embryonic stem cells ([Ochiai et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7299619))
- `ASE_CD8T_vanderVeeken.R` — B6 x SPRET virus-specific CD8 T cells ([van der Veeken et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31201093/))
- `ASE_CD8T_Pritykin.R` — B6 x SPRET CD8 T cells ([Pritykin et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33862018/))
- `ASE_liver_Weber.R` — B6 x seven-strain liver, per cross ([Weber et al. 2026](https://www.biorxiv.org/content/10.1101/2026.04.02.716195))
- `ASE_organoid_MedinaCano.R` — B6 x four-strain neocortical organoids, per cross ([Medina-Cano et al. 2025](https://pubmed.ncbi.nlm.nih.gov/40449493/))
- `ASE_summary_forest.R` — Figure 2C six-dataset forest (all genes vs. top 500 HVGs) and Supplementary Figure 8 (partial Spearman's correlation vs. genome-wide sequence divergence)

### `induced-damage-senescence/` Figure 2D-2E

Transcriptional dyscoordination under controlled perturbations of cellular integrity.

- `ionizing_radiation_ESCC_Wu.R` — dyscoordination in 0 vs. 12 Gy irradiated ESCC cells ([Wu et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6831193/))
- `senescent_mouse_MSC_deviation.R`, `senescent_mouse_MSC_analysis.R` — dyscoordination in mouse muscle satellite cells and fibro-adipogenic progenitors across control, doxorubicin, and senolytic (ABT-263) conditions ([Limbad et al. 2022](https://www.cell.com/iscience/fulltext/S2589-0042(22)00118-3))

### `tabula-muris-senis/` Figure 3A-3E

- `TMS_marrow_deviation.R`, `TMS_marrow_analysis.R` — dyscoordination across bone marrow cell types with age ([Tabula Muris Senis](https://www.nature.com/articles/s41586-020-2496-1)). Directly applicable to the other TMS tissues (lung, limb muscle, adipose).
- `TMS_marrow_benchmarks.R` — comparison scores: SenMayo senescence score ([Saul et al. 2022](https://www.nature.com/articles/s41467-022-32552-1)), SCENT signaling entropy, and CytoTRACE potency, each rank-correlated against cell-level dyscoordination. The distance-based (decibel) and Scallop membership benchmarks are computed in `decibel_scallop_analysis.ipynb`.

### `aging-human-BM/` Figure 3F

- `aging_human_BM_deviation.R`, `aging_human_BM_analysis.R` — dyscoordination in human hematopoietic stem and progenitor cells across age ([Oetjen et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6237339/))

### `disease-aging-T-cell/` Figure 4

Lineage-resolved dyscoordination in four human T cell datasets with paired TCR clonotypes.

- `glioma_T_cell_deviation.R`, `glioma_T_cell_analysis.R` — glioma-infiltrating CD8 T cells ([Wang et al. 2025](https://www.biorxiv.org/content/10.1101/2025.09.06.674490v1))
- `melanoma_T_cell_deviation.R` — ICB-treated melanoma CD8 T cells, per timepoint ([Wang et al. 2024](https://pubmed.ncbi.nlm.nih.gov/39214097/))
- `Terekhova_T_cell_deviation.R` — healthy-aging blood T cells ([Terekhova et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37963457/))
- `Wang_aging_T_cell_deviation.R` — healthy-aging PBMC T cells ([Wang et al. 2025](https://pubmed.ncbi.nlm.nih.gov/39881000/))
- `T_cell_common.R` — shared helpers (depth correction, within-donor partial Spearman, clone aggregation, plotting)
- `T_cell_clone_size.R` — clone size vs. clone-level dyscoordination (Figure 4A)
- `T_cell_subset_groups.R` — dyscoordination across coarse T cell subset groups (Figure 4B)
- `T_cell_gene_programs.R` — dyscoordination vs. six T cell gene-program scores (Figure 4C)
- `T_cell_tumor_grade.R` — glioma patient-level dyscoordination vs. tumor grade (Figure 4D)
- `T_cell_longitudinal_ICB.R` — within-clone dyscoordination change after immune checkpoint blockade (Figure 4E)

### `aging-rat-kidney/` Figure 5

Aging rat kidney proximal tubule Multiome analysis.

- `aging_rat_kidney_deviation.R`, `aging_rat_kidney_analysis.R` — proximal tubule dyscoordination by age and segment
- `aging_rat_kidney_epitrace.R` — EpiTrace mitotic age vs. dyscoordination (Figure 5C)
- `aging_rat_kidney_pathway.R` — segment markers, inflammation markers, and pathway enrichment (Figure 5D, 5E, 5F)

### `aging-human-kidney/` Figure 6

- `aging_human_kidney_deviation.R`, `aging_human_kidney_analysis.R` — dyscoordination in human tubular epithelial cells (PT, TAL, DCT) by age (Figure 6B, 6C)
- `aging_human_kidney_pathway.R` — pathway enrichment for genes gaining dyscoordination and genes activated in high-dyscoordination cells (Figure 6E)

### `aging-mouse-liver/` Figure 7

Aging mouse liver hepatocyte Multiome analysis.

- `aging_mouse_liver_deviation.R`, `aging_mouse_liver_analysis.R` — hepatocyte dyscoordination by age
- `aging_mouse_liver_epitrace.R` — EpiTrace mitotic age vs. dyscoordination (Figure 7B)
- `aging_mouse_liver_zonation.R` — zonation markers and per-zone dyscoordination (Figure 7C-D)
- `aging_mouse_liver_pathway.R` — per-zone pathway enrichment (Figure 7E-F)

### `in-silico-spike-in/` Supplementary Note 1

- `spike_in_simulation.R` — spike-in simulation on real substrates (rat kidney, mouse liver, mouse bone marrow): coordinated program, mean-preserving noise, and both, with prevalence and gene-count sweeps
- `spike_in_analysis.R` — figures of DE-gene count vs. dyscoordination ratio and the learnability sweeps

### `decibel_scallop_analysis.ipynb`

Python notebook computing the distance-based (decibel) and Scallop membership benchmark scores.
