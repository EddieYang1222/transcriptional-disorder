# Transcriptional Dyscoordination

This is the code repository for the [paper](https://www.biorxiv.org/content/10.64898/2026.01.24.701460v1):

**A Manifold-Based Measure of Transcriptional Entropy for Quantifying Aging in Single Cells**

## Standard Workflow

Standard workflow to compute transcriptional dyscoordination includes (1) loading and preparing data (the raw count matrix and desired cell type and age annotation vectors), (2) fitting the manifold of gene expression using SAVER, (3) selecting the best-fitting variance model for each gene, and (4) computing the gene-level and cell-level transcriptional dyscoordination measures across cell types and age groups.

## Repository Structure

### `Transcriptional_dyscoordination_functions.R`

Helper functions for computing transcriptional dyscoordination measures.

### `decibel_scallop_analysis.ipynb`

Workflow to compute mean Euclidean distance to cell type average, mean Euclidean distance to tissue average using invariant genes, and Scallop membership score from [Ibañez-Solé et al. (2022)](https://elifesciences.org/articles/80380) for (1) TMS marrow cells, (2) liver hepatocytes, and (3) kidney PT cells

### `allelic-specific-expression/`

Scripts for two allelic-specific expression (ASE) analyses:

- `ASE_mESC_Larsson.R` — ASE analysis on mESCs from [Larsson et al. (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7610481)
- `ASE_mESC_Ochiai.R` — ASE analysis on mESCs from [Ochiai et al. (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7299619)

### `ionizing-radiation/`

Scripts for analyzing ionizing radiated cells:

- `ionizing_radiation_ESCC_Wu.R` — Transcriptional dyscoordination by radiation dosage on ESCCs from [Wu et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6831193)

### `tabula-muris-senis/`

Scripts for analyzing cells across different tissues from Tabula Muris Senis:

- `TMS_marrow_deviation.R` — Transcriptional dyscoordination on bone marrow cells from [Tabula Muris Senis](https://www.nature.com/articles/s41586-020-2496-1). This script is for bone marrow but directly applicable to all other TMS tissues including lung, limb muscle, and adipose tissue.

- `TMS_marrow_benchmarks.R` — Comparisons between transcriptional dyscoordination and related aging scores: 
(1) SenMayo senescence score from [Saul et al. (2022)](https://www.nature.com/articles/s41467-022-32552-1)
(2) mean Euclidean distance to cell type average from [Ibañez-Solé et al. (2022)](https://elifesciences.org/articles/80380)
(3) mean Euclidean distance to tissue average using invariant genes from [Ibañez-Solé et al. (2022)](https://elifesciences.org/articles/80380)
(4) Scallop membership score from [Ibañez-Solé et al. (2022)](https://elifesciences.org/articles/80380)
