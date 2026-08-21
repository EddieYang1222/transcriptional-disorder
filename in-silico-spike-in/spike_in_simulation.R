# In-silico spike-in simulation of transcriptional dyscoordination (compute)
# Builds a baseline latent expression on three real substrates, simulates four
# arms (control, coordinated program, mean-preserving gene-wise noise, both),
# regenerates counts at the original library size, computes transcriptional
# dyscoordination per arm, and evaluates the dyscoordination ratio + DE-gene
# count per replicate. Also runs the prevalence and gene-count learnability sweeps.

# Online links
# https://www.biorxiv.org/content/10.64898/2026.01.24.701460v1

# Set up working directory
# setwd("path/to/working/directory")

library(Matrix)
library(SAVER)
library(doParallel)
library(pbapply)

# Load helper functions (get_var_model is used to score dyscoordination)
source('Transcriptional_dyscoordination_functions.R')

######################################################
# External substrate data. Each entry points to a raw-count source object.
# rat kidney PT (Parker/ParkerWilson DEFND), mouse liver hepatocytes (integrated
# Seurat), mouse bone-marrow HPCs (Tabula Muris Senis).
data_dir <- 'path/to/spike_in_substrates'
PT_RDATA      <- file.path(data_dir, 'kidney_aging_PT_only.RData')          # dataset.counts/age/celltype (+.levels)
META_TXT      <- file.path(data_dir, 'ParkerWilson_DEFND_metadata.txt')     # rat_id/library_id per barcode
LIVER_RDATA   <- file.path(data_dir, 'liver_data_integrated.RData')         # Seurat v5 object liver_data_integrated
TMS_RDATA     <- file.path(data_dir, 'TMS_marrow.RData')                    # dataset.counts/age/celltype (+.levels)

######################################################
# Design constants
SIM_SEED <- 20260729L

# Published QC thresholds
QC_MIN_GENES_PER_CELL <- 250L
QC_MIN_UMI_PER_CELL   <- 2500L
QC_MIN_CELLS_PER_GENE <- 5L

# Baseline (variant B: kNN-smoothed observed proportions)
SMOOTH_K   <- 20L    # neighbours for the kNN smoother
SMOOTH_PCS <- 30L    # PCs the neighbourhood is defined in
PSEUDO     <- 5      # pseudo-count shrinking toward the segment mean (keeps rates > 0)
N_HVG      <- 2000L
N_GENES_TARGET  <- 5000L
N_CELLS_TARGET  <- 8000L
MIN_DETECT_FRAC <- 0.10

# Spike-in perturbation design
BETA_MAG     <- 1.0
Z_RANGE      <- c(0.7, 1.3)
TAU          <- 0.8
GROUP_LEVELS <- c("control", "coord", "dyscoord", "both")

# Number of replicates: 10 for the primary design, 5 per sweep level (matches the
# reported runs).
N_REP_PRIMARY <- 10L
N_REP_SWEEP   <- 5L

# Label-permutation control resolution
N_PERM <- 1000L

# Gene classes excluded from the spike-in panel
MITO_PATTERN <- "^([Mm][Tt]-|MT-)"
RIBO_PATTERN <- "^([Rr][Pp][SsLl]|RP[SL])"
CC_GENES <- c(
  "Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7",
  "Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn",
  "Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45",
  "Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b",
  "Mrpl36","E2f8",
  "Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2",
  "Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb",
  "Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Cdc20",
  "Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2",
  "Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2",
  "G2e3","Gas2l3","Cbx5","Cenpa")

######################################################
# Helpers
# Timestamped log to stdout
say <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              paste0(..., collapse = "")))
  flush.console()
}

# Fit the variance model and report the label breakdown. Thin wrapper over
# get_var_model(); genes with no valid model carry name "" and drop out of delta.
var_model <- function(counts, mu, size.factor) {
  vm <- get_var_model(as.matrix(counts), mu, size.factor)
  lab <- names(vm)
  lab[is.na(lab) | lab == ""] <- "<none>"
  say("variance model: ", paste(sprintf("%s=%d", names(table(lab)),
      as.integer(table(lab))), collapse = " "))
  vm
}

# Cell- and gene-level dyscoordination. Vectorized delta = (est - mu)^2 / nu,
# nu per gene from the chosen variance model (cCV / cFF / cVar).
dyscoord <- function(est, mu, vm) {
  model <- names(vm)
  model[is.na(model)] <- ""
  nu <- matrix(NA_real_, nrow(est), ncol(est))
  i1 <- model == "cCV"
  i2 <- model == "cFF"
  i3 <- model == "cVar"
  if (any(i1)) nu[i1, ] <- (mu[i1, , drop = FALSE]^2) / vm[i1]
  if (any(i2)) nu[i2, ] <-  mu[i2, , drop = FALSE]    / vm[i2]
  if (any(i3)) nu[i3, ] <- vm[i3]
  delta <- (est - mu)^2 / nu
  delta[!is.finite(delta)] <- NA_real_
  list(delta   = delta,
       cell    = colMeans(delta, na.rm = TRUE),
       gene    = rowMeans(delta, na.rm = TRUE),
       n_genes = sum(rowSums(is.finite(delta)) > 0))
}

# Pull raw counts out of a Seurat object across the v4/v5 API split
seurat_counts <- function(obj, assay = "RNA") {
  m <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
                error = function(e) NULL)
  if (is.null(m) || !ncol(m))
    m <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, slot = "counts"),
                  error = function(e) NULL)
  if (is.null(m) || !ncol(m))
    m <- tryCatch(SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
                  error = function(e) NULL)
  if (is.null(m) || !ncol(m)) stop("could not extract a counts matrix")
  m
}

# Load exactly one named object out of an .RData without polluting globalenv
load_one <- function(path, object) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  get(object, envir = e)
}

# Decode integer-coded annotation to character labels
as_labels <- function(x, levels = NULL) {
  if (is.factor(x)) return(as.character(x))
  if (is.character(x)) return(x)
  levels[x]
}

# Recover rat_id/library_id by joining barcodes to the DEFND metadata
load_cell_metadata <- function(barcodes) {
  meta <- read.table(META_TXT, header = TRUE, sep = " ",
                     stringsAsFactors = FALSE, check.names = FALSE)
  keep <- intersect(c("rat_id", "library_id", "library_type"), colnames(meta))
  out <- meta[barcodes, keep, drop = FALSE]
  rownames(out) <- barcodes
  out
}

######################################################
# Substrate loaders. Each returns list(counts, cellmeta) with the five canonical
# columns rat_id / library_id / library_type / segment / age. "donor" is rats
# (kidney), donors (liver) or barcode channels (TMS).

# Rat kidney PT: young (16wk + 30wk) x segments PT-S1/S2/S3
load_kidney <- function() {
  e <- load_one_env(PT_RDATA)
  sel <- e$dataset.age %in% 1:2 & e$dataset.celltype %in% 1:3
  counts <- e$dataset.counts[, sel, drop = FALSE]
  cm <- load_cell_metadata(colnames(counts))
  cm$library_type <- "kidney10x"
  cm$segment <- e$dataset.celltype.levels[e$dataset.celltype[sel]]
  cm$age     <- e$dataset.age.levels[e$dataset.age[sel]]
  rownames(cm) <- colnames(counts)
  list(counts = counts, cellmeta = cm[, c("rat_id","library_id","library_type","segment","age")])
}

# Mouse liver: young hepatocyte singlets; donor = sample group
load_liver <- function() {
  obj <- load_one(LIVER_RDATA, "liver_data_integrated")
  md  <- obj@meta.data
  keep <- md$doublet == "Singlet" & md$annotation == "Hepatocyte" & md$age == "young"
  counts <- seurat_counts(obj)[, keep, drop = FALSE]
  md <- md[keep, , drop = FALSE]
  rm(obj)
  gc()
  donor <- as.character(md$group)
  cm <- data.frame(rat_id = donor, library_id = donor, library_type = "liver10x",
                   segment = "Hepatocyte", age = as.character(md$age),
                   stringsAsFactors = FALSE)
  rownames(cm) <- colnames(counts)
  list(counts = counts, cellmeta = cm)
}

# TMS marrow: hematopoietic precursor cells across all ages; donor = channel
load_tms <- function() {
  e <- load_one_env(TMS_RDATA)
  ct  <- as_labels(e$dataset.celltype, e$dataset.celltype.levels)
  sel <- ct == "hematopoietic precursor cell"
  counts <- as(as.matrix(e$dataset.counts[, sel, drop = FALSE]), "CsparseMatrix")
  counts <- counts[Matrix::rowSums(counts) > 0, , drop = FALSE]
  age <- as_labels(e$dataset.age, e$dataset.age.levels)[sel]
  age[age == "1m"] <- "01m"
  age[age == "3m"] <- "03m"
  bc <- colnames(counts)
  channel <- ifelse(grepl("^10X_", bc),
                    sub("^(10X_[^_]+_[^_]+)_.*$", "\\1", bc),
                    sub("^[ACGTN]+-1", "", bc))
  cm <- data.frame(rat_id = channel, library_id = channel, library_type = "TMS10x",
                   segment = "HPC", age = age, stringsAsFactors = FALSE)
  rownames(cm) <- bc
  list(counts = counts, cellmeta = cm)
}

# Load a whole .RData into a private env (dataset.counts style objects)
load_one_env <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  e
}

# Substrate registry: variant label (also the metrics-file token) and cell target.
# Sweeps are run on the rat kidney substrate only (variant B).
SUBSTRATES <- list(
  kidney = list(variant = "B",      load = load_kidney, n_cells_target = N_CELLS_TARGET, sweeps = TRUE),
  liver  = list(variant = "BLIVER", load = load_liver,  n_cells_target = N_CELLS_TARGET, sweeps = FALSE),
  tms    = list(variant = "BTMS",   load = load_tms,    n_cells_target = Inf,            sweeps = FALSE))

######################################################
# 1. Build the baseline substrate (published QC -> gene panel -> stratified
#    subsample -> variant-B latent expression). Mirrors 01/01b of the pipeline.

# Exact truncated PCA using base R only
pca_scores <- function(logmat, n_pc) {
  X <- t(as.matrix(logmat))
  X <- scale(X, center = TRUE, scale = TRUE)
  X[!is.finite(X)] <- 0
  cv <- crossprod(X) / (nrow(X) - 1)
  ev <- eigen(cv, symmetric = TRUE)
  X %*% ev$vectors[, seq_len(n_pc), drop = FALSE]
}

# Variant B baseline: kNN-smoothed observed proportions, independent of SAVER
build_baseline_B <- function(counts, lib_size, segment) {
  prop <- Matrix::t(Matrix::t(counts) / lib_size)         # genes x cells rates
  logn <- log1p(prop * 1e4)
  v    <- apply(as.matrix(logn), 1, var)
  hvg  <- names(sort(v, decreasing = TRUE))[seq_len(min(N_HVG, nrow(logn)))]
  sc   <- pca_scores(logn[hvg, , drop = FALSE], SMOOTH_PCS)

  n  <- nrow(sc)
  nn <- matrix(NA_integer_, n, SMOOTH_K)
  sq <- rowSums(sc^2)
  chunk <- 1000L
  for (start in seq(1L, n, by = chunk)) {
    ii <- start:min(start + chunk - 1L, n)
    d2 <- outer(sq[ii], sq, "+") - 2 * (sc[ii, , drop = FALSE] %*% t(sc))
    for (r in seq_along(ii)) {
      d2[r, ii[r]] <- Inf
      nn[ii[r], ] <- order(d2[r, ])[seq_len(SMOOTH_K)]
    }
  }
  propm <- as.matrix(prop)
  mu <- propm
  for (j in seq_len(n))
    mu[, j] <- (propm[, j] + rowSums(propm[, nn[j, ], drop = FALSE])) / (SMOOTH_K + 1)

  seg_mean <- sapply(split(seq_len(n), segment), function(ii)
    Matrix::rowSums(counts[, ii, drop = FALSE]) / sum(lib_size[ii]))
  if (is.null(dim(seg_mean)))                              # single segment level
    seg_mean <- matrix(seg_mean, ncol = 1, dimnames = list(NULL, unique(segment)))
  mu <- mu + PSEUDO * seg_mean[, segment] / mean(lib_size)
  mu[mu <= 0] <- .Machine$double.eps
  mu
}

build_baseline <- function(sub) {
  loaded <- sub$load()
  counts <- loaded$counts
  cellmeta <- loaded$cellmeta

  # Published QC: cells first, then genes
  cell_pass <- Matrix::colSums(counts > 0) >= QC_MIN_GENES_PER_CELL &
               Matrix::colSums(counts)     >= QC_MIN_UMI_PER_CELL
  counts <- counts[, cell_pass, drop = FALSE]
  cellmeta <- cellmeta[cell_pass, , drop = FALSE]
  counts <- counts[Matrix::rowSums(counts > 0) >= QC_MIN_CELLS_PER_GENE, , drop = FALSE]
  say("post-QC pool: ", ncol(counts), " cells x ", nrow(counts), " genes")

  # Gene panel: detected in >=10% of cells, non-mito/ribo/cell-cycle, below the
  # 99th expression percentile; keep the N_GENES_TARGET most-detected
  genes <- rownames(counts)
  detect_frac <- Matrix::rowSums(counts > 0) / ncol(counts)
  gene_mean   <- Matrix::rowMeans(counts)
  excluded <- grepl(MITO_PATTERN, genes) | grepl(RIBO_PATTERN, genes) | genes %in% CC_GENES
  eligible <- detect_frac >= MIN_DETECT_FRAC & !excluded
  eligible <- eligible & gene_mean < quantile(gene_mean[eligible], 0.99)
  panel <- names(sort(detect_frac[eligible], decreasing = TRUE))[
    seq_len(min(N_GENES_TARGET, sum(eligible)))]
  counts <- counts[panel, , drop = FALSE]
  say("gene panel: ", length(panel), " genes")

  # Stratified subsample (segment x age x library); drop = TRUE collapses constant
  # factors so the stratum reduces to donor/channel automatically
  set.seed(SIM_SEED)
  strata <- interaction(cellmeta$segment, cellmeta$age, cellmeta$library_id, drop = TRUE)
  n_target <- sub$n_cells_target
  if (ncol(counts) > n_target) {
    frac <- n_target / ncol(counts)
    idx <- unlist(lapply(split(seq_len(ncol(counts)), strata), function(ii) {
      take <- max(1L, round(length(ii) * frac))
      if (take >= length(ii)) ii else sample(ii, take)
    }), use.names = FALSE)
    idx <- sort(idx)
  } else {
    idx <- seq_len(ncol(counts))
  }
  counts <- counts[, idx, drop = FALSE]
  cellmeta <- cellmeta[idx, , drop = FALSE]
  lib_size <- Matrix::colSums(counts)
  say("final substrate: ", nrow(counts), " genes x ", ncol(counts), " cells")

  mu <- build_baseline_B(counts, lib_size, cellmeta$segment)
  list(counts = counts, mu = mu, lib_size = lib_size, cellmeta = cellmeta,
       genes = rownames(counts))
}

######################################################
# 2. Simulate one replicate.
# Arms:
#   control   z_i = 0, eps = 0
#   coord     z_i ~ U(0.7,1.3), eps = 0          coordinated rank-one program
#   dyscoord  z_i = 0, eps ~ N(-tau^2/2, tau^2)  mean-preserving gene-wise noise
#   both      both
# log mu*_ig = log mu_ig + z_i beta_g + eps_ig; Y*_i ~ Multinomial(L_i, p*_i) at
# the cell's ORIGINAL library size. The -tau^2/2 centering keeps E[exp(eps)] = 1,
# so the dyscoord arm perturbs covariance while leaving mean expression intact.

# Largest-remainder split of randomized indices into arms of the requested fractions
allocate_groups <- function(ii, fracs) {
  n <- length(ii)
  fracs <- fracs / sum(fracs)
  exact <- fracs * n
  base <- floor(exact)
  short <- as.integer(round(n - sum(base)))
  if (short > 0) {
    ord <- order(exact - base, runif(length(fracs)), decreasing = TRUE)
    base[ord[seq_len(short)]] <- base[ord[seq_len(short)]] + 1L
  }
  rep(names(fracs), times = base)
}

simulate_replicate <- function(base, variant, tag, rep_id, cfg) {
  # Deterministic seed as a function of (variant, tag, rep)
  seed <- as.integer(
    (sum(utf8ToInt(paste0(variant, tag))) * 1000L + rep_id * 7919L) %% .Machine$integer.max)
  set.seed(seed)

  mu <- base$mu
  lib_size <- base$lib_size
  cellmeta <- base$cellmeta
  genes <- rownames(mu)
  G <- nrow(mu)
  C <- ncol(mu)

  # Balanced group assignment within segment x age x library
  fracs <- c(control  = 1 - (cfg$coord_frac + cfg$dyscoord_frac + cfg$both_frac),
             coord    = cfg$coord_frac, dyscoord = cfg$dyscoord_frac, both = cfg$both_frac)
  if (any(fracs < 0)) stop("group fractions exceed 1")
  strata <- interaction(cellmeta$segment, cellmeta$age, cellmeta$library_id, drop = TRUE)
  group  <- character(C)
  for (lv in levels(strata)) {
    ii <- sample(which(strata == lv))
    group[ii] <- allocate_groups(ii, fracs)
  }
  group <- factor(group, levels = GROUP_LEVELS)

  # Spike-in design: same gene set for coordinated and dyscoordinated perturbations
  n_spike <- cfg$n_spike_up + cfg$n_spike_down
  spike    <- sample(genes, n_spike)
  spike_up <- spike[seq_len(cfg$n_spike_up)]
  spike_dn <- spike[cfg$n_spike_up + seq_len(cfg$n_spike_down)]
  beta <- setNames(numeric(G), genes)
  beta[spike_up] <-  BETA_MAG
  beta[spike_dn] <- -BETA_MAG

  # Latent perturbation: log mu* = log mu + z_i beta_g + eps_ig
  z <- numeric(C)
  coord_cells <- group %in% c("coord", "both")
  z[coord_cells] <- runif(sum(coord_cells), Z_RANGE[1], Z_RANGE[2])
  eps_cells <- group %in% c("dyscoord", "both")

  mu_star <- mu
  if (any(z != 0)) mu_star <- mu_star * exp(outer(beta, z))
  if (any(eps_cells)) {
    ii <- which(eps_cells)
    e  <- matrix(rnorm(n_spike * length(ii), mean = -TAU^2 / 2, sd = TAU),
                 nrow = n_spike, ncol = length(ii))
    mu_star[spike, ii] <- mu_star[spike, ii] * exp(e)
  }

  # Regenerate counts multinomially at the ORIGINAL library size
  p_star <- sweep(mu_star, 2, colSums(mu_star), "/")
  Y <- matrix(0L, G, C, dimnames = dimnames(mu))
  for (j in seq_len(C)) Y[, j] <- rmultinom(1, size = lib_size[j], prob = p_star[, j])
  Y <- Matrix(Y, sparse = TRUE)

  list(counts = Y, lib_size = lib_size, group = group, z = z, beta = beta,
       spike_up = spike_up, spike_dn = spike_dn, cellmeta = cellmeta,
       tag = tag, rep_id = rep_id, seed = seed, cfg = cfg)
}

######################################################
# 3. Fit the framework on one replicate (ONE global SAVER fit, blind to arm),
#    select the variance model, and score dyscoordination + program projections.

fit_and_score <- function(sim, ncores = 1L) {
  counts <- sim$counts
  group <- sim$group
  beta <- sim$beta
  z <- sim$z
  cm <- sim$cellmeta
  ncount <- Matrix::colSums(counts)
  sf <- ncount / mean(ncount)

  if (ncores > 1) registerDoParallel(cores = ncores)
  y   <- as.matrix(counts)
  fit <- saver(y, size.factor = sf, ncores = ncores)
  dimnames(fit$mu.out) <- dimnames(fit$estimate) <- dimnames(y)

  vm <- var_model(y, fit$mu.out, sf)
  sc <- dyscoord(fit$estimate, fit$mu.out, vm)

  # Program projections: P = predictable, R = residual component of the beta
  # program; gene-centred and unit-normalized by ||beta||
  bnorm <- sqrt(sum(beta^2))
  centre_rows <- function(M) M - rowMeans(M)
  log_mu  <- centre_rows(log(pmax(fit$mu.out,   .Machine$double.eps)))
  log_est <- centre_rows(log(pmax(fit$estimate, .Machine$double.eps)))
  P <- as.vector(crossprod(log_mu, beta)) / bnorm
  R <- as.vector(crossprod(log_est - log_mu, beta)) / bnorm

  scores <- data.frame(
    barcode = colnames(y), group = as.character(group), z = z,
    rat_id = cm$rat_id, library_id = cm$library_id, segment = cm$segment, age = cm$age,
    nCount = ncount, dyscoord = sc$cell, P = P, R = R, stringsAsFactors = FALSE)
  list(scores = scores, n_genes = sc$n_genes)
}

######################################################
# 4. Evaluate one replicate into the reviewer-facing metrics: the rat-level
#    dyscoordination ratio (arm vs control) and the DE-gene count of the spike-in
#    program. Inference unit is the rat, not the cell.

# Rat-level paired comparison of one arm against control
arm_vs_control <- function(sc, arm) {
  rat_arm <- tapply(sc$dyscoord, list(sc$rat_id, sc$group), mean, na.rm = TRUE)
  if (!all(c(arm, "control") %in% colnames(rat_arm))) return(NULL)
  d <- rat_arm[, arm] - rat_arm[, "control"]
  d <- d[is.finite(d)]
  if (length(d) < 3L) return(NULL)
  tt <- t.test(d)
  list(ratio = mean(rat_arm[, arm], na.rm = TRUE) / mean(rat_arm[, "control"], na.rm = TRUE),
       diff = unname(tt$estimate), p_two = tt$p.value,
       p_up = t.test(d, alternative = "greater")$p.value, n_rats = length(d))
}

# Within-rat label-permutation control; the null mean ratio should sit at ~1.0
perm_test <- function(sc, arm, n_perm = 1000L) {
  keep <- sc$group %in% c(arm, "control")
  s <- sc[keep, , drop = FALSE]
  if (!nrow(s)) return(NULL)
  obs_stat <- function(g) {
    ra <- tapply(s$dyscoord, list(s$rat_id, g), mean, na.rm = TRUE)
    if (!all(c(arm, "control") %in% colnames(ra))) return(NA_real_)
    mean(ra[, arm], na.rm = TRUE) / mean(ra[, "control"], na.rm = TRUE)
  }
  obs <- obs_stat(s$group)
  null <- replicate(n_perm, {
    g <- s$group
    for (r in unique(s$rat_id)) { ii <- which(s$rat_id == r)
    g[ii] <- sample(g[ii]) }
    obs_stat(g)
  })
  null <- null[is.finite(null)]
  list(p_perm = (1 + sum(null >= obs)) / (1 + length(null)),
       null_mean_ratio = mean(null))
}

# Differential expression of the spike-in program, arm vs control. Crude on
# purpose: depth-normalized log2 fold change + per-gene Welch t, same yardstick
# for both arms. Reports the DE-gene count and the mean-preserving rate ratio.
de_by_arm <- function(sim, arm) {
  grp <- sim$group
  ii <- which(grp == arm)
  jj <- which(grp == "control")
  if (length(ii) < 10L || length(jj) < 10L) return(NULL)
  prop <- Matrix::t(Matrix::t(sim$counts) / Matrix::colSums(sim$counts)) * 1e4
  prop_sub <- function(idx) as.matrix(prop[, idx, drop = FALSE])
  dens <- function(idx) log2(prop_sub(idx) + 1)
  c1 <- dens(ii)
  c0 <- dens(jj)
  m1 <- rowMeans(c1)
  m0 <- rowMeans(c0)
  lfc <- m1 - m0
  v1 <- apply(c1, 1, var)
  v0 <- apply(c0, 1, var)
  se <- sqrt(v1 / length(ii) + v0 / length(jj))
  tstat <- ifelse(se > 0, lfc / se, 0)
  df <- pmax(1, (v1/length(ii) + v0/length(jj))^2 /
              ((v1/length(ii))^2/(length(ii)-1) + (v0/length(jj))^2/(length(jj)-1)))
  fdr <- p.adjust(2 * pt(-abs(tstat), df), "BH")
  # Arithmetic-scale ratio: the mean-preservation claim is made on the rate scale
  ratio_up <- function(gs) mean(rowMeans(prop_sub(ii)[gs, , drop = FALSE])) /
                           mean(rowMeans(prop_sub(jj)[gs, , drop = FALSE]))
  up <- names(sim$beta)[sim$beta > 0]
  dn <- names(sim$beta)[sim$beta < 0]
  nl <- names(sim$beta)[sim$beta == 0]
  list(lfc_up = mean(lfc[up]), lfc_dn = mean(lfc[dn]), lfc_null = mean(lfc[nl]),
       rate_ratio_up = ratio_up(up), rate_ratio_dn = ratio_up(dn), rate_ratio_null = ratio_up(nl),
       n_de = sum(fdr < 0.05), n_de_spike = sum(fdr[c(up, dn)] < 0.05),
       n_de_null = sum(fdr[nl] < 0.05))
}

# One replicate -> one metrics row per non-control arm
evaluate_replicate <- function(sim, scores, rep_id) {
  rows <- list()
  for (arm in setdiff(GROUP_LEVELS, "control")) {
    a <- arm_vs_control(scores, arm)
    if (is.null(a)) next
    d  <- de_by_arm(sim, arm)
    pt <- perm_test(scores, arm, n_perm = N_PERM)
    ii <- scores$group == arm
    rows[[length(rows) + 1L]] <- data.frame(
      rep = rep_id, arm = arm, n_rats = a$n_rats, n_cells = sum(ii), tau = TAU,
      dys_ratio = a$ratio, dys_diff = a$diff, p_two = a$p_two, p_up = a$p_up,
      p_perm          = if (is.null(pt)) NA_real_ else pt$p_perm,
      perm_null_ratio = if (is.null(pt)) NA_real_ else pt$null_mean_ratio,
      cor_P_z = if (sd(scores$z[ii]) > 0) cor(scores$P[ii], scores$z[ii]) else NA_real_,
      cor_R_z = if (sd(scores$z[ii]) > 0) cor(scores$R[ii], scores$z[ii]) else NA_real_,
      lfc_up   = if (is.null(d)) NA_real_ else d$lfc_up,
      lfc_dn   = if (is.null(d)) NA_real_ else d$lfc_dn,
      lfc_null = if (is.null(d)) NA_real_ else d$lfc_null,
      rate_ratio_up   = if (is.null(d)) NA_real_ else d$rate_ratio_up,
      rate_ratio_dn   = if (is.null(d)) NA_real_ else d$rate_ratio_dn,
      rate_ratio_null = if (is.null(d)) NA_real_ else d$rate_ratio_null,
      n_de       = if (is.null(d)) NA_integer_ else d$n_de,
      n_de_spike = if (is.null(d)) NA_integer_ else d$n_de_spike,
      n_de_null  = if (is.null(d)) NA_integer_ else d$n_de_null,
      depth_cor = cor(scores$dyscoord[ii], scores$nCount[ii]),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

######################################################
# 5. Design table. Each tag sets per-arm prevalence and spike-in panel size.
#   primary        : all three arms at 20% prevalence, 100 up / 50 down
#   prevalence sweep: coordinated-arm prevalence varies (both arm off)
#   gene-count sweep: spike-in panel size varies (prev20 = 150-gene midpoint)
# NOTE: the sweep arms use both_frac = 0, so each sweep replicate carries the
# control / coord / dyscoord arms only (the "both" arm is the primary design).
mk_cfg <- function(tag, coord_frac, dyscoord_frac, both_frac, n_up, n_dn)
  list(tag = tag, coord_frac = coord_frac, dyscoord_frac = dyscoord_frac,
       both_frac = both_frac, n_spike_up = n_up, n_spike_down = n_dn)

primary_cfg <- mk_cfg("primary", 0.20, 0.20, 0.20, 100, 50)

# Prevalence sweep (coordinated-arm fraction on the x-axis)
prevalence_cfgs <- list(
  mk_cfg("prev05", 0.05, 0.20, 0, 100, 50),
  mk_cfg("prev10", 0.10, 0.20, 0, 100, 50),
  mk_cfg("prev20", 0.20, 0.20, 0, 100, 50),   # also the gene-sweep 150-gene midpoint
  mk_cfg("prev30", 0.30, 0.20, 0, 100, 50),
  mk_cfg("prev50", 0.50, 0.20, 0, 100, 50))

# Gene-count sweep (spike-in panel size on the x-axis; prev20 = 100 up / 50 down)
genecount_cfgs <- list(
  mk_cfg("ng20u10d",   0.20, 0.20, 0,  20,  10),
  mk_cfg("ng50u25d",   0.20, 0.20, 0,  50,  25),
  mk_cfg("ng200u100d", 0.20, 0.20, 0, 200, 100),
  mk_cfg("ng400u200d", 0.20, 0.20, 0, 400, 200))

######################################################
# 6. Driver. Build each substrate's baseline once, then run every tag x replicate
#    and write results/replicate_metrics_<variant>_<tag>.csv.
# Heavy compute (SAVER per replicate): run one substrate at a time on a cluster.
dir.create("results", showWarnings = FALSE)
ncores <- as.integer(Sys.getenv("NSLOTS", "1"))
if (is.na(ncores) || ncores < 1) ncores <- 1L

run_config <- function(base, variant, cfg, n_rep) {
  met <- do.call(rbind, lapply(seq_len(n_rep), function(rp) {
    say("  ", variant, "/", cfg$tag, " rep ", rp)
    sim    <- simulate_replicate(base, variant, cfg$tag, rp, cfg)
    scored <- fit_and_score(sim, ncores = ncores)
    evaluate_replicate(sim, scored$scores, rp)
  }))
  path <- file.path("results", sprintf("replicate_metrics_%s_%s.csv", variant, cfg$tag))
  write.csv(met, path, row.names = FALSE)
  say("  wrote ", path)
}

for (key in names(SUBSTRATES)) {
  sub <- SUBSTRATES[[key]]
  say("=== substrate ", key, " (variant ", sub$variant, ") ===")
  base <- build_baseline(sub)

  # Primary four-arm design (all substrates)
  run_config(base, sub$variant, primary_cfg, N_REP_PRIMARY)

  # Learnability sweeps (rat kidney only)
  if (isTRUE(sub$sweeps)) {
    for (cfg in c(prevalence_cfgs, genecount_cfgs))
      run_config(base, sub$variant, cfg, N_REP_SWEEP)
  }
}

say("done")

######################################################
