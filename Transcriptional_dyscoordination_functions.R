# Get the variance model that maximizes the likelihood
get_var_model <- function(y_mat, mu_mat, sf) {
  
  # Define a helper function that processes each single gene
  get_var_model_gene <- function(g) {
    y <- y_mat[g, ]
    mu <- mu_mat[g, ]
    
    if (var(mu) == 0) {
      a <- c(0, Inf)
      b <- c(0, Inf)
      k <- c(0, Inf)
    } else{
      # Calculate parameters for the three models
      a <- tryCatch(calc.a(y, mu, sf), error = function(cond) return(c(0, Inf)))
      b <- tryCatch(calc.b(y, mu, sf), error = function(cond) return(c(0, Inf)))
      k <- tryCatch(calc.k(y, mu, sf), error = function(cond) return(c(0, Inf)))
    }
    
    # Compare log-likelihoods and choose the model with the maximum likelihood
    if (min(a[2], b[2], k[2]) == Inf) {
      return(NA) # No valid model
    }
      
    var.method <- which.min(c(a[2], b[2], k[2]))
    
    # Return the parameter value with name as the model
    if (var.method == 1) {
      result <- a[1]
      names(result) <- "cCV"
      return(result)
    } else if (var.method == 2) {
      result <- b[1]
      names(result) <- "cFF"
      return(result)
    } else {
      result <- k[1]
      names(result) <- "cVar"
      return(result)
    }
  }
  
  best_models <- pblapply(1:nrow(y_mat), FUN = get_var_model_gene)
  best_models <- unlist(best_models)
  
  return(best_models)
}

# Get the variance model that maximizes the likelihood
get_var_model_old <- function(y_mat, mu_mat, sf) {
  best_models <- rep(NA, nrow(y_mat))
  
  for (g in 1:nrow(y_mat)) {
    y <- y_mat[g, ]
    mu <- mu_mat[g, ]
    
    # Calculate parameters for the three models
    a <- tryCatch(calc.a(y, mu, sf), error = function(cond) return(c(0, Inf)))
    b <- tryCatch(calc.b(y, mu, sf), error = function(cond) return(c(0, Inf)))
    k <- tryCatch(calc.k(y, mu, sf), error = function(cond) return(c(0, Inf)))
    
    # Compare log-likelihoods and choose the model with the maximum likelihood
    if (min(a[2], b[2], k[2]) == Inf) {
      var.method = 0
    } else {
      var.method <- which.min(c(a[2], b[2], k[2]))
    }
    
    # Store the model name based on the selected index
    if (var.method == 0) {
      best_models[g] <- NA
      message(g, ": ", "No model")
    } else if (var.method == 1) {
      best_models[g] <- "cCV"
      message(g, ": ", "cCV")
    } else if (var.method == 2) {
      best_models[g] <- "cFF"
      message(g, ": ", "cFF")
    } else {
      best_models[g] <- "cVar"
      message(g, ": ", "cVar")
    }
  }
  return(best_models)
}

# Clean the dispersion tables
remove_outliers <- function(x) {
  x[is.infinite(x) | x==0] <- NA
  return(x)
}

# Fit the manifold of young cells to the old cells
compute_mu <- function(y, sf, glm.param, num_cores = 4) {
  n_genes <- nrow(y)
  n_cells <- ncol(y)
  
  mu <- matrix(0, nrow = n_genes, ncol = n_cells)
  dimnames(mu) <- dimnames(y)
  
  log_sf <- log(sf)
  log_y <- log(y + 1)
  
  for (g in 1:n_genes) {
    message(paste0("Estimate the prior mean for ", rownames(y)[g]))
    
    glm_vector <- glm.param[[g]]
    
    if (is.na(glm_vector[1]) || glm_vector[1] == "(No Model)") {
      mu[g, ] <- rowMeans(y[g, , drop = FALSE])
    } else {
      intercept <- glm_vector[1]
      gene_names <- names(glm_vector)[-1]
      betas <- glm_vector[-1]
      
      gene_indices <- match(gene_names, rownames(y))
      func1 <- log_y[gene_indices, , drop = FALSE] - log_sf
      log_mu <- intercept + colSums(betas * func1)
      mu[g, ] <- exp(log_mu)
    }
  }
  return(mu)
}

# Plot cell types by development from hematopoiesis 
lineages <- list(
  "Erythroid" = c("hematopoietic precursor cell", "megakaryocyte-erythroid progenitor cell", 
                  "proerythroblast", "erythroblast"),
  "Myeloid" = c("hematopoietic precursor cell", "promonocyte", "monocyte", "macrophage", 
                "granulocytopoietic cell", "granulocyte"),
  "B cell" = c("hematopoietic precursor cell", "late pro-B cell", 
               "precursor B cell", "immature B cell", "naive B cell", "plasma cell"),
  "T and NK cell" = c("hematopoietic precursor cell", "naive T cell", "NK cell")
)

plot_lineage <- function(lineage_name, lineage_cells, combined_data, pairwise_gene, pairwise_cell) {
  
  lineage_data <- combined_data %>%
    filter(CellType %in% lineage_cells) %>%
    mutate(CellType = factor(CellType, levels = lineage_cells))
  
  p <- ggplot(lineage_data, aes(x = Age, y = Dispersion_Value, fill = Age)) +
    geom_violin(trim = TRUE) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2) +
    scale_fill_brewer(palette = "Blues") +
    facet_grid(Dispersion_Type ~ CellType, scales = "free") +
    theme_minimal() +
    ggtitle(paste("Dispersions in", lineage_name, "Lineage"))
  
  comparisons_list <- list()
  for (cell_type in lineage_cells) {
    temp_data <- subset(lineage_data, CellType == cell_type)
    age_groups <- unique(temp_data$Age)
    if (length(age_groups) >= 2) {
      comparisons_list[[cell_type]] <- combn(age_groups, 2, simplify = FALSE)
    } else {
      comparisons_list[[cell_type]] <- NULL  # No comparisons if fewer than 2 age groups
    }
  }
  
  for (cell_type in lineage_cells) {
    # Gene-level significance
    gene_data <- filter(pairwise_gene, CellType == cell_type)
    if (!is.null(comparisons_list[[cell_type]]) && nrow(gene_data) > 0) {
      p <- p + geom_signif(
        comparisons = comparisons_list[[cell_type]],
        map_signif_level = TRUE,
        step_increase = 0.1,
        data = subset(lineage_data, CellType == cell_type & Dispersion_Type == "Gene-level")
      )
    }
    
    # Cell-level significance
    cell_data <- filter(pairwise_cell, CellType == cell_type)
    if (!is.null(comparisons_list[[cell_type]]) && nrow(cell_data) > 0) {
      p <- p + geom_signif(
        comparisons = comparisons_list[[cell_type]],
        map_signif_level = TRUE,
        step_increase = 0.1,
        data = subset(lineage_data, CellType == cell_type & Dispersion_Type == "Cell-level")
      )
    }
  }
  p
}