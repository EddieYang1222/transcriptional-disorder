# Helper Functions for Computing Transcriptional Dyscoordination

######################################################
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

######################################################
# Clean the resulting dispersion tables
remove_outliers <- function(x) {
  x[is.infinite(x) | x==0] <- NA
  return(x)
}

######################################################
# Fit the manifold of young cells to old cells
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