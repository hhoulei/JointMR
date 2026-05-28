#' Meta-analyze exposure data for MR-Egger
#'
#' @description Performs fixed-effect inverse-variance weighting (IVW) meta-analysis 
#' to collapse multiple exposure datasets into a single set of exposure estimates. 
#' This is a prerequisite for calculating the MR-Egger intercept across datasets.
#'
#' @param data_x_beta A numeric matrix of SNP-exposure effect sizes.
#' @param data_x_se A numeric matrix of SNP-exposure standard errors.
#'
#' @return A list containing the meta-analyzed exposure effects (`bx`) and their standard errors (`bxse`).
#' @noRd
meta_exposure_for_egger <- function(data_x_beta, data_x_se) {
  w <- 1 / data_x_se^2
  bx_meta <- rowSums(w * data_x_beta) / rowSums(w)
  bxse_meta <- sqrt(1 / rowSums(w))
  list(bx = bx_meta, bxse = bxse_meta)
}

#' Calculate weighted MR-Egger intercept
#'
#' @description Fits a weighted linear regression (MR-Egger) to estimate directional 
#' pleiotropy (the intercept) for a single outcome dataset.
#'
#' @param bx Numeric vector of meta-analyzed SNP-exposure effects.
#' @param by Numeric vector of SNP-outcome effects.
#' @param byse Numeric vector of SNP-outcome standard errors.
#'
#' @return A numeric vector containing the estimated intercept (`alpha`) and its standard error (`se`).
#' @importFrom MASS ginv
#' @noRd
weighted_egger_intercept <- function(bx, by, byse) {
  keep <- is.finite(bx) & is.finite(by) & is.finite(byse) & byse > 0
  bx <- bx[keep]
  by <- by[keep]
  byse <- byse[keep]
  
  if (length(by) < 3 || length(unique(bx)) < 2) {
    return(c(alpha = NA_real_, se = NA_real_))
  }
  
  w <- 1 / byse^2
  X <- cbind(1, bx)
  XtWX <- crossprod(X, X * w)
  XtWy <- crossprod(X, by * w)
  
  # Ensure non-singular matrix inversion
  XtWX_inv <- MASS::ginv(XtWX)
  coef_hat <- as.numeric(XtWX_inv %*% XtWy)
  
  resid <- by - as.numeric(X %*% coef_hat)
  sigma2 <- sum(w * resid^2) / max(length(by) - 2, 1)
  se_alpha <- sqrt(max(0, sigma2 * XtWX_inv[1, 1]))
  
  c(alpha = coef_hat[1], se = se_alpha)
}

#' Estimate MR-Egger intercepts across all outcomes
#'
#' @description Computes outcome-specific MR-Egger intercepts and constructs the 
#' covariance matrix for these intercepts based on outcome genetic correlations.
#'
#' @param data_x_beta Matrix of SNP-exposure effects.
#' @param data_x_se Matrix of SNP-exposure standard errors.
#' @param data_y_beta Matrix of SNP-outcome effects.
#' @param data_y_se Matrix of SNP-outcome standard errors.
#' @param cor_matrix Overall ratio correlation matrix.
#' @param nx Integer, number of exposure datasets.
#' @param ny Integer, number of outcome datasets.
#' @param snp_indices Integer vector specifying which SNPs to use for estimation.
#'
#' @return A list containing the estimated intercepts (`alpha_hat`), their standard errors (`alpha_se`), 
#' and their full covariance matrix (`V_alpha`).
#' @noRd
estimate_egger_intercepts <- function(data_x_beta, data_x_se, data_y_beta,
                                      data_y_se, cor_matrix, nx, ny,
                                      snp_indices = NULL) {
  if (is.null(snp_indices)) {
    snp_indices <- seq_len(nrow(data_x_beta))
  }
  
  exposure_meta <- meta_exposure_for_egger(data_x_beta, data_x_se)
  bx <- exposure_meta$bx[snp_indices]
  
  alpha_hat <- rep(NA_real_, ny)
  alpha_se <- rep(NA_real_, ny)
  
  for (m in seq_len(ny)) {
    by <- data_y_beta[snp_indices, m]
    byse <- data_y_se[snp_indices, m]
    fit <- weighted_egger_intercept(bx, by, byse)
    alpha_hat[m] <- fit[["alpha"]]
    alpha_se[m] <- fit[["se"]]
  }
  
  alpha_hat[!is.finite(alpha_hat)] <- 0
  alpha_se[!is.finite(alpha_se)] <- 0
  
  # Reconstruct the outcome-outcome correlation structure for the intercepts
  rho_y <- matrix(1, nrow = ny, ncol = ny)
  for (m in seq_len(ny)) {
    for (q in seq_len(ny)) {
      rho_y[m, q] <- cor_matrix[(m - 1) * nx + 1, (q - 1) * nx + 1]
    }
  }
  
  V_alpha <- rho_y * outer(alpha_se, alpha_se)
  V_alpha[!is.finite(V_alpha)] <- 0
  
  list(alpha_hat = alpha_hat, alpha_se = alpha_se, V_alpha = V_alpha)
}

#' Build pleiotropy-adjusted Wald ratio inputs
#'
#' @description Calculates Wald ratios adjusted for the estimated MR-Egger intercepts 
#' and prepares the base covariance structures and design matrices needed for the 
#' plug-in GLS estimator.
#'
#' @param data_x_beta Matrix of SNP-exposure effects.
#' @param data_x_se Matrix of SNP-exposure standard errors.
#' @param data_y_beta Matrix of SNP-outcome effects.
#' @param data_y_se Matrix of SNP-outcome standard errors.
#' @param cor_matrix Overall ratio correlation matrix.
#' @param nx Integer, number of exposure datasets.
#' @param ny Integer, number of outcome datasets.
#' @param valid_indices Integer vector of valid SNP indices to process.
#' @param egger_indices Integer vector of SNP indices to use for estimating pleiotropy.
#'
#' @return A list containing adjusted Wald ratios (`WR_adj`), standard errors (`seWR`), 
#' base covariance matrices (`base_cov_list`), intercept design matrices (`BA_list`), 
#' and pleiotropy estimates (`alpha_hat`, `alpha_se`, `V_alpha`).
#' @noRd
build_plugin_adjusted_inputs <- function(data_x_beta, data_x_se, data_y_beta,
                                         data_y_se, cor_matrix, nx, ny,
                                         valid_indices,
                                         egger_indices = NULL) {
  n_ratio <- nx * ny
  
  egger_fit <- estimate_egger_intercepts(
    data_x_beta = data_x_beta,
    data_x_se = data_x_se,
    data_y_beta = data_y_beta,
    data_y_se = data_y_se,
    cor_matrix = cor_matrix,
    nx = nx,
    ny = ny,
    snp_indices = egger_indices
  )
  
  g <- nrow(data_x_beta)
  WR_adj_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
  seWR_matrix_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
  
  for (m in seq_len(ny)) {
    for (n in seq_len(nx)) {
      row_id <- (m - 1) * nx + n
      # Subtract pleiotropy intercept from the outcome effect
      WR_adj_full[row_id, ] <- (data_y_beta[, m] - egger_fit$alpha_hat[m]) / data_x_beta[, n]
      seWR_matrix_full[row_id, ] <- sqrt(data_y_se[, m]^2 / data_x_beta[, n]^2)
    }
  }
  
  A <- matrix(0, nrow = n_ratio, ncol = ny)
  for (m in seq_len(ny)) {
    row_ids <- ((m - 1) * nx + 1):(m * nx)
    A[row_ids, m] <- 1
  }
  
  WR_adj <- WR_adj_full[, valid_indices, drop = FALSE]
  seWR <- seWR_matrix_full[, valid_indices, drop = FALSE]
  beta_x_selected <- data_x_beta[valid_indices, , drop = FALSE]
  
  base_cov_list <- vector("list", length(valid_indices))
  BA_list <- vector("list", length(valid_indices))
  
  for (idx in seq_along(valid_indices)) {
    se_vector <- seWR[, idx]
    Omega_j <- cor_matrix * outer(se_vector, se_vector)
    
    beta_x_vec <- as.numeric(beta_x_selected[idx, ])
    B_j <- diag(rep(1 / beta_x_vec, times = ny), nrow = n_ratio)
    BA_j <- B_j %*% A
    
    dimnames(Omega_j) <- list(rownames(WR_adj), rownames(WR_adj))
    base_cov_list[[idx]] <- Omega_j
    BA_list[[idx]] <- BA_j
  }
  
  list(
    WR_adj = WR_adj,
    seWR = seWR,
    base_cov_list = base_cov_list,
    BA_list = BA_list,
    alpha_hat = egger_fit$alpha_hat,
    alpha_se = egger_fit$alpha_se,
    V_alpha = egger_fit$V_alpha
  )
}

#' Build baseline Wald ratio inputs
#'
#' @description Calculates unadjusted (or globally adjusted via fixed alpha/tau2) Wald ratios 
#' and prepares the base covariance structures needed for the standard GLS estimator.
#'
#' @param data_x_beta Matrix of SNP-exposure effects.
#' @param data_y_beta Matrix of SNP-outcome effects.
#' @param data_y_se Matrix of SNP-outcome standard errors.
#' @param cor_matrix Overall ratio correlation matrix.
#' @param nx Integer, number of exposure datasets.
#' @param ny Integer, number of outcome datasets.
#' @param valid_indices Integer vector of valid SNP indices.
#' @param alpha_hat Optional numeric vector of fixed intercept adjustments.
#' @param tau2 Optional numeric scalar for baseline heterogeneity adjustment.
#'
#' @return A list containing Wald ratios (`WR`), standard errors (`seWR`), selected exposure 
#' betas (`beta_x_selected`), and base covariance matrices (`base_cov_list`).
#' @noRd
build_wald_ratio_inputs <- function(data_x_beta, data_y_beta, data_y_se,
                                    cor_matrix, nx, ny, valid_indices,
                                    alpha_hat = rep(0, ny),
                                    tau2 = 0) {
  n_ratio <- nx * ny
  g <- nrow(data_x_beta)
  WR_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
  seWR_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
  
  for (m in seq_len(ny)) {
    for (n in seq_len(nx)) {
      row_id <- (m - 1) * nx + n
      WR_full[row_id, ] <- (data_y_beta[, m] - alpha_hat[m]) / data_x_beta[, n]
      seWR_full[row_id, ] <- sqrt(data_y_se[, m]^2 / data_x_beta[, n]^2)
    }
  }
  
  WR <- WR_full[, valid_indices, drop = FALSE]
  seWR <- seWR_full[, valid_indices, drop = FALSE]
  beta_x_selected <- data_x_beta[valid_indices, , drop = FALSE]
  base_cov_list <- vector("list", length(valid_indices))
  
  for (idx in seq_along(valid_indices)) {
    se_vector <- seWR[, idx]
    Omega_j <- cor_matrix * outer(se_vector, se_vector) + tau2
    dimnames(Omega_j) <- list(rownames(WR), rownames(WR))
    base_cov_list[[idx]] <- Omega_j
  }
  
  list(
    WR = WR,
    seWR = seWR,
    beta_x_selected = beta_x_selected,
    base_cov_list = base_cov_list
  )
}