#' Safe matrix inversion and linear system solution
#'
#' @description Solves a linear system using the generalized inverse from the MASS package.
#' Forces the matrix to be symmetric before inversion.
#'
#' @param Sigma A numeric matrix (typically a covariance matrix).
#' @param rhs A numeric matrix or vector representing the right-hand side of the system.
#' @param tol Scalar, tolerance for the generalized inverse.
#'
#' @return A matrix representing the solution to the linear system.
#' @importFrom MASS ginv
#' @noRd
ginv_solve <- function(Sigma, rhs, tol = sqrt(.Machine$double.eps)) {
  Sigma <- as.matrix(Sigma)
  Sigma <- (Sigma + t(Sigma)) / 2
  rhs <- as.matrix(rhs)
  MASS::ginv(Sigma, tol = tol) %*% rhs
}

#' Make a matrix positive semi-definite (PSD)
#'
#' @description Adjusts a matrix by replacing negative eigenvalues with a small positive threshold
#' to ensure numerical stability in GLS operations.
#'
#' @param Sigma A numeric matrix to be adjusted.
#' @param eps Scalar, minimum threshold for eigenvalues.
#'
#' @return A symmetric, positive semi-definite numeric matrix.
#' @noRd
make_psd <- function(Sigma, eps = sqrt(.Machine$double.eps)) {
  Sigma <- as.matrix(Sigma)
  Sigma[!is.finite(Sigma)] <- 0
  Sigma <- (Sigma + t(Sigma)) / 2
  eig <- eigen(Sigma, symmetric = TRUE)
  values <- pmax(eig$values, eps)
  eig$vectors %*% diag(values, length(values)) %*% t(eig$vectors)
}

#' Solve generalized least squares with shared random effects
#'
#' @description Internal optimization wrapper to calculate the inverse operation under
#' a shared random effect setting.
#'
#' @param Omega Base covariance matrix.
#' @param rhs Right-hand side matrix/vector.
#' @param tau2 Variance component for the random effects.
#' @param tol Tolerance parameter passed to ginv_solve.
#'
#' @return Solved matrix block.
#' @noRd
solve_shared_re <- function(Omega, rhs, tau2 = 0, tol = sqrt(.Machine$double.eps)) {
  rhs <- as.matrix(rhs)
  if (tau2 <= 0) {
    return(ginv_solve(Omega, rhs, tol = tol))
  }
  
  n <- nrow(Omega)
  z <- matrix(1, nrow = n, ncol = 1)
  Omega_inv_rhs <- ginv_solve(Omega, rhs, tol = tol)
  Omega_inv_z <- ginv_solve(Omega, z, tol = tol)
  middle <- as.numeric(1 / tau2 + crossprod(z, Omega_inv_z))
  
  Omega_inv_rhs - Omega_inv_z %*% (crossprod(z, Omega_inv_rhs) / middle)
}

#' Generalized Least Squares estimator for pleiotropy plug-in model
#'
#' @description Computes the adjusted causal effect estimate (theta) and its standard error
#' after incorporating estimated horizontal pleiotropy components.
#'
#' @param base_cov_list A list of baseline covariance matrices for each SNP.
#' @param WR_matrix Matrix of adjusted Wald ratios.
#' @param BA_list Design matrices capturing relationship with estimated intercepts.
#' @param V_alpha Covariance matrix of the estimated pleiotropy intercepts.
#' @param tau2 Heterogeneity variance component.
#'
#' @return A named numeric vector with `theta_hat`, `theta_se`, and `theta_p_value`.
#' @importFrom stats pnorm
#' @importFrom MASS ginv
#' @noRd
gls_theta_plugin <- function(base_cov_list, WR_matrix, BA_list, V_alpha, tau2 = 0) {
  n_ratio <- nrow(WR_matrix)
  n_snp <- ncol(WR_matrix)
  n_total <- n_ratio * n_snp
  one <- matrix(1, nrow = n_total, ncol = 1)
  y <- matrix(as.vector(WR_matrix), nrow = n_total, ncol = 1)
  U <- do.call(rbind, BA_list)
  
  D_inv_apply <- function(X) {
    X <- as.matrix(X)
    out <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X))
    for (j in seq_len(n_snp)) {
      row_ids <- ((j - 1) * n_ratio + 1):(j * n_ratio)
      out[row_ids, ] <- solve_shared_re(
        base_cov_list[[j]],
        X[row_ids, , drop = FALSE],
        tau2 = tau2
      )
    }
    out
  }
  
  D_inv_one <- D_inv_apply(one)
  D_inv_y <- D_inv_apply(y)
  D_inv_U <- D_inv_apply(U)
  
  middle <- MASS::ginv(V_alpha) + t(U) %*% D_inv_U
  middle_inv <- MASS::ginv((middle + t(middle)) / 2)
  
  Sigma_inv_one <- D_inv_one - D_inv_U %*% middle_inv %*% t(U) %*% D_inv_one
  Sigma_inv_y <- D_inv_y - D_inv_U %*% middle_inv %*% t(U) %*% D_inv_y
  
  numerator <- as.numeric(crossprod(one, Sigma_inv_y))
  denominator <- as.numeric(crossprod(one, Sigma_inv_one))
  theta_hat <- numerator / denominator
  theta_se <- sqrt(1 / denominator)
  theta_p_value <- 2 * stats::pnorm(-abs(theta_hat / theta_se))
  
  c(theta_hat = theta_hat, theta_se = theta_se, theta_p_value = theta_p_value)
}

#' Structured Covariance Maximum Likelihood / GLS Engine (MLE_S)
#'
#' @description Jointly analyzes correlated Wald-ratio estimates using a structured
#' covariance model under the fixed-intercept assumption.
#'
#' @param nx Integer, number of exposure datasets.
#' @param ny Integer, number of outcome datasets.
#' @param cov_list List of SNP-specific covariance matrices.
#' @param WR_matrix Matrix of calculated Wald ratios.
#' @param SNPnew Integer, number of valid independent instrumental variables.
#'
#' @return A scalar numeric value representing the estimated theta, with a `se` attribute.
#' @importFrom MASS ginv
#' @noRd
MLE_S <- function(nx, ny, cov_list, WR_matrix, SNPnew) {
  n <- nx * ny
  I <- matrix(1, n, 1)
  precision_sum <- 0
  weighted_sum <- 0
  
  for (j in seq_len(SNPnew)) {
    Sigma <- cov_list[[j]]
    Sigma <- as.matrix(Sigma)
    Sigma <- (Sigma + t(Sigma)) / 2
    eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
    min_eig <- min(eig, na.rm = TRUE)
    if (!is.finite(min_eig) || min_eig <= sqrt(.Machine$double.eps)) {
      jitter <- if (is.finite(min_eig)) abs(min_eig) + sqrt(.Machine$double.eps) else sqrt(.Machine$double.eps)
      Sigma <- Sigma + diag(jitter, nrow(Sigma))
    }
    Sigma_inv <- MASS::ginv(Sigma)
    precision_sum <- precision_sum + as.numeric(crossprod(I, Sigma_inv %*% I))
    weighted_sum <- weighted_sum + as.numeric(crossprod(I, Sigma_inv %*% WR_matrix[, j, drop = FALSE]))
  }
  
  theta <- weighted_sum / precision_sum
  attr(theta, "se") <- sqrt(1 / precision_sum)
  theta
}

#' Profile Q function for Relaxed-NOME model
#'
#' @description Calculates the global joint residual Q-statistic profile conditional on a specific theta.
#' Used as the objective function for optimization in Relaxed-NOME.
#'
#' @param theta Numeric value of the causal effect parameter.
#' @param beta_x Matrix of SNP-exposure effects.
#' @param se_x Matrix of SNP-exposure standard errors.
#' @param beta_y Matrix of SNP-outcome effects.
#' @param se_y Matrix of SNP-outcome standard errors.
#' @param alpha_hat Vector of estimated multi-pleiotropy intercepts.
#' @param V_alpha Covariance matrix of the estimated intercepts.
#' @param rho_y Matrix of outcome database genetic correlations.
#' @param rho_x Matrix of exposure database genetic correlations.
#' @param xy_rho Matrix of cross-database overlap intercepts.
#' @param tau2 Residual heterogeneity component.
#'
#' @return A single numeric value representing the aggregated Q-statistic.
#' @noRd
jointmr_profile_q <- function(theta, beta_x, se_x, beta_y, se_y,
                              alpha_hat, V_alpha, rho_y, rho_x,
                              xy_rho, tau2 = 0) {
  q_value <- 0
  n_snp <- nrow(beta_x)
  nx_local <- ncol(beta_x)
  ny_local <- ncol(beta_y)
  
  for (idx in seq_len(n_snp)) {
    beta_x_vec <- as.numeric(beta_x[idx, ])
    se_x_vec <- as.numeric(se_x[idx, ])
    beta_y_vec <- as.numeric(beta_y[idx, ] - alpha_hat)
    se_y_vec <- as.numeric(se_y[idx, ])
    
    cov_x <- rho_x * outer(se_x_vec, se_x_vec)
    cov_y <- rho_y * outer(se_y_vec, se_y_vec)
    if (!is.null(V_alpha) && all(is.finite(V_alpha)) && max(abs(V_alpha)) > .Machine$double.eps) {
      cov_y <- cov_y + V_alpha
    }
    cov_xy <- xy_rho * outer(se_x_vec, se_y_vec)
    
    Sigma <- rbind(
      cbind(cov_x, cov_xy),
      cbind(t(cov_xy), cov_y)
    )
    Sigma <- make_psd(Sigma)
    
    observed <- matrix(c(beta_x_vec, beta_y_vec), ncol = 1)
    design <- matrix(c(rep(1, nx_local), rep(theta, ny_local)), ncol = 1)
    denom <- as.numeric(t(design) %*% ginv_solve(Sigma, design))
    if (!is.finite(denom) || denom <= 0) {
      return(Inf)
    }
    gamma_hat <- as.numeric(t(design) %*% ginv_solve(Sigma, observed) / denom)
    
    if (tau2 > 0) {
      cov_y_re <- cov_y + tau2 * gamma_hat^2 * matrix(1, nrow = ny_local, ncol = ny_local)
      Sigma <- rbind(
        cbind(cov_x, cov_xy),
        cbind(t(cov_xy), cov_y_re)
      )
      Sigma <- make_psd(Sigma)
      denom <- as.numeric(t(design) %*% ginv_solve(Sigma, design))
      if (!is.finite(denom) || denom <= 0) {
        return(Inf)
      }
      gamma_hat <- as.numeric(t(design) %*% ginv_solve(Sigma, observed) / denom)
    }
    
    residual <- observed - design * gamma_hat
    q_j <- as.numeric(t(residual) %*% ginv_solve(Sigma, residual))
    if (!is.finite(q_j)) {
      return(Inf)
    }
    q_value <- q_value + q_j
  }
  
  q_value
}

#' Relaxed-NOME parameter optimizer
#'
#' @description Implements the relaxed-NOME extension by accounting for measurement errors
#' in the SNP-exposure associations. It profile-optimizes the integrated Q-statistic.
#'
#' @param beta_x Matrix of SNP-exposure effects.
#' @param se_x Matrix of SNP-exposure standard errors.
#' @param beta_y Matrix of SNP-outcome effects.
#' @param se_y Matrix of SNP-outcome standard errors.
#' @param alpha_hat Vector of estimated pleiotropy intercepts.
#' @param V_alpha Covariance matrix of the intercepts.
#' @param rho_y Matrix of outcome genetic correlations.
#' @param rho_x Matrix of exposure genetic correlations.
#' @param xy_rho Matrix of cross-study overlap intercepts.
#' @param tau2 Heterogeneity variance component.
#'
#' @return A named numeric vector with `theta_hat`, `theta_se`, and `theta_p_value`.
#' @importFrom stats quantile optimize pnorm
#' @noRd
estimate_relaxed_nome <- function(beta_x, se_x, beta_y, se_y,
                                  alpha_hat, V_alpha, rho_y, rho_x,
                                  xy_rho, tau2 = 0) {
  theta_start_values <- c()
  ratio_values <- c()
  nx_local <- ncol(beta_x)
  ny_local <- ncol(beta_y)
  
  for (m in seq_len(ny_local)) {
    for (n in seq_len(nx_local)) {
      bx <- beta_x[, n]
      by <- beta_y[, m] - alpha_hat[m]
      keep <- is.finite(bx) & is.finite(by) & abs(bx) > .Machine$double.eps
      if (any(keep)) {
        theta_start_values <- c(theta_start_values, sum(bx[keep] * by[keep]) / sum(bx[keep]^2))
        ratio_values <- c(ratio_values, by[keep] / bx[keep])
      }
    }
  }
  
  theta_start <- mean(theta_start_values, na.rm = TRUE)
  if (!is.finite(theta_start)) theta_start <- 0
  
  ratio_values <- ratio_values[is.finite(ratio_values)]
  if (length(ratio_values) >= 5) {
    ratio_range <- as.numeric(stats::quantile(ratio_values, c(0.01, 0.99), na.rm = TRUE))
    lower <- min(-1, theta_start - 1, ratio_range[1] - 0.5)
    upper <- max(1, theta_start + 1, ratio_range[2] + 0.5)
  } else {
    lower <- theta_start - 1
    upper <- theta_start + 1
  }
  
  q_fun <- function(theta) {
    jointmr_profile_q(theta, beta_x, se_x, beta_y, se_y,
                      alpha_hat, V_alpha, rho_y, rho_x, xy_rho,
                      tau2 = tau2)
  }
  
  grid <- seq(lower, upper, length.out = 201)
  grid_q <- vapply(grid, q_fun, numeric(1))
  best <- which.min(grid_q)
  if (length(best) == 0 || !is.finite(grid_q[best])) {
    return(c(theta_hat = NA_real_, theta_se = NA_real_, theta_p_value = NA_real_))
  }
  
  opt_interval <- if (best == 1 || best == length(grid)) {
    c(lower, upper)
  } else {
    c(grid[best - 1], grid[best + 1])
  }
  
  opt <- stats::optimize(q_fun, interval = opt_interval)
  theta_hat <- opt$minimum
  q0 <- opt$objective
  
  h <- max(1e-4, abs(theta_hat) * 1e-4)
  q_minus <- q_fun(theta_hat - h)
  q_plus <- q_fun(theta_hat + h)
  hessian <- (q_plus - 2 * q0 + q_minus) / h^2
  if (!is.finite(hessian) || hessian <= 0) {
    h <- h * 10
    q_minus <- q_fun(theta_hat - h)
    q_plus <- q_fun(theta_hat + h)
    hessian <- (q_plus - 2 * q0 + q_minus) / h^2
  }
  
  theta_se <- if (is.finite(hessian) && hessian > 0) sqrt(2 / hessian) else NA_real_
  theta_p_value <- if (is.finite(theta_se) && theta_se > 0) {
    2 * stats::pnorm(-abs(theta_hat / theta_se))
  } else {
    NA_real_
  }
  
  c(theta_hat = theta_hat, theta_se = theta_se, theta_p_value = theta_p_value)
}