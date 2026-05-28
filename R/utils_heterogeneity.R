#' Collapse correlated Wald ratios for a single SNP using GLS
#'
#' @description Combines repeated Wald-ratio estimates for a single SNP across multiple 
#' outcome datasets into a single SNP-level causal estimate using Generalized Least Squares (GLS).
#'
#' @param theta_vec Numeric vector of Wald ratio estimates for a single SNP across outcomes.
#' @param Sigma Covariance matrix of the Wald ratio estimates for this SNP.
#'
#' @return A named numeric vector containing the collapsed estimate (`theta_hat`) and its standard error (`se_hat`).
#' @importFrom MASS ginv
#' @noRd
collapse_gls_snp <- function(theta_vec, Sigma) {
  # Note: make_psd() is defined in utils_gls.R and is available internally within the package
  Sigma <- make_psd(Sigma)
  one <- rep(1, length(theta_vec))
  Sigma_inv <- MASS::ginv(Sigma)
  
  precision <- as.numeric(crossprod(one, Sigma_inv %*% one))
  
  if (!is.finite(precision) || precision <= 0) {
    stop("Invalid GLS precision in collapse_gls_snp: ", precision)
  }
  
  theta_hat <- as.numeric(crossprod(one, Sigma_inv %*% theta_vec) / precision)
  se_hat <- sqrt(1 / precision)
  
  if (!is.finite(theta_hat) || !is.finite(se_hat) || se_hat <= 0) {
    stop("Invalid collapsed SNP theta/se in collapse_gls_snp: theta=",
         theta_hat, ", se=", se_hat)
  }
  
  c(theta_hat = theta_hat, se_hat = se_hat)
}

#' Calculate DerSimonian-Laird heterogeneity variance (tau^2)
#'
#' @description Estimates the residual heterogeneity across SNP-specific causal estimates 
#' using the DerSimonian-Laird (DL) method. This variance component represents variation 
#' among instrumental variables rather than differences among GWAS databases.
#'
#' @param theta_j Numeric vector of SNP-specific causal estimates.
#' @param se_j Numeric vector of standard errors for the SNP-specific estimates.
#'
#' @return A list containing the fixed-effects meta-analysis estimate (`theta_fixed`), 
#' estimated heterogeneity variance (`tau2`), Cochran's Q statistic (`Q`), degrees of freedom (`df`), 
#' p-value (`p_value`), and I-squared statistic (`I2`).
#' @importFrom stats pchisq
#' @noRd
calculate_tau2_DL <- function(theta_j, se_j) {
  # Filter out invalid estimates before calculation
  keep <- is.finite(theta_j) & is.finite(se_j) & se_j > 0
  if (sum(keep) < length(theta_j)) {
    warning("Dropping ", length(theta_j) - sum(keep),
            " SNPs with non-finite/non-positive collapsed theta or se for tau2 DL.")
  }
  
  theta_j <- theta_j[keep]
  se_j <- se_j[keep]
  
  if (length(theta_j) < 2) {
    stop("Fewer than 2 valid SNP-level estimates remain for tau2 DL.")
  }
  
  # Inverse-variance weights
  w_j <- 1 / (se_j^2)
  
  # Fixed-effects estimate
  theta_fixed <- sum(w_j * theta_j) / sum(w_j)
  
  # Cochran's Q statistic
  Q <- sum(w_j * (theta_j - theta_fixed)^2)
  k <- length(theta_j)
  df <- k - 1
  
  # C scaling factor for DL method
  C <- sum(w_j) - sum(w_j^2) / sum(w_j)
  
  if (!is.finite(Q) || !is.finite(C) || C <= 0 || df <= 0) {
    stop("Invalid tau2 DL components: Q=", Q,
         ", C=", C, ", df=", df)
  }
  
  # tau^2 estimation (truncated at 0)
  tau2 <- max(0, (Q - df) / C)
  
  # I^2 estimation
  I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
  
  list(
    theta_fixed = theta_fixed,
    tau2 = tau2,
    Q = Q,
    df = df,
    p_value = 1 - stats::pchisq(Q, df),
    I2 = I2
  )
}