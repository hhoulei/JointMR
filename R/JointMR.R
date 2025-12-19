#'
#' @title JointMR
#'
#' @description A joint likelihood-based approach designed to integrate multiple GWAS summary databases
#' while explicitly accounting for the covariance matrix of the Wald ratio estimates. Details refer to
#' JointMR: A joint likelihood-based approach for causal effect estimation in overlapping Mendelian Randomization studies
#' Summary Datasets by  Sijia Wu, Lei Hou, Zhongshang Yuan, Xiaoru Sun, Yuanyuan Yu, Hao Chen, Lin Huang,Hongkai Li, Fuzhong Xue
#'
#' @param exposure_list A list of exposure GWAS summary statistics datasets that have
#'   been pre-processed for SNP alignment. Each data frame should contain the following
#'   columns: SNP (rsID), effect_allele, other_allele, beta, se, pval.
#' @param outcome_list A list of outcome GWAS summary statistics datasets that have
#'   been pre-processed to match the SNPs in the exposure_list. Each data frame should
#'   contain the same columns as exposure_list.
#' @param original_outcome_list A list of complete outcome GWAS summary statistics datasets
#'   without SNP filtering. Used for LDSC analysis. Each data frame should contain:
#'   SNP (rsID), effect_allele, other_allele, beta, se, pval.
#' @param N_list A named list of sample sizes for each outcome dataset. Names should
#'   correspond to the names in original_outcome_list.
#' @param ancestry Character string specifying the ancestry of the study population.
#'   Used for LDSC reference panel selection (e.g., "EUR", "EAS", "AFR", "SAS", "mixed").
#' @param bootstrap_time Integer specifying the number of bootstrap iterations for
#'   standard error estimation. Typically 1000-10000 iterations.
#'
#' @return A list containing:
#'   \item{theta}{The joint MR estimate (causal effect).}
#'   \item{se}{Bootstrap standard error of the estimate.}
#'   \item{pvalue}{Two-sided p-value for the estimate.}
#'   \item{nSNP}{Number of SNPs used in the analysis.}
#'   \item{cor_matrix}{Genetic correlation matrix estimated from LDSC.}
#'
#' @details
#' The function performs the following steps:
#' 1. Applies MR-Egger regression to each exposure-outcome pair to estimate
#'    and correct for horizontal pleiotropy (intercept).
#' 2. Computes Wald ratios for each SNP, adjusted for pleiotropy.
#' 3. Uses LDSC to estimate genetic correlations between outcomes.
#' 4. Constructs a covariance matrix incorporating SNP-specific variances,
#'    between-outcome genetic correlations, and between-study heterogeneity (tau²).
#' 5. Computes a maximum likelihood estimate (MLE) of the joint causal effect.
#' 6. Estimates standard errors using bootstrap resampling of SNPs.
#'
#' @note
#' Required R packages: MASS (for ginv), TwoSampleMR (for mr_egger_regression),
#' and ldscr (for ldsc_rg, available from GitHub: jean997/ldsc). Ensure these
#' packages are installed and loaded before running the function.
#'
#' @examples
#' #load("exposure_list_trimmed.Rdata")
#' #load("outcome_list_example.Rdata")
#' #load("outcome_list_all.Rdata")
#'
#' #N_list_example <- list(
#' #finn    = 486367,
#' #MVP       = 432648,
#' #EPIC      = 22326)
#'
#' #JointMR(exposure_list = exposure_list_trimmed,
#' #outcome_list = outcome_list_example,
#' #original_outcome_list = outcome_list_all,
#' #N_list = N_list_example,
#' #ancestry = "EUR",
#' #bootstrap_time = 1000)


JointMR <- function(exposure_list,
                    outcome_list,
                    original_outcome_list,
                    N_list,
                    ancestry,
                    bootstrap_time) {

  MLE_S <- function(nx,ny,cov_list,WR_matrix,SNPnew){
    n <- nx * ny
    I <- matrix(1, n, 1)

    fenmu_vec <- numeric(SNPnew)
    fenzi_vec <- numeric(SNPnew)

    for (i in 1:SNPnew) {
      Sigma <- cov_list[[i]]

      Sigma_inv <- MASS::ginv(Sigma)

      fenmu_vec[i] <- as.numeric(crossprod(I, Sigma_inv) %*% I)
      fenzi_vec[i] <- as.numeric(crossprod(I, Sigma_inv) %*% WR_matrix[, i])
    }

    fenzi_all <- sum(fenzi_vec)
    fenmu_all <- sum(fenmu_vec)

    return(as.numeric(fenzi_all / fenmu_all))
  }

  Wald_Ratio <- function(betaXG,betaYG,sebetaXG,sebetaYG){
    WR <- betaYG/betaXG
    varWR_2 <- (sebetaYG^2)/(betaXG^2)
    seWR <- sqrt(varWR_2)

    return(data.frame(WR=WR,
                      varWR_2=varWR_2,
                      seWR=seWR))
  }

  calculate_tau2_DL <- function(theta_j, se_j) {

    w_j <- 1 / (se_j^2)

    theta_fixed <- sum(w_j * theta_j) / sum(w_j)

    Q <- sum(w_j * (theta_j - theta_fixed)^2)

    k <- length(theta_j)  # 研究数量
    df <- k - 1

    C <- sum(w_j) - sum(w_j^2) / sum(w_j)

    tau2 <- max(0, (Q - df) / C)

    return(list(
      tau2 = tau2,
      Q = Q,
      df = df,
      p_value = 1 - pchisq(Q, df),
      I2 = max(0, (Q - df) / Q) * 100  # I² 统计量
    ))
  }

  filtered_exposure_list <- exposure_list
  filtered_outcome_list <- outcome_list

  nx <- length(filtered_exposure_list)
  ny <- length(filtered_outcome_list)

  gamma_all <- NULL
  for (i in 1:ny) {
    for (j in 1:nx) {
      exposure_df <- filtered_exposure_list[[j]]
      outcome_df <- filtered_outcome_list[[i]]

      res_egger <- mr_egger_regression(exposure_df$beta,
                                       outcome_df$beta,
                                       exposure_df$se,
                                       outcome_df$se)
      if(res_egger$pval_i<0.05){
        gamma1 <- res_egger$b_i
      }else{
        gamma1 <- 0
      }
      gamma_all <- c(gamma_all,gamma1)
    }
  }

  Wald_list <- list()

  for (n in 1:ny) {
    for (m in 1:nx) {
      idx <- (n-1)*nx + m
      Wald_list[[idx]] <- Wald_Ratio(
        filtered_exposure_list[[m]]$beta-gamma_all[idx],
        filtered_outcome_list[[n]]$beta-gamma_all[idx],
        filtered_exposure_list[[m]]$se,
        filtered_outcome_list[[n]]$se
      )
    }
  }

  WR_vecs <- lapply(Wald_list, function(df) df$WR)
  se_vecs <- lapply(Wald_list, function(df) df$seWR)

  WR_matrix <- do.call(rbind, WR_vecs)
  seWR_matrix <- do.call(rbind, se_vecs)

  colname <- NULL
  for (n in 1:ny) {
    for (m in 1:nx) {
      colname1 <- paste0("X", m, "-Y", n)
      colname <- c(colname, colname1)
    }
  }
  rownames(WR_matrix) <- colname
  rownames(seWR_matrix) <- colname

  if (nx > 0 && !is.null(filtered_exposure_list[[1]]$SNP)) {
    colnames(WR_matrix) <- filtered_exposure_list[[1]]$SNP
    colnames(seWR_matrix) <- filtered_exposure_list[[1]]$SNP
  }

  bad_col <- unique(c(
    which(apply(WR_matrix, 2, function(x) any(!is.finite(x)))),
    which(apply(seWR_matrix, 2, function(x) any(!is.finite(x))))
  ))

  if (length(bad_col) > 0) {
    WR_matrix <- WR_matrix[, -bad_col, drop = FALSE]
    seWR_matrix <- seWR_matrix[, -bad_col, drop = FALSE]

    outcome_list_1 <- lapply(outcome_list, function(df) {
      df[-bad_col, , drop = FALSE]
    })
  } else {
    outcome_list_1 <- outcome_list
  }


  outcome_names <- names(original_outcome_list)

  munged_sumstats_list <- list()
  for (name in outcome_names) {
    df <- as.data.frame(original_outcome_list[[name]])
    N_val <- as.numeric(N_list[[name]][1])
    munged_df <- data.frame(
      SNP = df$SNP, A1 = df$effect_allele, A2 = df$other_allele,
      Z = df$beta / df$se, N = N_val
    )
    munged_df <- munged_df[is.finite(munged_df$Z), ]
    munged_sumstats_list[[name]] <- munged_df
  }


  rg_matrix <- matrix(1, nrow = ny, ncol = ny,
                      dimnames = list(outcome_names, outcome_names))

  for (i in 1:(ny - 1)) {
    for (j in (i + 1):ny) {

      name_i <- outcome_names[i]
      name_j <- outcome_names[j]

      munge_pair <- list(munged_sumstats_list[[name_i]], munged_sumstats_list[[name_j]])
      names(munge_pair) <- c(name_i, name_j)

      ldsc_results_pair <- tryCatch({
        ldscr::ldsc_rg(
          munged_sumstats = munge_pair,
          ancestry = ancestry
        )
      }, error = function(e) {
        warning(paste("ldsc_rg failed for pair", name_i, "vs", name_j, ". Error:", e$message))
        return(NULL)
      })

      if (!is.null(ldsc_results_pair)) {
        rg_val <- ldsc_results_pair$rg$rg
        rg_matrix[name_i, name_j] <- rg_val
        rg_matrix[name_j, name_i] <- rg_val
      } else {
        rg_matrix[name_i, name_j] <- 0 # 或 NA (如果失败，则假定不相关)
        rg_matrix[name_j, name_i] <- 0
      }
    }
  }
  rg_matrix[rg_matrix > 1] <- 0.99

  n_total <- nx * ny
  cor_matrix <- matrix(0, nrow = n_total, ncol = n_total)

  for (i in 1:ny) {
    for (j in 1:ny) {
      start_row_i <- (i - 1) * nx + 1; end_row_i <- i * nx
      start_col_j <- (j - 1) * nx + 1; end_col_j <- j * nx

      if (i == j) {
        cor_matrix[start_row_i:end_row_i, start_col_j:end_col_j] <- 1
      } else {
        name_i <- outcome_names[i]
        name_j <- outcome_names[j]
        rho_value <- rg_matrix[name_i, name_j]
        cor_matrix[start_row_i:end_row_i, start_col_j:end_col_j] <- rho_value
      }
    }
  }

  samples <- rownames(WR_matrix)
  SNPnew <- ncol(WR_matrix)

  y1 = as.vector(WR_matrix)
  se1 = as.vector(seWR_matrix)

  meta_res <- tryCatch(
    calculate_tau2_DL(y1, se1),
    error = function(e) list(tau2 = 0, I2 = 0)
  )

  tau_sq <- if (meta_res$I2 > 50 && meta_res$p_value < 0.05) meta_res$tau2 else 0

  cov_list <- lapply(1:SNPnew, function(j) {
    y <- WR_matrix[, j]
    s <- seWR_matrix[, j]

    se_outer <- outer(s, s)
    cov_matrix <- cor_matrix * se_outer + tau_sq

    dimnames(cov_matrix) <- list(samples, samples)
    return(cov_matrix)
  })

  theta_hat <- MLE_S(nx, ny, cov_list, WR_matrix, SNPnew)

  boot_effects <- numeric(bootstrap_time)
  for (k in 1:bootstrap_time) {
    boot_sam <- sample(SNPnew, replace = TRUE)
    WR_matrix_sample <- WR_matrix[, boot_sam, drop = FALSE]
    cov_list_sample <- cov_list[boot_sam]

    boot_effects[k] <- tryCatch(
      MLE_S(nx, ny, cov_list_sample, WR_matrix_sample, SNPnew),
      error = function(e) NA
    )
  }

  boot_effects <- boot_effects[!is.na(boot_effects)]

  theta_se <- sd(boot_effects)
  z_score <- theta_hat / theta_se
  theta_p_value <- 2 * pnorm(-abs(z_score))

  return(list(
    theta = theta_hat,
    se = theta_se,
    pvalue = theta_p_value,
    nSNP = SNPnew,
    cor_matrix =  rg_matrix
  ))
}
