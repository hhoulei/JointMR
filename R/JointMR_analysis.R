
#' Estimate the causal effect using JointMR
#'
#' @description This function integrates repeated MR information obtained from all
#' exposure-outcome GWAS database pairings.
#'
#' @param exposure_filter_list A list of filtered exposure GWAS data frames.
#' @param outcome_filter_list A list of filtered outcome GWAS data frames.
#' @param original_outcome_list A list of original outcome GWAS data frames.
#' @param N_list A list containing sample sizes for LDSC estimation.
#' @param ancestry Character, e.g., "EUR".
#' @param bootstrap_time Integer, number of bootstrap iterations.
#' @param relaxed_nome_bootstrap_time Integer, number of relaxed NOME bootstraps.
#' @param original_exposure_list Optional list of original exposure GWAS data.
#' @param ldsc_cache Optional cached LDSC correlation matrices.
#' @param run_relaxed_nome Logical, whether to run the relaxed-NOME extension.
#' @param ldsc_scope Character, scope for LDSC correlation estimation ("all" or "yy").
#'
#' @export
JointMR_analysis <- function(exposure_filter_list, outcome_filter_list, original_outcome_list,
                                N_list, ancestry = "EUR", bootstrap_time = 500,
                                relaxed_nome_bootstrap_time = NULL,
                                original_exposure_list = NULL,
                                ldsc_cache = NULL,
                                run_relaxed_nome = TRUE,
                                ldsc_scope = c("all", "yy")) {

  ldsc_scope <- match.arg(ldsc_scope)
  if (isTRUE(run_relaxed_nome)) {
    ldsc_scope <- "all"
  }

  # --- 1. 数据对齐与准备 ---
  nx <- length(exposure_filter_list)
  ny <- length(outcome_filter_list)

  data_x_beta <- do.call(cbind, lapply(exposure_filter_list, function(df) df$beta))
  data_x_se <- do.call(cbind, lapply(exposure_filter_list, function(df) df$se))
  data_y_beta <- do.call(cbind, lapply(outcome_filter_list, function(df) df$beta))
  data_y_se <- do.call(cbind, lapply(outcome_filter_list, function(df) df$se))

  # --- 2. 估计 GWAS 数据库间的相关性结构 (LDSC) ---
  # 调用内部函数 compute_ldsc_matrices
  if (is.null(ldsc_cache)) {
    ldsc_cache <- compute_ldsc_matrices(original_exposure_list, original_outcome_list, scope = ldsc_scope)
  }

  cor_matrix <- build_ratio_cor_matrix(ldsc_cache$rg_y, nx, ny, names(original_outcome_list))
  valid_indices <- which(apply(data_x_beta, 1, function(x) all(is.finite(x) & x != 0))) # 简化的合法性检查

  # --- 3. 构建 Wald-ratio 并进行初始联合估计 ---
  #
  original_inputs_tau0 <- build_wald_ratio_inputs(
    data_x_beta, data_y_beta, data_y_se, cor_matrix, nx, ny, valid_indices
  )

  # 计算初始异质性 (DerSimonian-Laird method) [cite: 20, 21, 23]
  original_gls_summary <- t(sapply(seq_len(ncol(original_inputs_tau0$WR)), function(j) {
    collapse_gls_snp(as.numeric(original_inputs_tau0$WR[, j]), original_inputs_tau0$base_cov_list[[j]])
  }))
  original_dl_res <- calculate_tau2_DL(original_gls_summary[, "theta_hat"], original_gls_summary[, "se_hat"])
  original_tau_sq <- original_dl_res$tau2

  # 使用 MLE_S 评估原始因果效应 [cite: 19]
  theta_original_fit <- MLE_S(nx, ny, original_inputs_tau0$base_cov_list, original_inputs_tau0$WR, ncol(original_inputs_tau0$WR))

  # --- 4. 评估并校正水平多效性 (MR-Egger) ---
  # [cite: 24, 25, 26]
  plugin_inputs <- build_plugin_adjusted_inputs(
    data_x_beta, data_x_se, data_y_beta, data_y_se, cor_matrix, nx, ny, valid_indices, valid_indices
  )

  # 重新评估异质性 [cite: 23]
  gls_summary <- t(sapply(seq_len(ncol(plugin_inputs$WR_adj)), function(j) {
    collapse_gls_snp(as.numeric(plugin_inputs$WR_adj[, j]), plugin_inputs$base_cov_list[[j]])
  }))
  dl_res <- calculate_tau2_DL(gls_summary[, "theta_hat"], gls_summary[, "se_hat"])
  tau_sq <- dl_res$tau2

  # 多效性校正后的估计
  est <- gls_theta_plugin(plugin_inputs$base_cov_list, plugin_inputs$WR_adj, plugin_inputs$BA_list, plugin_inputs$V_alpha, tau_sq)

  # --- 5. Relaxed-NOME Extension ---
  #
  relaxed_nome_est <- list(theta_hat = NA, theta_se = NA, theta_p_value = NA)
  if (isTRUE(run_relaxed_nome)) {
    relaxed_nome_est <- estimate_relaxed_nome(
      beta_x = data_x_beta[valid_indices, , drop = FALSE],
      se_x = data_x_se[valid_indices, , drop = FALSE],
      beta_y = data_y_beta[valid_indices, , drop = FALSE],
      se_y = data_y_se[valid_indices, , drop = FALSE],
      alpha_hat = plugin_inputs$alpha_hat,
      V_alpha = plugin_inputs$V_alpha,
      rho_y = as.matrix(ldsc_cache$rg_y),
      rho_x = as.matrix(ldsc_cache$rg_x),
      xy_rho = as.matrix(ldsc_cache$xy_rho),
      tau2 = tau_sq
    )
  }

  # --- 6. 返回结果 ---
  # [cite: 29]
  return(list(
    theta_original = as.numeric(theta_original_fit),
    theta_original_se = attr(theta_original_fit, "se"),
    theta_hat = unname(est[["theta_hat"]]),
    theta_se = unname(est[["theta_se"]]),
    theta_relaxed_nome = unname(relaxed_nome_est[["theta_hat"]]),
    theta_relaxed_nome_se = unname(relaxed_nome_est[["theta_se"]]),
    tau_sq = tau_sq,
    cor_matrix = ldsc_cache$rg_y
    # ... 其他需要的返回指标
  ))
}
