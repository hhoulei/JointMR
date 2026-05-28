#' Extract sample size from a GWAS data frame
#'
#' @description Safely attempts to extract the median sample size from common sample size 
#' column names in a GWAS data frame, or uses a provided list of sample sizes.
#'
#' @param df A data frame containing GWAS summary statistics.
#' @param name Character, the name of the dataset used to look up `N_list`.
#'
#' @return A numeric scalar representing the sample size, or `NA_real_` if not found.
#' @importFrom stats median
#' @noRd
get_sample_size <- function(df, name) {
  # N_list must be available in the parent environment or explicitly passed.
  # Assuming N_list is handled upstream or passed as an implicit global/closure.
  # For safety in package development, it's better to explicitly pass N_list if possible,
  # but keeping original logic structure here.
  if (exists("N_list") && !is.null(N_list[[name]]) && is.finite(as.numeric(N_list[[name]][1]))) {
    return(as.numeric(N_list[[name]][1]))
  }
  
  candidate_cols <- c("N", "n", "samplesize", "sample_size", "N_total",
                      "n_total", "samplesize.exposure", "samplesize.outcome")
  hit <- candidate_cols[candidate_cols %in% names(df)]
  
  if (length(hit) > 0) {
    vals <- suppressWarnings(as.numeric(df[[hit[1]]]))
    vals <- vals[is.finite(vals) & vals > 0]
    if (length(vals) > 0) {
      return(stats::median(vals))
    }
  }
  
  NA_real_
}

#' Prepare GWAS summary statistics for LDSC
#'
#' @description Formats a list of GWAS data frames into the specific munged summary 
#' statistics format required by the `ldscr` package.
#'
#' @param data_list A list of GWAS data frames.
#' @param label Character string indicating the type of dataset (e.g., "exposure" or "outcome") for logging.
#'
#' @return A list of munged data frames ready for LDSC.
#' @noRd
prepare_ldsc_sumstats <- function(data_list, label) {
  if (is.null(data_list)) {
    return(list())
  }
  
  out <- list()
  for (name in names(data_list)) {
    df <- as.data.frame(data_list[[name]])
    N_val <- get_sample_size(df, name)
    
    if (!is.finite(N_val) || N_val <= 0) {
      warning("No valid sample size found for ", label, " dataset ", name,
              "; LDSC pairwise correlation involving this dataset will be set to 0.")
      next
    }
    
    munged_df <- data.frame(
      SNP = df$SNP,
      A1 = df$effect_allele,
      A2 = df$other_allele,
      Z = df$beta / df$se,
      N = N_val
    )
    
    munged_df <- munged_df[is.finite(munged_df$Z), , drop = FALSE]
    if (nrow(munged_df) > 0) {
      out[[name]] <- munged_df
    }
  }
  out
}

#' Extract specific value from an LDSC result object
#'
#' @description Safely extracts either the genetic correlation (rg) or the cross-trait 
#' sample overlap intercept from an `ldscr` output object.
#'
#' @param ldsc_res The result object returned by `ldscr::ldsc_rg()`.
#' @param value_type Character, either "rg" or "intercept".
#'
#' @return A numeric scalar representing the requested value, or `NA_real_` if missing/failed.
#' @noRd
extract_ldsc_value <- function(ldsc_res, value_type = c("rg", "intercept")) {
  value_type <- match.arg(value_type)
  
  if (is.null(ldsc_res)) {
    return(NA_real_)
  }
  
  if (value_type == "rg") {
    if (is.null(ldsc_res$rg)) {
      return(NA_real_)
    }
    rg_tab <- as.data.frame(ldsc_res$rg)
    if ("rg" %in% names(rg_tab)) {
      return(as.numeric(rg_tab$rg[1]))
    }
    return(NA_real_)
  }
  
  # For XY sample-overlap-related covariance, use the off-diagonal element
  if (!is.null(ldsc_res$raw) && !is.null(ldsc_res$raw$I)) {
    I_mat <- as.matrix(ldsc_res$raw$I)
    if (nrow(I_mat) >= 2 && ncol(I_mat) >= 2 && is.finite(as.numeric(I_mat[1, 2]))) {
      return(as.numeric(I_mat[1, 2]))
    }
  }
  
  # Fallback for LDSC wrapper versions exposing intercept directly
  if (!is.null(ldsc_res$rg)) {
    rg_tab <- as.data.frame(ldsc_res$rg)
    intercept_cols <- grep("intercept|gcov.*int|int.*gcov|cov.*int|int.*cov",
                           names(rg_tab), ignore.case = TRUE, value = TRUE)
    if (length(intercept_cols) > 0) {
      return(as.numeric(rg_tab[[intercept_cols[1]]][1]))
    }
  }
  
  NA_real_
}

#' Run pairwise LDSC regression
#'
#' @description Wrapper around `ldscr::ldsc_rg()` to compute genetic correlation 
#' with safe error catching.
#'
#' @param sumstats_i Munged summary statistics for dataset 1.
#' @param sumstats_j Munged summary statistics for dataset 2.
#' @param name_i Character, name of dataset 1.
#' @param name_j Character, name of dataset 2.
#' @param ancestry Character, e.g., "EUR", passed to ldscr.
#'
#' @return The LDSC result object, or NULL if it fails.
#' @noRd
run_ldsc_pair <- function(sumstats_i, sumstats_j, name_i, name_j, ancestry = "EUR") {
  munge_pair <- list(sumstats_i, sumstats_j)
  names(munge_pair) <- c(name_i, name_j)
  
  tryCatch(
    ldscr::ldsc_rg(
      munged_sumstats = munge_pair,
      ancestry = ancestry
    ),
    error = function(e) {
      warning(paste("ldsc_rg failed for pair", name_i, "vs", name_j, ". Error:", e$message))
      NULL
    }
  )
}

#' Infer dataset source cohort from name
#'
#' @description Uses regex matching to infer the originating cohort (e.g., UKB, FINNGEN) 
#' from the dataset name to determine potential sample overlap.
#'
#' @param dataset_name Character string of the dataset name.
#'
#' @return A character string representing the inferred source, or `NA_character_`.
#' @noRd
infer_dataset_source <- function(dataset_name) {
  dataset_name <- toupper(gsub("[^A-Z0-9]+", "_", dataset_name))
  if (grepl("UKB|UK_BIOBANK|LOH", dataset_name)) return("UKB")
  if (grepl("MVP|MILLION", dataset_name)) return("MVP")
  if (grepl("FINN|FINNGEN", dataset_name)) return("FINNGEN")
  if (grepl("EPIC", dataset_name)) return("EPIC")
  if (grepl("SWED", dataset_name)) return("SWEDISH")
  if (grepl("GLGC|WILLER", dataset_name)) return("GLGC")
  NA_character_
}

#' Check for XY sample overlap
#'
#' @description Checks if an exposure dataset and an outcome dataset likely share sample 
#' overlap based on their inferred source cohorts.
#'
#' @param exposure_name Character, name of the exposure dataset.
#' @param outcome_name Character, name of the outcome dataset.
#'
#' @return Logical, `TRUE` if sample overlap is suspected.
#' @noRd
has_xy_sample_overlap <- function(exposure_name, outcome_name) {
  exposure_source <- infer_dataset_source(exposure_name)
  outcome_source <- infer_dataset_source(outcome_name)
  !is.na(exposure_source) && !is.na(outcome_source) && (exposure_source == outcome_source)
}

#' Compute all LDSC correlation matrices
#'
#' @description Systematically runs pairwise LDSC regression across datasets to build the 
#' outcome-outcome (yy), exposure-exposure (xx), and cross-trait (xy) correlation matrices.
#'
#' @param original_exposure_list List of original exposure GWAS datasets.
#' @param original_outcome_list List of original outcome GWAS datasets.
#' @param scope Character, whether to compute "all" correlations or just outcome-outcome ("yy").
#' @param ancestry Character, passed to the LDSC engine.
#'
#' @return A list containing the matrices: `rg_x`, `rg_y`, `xy_rho`, `xy_overlap`, and `xy_rho_available`.
#' @noRd
compute_ldsc_matrices <- function(original_exposure_list, original_outcome_list,
                                  scope = c("all", "yy"), ancestry = "EUR") {
  scope <- match.arg(scope)
  
  exposure_names <- names(original_exposure_list)
  outcome_names <- names(original_outcome_list)
  nx0 <- length(exposure_names)
  ny0 <- length(outcome_names)
  
  exposure_sumstats <- if (scope == "all") {
    prepare_ldsc_sumstats(original_exposure_list, "exposure")
  } else {
    list()
  }
  outcome_sumstats <- prepare_ldsc_sumstats(original_outcome_list, "outcome")
  
  # Initialize XX matrix
  rg_x <- matrix(0, nrow = nx0, ncol = nx0, dimnames = list(exposure_names, exposure_names))
  diag(rg_x) <- 1
  if (scope == "all" && nx0 > 1) {
    for (i in seq_len(nx0 - 1)) {
      for (j in (i + 1):nx0) {
        name_i <- exposure_names[i]
        name_j <- exposure_names[j]
        if (!is.null(exposure_sumstats[[name_i]]) && !is.null(exposure_sumstats[[name_j]])) {
          cat("  Calculating XX rg between:", name_i, "and", name_j, "\n")
          ldsc_pair <- run_ldsc_pair(exposure_sumstats[[name_i]], exposure_sumstats[[name_j]], name_i, name_j, ancestry)
          rg_val <- extract_ldsc_value(ldsc_pair, "rg")
          if (is.finite(rg_val)) {
            rg_x[name_i, name_j] <- rg_val
            rg_x[name_j, name_i] <- rg_val
          }
        }
      }
    }
  }
  
  # Initialize YY matrix
  rg_y <- matrix(0, nrow = ny0, ncol = ny0, dimnames = list(outcome_names, outcome_names))
  diag(rg_y) <- 1
  if (ny0 > 1) {
    for (i in seq_len(ny0 - 1)) {
      for (j in (i + 1):ny0) {
        name_i <- outcome_names[i]
        name_j <- outcome_names[j]
        if (!is.null(outcome_sumstats[[name_i]]) && !is.null(outcome_sumstats[[name_j]])) {
          cat("  Calculating YY rg between:", name_i, "and", name_j, "\n")
          ldsc_pair <- run_ldsc_pair(outcome_sumstats[[name_i]], outcome_sumstats[[name_j]], name_i, name_j, ancestry)
          rg_val <- extract_ldsc_value(ldsc_pair, "rg")
          if (is.finite(rg_val)) {
            rg_y[name_i, name_j] <- rg_val
            rg_y[name_j, name_i] <- rg_val
          }
        }
      }
    }
  }
  
  # Initialize XY matrix
  xy_rho <- matrix(0, nrow = nx0, ncol = ny0, dimnames = list(exposure_names, outcome_names))
  xy_overlap <- matrix(FALSE, nrow = nx0, ncol = ny0, dimnames = list(exposure_names, outcome_names))
  xy_rho_available <- matrix(FALSE, nrow = nx0, ncol = ny0, dimnames = list(exposure_names, outcome_names))
  
  if (scope == "all") {
    for (i in seq_len(nx0)) {
      for (j in seq_len(ny0)) {
        name_x <- exposure_names[i]
        name_y <- outcome_names[j]
        if (!has_xy_sample_overlap(name_x, name_y)) {
          next
        }
        xy_overlap[name_x, name_y] <- TRUE
        if (!is.null(exposure_sumstats[[name_x]]) && !is.null(outcome_sumstats[[name_y]])) {
          cat("  Calculating XY LDSC intercept between:", name_x, "and", name_y, "\n")
          ldsc_pair <- run_ldsc_pair(exposure_sumstats[[name_x]], outcome_sumstats[[name_y]], name_x, name_y, ancestry)
          intercept_val <- extract_ldsc_value(ldsc_pair, "intercept")
          if (is.finite(intercept_val)) {
            xy_rho[name_x, name_y] <- intercept_val
            xy_rho_available[name_x, name_y] <- TRUE
          } else {
            warning("XY LDSC overlap intercept was not available for ", name_x, " vs ", name_y)
          }
        }
      }
    }
  }
  
  # Constrain bounds to prevent matrix instability
  rg_x[!is.finite(rg_x)] <- 0
  rg_y[!is.finite(rg_y)] <- 0
  xy_rho[!is.finite(xy_rho)] <- 0
  
  rg_x[rg_x > 0.99] <- 0.99
  rg_x[rg_x < -0.99] <- -0.99
  rg_y[rg_y > 0.99] <- 0.99
  rg_y[rg_y < -0.99] <- -0.99
  xy_rho[xy_rho > 0.99] <- 0.99
  xy_rho[xy_rho < -0.99] <- -0.99
  
  diag(rg_x) <- 1
  diag(rg_y) <- 1
  
  list(rg_x = rg_x, rg_y = rg_y, xy_rho = xy_rho,
       xy_overlap = xy_overlap, xy_rho_available = xy_rho_available)
}

#' Build a block correlation matrix for Wald ratios
#'
#' @description Expands the outcome-outcome genetic correlation matrix into the full block 
#' ratio correlation matrix required for the GLS covariance structure.
#'
#' @param rg_y Numeric matrix of outcome-outcome genetic correlations.
#' @param nx Integer, number of exposure datasets.
#' @param ny Integer, number of outcome datasets.
#' @param outcome_names Character vector of outcome dataset names.
#'
#' @return A numeric matrix representing the expanded correlation structure.
#' @noRd
build_ratio_cor_matrix <- function(rg_y, nx, ny, outcome_names) {
  n_total <- nx * ny
  out <- matrix(0, nrow = n_total, ncol = n_total)
  
  for (i in seq_len(ny)) {
    for (j in seq_len(ny)) {
      start_row_i <- (i - 1) * nx + 1
      end_row_i <- i * nx
      start_col_j <- (j - 1) * nx + 1
      end_col_j <- j * nx
      
      if (i == j) {
        # Block diagonal: perfectly correlated within the same outcome
        out[start_row_i:end_row_i, start_col_j:end_col_j] <- 1
      } else {
        # Off-diagonal blocks: bounded by the genetic correlation between outcomes
        out[start_row_i:end_row_i, start_col_j:end_col_j] <- rg_y[outcome_names[i], outcome_names[j]]
      }
    }
  }
  out
}