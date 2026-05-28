# Estimate the causal effect using JointMR

This function integrates repeated MR information obtained from all
exposure-outcome GWAS database pairings.

## Usage

``` r
JointMR_analysis(
  exposure_filter_list,
  outcome_filter_list,
  original_outcome_list,
  N_list,
  ancestry = "EUR",
  bootstrap_time = 500,
  relaxed_nome_bootstrap_time = NULL,
  original_exposure_list = NULL,
  ldsc_cache = NULL,
  run_relaxed_nome = TRUE,
  ldsc_scope = c("all", "yy")
)
```

## Arguments

- exposure_filter_list:

  A list of filtered exposure GWAS data frames.

- outcome_filter_list:

  A list of filtered outcome GWAS data frames.

- original_outcome_list:

  A list of original outcome GWAS data frames.

- N_list:

  A list containing sample sizes for LDSC estimation.

- ancestry:

  Character, e.g., "EUR".

- bootstrap_time:

  Integer, number of bootstrap iterations.

- relaxed_nome_bootstrap_time:

  Integer, number of relaxed NOME bootstraps.

- original_exposure_list:

  Optional list of original exposure GWAS data.

- ldsc_cache:

  Optional cached LDSC correlation matrices.

- run_relaxed_nome:

  Logical, whether to run the relaxed-NOME extension.

- ldsc_scope:

  Character, scope for LDSC correlation estimation ("all" or "yy").
