# JointMR: Joint Mendelian Randomization for Estimating Causal Effects

**JointMR** is a novel analytical framework designed to estimate causal effects by integrating repeated Mendelian Randomization (MR) information obtained from all exposure-outcome GWAS database pairings. 

Traditional MR relies on single GWAS datasets. In contrast, JointMR leverages multiple GWAS summary datasets for the same exposure and multiple GWAS summary datasets for the same outcome to improve statistical power and robustness.

---

## 📦 Part 1: Installation

You can install the development version of JointMR from GitHub. 

Since JointMR depends on `ldscr` for genetic correlation estimation, please ensure it is installed first. If you don't have `devtools` installed, please run `install.packages("devtools")` first.

```r
# Install the ldscr dependency
devtools::install_github("Bulik/ldsc") # Replace with your specific ldscr source if different

# Install JointMR
devtools::install_github("hhoulei/JointMR")
```

## 🚀 Part 2: Quick Start (Using Built-in Sample Data)

For a quick demonstration of the JointMR package, you can run the core model using the built-in, pre-harmonized sample datasets. This allows you to test the model fitting and output structure instantly.

```r
library(JointMR)

# 1. Load the built-in sample data provided by the package
data("sample_exposure")
data("sample_outcome")
data("sample_N_list")

# 2. Run the core JointMR method
results <- JointMR_analysis(
  exposure_filter_list = sample_exposure,
  outcome_filter_list = sample_outcome,
  original_outcome_list = sample_outcome,
  N_list = sample_N_list,
  ancestry = "EUR",
  bootstrap_time = 500,
  run_relaxed_nome = TRUE,
  ldsc_scope = "all"
)

# 3. View the key outputs
print(paste("Original JointMR Estimate:", round(results$theta_original, 4)))
print(paste("Pleiotropy-adjusted Estimate:", round(results$theta_hat, 4)))
print(paste("Relaxed-NOME Estimate:", round(results$theta_relaxed_nome, 4)))
print(paste("Estimated Heterogeneity (tau^2):", round(results$tau_sq, 4)))
```

## 📊 Part 3: Real-data Analysis Workflow

For actual data analysis, preparing the instrumental variables requires rigorous quality control and harmonization. Below is the standard workflow explaining how raw GWAS data is processed into the final JointMR estimates. The R function in Step 1 and 2 can be download from https://github.com/hhoulei/JointMR_Simul/application-final.R, and Step 3 can be implemented using functions `JointMR_analysis()` from above R package JointMR.  

### Step 1: Input GWAS summary datasets

JointMR requires multiple GWAS summary datasets for the same exposure and multiple GWAS summary datasets for the same outcome. In the lipid-T2D example, separate exposure analyses were conducted for HDL-C, LDL-C, total cholesterol and triglycerides, while five T2D GWAS datasets were used as outcome sources.

Each GWAS dataset should be provided as a data frame containing at least the following columns:

- `SNP`: SNP identifier
- `beta`: estimated SNP-trait association
- `se`: standard error of the association estimate
- `pval`: association P-value
- `effect_allele`: effect allele
- `other_allele`: other allele

For analyses using LDSC to estimate correlations among GWAS databases, genome-wide summary statistics and corresponding sample sizes are additionally required. Exposure and outcome datasets should be derived from compatible ancestry groups and use harmonizable allele coding.

### Step 2: Prepare the instrument set

The instrument preparation step is implemented using the custom function `full_dataset_processing()`.

1. **Initial Identification:** For each exposure GWAS dataset, candidate SNPs are first identified using a genome-wide significance threshold of 5e-8 followed by linkage disequilibrium clumping. 
2. **Instrument Strength Filtering:** The primary JointMR analysis retains SNPs with an F-statistic greater than 10 in all exposure GWAS datasets. A relaxed-F sensitivity analysis instead retains SNPs with an F-statistic greater than 10 in at least one exposure GWAS dataset. 
3. **Global Clumping:** The SNPs selected from the exposure datasets are then combined into a union set and subjected to a final LD clumping procedure. For this global clumping step, the exposure evidence is summarized using the custom function `Correlated_Stouffer()`, which combines signed exposure Z statistics while accounting for correlation among exposure GWAS databases. LD clumping is performed using `ld_clump_local()`. 
4. **Harmonization:** Finally, the retained SNPs are matched to all outcome GWAS datasets and aligned to a common effect allele using the custom function `harmonize_to_reference()`. The resulting aligned exposure and outcome datasets form the input for JointMR estimation.

### Step 3: Estimate the causal effect using JointMR

The main analysis is performed using our R package's core function `JointMR_analysis()`. This function integrates repeated MR information obtained from all exposure-outcome GWAS database pairings and includes several key components.

* **Correlation Estimation:** First, JointMR estimates the correlation structure among GWAS databases. Correlations among outcome GWAS datasets are estimated from genome-wide summary statistics using `ldscr::ldsc_rg()`. For the relaxed-NOME analysis, the function additionally evaluates correlations among exposure GWAS datasets and exposure-outcome covariance components when sample overlap information is available.
* **Wald Ratio Construction:** Second, the function constructs Wald-ratio estimates for every selected SNP and every exposure-outcome database pairing. Each Wald ratio is obtained by comparing the SNP-outcome association estimate with the corresponding SNP-exposure association estimate. These correlated Wald-ratio estimates are then jointly analyzed using the structured covariance model implemented in `MLE_S()`.
* **Heterogeneity Assessment:** Third, JointMR estimates IV-level heterogeneity. For each SNP, the repeated correlated Wald-ratio estimates are first combined into one SNP-level estimate using the custom function `collapse_gls_snp()`. The residual heterogeneity across SNP-specific causal estimates is then estimated using the DerSimonian-Laird method implemented in `calculate_tau2_DL()`. This heterogeneity component represents variation among instrumental variables rather than differences among GWAS databases. 
* **Pleiotropy Correction:** Fourth, horizontal pleiotropy is evaluated using outcome-specific MR-Egger intercepts estimated by the custom function `estimate_egger_intercepts()`. When directional horizontal pleiotropy is considered, `build_plugin_adjusted_inputs()` constructs pleiotropy-adjusted Wald ratios and propagates the uncertainty of the estimated intercepts into the covariance structure. JointMR then returns the corresponding pleiotropy-adjusted causal estimate.
* **Relaxed-NOME Extension:** Finally, the relaxed-NOME extension is implemented using `estimate_relaxed_nome()`. Unlike the standard model, this extension does not treat SNP-exposure associations as measured without error. It jointly models SNP-exposure and SNP-outcome association estimates while incorporating exposure-exposure, outcome-outcome and available exposure-outcome covariance components.

The main results returned include the original JointMR estimate, the pleiotropy-adjusted estimate, the relaxed-NOME estimate, the estimated between-SNP heterogeneity variance, the horizontal pleiotropy test result and the estimated database correlation matrices.

---

## 📝 Citation

If you use JointMR in your research, please cite our paper:

> **JointMR: A joint likelihood-based approach for causal effect estimation using multiple GWAS summary databases in Mendelian randomization.**
> [Sijia Wu, Lei Hou, Hongkai Li, Fuzhong Xue]  
> 2026.  

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details.
