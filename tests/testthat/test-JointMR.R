# tests/testthat/test-JointMR.R

test_that("JointMR_analysis runs successfully on mock data", {

  # 1. 生成模拟的 GWAS 数据 (Mock Data)
  set.seed(2026)

  generate_mock_gwas <- function(snps) {
    n_snps <- length(snps)
    data.frame(
      SNP = snps,
      beta = rnorm(n_snps, mean = 0.05, sd = 0.1),
      se = runif(n_snps, min = 0.01, max = 0.05),
      pval = runif(n_snps, min = 1e-12, max = 1e-6),
      effect_allele = rep("A", n_snps),
      other_allele = rep("G", n_snps),
      stringsAsFactors = FALSE
    )
  }

  shared_snps <- paste0("rs", 1:150)
  exposure_names <- c("Mock_HDL", "Mock_LDL")
  outcome_names <- c("Mock_T2D_1", "Mock_T2D_2")

  exposure_filter_list <- list(
    Mock_HDL = generate_mock_gwas(shared_snps),
    Mock_LDL = generate_mock_gwas(shared_snps)
  )

  outcome_filter_list <- list(
    Mock_T2D_1 = generate_mock_gwas(shared_snps),
    Mock_T2D_2 = generate_mock_gwas(shared_snps)
  )

  N_list <- list(
    Mock_HDL = 100000, Mock_LDL = 120000,
    Mock_T2D_1 = 150000, Mock_T2D_2 = 180000
  )

  # 2. 构造虚拟的 LDSC 缓存
  mock_ldsc_cache <- list(
    rg_x = matrix(c(1.0, -0.1, -0.1,  1.0), nrow = 2,
                  dimnames = list(exposure_names, exposure_names)),
    rg_y = matrix(c(1.0, 0.8, 0.8, 1.0), nrow = 2,
                  dimnames = list(outcome_names, outcome_names)),
    xy_rho = matrix(0, nrow = 2, ncol = 2,
                    dimnames = list(exposure_names, outcome_names)),
    xy_overlap = matrix(FALSE, nrow = 2, ncol = 2,
                        dimnames = list(exposure_names, outcome_names)),
    xy_rho_available = matrix(TRUE, nrow = 2, ncol = 2,
                              dimnames = list(exposure_names, outcome_names))
  )

  # 3. 运行主函数
  # 注意：这里调用的是你改名后的 JointMR_analysis
  test_results <- JointMR_analysis(
    exposure_filter_list = exposure_filter_list,
    outcome_filter_list = outcome_filter_list,
    original_outcome_list = outcome_filter_list,
    N_list = N_list,
    bootstrap_time = 5,                 # 测试环境设小一点，加快测试速度
    relaxed_nome_bootstrap_time = 5,
    ldsc_cache = mock_ldsc_cache,
    run_relaxed_nome = TRUE,
    ldsc_scope = "all"
  )

  # 4. 核心断言 (Assertions)
  # 检查返回的是否为列表
  expect_type(test_results, "list")

  # 检查关键结果是否存在且为数值
  expect_true(is.numeric(test_results$theta_original))
  expect_true(is.numeric(test_results$theta_hat))
  expect_true(is.numeric(test_results$theta_relaxed_nome))
  expect_true(is.numeric(test_results$tau_sq))

  # 检查标准误是否大于0
  expect_true(test_results$theta_se > 0)
})
