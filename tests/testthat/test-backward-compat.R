# Backward compatibility tests
#
# Verifies that all functionality is unchanged:
# - run_single_iteration() with gamma_G = 0.6 (default) produces
# identical output to the default
# - analyze_methods_robust() / analyze_mediation_robust() unchanged
# - run_simulation(), sweep_param(), run_mediation_sim(),
# sweep_mediation_param() unchanged
# - gan_sensitivity(), gan_mediation_sensitivity(),
# recommend_estimator() unchanged

# ═══════════════════════════════════════════════════════════════
# gamma_G backward compatibility
# ═══════════════════════════════════════════════════════════════

test_that("gamma_G default 0.6 produces identical output to explicit 0.6", {
  set.seed(42)
  d1 <- run_single_iteration(n_synthetic_samples = 50, n_features = 3, seed = 42)
  set.seed(42)
  d2 <- run_single_iteration(n_synthetic_samples = 50, n_features = 3,
                             gamma_G = 0.6, seed = 42)
  expect_identical(d1$Z, d2$Z)
  expect_identical(d1$Y, d2$Y)
  expect_identical(d1$M, d2$M)
  expect_identical(d1$W, d2$W)
  expect_identical(d1$G, d2$G)
})

test_that("gamma_G default produces identical output to no gamma_G arg", {
  # Calling without gamma_G should be the same as gamma_G = 0.6
  set.seed(100)
  d1 <- run_single_iteration(n_synthetic_samples = 80, n_features = 4,
                             mo_confounding = 0.8, phi = 0.8, seed = 100)
  set.seed(100)
  d2 <- run_single_iteration(n_synthetic_samples = 80, n_features = 4,
                             mo_confounding = 0.8, phi = 0.8,
                             gamma_G = 0.6, seed = 100)
  expect_identical(d1$Z, d2$Z)
  expect_identical(d1$Y, d2$Y)
})

test_that("gamma_G changes instrument strength", {
  set.seed(42)
  dw <- run_single_iteration(n_synthetic_samples = 100, n_features = 3,
                             mo_confounding = 0.8, phi = 0.8,
                             gamma_G = 0.2, seed = 42)
  set.seed(42)
  ds <- run_single_iteration(n_synthetic_samples = 100, n_features = 3,
                             mo_confounding = 0.8, phi = 0.8,
                             gamma_G = 1.0, seed = 42)
  # Stronger gamma_G should produce stronger first-stage relationship
  cor_w <- cor(dw$G[, 1], dw$Z)
  cor_s <- cor(ds$G[, 1], ds$Z)
  expect_lt(abs(cor_w), abs(cor_s))
})

# ═══════════════════════════════════════════════════════════════
# analyze_methods_robust unchanged
# ═══════════════════════════════════════════════════════════════

test_that("analyze_methods_robust produces 5 methods on total-effect data", {
  set.seed(42)
  dat <- run_single_iteration(n_synthetic_samples = 100, n_features = 5, seed = 42)
  res <- analyze_methods_robust(dat)
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC") %in%
                  res$method))
  expect_true("significant" %in% names(res))
  expect_true("beta" %in% names(res))
  expect_true("se" %in% names(res))
  expect_true("pvalue" %in% names(res))
})

test_that("analyze_methods_robust respects test_features", {
  set.seed(42)
  dat <- run_single_iteration(n_synthetic_samples = 100, n_features = 5, seed = 42)
  res <- analyze_methods_robust(dat, test_features = 1:2)
  expect_true(all(res$feature %in% 1:2))
})

# ═══════════════════════════════════════════════════════════════
# analyze_mediation_robust unchanged
# ═══════════════════════════════════════════════════════════════

test_that("analyze_mediation_robust produces 8 methods with phi and", {
  set.seed(42)
  dat <- run_single_iteration(n_synthetic_samples = 200, n_features = 5,
                              mo_confounding = 0.8, phi = 0.8,
                              separate_U = TRUE, omega_1 = 0.7, omega_2 = 0.7,
                              seed = 42)
  res <- analyze_mediation_robust(dat)
  methods_present <- unique(res$method)
  # With phi > 0 and separate_U, all 8 methods should be present
  for (m in c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
              "IV2SLS2", "PGC2", "PGC2Gm")) {
    expect_true(m %in% methods_present, label = paste("method", m))
  }
  expect_true("NDE" %in% names(res))
  expect_true("NIE" %in% names(res))
  expect_true("NDE_significant" %in% names(res))
  expect_true("NIE_significant" %in% names(res))
})

# ═══════════════════════════════════════════════════════════════
# run_simulation unchanged
# ═══════════════════════════════════════════════════════════════

test_that("run_simulation returns summary with 5 methods", {
  sim <- run_simulation(n_iter = 5, n_samples = 100, n_features = 5,
                        base_seed = 42)
  expect_true(is.data.frame(sim$summary))
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC") %in%
                  sim$summary$method))
  expect_true("bias" %in% names(sim$summary))
  expect_true("rmse" %in% names(sim$summary))
})

# ═══════════════════════════════════════════════════════════════
# sweep_param unchanged
# ═══════════════════════════════════════════════════════════════

test_that("sweep_param returns summary across grid", {
  sw <- sweep_param("conf_str", c(0.2, 0.8),
                    n_iter = 3, n_samples = 100, n_features = 5, base_seed = 42)
  expect_true(is.data.frame(sw$summary))
  expect_true("param" %in% names(sw$summary))
  expect_true("param_value" %in% names(sw$summary))
  expect_true("method" %in% names(sw$summary))
  expect_equal(length(unique(sw$summary$param_value)), 2)
})

# ═══════════════════════════════════════════════════════════════
# run_mediation_sim unchanged
# ═══════════════════════════════════════════════════════════════

test_that("run_mediation_sim returns summary with NDE/NIE bias", {
  sim <- run_mediation_sim(n_iter = 5, n_samples = 200, n_features = 5,
                           mo_confounding = 0.8, phi = 0.8, base_seed = 42)
  expect_true(is.data.frame(sim$summary))
  expect_true("NDE_bias" %in% names(sim$summary))
  expect_true("NIE_bias" %in% names(sim$summary))
  expect_true("NDE_rmse" %in% names(sim$summary))
  expect_true("NIE_rmse" %in% names(sim$summary))
  expect_true("NIE_type1" %in% names(sim$summary))
  # With phi > 0, IV2SLS2 should be present
  expect_true("IV2SLS2" %in% sim$summary$method)
})

# ═══════════════════════════════════════════════════════════════
# sweep_mediation_param unchanged
# ═══════════════════════════════════════════════════════════════

test_that("sweep_mediation_param returns summary across grid", {
  sw <- sweep_mediation_param("conf_str", c(0.2, 0.8),
                              n_iter = 3, n_samples = 200, n_features = 5,
                              mo_confounding = 0.8, phi = 0.8, base_seed = 42)
  expect_true(is.data.frame(sw$summary))
  expect_true("param" %in% names(sw$summary))
  expect_true("param_value" %in% names(sw$summary))
  expect_true("NDE_bias" %in% names(sw$summary))
  expect_equal(length(unique(sw$summary$param_value)), 2)
})

# ═══════════════════════════════════════════════════════════════
# gan_sensitivity unchanged
# ═══════════════════════════════════════════════════════════════

test_that("gan_sensitivity returns summary with 5 methods", {
  set.seed(42)
  sens <- gan_sensitivity(NULL,
                          conf_grid = 0.8, coverage_grid = 0.7, k_grid = 1,
                          n_iter = 3, n_samples = 100, n_features = 5,
                          base_seed = 42)
  expect_true(is.data.frame(sens$summary))
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC") %in%
                  sens$summary$method))
  expect_true("conf_strength" %in% names(sens$summary))
  expect_true("coverage" %in% names(sens$summary))
})

# ═══════════════════════════════════════════════════════════════
# gan_mediation_sensitivity unchanged
# ═══════════════════════════════════════════════════════════════

test_that("gan_mediation_sensitivity returns summary with NDE/NIE bias", {
  set.seed(42)
  sens <- gan_mediation_sensitivity(NULL,
                                    conf_grid = 0.8, coverage_grid = 0.7,
                                    k_grid = 1, mo_confounding = 0.8, phi = 0.8,
                                    n_iter = 3, n_samples = 200, n_features = 5,
                                    base_seed = 42)
  expect_true(is.data.frame(sens$summary))
  expect_true("NDE_bias" %in% names(sens$summary))
  expect_true("NIE_bias" %in% names(sens$summary))
  expect_true("IV2SLS2" %in% sens$summary$method)
})

# ═══════════════════════════════════════════════════════════════
# recommend_estimator unchanged
# ═══════════════════════════════════════════════════════════════

test_that("recommend_estimator returns overall and per-scenario", {
  sens <- gan_sensitivity(NULL,
                          conf_grid = c(0.2, 0.8), coverage_grid = 0.7,
                          k_grid = 1, n_iter = 3, n_samples = 100,
                          n_features = 5, base_seed = 42)
  rec <- recommend_estimator(sens)
  expect_true(!is.null(rec$overall))
  expect_true(is.data.frame(rec$per_scenario))
  expect_true("best_method" %in% names(rec$per_scenario))
})

# ═══════════════════════════════════════════════════════════════
# Null simulation functions unchanged
# ═══════════════════════════════════════════════════════════════

test_that("run_null_sim returns Type I error rates", {
  null <- run_null_sim(n_iter = 5, n_samples = 100, n_features = 5,
                       base_seed = 42)
  expect_true(is.data.frame(null$rates))
  expect_true("method" %in% names(null$rates))
  expect_true("type1_error" %in% names(null$rates))
})

test_that("run_null_mediation_sim returns NDE/NIE Type I rates", {
  null <- run_null_mediation_sim(n_iter = 5, n_samples = 200, n_features = 5,
                                 mo_confounding = 0.8, phi = 0.8, base_seed = 42)
  expect_true(is.data.frame(null$rates))
  expect_true("NDE_type1" %in% names(null$rates))
  expect_true("NIE_type1" %in% names(null$rates))
})

# ═══════════════════════════════════════════════════════════════
# generate_toy_data unchanged (gamma_G default)
# ═══════════════════════════════════════════════════════════════

test_that("generate_toy_data with default gamma_G produces valid data", {
  set.seed(42)
  dat <- generate_toy_data(n = 100, n_features = 5, mo_confounding = 0.8,
                           phi = 0.8, seed = 42)
  expect_equal(nrow(dat$Y), 100)
  expect_equal(ncol(dat$Y), 5)
  expect_true(!is.null(dat$G))
  expect_true(!is.null(dat$Gm))
  expect_true(!is.null(dat$W))
})

test_that("generate_toy_data gamma_G=0.6 identical to default", {
  set.seed(42)
  d1 <- generate_toy_data(n = 50, n_features = 3, mo_confounding = 0.8,
                          phi = 0.8, seed = 42)
  set.seed(42)
  d2 <- generate_toy_data(n = 50, n_features = 3, mo_confounding = 0.8,
                          phi = 0.8, gamma_G = 0.6, seed = 42)
  expect_identical(d1$Z, d2$Z)
  expect_identical(d1$Y, d2$Y)
})
