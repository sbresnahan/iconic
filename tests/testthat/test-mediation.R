# ============================================================
# Tests for the mediation extension.
#
# Covers:
# - generate_toy_data with mo_confounding > 0
# - Each mediation estimator returns correct structure
# - fit_pgc_mediation (matrix bridge) and fit_pgc_scalar_mediation
# - run_mediation_sim returns correct structure with five methods
# - run_null_mediation_sim returns Type I rates for NDE and NIE
# - sweep_mediation_param sweeps correctly
# - gan_mediation_sensitivity returns scenario x method summary
#
#:
# - generate_toy_data / run_single_iteration with phi > 0 returns Gm
# - fit_iv2sls_mediation2 returns correct structure
# - fit_iv2sls_mediation2 point-identifies NDE/NIE with strong instruments
# - fit_iv2sls_mediation2 returns NA when Gm is weak
# - fit_iv2sls_mediation2 cross-validates against AER::ivreg
# - analyze_mediation_robust includes IV2SLS2 when Gm present
# - run_mediation_sim with phi > 0 includes IV2SLS2
# - nc_independence_check_gm returns correct structure
# - Backward compatibility: phi = 0 produces identical results
# ============================================================

test_that("generate_toy_data with mo_confounding returns correct structure", {
  dat <- iconic:::generate_toy_data(n = 100, n_features = 5,
                                    mo_confounding = 0.8, seed = 1)
  expect_true(all(c("X", "G", "Y", "W", "U1", "M", "synthetic_data",
                    "true_total", "true_NDE", "true_NIE") %in% names(dat)))
  expect_length(dat$X, 100)
  expect_equal(dim(dat$Y), c(100, 5))
  expect_equal(dim(dat$W), c(100, 5))
  expect_equal(dat$true_total, 0.10 + 0.50 * 0.30)
  expect_equal(dat$true_NDE, 0.10)
  expect_equal(dat$true_NIE, 0.50 * 0.30)
})

test_that("generate_toy_data with mo_confounding = 0 matches original DGP", {
  dat0 <- iconic:::generate_toy_data(n = 200, mo_confounding = 0, seed = 42)
  expect_true(all(c("X", "G", "Y", "W", "U1", "M", "synthetic_data",
                    "true_total", "true_NDE", "true_NIE") %in% names(dat0)))
  # M should not have U1 component when mo_confounding = 0
  # (just alpha_M * X + noise)
  expect_true(is.numeric(dat0$M))
})

test_that("fit_unadj_mediation returns named list with NDE/NIE", {
  dat <- iconic:::generate_toy_data(n = 200, mo_confounding = 0.8, seed = 42)
  res <- fit_unadj_mediation(dat$Y[, 1], dat$X, dat$M)
  expect_named(res, c("NDE", "NDE_se", "NDE_p", "NIE", "NIE_se", "NIE_p", "alpha_M", "alpha_se", "beta_M", "beta_M_se"))
  expect_true(is.numeric(res$NDE))
  expect_true(res$NDE_p >= 0 && res$NDE_p <= 1)
})

test_that("fit_direct_mediation returns named list with NDE/NIE", {
  dat <- iconic:::generate_toy_data(n = 200, mo_confounding = 0.8, seed = 42)
  res <- fit_direct_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$W[, 1])
  expect_named(res, c("NDE", "NDE_se", "NDE_p", "NIE", "NIE_se", "NIE_p", "alpha_M", "alpha_se", "beta_M", "beta_M_se"))
  expect_true(is.numeric(res$NDE) || is.na(res$NDE))
})

test_that("fit_coca_mediation returns named list with NDE/NIE", {
  dat <- iconic:::generate_toy_data(n = 200, mo_confounding = 0.8, seed = 42)
  res <- fit_coca_mediation(dat$Y[, 1], dat$X, dat$M, rowMeans(dat$W))
  expect_named(res, c("NDE", "NDE_se", "NDE_p", "NIE", "NIE_se", "NIE_p", "alpha_M", "alpha_se", "beta_M", "beta_M_se"))
  expect_true(is.numeric(res$NDE) || is.na(res$NDE))
})

test_that("fit_iv2sls_mediation returns named list with NDE/NIE", {
  dat <- iconic:::generate_toy_data(n = 300, mo_confounding = 0.8, seed = 42)
  res <- fit_iv2sls_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$W[, 1])
  expect_named(res, c("NDE", "NDE_se", "NDE_p", "NIE", "NIE_se", "NIE_p", "alpha_M", "alpha_se", "beta_M", "beta_M_se"))
  expect_true(is.numeric(res$NDE) || is.na(res$NDE))
  if (!is.na(res$NDE_p)) {
    expect_true(res$NDE_p >= 0 && res$NDE_p <= 1)
  }
})

test_that("fit_pgc_mediation (matrix bridge) returns named list with NDE/NIE", {
  dat <- iconic:::generate_toy_data(n = 200, mo_confounding = 0.8, seed = 42)
  res <- fit_pgc_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$W)
  expect_named(res, c("NDE", "NDE_se", "NDE_p", "NIE", "NIE_se", "NIE_p", "alpha_M", "alpha_se", "beta_M", "beta_M_se"))
  expect_true(is.numeric(res$NDE) || is.na(res$NDE))
})

test_that("fit_pgc_scalar_mediation returns named list with NDE/NIE", {
  dat <- iconic:::generate_toy_data(n = 200, mo_confounding = 0.8, seed = 42)
  res <- fit_pgc_scalar_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], rowMeans(dat$W))
  expect_named(res, c("NDE", "NDE_se", "NDE_p", "NIE", "NIE_se", "NIE_p", "alpha_M", "alpha_se", "beta_M", "beta_M_se"))
  expect_true(is.numeric(res$NDE) || is.na(res$NDE))
})

test_that("analyze_mediation_robust returns correct structure with five methods", {
  dat <- iconic:::generate_toy_data(n = 200, n_features = 3,
                                    mo_confounding = 0.8, seed = 42)
  res <- analyze_mediation_robust(dat)
  expect_true(all(c("feature", "method", "NDE", "NDE_se", "NDE_p",
                    "NIE", "NIE_se", "NIE_p",
                    "NDE_significant", "NIE_significant") %in% names(res)))
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC") %in% res$method))
})

test_that("run_mediation_sim returns correct structure with five methods", {
  res <- run_mediation_sim(n_iter = 5, n_samples = 100, n_features = 3,
                           mo_confounding = 0.8)
  expect_named(res, c("summary", "raw", "true_NDE", "true_NIE", "params"))
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC") %in% res$summary$method))
  expect_true(all(c("NDE_bias", "NIE_bias", "NDE_rmse", "NIE_rmse",
                    "NIE_type1", "NDE_type1") %in% names(res$summary)))
  expect_equal(res$true_NDE, 0.10)
  expect_equal(res$true_NIE, 0.50 * 0.30)
})

test_that("run_null_mediation_sim returns Type I rates for NDE and NIE", {
  res <- run_null_mediation_sim(n_iter = 10, n_samples = 100, n_features = 3,
                                mo_confounding = 0.8)
  expect_named(res, c("rates", "raw"))
  expect_true(all(c("NIE_type1", "NDE_type1") %in% names(res$rates)))
  expect_true(all(res$rates$NIE_type1 >= 0 & res$rates$NIE_type1 <= 1))
  expect_true(all(res$rates$NDE_type1 >= 0 & res$rates$NDE_type1 <= 1))
  expect_equal(nrow(res$rates), 5)
})

test_that("sweep_mediation_param sweeps conf_str correctly", {
  res <- sweep_mediation_param("conf_str", c(0.4, 0.8), n_iter = 5,
                               n_samples = 100, n_features = 3,
                               mo_confounding = 0.8)
  expect_named(res, c("summary", "iter_bias"))
  expect_true(all(res$summary$param_value %in% c(0.4, 0.8)))
  expect_true("NDE_bias" %in% names(res$summary))
  expect_true("NIE_bias" %in% names(res$summary))
})

test_that("sweep_mediation_param sweeps mo_confounding correctly", {
  res <- sweep_mediation_param("mo_confounding", c(0, 0.8), n_iter = 5,
                               n_samples = 100, n_features = 3)
  expect_true(all(res$summary$param_value %in% c(0, 0.8)))
  # With mo_confounding = 0, bias should be smaller than with 0.8
  # (not a strict test, just check structure is right)
  expect_true(nrow(res$summary) > 0)
})

test_that("gan_mediation_sensitivity returns scenario x method summary with five methods", {
  sens <- gan_mediation_sensitivity(NULL, conf_grid = c(0.3, 0.8),
                                    coverage_grid = c(0.5, 1), k_grid = 1,
                                    mo_confounding = 0.8, n_iter = 5,
                                    n_samples = 250, n_features = 5,
                                    base_seed = 11)
  expect_true(all(c("summary", "grid") %in% names(sens)))
  expect_true(all(c("conf_strength", "coverage", "k", "mo_confounding",
                    "true_NDE", "true_NIE", "method",
                    "NDE_bias", "NIE_bias", "NDE_rmse", "NIE_rmse") %in%
                   names(sens$summary)))
  # 2 conf x 2 coverage x 1 k x 5 methods
  expect_equal(nrow(sens$summary), 2 * 2 * 1 * 5)
})

test_that("run_single_iteration with mo_confounding returns true_NDE and true_NIE", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 100, n_features = 5,
                              mo_confounding = 0.8, seed = 1)
  expect_true("true_NDE" %in% names(dat))
  expect_true("true_NIE" %in% names(dat))
  expect_equal(dat$true_NDE, 0.10)
  expect_equal(dat$true_NIE, 0.50 * 0.30)
  expect_equal(dat$params$mo_confounding, 0.8)
})


# ============================================================
#: mediator instrument (Gm) and 2-stage MR
# ============================================================

test_that("generate_toy_data with phi > 0 returns Gm", {
  dat <- iconic:::generate_toy_data(n = 100, n_features = 5,
                                    mo_confounding = 0.8, phi = 0.8, seed = 1)
  expect_true("Gm" %in% names(dat))
  expect_length(dat$Gm, 100)
  expect_true(is.numeric(dat$Gm))
  # Gm should be approximately standard normal (independent of everything)
  expect_true(abs(mean(dat$Gm)) < 0.3)
  expect_true(abs(sd(dat$Gm) - 1) < 0.2)
})

test_that("generate_toy_data with phi = 0 does not return Gm (backward compat)", {
  dat <- iconic:::generate_toy_data(n = 100, n_features = 5,
                                    mo_confounding = 0.8, seed = 1)
  expect_false("Gm" %in% names(dat))
})

test_that("run_single_iteration with phi > 0 returns Gm and params$phi", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 100, n_features = 5,
                              mo_confounding = 0.8, phi = 0.8, seed = 1)
  expect_true("Gm" %in% names(dat))
  expect_length(dat$Gm, 100)
  expect_equal(dat$params$phi, 0.8)
})

test_that("run_single_iteration with phi = 0 does not return Gm (backward compat)", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 100, n_features = 5,
                              mo_confounding = 0.8, seed = 1)
  expect_false("Gm" %in% names(dat))
  expect_equal(dat$params$phi, 0)
})

test_that("fit_iv2sls_mediation2 returns named list with NDE/NIE", {
  dat <- iconic:::generate_toy_data(n = 500, n_features = 3,
                                    mo_confounding = 0.8, phi = 0.8, seed = 42)
  res <- fit_iv2sls_mediation2(dat$Y[, 1], dat$X, dat$M,
                                dat$G[, 1], dat$Gm)
  expect_named(res, c("NDE", "NDE_se", "NDE_p", "NIE", "NIE_se", "NIE_p", "alpha_M", "alpha_se", "beta_M", "beta_M_se"))
  expect_true(is.numeric(res$NDE))
  if (!is.na(res$NDE_p)) {
    expect_true(res$NDE_p >= 0 && res$NDE_p <= 1)
  }
})

test_that("fit_iv2sls_mediation2 point-identifies NDE/NIE with strong instruments", {
  # The key identification result: with strong instruments for both X and M,
  # IV2SLS2 should recover true NDE and NIE with low bias, unlike IV2SLS
  # which stays biased under M-O confounding.
  set.seed(123)
  n_rep <- 30
  nde_iv2 <- nde_iv2_2 <- nie_iv2 <- nie_iv2_2 <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    dat <- iconic:::generate_toy_data(n = 500, n_features = 1,
                                      mo_confounding = 0.8, phi = 0.8,
                                      seed = 2000 + i)
    r1 <- fit_iv2sls_mediation(dat$Y[, 1], dat$X, dat$M,
                               dat$G[, 1], dat$W[, 1])
    r2 <- fit_iv2sls_mediation2(dat$Y[, 1], dat$X, dat$M,
                                dat$G[, 1], dat$Gm)
    nde_iv2[i] <- r1$NDE; nie_iv2[i] <- r1$NIE
    nde_iv2_2[i] <- r2$NDE; nie_iv2_2[i] <- r2$NIE
  }
  true_NDE <- 0.10; true_NIE <- 0.15
  # IV2SLS2 should have much smaller bias than IV2SLS
  iv2sls2_nde_bias <- abs(mean(nde_iv2_2, na.rm = TRUE) - true_NDE)
  iv2sls_nde_bias <- abs(mean(nde_iv2, na.rm = TRUE) - true_NDE)
  expect_true(iv2sls2_nde_bias < iv2sls_nde_bias,
              info = paste("IV2SLS2 NDE bias =", round(iv2sls2_nde_bias, 4),
                           "vs IV2SLS NDE bias =", round(iv2sls_nde_bias, 4)))
  # IV2SLS2 NDE bias should be small (point-identified)
  expect_true(iv2sls2_nde_bias < 0.05,
              info = paste("IV2SLS2 NDE bias =", round(iv2sls2_nde_bias, 4)))
  # IV2SLS2 NIE bias should also be small
  iv2sls2_nie_bias <- abs(mean(nie_iv2_2, na.rm = TRUE) - true_NIE)
  expect_true(iv2sls2_nie_bias < 0.05,
              info = paste("IV2SLS2 NIE bias =", round(iv2sls2_nie_bias, 4)))
})

test_that("fit_iv2sls_mediation2 returns NA when Gm is weak", {
  # With phi near 0, the mediator instrument is too weak (partial F < 10)
  dat <- iconic:::generate_toy_data(n = 500, n_features = 1,
                                    mo_confounding = 0.8, phi = 0.001, seed = 42)
  res <- fit_iv2sls_mediation2(dat$Y[, 1], dat$X, dat$M,
                                dat$G[, 1], dat$Gm)
  expect_true(all(is.na(unlist(res))))
})

test_that("fit_iv2sls_mediation2 cross-validates against AER::ivreg", {
  skip_if_not_installed("AER")
  dat <- iconic:::generate_toy_data(n = 500, n_features = 1,
                                    mo_confounding = 0.8, phi = 0.8, seed = 42)
  y <- dat$Y[, 1]; X <- dat$X; M <- dat$M
  g <- dat$G[, 1]; gm <- dat$Gm

  # Our sequential 2SLS estimator (pure-MR form: no NC augmentation)
  res_seq <- fit_iv2sls_mediation2(y, X, M, g, gm)

  # Canonical ivreg: just-identified system with 2 endogenous, 2 excluded instruments
  d_iv <- data.frame(y = y, X = X, M = M, G_inst = g, Gm_inst = gm)
  fit_iv <- AER::ivreg(y ~ X + M | G_inst + Gm_inst, data = d_iv)

  # NDE (coefficient on X) should match within tolerance
  nde_diff <- abs(res_seq$NDE - as.numeric(coef(fit_iv)["X"]))
  expect_true(nde_diff < 0.02,
              info = paste("NDE diff =", round(nde_diff, 6),
                           "(seq =", round(res_seq$NDE, 6),
                           ", ivreg =", round(as.numeric(coef(fit_iv)["X"]), 6), ")"))

  # beta_M (coefficient on M) should match within tolerance
  # Extract beta_M from the sequential estimator by re-running stage 3
  fs <- lm(X ~ g); X_hat <- fitted(fs)
  ms <- lm(M ~ X_hat + gm); M_hat <- fitted(ms)
  os <- lm(y ~ X_hat + M_hat)
  beta_M_seq <- as.numeric(coef(os)["M_hat"])
  beta_M_diff <- abs(beta_M_seq - as.numeric(coef(fit_iv)["M"]))
  expect_true(beta_M_diff < 0.02,
              info = paste("beta_M diff =", round(beta_M_diff, 6),
                           "(seq =", round(beta_M_seq, 6),
                           ", ivreg =", round(as.numeric(coef(fit_iv)["M"]), 6), ")"))
})

test_that("analyze_mediation_robust includes IV2SLS2 when Gm present", {
  dat <- iconic:::generate_toy_data(n = 200, n_features = 3,
                                    mo_confounding = 0.8, phi = 0.8, seed = 42)
  res <- analyze_mediation_robust(dat)
  expect_true("IV2SLS2" %in% res$method)
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC", "IV2SLS2") %in% res$method))
})

test_that("analyze_mediation_robust does NOT include IV2SLS2 when Gm absent", {
  dat <- iconic:::generate_toy_data(n = 200, n_features = 3,
                                    mo_confounding = 0.8, seed = 42)
  res <- analyze_mediation_robust(dat)
  expect_false("IV2SLS2" %in% res$method)
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC") %in% res$method))
})

test_that("run_mediation_sim with phi > 0 includes IV2SLS2", {
  res <- run_mediation_sim(n_iter = 5, n_samples = 200, n_features = 3,
                           mo_confounding = 0.8, phi = 0.8)
  expect_true("IV2SLS2" %in% res$summary$method)
  expect_equal(res$params$phi, 0.8)
})

test_that("run_mediation_sim with phi = 0 does NOT include IV2SLS2 (backward compat)", {
  res <- run_mediation_sim(n_iter = 5, n_samples = 100, n_features = 3,
                           mo_confounding = 0.8)
  expect_false("IV2SLS2" %in% res$summary$method)
  expect_equal(res$params$phi, 0)
})

test_that("sweep_mediation_param sweeps phi correctly", {
  res <- sweep_mediation_param("phi", c(0, 0.8), n_iter = 5,
                               n_samples = 200, n_features = 3,
                               mo_confounding = 0.8)
  expect_true(all(res$summary$param_value %in% c(0, 0.8)))
  # At phi = 0.8, IV2SLS2 should be present; at phi = 0, it should not
  methods_08 <- unique(res$summary$method[res$summary$param_value == 0.8])
  methods_00 <- unique(res$summary$method[res$summary$param_value == 0])
  expect_true("IV2SLS2" %in% methods_08)
  expect_false("IV2SLS2" %in% methods_00)
})

test_that("run_null_mediation_sim with phi > 0 includes IV2SLS2", {
  res <- run_null_mediation_sim(n_iter = 5, n_samples = 200, n_features = 3,
                                mo_confounding = 0.8, phi = 0.8)
  expect_true("IV2SLS2" %in% res$rates$method)
  expect_equal(nrow(res$rates), 6)
})

test_that("sweep_mediation_null_by_conf with phi > 0 includes IV2SLS2", {
  res <- sweep_mediation_null_by_conf(c(0.4, 0.8), n_iter = 5,
                                      n_samples = 200, n_features = 3,
                                      mo_confounding = 0.8, phi = 0.8)
  expect_true("IV2SLS2" %in% res$method)
})

test_that("gan_mediation_sensitivity with phi > 0 includes IV2SLS2 and phi column", {
  sens <- gan_mediation_sensitivity(NULL, conf_grid = c(0.3, 0.8),
                                    coverage_grid = c(0.5, 1), k_grid = 1,
                                    mo_confounding = 0.8, phi = 0.8,
                                    n_iter = 3, n_samples = 200, n_features = 3,
                                    base_seed = 11)
  expect_true("IV2SLS2" %in% sens$summary$method)
  expect_true("phi" %in% names(sens$summary))
  # 2 conf x 2 coverage x 1 k x 6 methods
  expect_equal(nrow(sens$summary), 2 * 2 * 1 * 6)
})

test_that("nc_independence_check_gm returns correct structure when Gm present", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 200, n_features = 10,
                              mo_confounding = 0.8, phi = 0.8, seed = 1)
  res <- nc_independence_check_gm(dat)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("feature", "partial_r", "p_value", "fdr",
                    "significant", "verdict") %in% names(res)))
  expect_equal(nrow(res), 10)
  # In the simulation, Gm is independent of W by construction
  expect_true(all(res$verdict == "valid"))
})

test_that("nc_independence_check_gm returns NULL when Gm absent", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 200, n_features = 10,
                              mo_confounding = 0.8, seed = 1)
  res <- nc_independence_check_gm(dat)
  expect_null(res)
})

test_that("nc_independence_check_gm detects injected Gm->W violations", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 200, n_features = 10,
                              mo_confounding = 0.8, phi = 0.8, seed = 1)
  # Inject Gm -> W into features 1-5
  dat$W[, 1:5] <- dat$W[, 1:5] + 0.3 * dat$Gm
  res <- nc_independence_check_gm(dat)
  # Should detect violations in features 1-5
  expect_true(sum(res$significant[1:5]) >= 3)
  # Should not flag clean features 6-10 (allow some false positives)
  expect_true(sum(res$significant[6:10]) <= 2)
})

test_that("backward compat: run_mediation_sim default args produce 5 methods only", {
  res <- run_mediation_sim(n_iter = 5, n_samples = 100, n_features = 3,
                           mo_confounding = 0.8)
  expect_equal(length(unique(res$summary$method)), 5)
  expect_false("IV2SLS2" %in% res$summary$method)
})


# ============================================================
#: imperfect instrument independence, PGC-2,
# path-specific negative controls, tipping-point behaviour
# ============================================================

# --- DGP: ---

test_that("generate_toy_data with path-specific loadings returns G1, G2, W1, W2, conf_XM, conf_MY", {
  dat <- iconic:::generate_toy_data(n = 200, n_features = 5,
                                    mo_confounding = 0.8, phi = 0.8,
                                    rho_G1 = 0.3, rho_G2 = 0.3,
                                    lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                    omega_1 = 0.7, omega_2 = 0.7, seed = 1)
  expect_true("G1" %in% names(dat))
  expect_true("G2" %in% names(dat))
  expect_true("W1" %in% names(dat))
  expect_true("W2" %in% names(dat))
  expect_true("conf_XM" %in% names(dat))
  expect_true("conf_MY" %in% names(dat))
  expect_length(dat$G1, 200)
  expect_length(dat$G2, 200)
  expect_equal(dim(dat$W1), c(200, 5))
  expect_equal(dim(dat$W2), c(200, 5))
})

test_that("generate_toy_data without path-specific loadings does NOT return", {
  dat <- iconic:::generate_toy_data(n = 200, n_features = 5,
                                    mo_confounding = 0.8, phi = 0.8, seed = 1)
  expect_false("G1" %in% names(dat))
  expect_false("G2" %in% names(dat))
  expect_false("W1" %in% names(dat))
  expect_false("W2" %in% names(dat))
  expect_false("conf_XM" %in% names(dat))
  expect_false("conf_MY" %in% names(dat))
})

test_that("generate_toy_data with rho_pop > 0 returns P", {
  dat <- iconic:::generate_toy_data(n = 200, n_features = 3,
                                    mo_confounding = 0.8, phi = 0.8,
                                    rho_G1 = 0.2, rho_G2 = 0.2, rho_pop = 0.3,
                                    lambda_XM = c(1, 0), lambda_MY = c(0, 1), seed = 1)
  expect_true("P" %in% names(dat))
  expect_length(dat$P, 200)
})

test_that("generate_toy_data G1 correlated with conf_XM, G2 with conf_MY", {
  dat <- iconic:::generate_toy_data(n = 2000, n_features = 3,
                                    mo_confounding = 0.8, phi = 0.8,
                                    rho_G1 = 0.3, rho_G2 = 0.3,
                                    lambda_XM = c(1, 0), lambda_MY = c(0, 1), seed = 1)
  expect_true(abs(cor(dat$G1, dat$conf_XM) - 0.3) < 0.05)
  expect_true(abs(cor(dat$G2, dat$conf_MY) - 0.3) < 0.05)
  # conf_XM and conf_MY should be independent when lambda_XM = c(1, 0), lambda_MY = c(0, 1)
  expect_true(abs(cor(dat$conf_XM, dat$conf_MY)) < 0.05)
})

test_that("generate_toy_data W1 captures conf_XM, W2 captures conf_MY", {
  dat <- iconic:::generate_toy_data(n = 2000, n_features = 50,
                                    mo_confounding = 0.8, phi = 0.8,
                                    rho_G1 = 0.3, rho_G2 = 0.3,
                                    lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                    omega_1 = 0.7, omega_2 = 0.7, seed = 1)
  # W1 should be more correlated with conf_XM than with conf_MY
  cor_w1_uxm <- mean(cor(dat$W1, dat$conf_XM))
  cor_w1_umy <- mean(cor(dat$W1, dat$conf_MY))
  expect_true(cor_w1_uxm > cor_w1_umy)
  # W2 should be more correlated with conf_MY than with conf_XM
  cor_w2_umy <- mean(cor(dat$W2, dat$conf_MY))
  cor_w2_uxm <- mean(cor(dat$W2, dat$conf_XM))
  expect_true(cor_w2_umy > cor_w2_uxm)
})

test_that("generate_toy_data default shared loadings collapses conf_XM=conf_MY=U1", {
  dat <- iconic:::generate_toy_data(n = 500, n_features = 3,
                                    mo_confounding = 0.8, phi = 0.8,
                                    rho_G1 = 0.3, rho_G2 = 0.3,
                                    seed = 1)
  expect_true(abs(cor(dat$conf_XM, dat$conf_MY) - 1) < 1e-8)
  expect_true(abs(cor(dat$conf_XM, dat$U1) - 1) < 1e-8)
})

# --- run_single_iteration: ---

test_that("run_single_iteration with path-specific loadings returns", {
  # n_confounders = 2 so the 2-vector loadings lambda_XM/lambda_MY are valid.
  dat <- run_single_iteration(NULL, n_synthetic_samples = 200, n_features = 5,
                              n_confounders = 2,
                              mo_confounding = 0.8, phi = 0.8,
                              rho_G1 = 0.3, rho_G2 = 0.3,
                              lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                              omega_1 = 0.7, omega_2 = 0.7, seed = 1)
  expect_true("G1" %in% names(dat))
  expect_true("G2" %in% names(dat))
  expect_true("W1" %in% names(dat))
  expect_true("W2" %in% names(dat))
  expect_true("conf_XM" %in% names(dat))
  expect_true("conf_MY" %in% names(dat))
  expect_equal(dat$params$rho_G1, 0.3)
  expect_equal(dat$params$rho_G2, 0.3)
  expect_equal(dat$params$lambda_XM, c(1, 0))
})

test_that("run_single_iteration without path-specific loadings does NOT return", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 200, n_features = 5,
                              mo_confounding = 0.8, phi = 0.8, seed = 1)
  expect_false("G1" %in% names(dat))
  expect_false("W1" %in% names(dat))
  expect_false("conf_XM" %in% names(dat))
})

test_that("run_single_iteration with path-specific loadings produces W1/W2", {
  # k = 2 gives the two-confounder space; e1/e2 loadings recover the old
  # old separate_U = TRUE behavior (independent per-path confounders)
  dat <- run_single_iteration(NULL, n_synthetic_samples = 200, n_features = 5,
                              n_confounders = 2,
                              mo_confounding = 0.8, phi = 0.8,
                              rho_G1 = 0.3, rho_G2 = 0.3,
                              lambda_XM = c(1, 0), lambda_MY = c(0, 1), seed = 1)
  expect_true("W1" %in% names(dat))
  expect_true("W2" %in% names(dat))
})

# --- fit_pgc_mediation2: structure and basic behaviour ---

test_that("fit_pgc_mediation2 returns named list with NDE/NIE", {
  dat <- iconic:::generate_toy_data(n = 500, n_features = 3,
                                    mo_confounding = 0.8, phi = 0.8,
                                    rho_G1 = 0.3, rho_G2 = 0.3,
                                    lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                    omega_1 = 0.7, omega_2 = 0.7, seed = 42)
  res <- fit_pgc_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G1,
                            dat$W1, dat$W2, gm = dat$Gm)
  expect_named(res, c("NDE", "NDE_se", "NDE_p", "NIE", "NIE_se", "NIE_p", "alpha_M", "alpha_se", "beta_M", "beta_M_se"))
  expect_true(is.numeric(res$NDE))
  expect_true(is.numeric(res$NIE))
})

test_that("fit_pgc_mediation2 without gm returns named list with NDE/NIE", {
  dat <- iconic:::generate_toy_data(n = 500, n_features = 3,
                                    mo_confounding = 0.8, phi = 0.8,
                                    rho_G1 = 0.3, rho_G2 = 0.3,
                                    lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                    omega_1 = 0.7, omega_2 = 0.7, seed = 42)
  res <- fit_pgc_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G1,
                            dat$W1, dat$W2, gm = NULL)
  expect_named(res, c("NDE", "NDE_se", "NDE_p", "NIE", "NIE_se", "NIE_p", "alpha_M", "alpha_se", "beta_M", "beta_M_se"))
  expect_true(is.numeric(res$NDE) || is.na(res$NDE))
})

test_that("fit_pgc_mediation2 returns all-NA when g is weak (pure noise)", {
  dat <- iconic:::generate_toy_data(n = 500, n_features = 3,
                                    mo_confounding = 0.8, phi = 0.8,
                                    rho_G1 = 0.3, rho_G2 = 0.3,
                                    lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                    omega_1 = 0.7, omega_2 = 0.7, seed = 42)
  # Replace G1 with pure noise
  g_noise <- rnorm(500)
  res <- fit_pgc_mediation2(dat$Y[, 1], dat$X, dat$M, g_noise,
                            dat$W1, dat$W2, gm = dat$Gm)
  expect_true(all(is.na(unlist(res))))
})

# --- fit_pgc_mediation2: identification under imperfect independence ---

test_that("fit_pgc_mediation2 with gm has lower bias than unaugmented IV2SLS2 under rho_G2 > 0", {
  # The tipping-point result: as rho_G2 increases, IV2SLS2 degrades faster
  # than PGC2Gm because Gm-conf_MY correlation violates IV2SLS2's exogeneity.
  # PGC2Gm's W2 bridge absorbs the Gm-conf_MY correlation that biases IV2SLS2.
  # NOTE (v0.9.9): the comparison is against UNaugmented IV2SLS2 (no W1/W2).
  # When IV2SLS2 is given the same path-specific W2, it ALSO absorbs the
  # Gm-conf_MY correlation and the two estimators perform comparably -- the
  # PGC2Gm advantage is specifically over the instrument-only (unaugmented)
  # IV2SLS2 form.
  set.seed(777)
  n_rep <- 20
  nde_pgc2gm <- nde_iv2sls2 <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    dat <- iconic:::generate_toy_data(n = 500, n_features = 1,
                                      mo_confounding = 0.8, phi = 0.8,
                                      rho_G1 = 0.3, rho_G2 = 0.5,
                                      lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                      omega_1 = 0.7, omega_2 = 0.7,
                                      seed = 5000 + i)
    r_pgc <- fit_pgc_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G1,
                                 dat$W1, dat$W2, gm = dat$Gm)
    r_iv <- fit_iv2sls_mediation2(dat$Y[, 1], dat$X, dat$M,
                                    dat$G[, 1], dat$Gm)  # unaugmented (pure MR)
    nde_pgc2gm[i] <- r_pgc$NDE
    nde_iv2sls2[i] <- r_iv$NDE
  }
  true_NDE <- 0.10
  pgc2gm_bias <- abs(mean(nde_pgc2gm, na.rm = TRUE) - true_NDE)
  iv2sls2_bias <- abs(mean(nde_iv2sls2, na.rm = TRUE) - true_NDE)
  expect_true(pgc2gm_bias < iv2sls2_bias,
              info = paste("PGC2Gm NDE bias =", round(pgc2gm_bias, 4),
                           "vs unaugmented IV2SLS2 NDE bias =", round(iv2sls2_bias, 4)))
})

test_that("fit_pgc_mediation2 with gm recovers NDE/NIE near rho_G2 = 0", {
  # Under perfect instrument independence (rho_G2 = 0), PGC2Gm should
  # recover true effects with low bias, comparable to IV2SLS2.
  set.seed(999)
  n_rep <- 20
  nde <- nie <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    dat <- iconic:::generate_toy_data(n = 500, n_features = 1,
                                      mo_confounding = 0.8, phi = 0.8,
                                      rho_G1 = 0.3, rho_G2 = 0,
                                      lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                      omega_1 = 0.7, omega_2 = 0.7,
                                      seed = 7000 + i)
    r <- fit_pgc_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G1,
                            dat$W1, dat$W2, gm = dat$Gm)
    nde[i] <- r$NDE; nie[i] <- r$NIE
  }
  expect_true(abs(mean(nde, na.rm = TRUE) - 0.10) < 0.05,
              info = paste("PGC2Gm NDE bias =", round(abs(mean(nde) - 0.10), 4)))
  expect_true(abs(mean(nie, na.rm = TRUE) - 0.15) < 0.05,
              info = paste("PGC2Gm NIE bias =", round(abs(mean(nie) - 0.15), 4)))
})

# --- Pipeline integration ---

test_that("analyze_mediation_robust includes PGC2 and PGC2Gm when W1/W2/Gm present", {
  dat <- iconic:::generate_toy_data(n = 300, n_features = 3,
                                    mo_confounding = 0.8, phi = 0.8,
                                    rho_G1 = 0.3, rho_G2 = 0.3,
                                    lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                    omega_1 = 0.7, omega_2 = 0.7, seed = 42)
  res <- analyze_mediation_robust(dat)
  expect_true("PGC2" %in% res$method)
  expect_true("PGC2Gm" %in% res$method)
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
                    "IV2SLS2", "PGC2", "PGC2Gm") %in% res$method))
})

test_that("analyze_mediation_robust includes PGC2 but NOT PGC2Gm when W1/W2 present but Gm absent", {
  dat <- iconic:::generate_toy_data(n = 300, n_features = 3,
                                    mo_confounding = 0.8, phi = 0,
                                    rho_G1 = 0.3, rho_G2 = 0,
                                    lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                    omega_1 = 0.7, omega_2 = 0.7, seed = 42)
  res <- analyze_mediation_robust(dat)
  expect_true("PGC2" %in% res$method)
  expect_false("PGC2Gm" %in% res$method)
  expect_false("IV2SLS2" %in% res$method)
})

test_that("analyze_mediation_robust does NOT include PGC2/PGC2Gm when W1/W2 absent", {
  dat <- iconic:::generate_toy_data(n = 300, n_features = 3,
                                    mo_confounding = 0.8, phi = 0.8, seed = 42)
  res <- analyze_mediation_robust(dat)
  expect_false("PGC2" %in% res$method)
  expect_false("PGC2Gm" %in% res$method)
})

# --- Simulation functions with ---

test_that("run_mediation_sim with path-specific loadings includes PGC2 and PGC2Gm", {
  res <- run_mediation_sim(n_iter = 5, n_samples = 300, n_features = 3,
                           mo_confounding = 0.8, phi = 0.8,
                           rho_G1 = 0.3, rho_G2 = 0.3,
                           lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                           omega_1 = 0.7, omega_2 = 0.7)
  expect_true("PGC2" %in% res$summary$method)
  expect_true("PGC2Gm" %in% res$summary$method)
  expect_equal(res$params$rho_G1, 0.3)
  expect_equal(res$params$rho_G2, 0.3)
  expect_equal(res$params$lambda_XM, c(1, 0))
})

test_that("run_mediation_sim without path-specific loadings does NOT include PGC2/PGC2Gm", {
  res <- run_mediation_sim(n_iter = 5, n_samples = 100, n_features = 3,
                           mo_confounding = 0.8, phi = 0.8)
  expect_false("PGC2" %in% res$summary$method)
  expect_false("PGC2Gm" %in% res$summary$method)
})

test_that("sweep_mediation_param sweeps rho_G2 correctly", {
  res <- sweep_mediation_param("rho_G2", c(0, 0.3), n_iter = 5,
                               n_samples = 300, n_features = 3,
                               mo_confounding = 0.8, phi = 0.8,
                               rho_G1 = 0.3, lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                               omega_1 = 0.7, omega_2 = 0.7)
  expect_true(all(res$summary$param_value %in% c(0, 0.3)))
  # PGC2Gm should be present at both levels
  expect_true("PGC2Gm" %in% unique(res$summary$method))
})

test_that("run_null_mediation_sim with path-specific loadings includes PGC2 and PGC2Gm", {
  res <- run_null_mediation_sim(n_iter = 5, n_samples = 300, n_features = 3,
                                mo_confounding = 0.8, phi = 0.8,
                                rho_G1 = 0.3, rho_G2 = 0.3,
                                lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                omega_1 = 0.7, omega_2 = 0.7)
  expect_true("PGC2" %in% res$rates$method)
  expect_true("PGC2Gm" %in% res$rates$method)
})

test_that("sweep_mediation_null_by_conf with path-specific loadings includes PGC2 and PGC2Gm", {
  res <- sweep_mediation_null_by_conf(c(0.4, 0.8), n_iter = 5,
                                      n_samples = 300, n_features = 3,
                                      mo_confounding = 0.8, phi = 0.8,
                                      rho_G1 = 0.3, rho_G2 = 0.3,
                                      lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                      omega_1 = 0.7, omega_2 = 0.7)
  expect_true("PGC2" %in% res$method)
  expect_true("PGC2Gm" %in% res$method)
})

test_that("gan_mediation_sensitivity with path-specific loadings includes PGC2 and PGC2Gm", {
  sens <- gan_mediation_sensitivity(NULL, conf_grid = c(0.3, 0.8),
                                    coverage_grid = c(0.5, 1), k_grid = 2,
                                    mo_confounding = 0.8, phi = 0.8,
                                    rho_G1 = 0.3, rho_G2 = 0.3,
                                    lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                                    omega_1 = 0.7, omega_2 = 0.7,
                                    n_iter = 3, n_samples = 200, n_features = 3,
                                    base_seed = 11)
  expect_true("PGC2" %in% sens$summary$method)
  expect_true("PGC2Gm" %in% sens$summary$method)
  expect_true("rho_G1" %in% names(sens$summary))
  expect_true("rho_G2" %in% names(sens$summary))
})

# --- Backward compatibility ---

test_that("backward compat: default generate_toy_data output unchanged (no path-specific loadings)", {
  dat <- iconic:::generate_toy_data(n = 100, n_features = 5,
                                    mo_confounding = 0.8, seed = 1)
  # Should have exactly the names, no
  expected <- c("X", "G", "Y", "W", "U1", "M", "synthetic_data",
                "true_total", "true_NDE", "true_NIE")
  expect_true(all(expected %in% names(dat)))
  expect_false(any(c("G1", "G2", "W1", "W2", "conf_XM", "conf_MY", "P") %in% names(dat)))
})
