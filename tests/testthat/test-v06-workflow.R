# Tests for the model selection workflow:
# diagnose, estimate, recommend, sensitivity, prospect

# ── Helper: create a small iconic_data from run_single_iteration ──
.make_test_data <- function(n = 100, n_features = 5, mediation = TRUE,
                            phi = 0.8, seed = 42) {
  dat <- run_single_iteration(
    n_synthetic_samples = n, n_features = n_features,
    mo_confounding = if (mediation) 0.8 else 0,
    phi = if (mediation) phi else 0,
    seed = seed)
  if (mediation) {
    iconic_data(X = dat$X, Y = t(dat$Y), M = dat$M, W = t(dat$W),
                G = dat$G[, 1], Gm = dat$Gm, covariates = dat$synthetic_data)
  } else {
    iconic_data(X = dat$X, Y = t(dat$Y), W = t(dat$W),
                G = dat$G[, 1], covariates = dat$synthetic_data)
  }
}

# ═══════════════════════════════════════════════════════════════
# iconic_estimate
# ═══════════════════════════════════════════════════════════════

test_that("iconic_estimate total-effect matches simulation pipeline", {
  set.seed(42)
  dat <- run_single_iteration(n_synthetic_samples = 100, n_features = 5, seed = 42)
  # scale = FALSE so idata is on the same (raw) scale as the reference
  # simulation pipeline analyze_methods_robust(dat), which runs unscaled.
  idata <- iconic_data(X = dat$X, Y = t(dat$Y), W = t(dat$W),
                       G = dat$G[, 1], covariates = dat$synthetic_data,
                       scale = FALSE)
  sim <- analyze_methods_robust(dat)
  real <- iconic_estimate(idata)
  for (m in unique(sim$method)) {
    s <- sim$beta[sim$method == m]
    r <- real$beta[real$method == m]
    expect_lt(max(abs(s - r)), 1e-10)
  }
})

test_that("iconic_estimate mediation matches simulation pipeline", {
  set.seed(42)
  dat <- run_single_iteration(n_synthetic_samples = 100, n_features = 5,
                              mo_confounding = 0.8, phi = 0.8, seed = 42)
  idata <- iconic_data(X = dat$X, Y = t(dat$Y), M = dat$M, W = t(dat$W),
                       G = dat$G[, 1], Gm = dat$Gm,
                       covariates = dat$synthetic_data, scale = FALSE)
  sim <- analyze_mediation_robust(dat)
  real <- iconic_estimate(idata)
  for (m in unique(sim$method)) {
    s <- sim$NDE[sim$method == m]
    r <- real$NDE[real$method == m]
    expect_lt(max(abs(s - r)), 1e-10)
    s <- sim$NIE[sim$method == m]
    r <- real$NIE[real$method == m]
    expect_lt(max(abs(s - r)), 1e-10)
  }
})

test_that("iconic_estimate with diagnosis gives same as without", {
  idata <- .make_test_data()
  diag <- iconic_diagnose(idata)
  est1 <- iconic_estimate(idata)
  est2 <- iconic_estimate(idata, diagnosis = diag)
  expect_true(setequal(est1$method, est2$method))
  expect_lt(max(abs(est1$beta - est2$beta), na.rm = TRUE), 1e-10)})

test_that("iconic_estimate respects methods override", {
  idata <- .make_test_data(mediation = FALSE)
  est <- iconic_estimate(idata, methods = c("UNADJ", "IV2SLS"))
  expect_true(all(est$method %in% c("UNADJ", "IV2SLS")))
})

test_that("iconic_estimate rejects non-iconic_data", {
  expect_error(iconic_estimate(list()), "iconic_data")
})

# ═══════════════════════════════════════════════════════════════
# iconic_diagnose
# ═══════════════════════════════════════════════════════════════

test_that("iconic_diagnose returns eligibility table", {
  idata <- .make_test_data()
  diag <- iconic_diagnose(idata)
  expect_s3_class(diag, "iconic_diagnosis")
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
                     "IV2SLS2", "PGC2", "PGC2Gm") %in%
                   diag$eligibility$estimator))
  expect_true(all(c("eligible", "reason") %in% names(diag$eligibility)))
})

test_that("iconic_diagnose bare data only UNADJ eligible", {
  bare <- iconic_data(X = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50))
  diag <- iconic_diagnose(bare)
  elig <- diag$eligibility
  expect_true(elig$eligible[elig$estimator == "UNADJ"])
  expect_false(any(elig$eligible[elig$estimator != "UNADJ"]))
})

test_that("iconic_diagnose reports instrument strength", {
  idata <- .make_test_data()
  diag <- iconic_diagnose(idata)
  expect_true(!is.na(diag$instrument_strength$F_G))
  expect_true(!is.na(diag$instrument_strength$F_Gm))
  expect_false(diag$instrument_strength$weak_G)
})

test_that("iconic_diagnose total-effect excludes mediation estimators", {
  idata <- .make_test_data(mediation = FALSE)
  diag <- iconic_diagnose(idata)
  elig <- diag$eligibility
  expect_false(elig$eligible[elig$estimator == "IV2SLS2"])
  expect_false(elig$eligible[elig$estimator == "PGC2"])
  expect_false(elig$eligible[elig$estimator == "PGC2Gm"])
})

test_that("print.iconic_diagnosis produces output", {
  idata <- .make_test_data()
  diag <- iconic_diagnose(idata)
  out <- capture.output(print(diag))
  expect_true(any(grepl("iconic_diagnosis", out)))
})

# ═══════════════════════════════════════════════════════════════
# iconic_recommend
# ═══════════════════════════════════════════════════════════════

test_that("iconic_recommend returns ranking with tiers", {
  idata <- .make_test_data()
  diag <- iconic_diagnose(idata)
  est <- iconic_estimate(idata, diagnosis = diag)
  # auto_sensitivity = FALSE keeps this eligibility-ranking test fast and
  # independent of the torch backend (no GAN auto-run).
  rec <- iconic_recommend(idata, diagnosis = diag, estimate = est,
                          auto_sensitivity = FALSE)
  expect_s3_class(rec, "iconic_recommendation")
  expect_true(!is.na(rec$recommended))
})

test_that("iconic_recommend bare data recommends UNADJ", {
  bare <- iconic_data(X = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50))
  rec <- iconic_recommend(bare, auto_sensitivity = FALSE)
  expect_equal(rec$recommended, "UNADJ")
})

test_that("iconic_recommend mediation recommends an eligible estimator", {
  idata <- .make_test_data()
  diag <- iconic_diagnose(idata)
  est <- iconic_estimate(idata, diagnosis = diag)
  rec <- iconic_recommend(idata, diagnosis = diag, estimate = est,
                          auto_sensitivity = FALSE)
  # no tier system: recommendation is the top eligible estimator
  expect_true(rec$recommended %in% diag$eligibility$estimator[diag$eligibility$eligible])
})

test_that("iconic_recommend total-effect recommends an eligible estimator", {
  idata <- .make_test_data(mediation = FALSE)
  diag <- iconic_diagnose(idata)
  est <- iconic_estimate(idata, diagnosis = diag)
  rec <- iconic_recommend(idata, diagnosis = diag, estimate = est,
                          auto_sensitivity = FALSE)
  expect_true(rec$recommended %in% diag$eligibility$estimator[diag$eligibility$eligible])
})

test_that("print.iconic_recommendation produces output", {
  idata <- .make_test_data()
  rec <- iconic_recommend(idata, auto_sensitivity = FALSE)
  out <- capture.output(print(rec))
  expect_true(any(grepl("iconic_recommendation", out)))
})

test_that("iconic_recommend auto-runs sensitivity when not supplied", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  diag <- iconic_diagnose(idata)
  # auto_sensitivity = TRUE (default) with a tiny grid: the recommendation
  # should be robustness-based (robustness columns populated), not
  # eligibility-only.
  rec <- iconic_recommend(idata, diagnosis = diag,
                          rho_G1_grid = c(0, 0.5), rho_G2_grid = c(0, 0.5),
                          omega_1 = c(0.3, 1.0), omega_2 = c(0.3, 1.0),
                          n_iter_sens = 2, gan_epochs = 5, n_cores = 1)
  expect_s3_class(rec, "iconic_recommendation")
  expect_true("robustness_NDE" %in% names(rec$ranking))
  expect_true(!is.na(rec$recommended))
})

# ═══════════════════════════════════════════════════════════════
# iconic_sensitivity
# ═══════════════════════════════════════════════════════════════

test_that("iconic_sensitivity returns surface and tipping points", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0, 0.3),
                             rho_G2_grid = c(0, 0.3),
                             n_iter = 3, n_features = 3)
  expect_s3_class(sens, "iconic_sensitivity")
  expect_true("NDE_bias" %in% names(sens$surface))
  expect_true("tipped" %in% names(sens$surface))
  expect_true(nrow(sens$tipping_points) > 0)
})

test_that("iconic_sensitivity rejects non-mediation data", {
  idata <- .make_test_data(mediation = FALSE)
  expect_error(iconic_sensitivity(idata), "mediation data")
})

test_that("iconic_sensitivity PGC2Gm robust at origin", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0),
                             rho_G2_grid = c(0),
                             omega_1 = 0.7, omega_2 = 0.7,
                             n_iter = 10, n_features = 3)
  pgc2gm <- sens$surface[sens$surface$method == "PGC2Gm", ]
  expect_lt(abs(pgc2gm$NDE_bias), 0.1)
})

test_that("print.iconic_sensitivity produces output", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3)
  out <- capture.output(print(sens))
  expect_true(any(grepl("iconic_sensitivity", out)))
})

# ═══════════════════════════════════════════════════════════════
# iconic_prospect
# ═══════════════════════════════════════════════════════════════

test_that("iconic_prospect returns strength surface and prospective", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  bare <- iconic_data(X = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50),
                      M = rnorm(50))
  result <- iconic_prospect(bare, gamma_G_grid = c(0.3, 0.6),
                            n_iter = 3, n_features = 3)
  expect_s3_class(result, "iconic_prospect")
  expect_true("gamma_G" %in% names(result$strength_surface))
  expect_true("NDE_bias" %in% names(result$prospective))
})

test_that("iconic_prospect rejects non-mediation data", {
  bare <- iconic_data(X = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50))
  expect_error(iconic_prospect(bare), "mediation data")
})

test_that("iconic_prospect UNADJ bias exceeds IV2SLS2 bias", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  bare <- iconic_data(X = rnorm(80), Y = matrix(rnorm(80 * 3), 3, 80),
                      M = rnorm(80))
  result <- iconic_prospect(bare, gamma_G_grid = c(0.6),
                            n_iter = 5, n_features = 3)
  unadj <- result$prospective$NDE_bias[result$prospective$method == "UNADJ"]
  iv2sls2 <- result$prospective$NDE_bias[result$prospective$method == "IV2SLS2"]
  expect_gt(abs(unadj), abs(iv2sls2))
})

test_that("print.iconic_prospect produces output", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  bare <- iconic_data(X = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50),
                      M = rnorm(50))
  result <- iconic_prospect(bare, gamma_G_grid = c(0.6),
                            n_iter = 2, n_features = 3)
  out <- capture.output(print(result))
  expect_true(any(grepl("iconic_prospect", out)))
})

# ═══════════════════════════════════════════════════════════════
# gamma_G backward compatibility
# ═══════════════════════════════════════════════════════════════

test_that("generate_toy_data is exported and works without ::: prefix", {
  dat <- generate_toy_data(n = 100, n_features = 5, seed = 42)
  expect_true(is.matrix(dat$G))
  expect_equal(nrow(dat$G), 100)
  expect_equal(ncol(dat$G), 5)
})

test_that("generate_toy_data mediation data works with iconic_data (no W workaround)", {
  dat <- generate_toy_data(n = 100, n_features = 5, phi = 0.8,
                           mo_confounding = 0.8,
                           omega_1 = 0.7, omega_2 = 0.7, seed = 42)
  # Pass only W1/W2 — iconic_data should derive W internally
  idata <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                       G = dat$G[, 1], Gm = dat$Gm,
                       W1 = dat$W1, W2 = dat$W2)
  expect_true(idata$has_nc)
  expect_true(idata$has_path_nc)
  expect_false(is.null(idata$W))
  # Diagnose should work without crashing
  diag <- iconic_diagnose(idata, k = 1)
  expect_s3_class(diag, "iconic_diagnosis")
})

test_that("generate_toy_data mediation data works with iconic_data (matrix G)", {
  dat <- generate_toy_data(n = 100, n_features = 5, phi = 0.8,
                           mo_confounding = 0.8,
                           omega_1 = 0.7, omega_2 = 0.7, seed = 42)
  # Pass matrix G directly — iconic_data should extract first column
  idata <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                       G = dat$G, Gm = dat$Gm,
                       W1 = dat$W1, W2 = dat$W2)
  expect_true(idata$has_instrument)
  expect_length(idata$G, 100)
  expect_equal(idata$G, as.numeric(dat$G[, 1]))
})

# ═══════════════════════════════════════════════════════════════
# gamma_G backward compatibility
# ═══════════════════════════════════════════════════════════════

test_that("gamma_G default 0.6 is backward-compatible", {
  set.seed(42)
  d1 <- run_single_iteration(n_synthetic_samples = 50, n_features = 3, seed = 42)
  set.seed(42)
  d2 <- run_single_iteration(n_synthetic_samples = 50, n_features = 3,
                             gamma_G = 0.6, seed = 42)
  expect_identical(d1$X, d2$X)
  expect_identical(d1$Y, d2$Y)
})

test_that("gamma_G affects instrument strength", {
  set.seed(42)
  dw <- run_single_iteration(n_synthetic_samples = 100, n_features = 3,
                             mo_confounding = 0.8, phi = 0.8,
                             gamma_G = 0.2, seed = 42)
  set.seed(42)
  ds <- run_single_iteration(n_synthetic_samples = 100, n_features = 3,
                             mo_confounding = 0.8, phi = 0.8,
                             gamma_G = 1.0, seed = 42)
  iw <- iconic_data(X = dw$X, Y = t(dw$Y), M = dw$M, W = t(dw$W),
                    G = dw$G[, 1], Gm = dw$Gm, covariates = dw$synthetic_data)
  is_ <- iconic_data(X = ds$X, Y = t(ds$Y), M = ds$M, W = t(ds$W),
                     G = ds$G[, 1], Gm = ds$Gm, covariates = ds$synthetic_data)
  dw_diag <- iconic_diagnose(iw)
  ds_diag <- iconic_diagnose(is_)
  expect_lt(dw_diag$instrument_strength$F_G, ds_diag$instrument_strength$F_G)
})

# ═══════════════════════════════════════════════════════════════
# GAN auto-training and texture_source
# ═══════════════════════════════════════════════════════════════

test_that("iconic_sensitivity reports texture_source field", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3)
  expect_true("texture_source" %in% names(sens))
  expect_true(sens$texture_source %in% c("auto-trained from data",
                                          "data-attached GAN",
                                          "user-supplied GAN"))
})

test_that("iconic_sensitivity auto-trains GAN when no trained_gan supplied", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3,
                             gan_epochs = 10)
  expect_equal(sens$texture_source, "auto-trained from data")
})

test_that("iconic_sensitivity uses attached GAN from iconic_data", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  # Attach a GAN by training one (light)
  trained <- iconic:::.auto_train_gan(idata, epochs = 5)
  idata2 <- iconic_data(X = idata$X, Y = idata$Y, M = idata$M,
                        W = idata$W, G = idata$G, Gm = idata$Gm,
                        covariates = idata$covariates,
                        trained_gan = trained)
  sens <- iconic_sensitivity(idata2,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3)
  expect_equal(sens$texture_source, "data-attached GAN")
})

test_that("iconic_sensitivity uses supplied trained_gan argument", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  trained <- iconic:::.auto_train_gan(idata, epochs = 5)
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3,
                             trained_gan = trained)
  expect_equal(sens$texture_source, "user-supplied GAN")
})

test_that("iconic_sensitivity confounding=default is backward compatible", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3,
                             confounding = "default")
  expect_null(sens$inferred_confounding)
  # Default omega values are swept over c(0.3, 0.7, 1.0)
  expect_equal(sens$omega_1, c(0.3, 0.7, 1.0))
  expect_equal(sens$omega_2, c(0.3, 0.7, 1.0))
})

test_that("iconic_sensitivity confounding=inferred populates inferred_confounding", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data(n = 80, n_features = 5)
  diag <- iconic_diagnose(idata)
  sens <- iconic_sensitivity(idata, diagnosis = diag,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3,
                             confounding = "inferred",
                             gan_epochs = 10)
  expect_true(!is.null(sens$inferred_confounding))
  expect_s3_class(sens$inferred_confounding, "iconic_confounding")
})

test_that("iconic_sensitivity confounding=manual uses provided values", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3,
                             confounding = "manual",
                             mo_confounding = 0.5,
                             omega_1 = 0.3, omega_2 = 0.4)
  expect_equal(sens$omega_1, 0.3)
  expect_equal(sens$omega_2, 0.4)
})

# ═══════════════════════════════════════════════════════════════
# iconic_prospect GAN and confounding
# ═══════════════════════════════════════════════════════════════

test_that("iconic_prospect reports texture_source field", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  bare <- iconic_data(X = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50),
                      M = rnorm(50))
  result <- iconic_prospect(bare, gamma_G_grid = c(0.6),
                            n_iter = 2, n_features = 3,
                            gan_epochs = 10)
  expect_true("texture_source" %in% names(result))
})

test_that("iconic_prospect confounding=inferred populates inferred_confounding", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data(n = 80, n_features = 5)
  result <- iconic_prospect(idata,
                            gamma_G_grid = c(0.6),
                            n_iter = 2, n_features = 3,
                            confounding = "inferred",
                            gan_epochs = 10)
  expect_true(!is.null(result$inferred_confounding))
  expect_s3_class(result$inferred_confounding, "iconic_confounding")
})

# ═══════════════════════════════════════════════════════════════
# iconic_data trained_gan slot
# ═══════════════════════════════════════════════════════════════

test_that("iconic_data stores trained_gan slot", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  trained <- iconic:::.auto_train_gan(idata, epochs = 5)
  idata2 <- iconic_data(X = idata$X, Y = idata$Y, M = idata$M,
                        W = idata$W, G = idata$G, Gm = idata$Gm,
                        covariates = idata$covariates,
                        trained_gan = trained)
  expect_false(is.null(idata2$trained_gan))
})

test_that("iconic_data print reports GAN attachment", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_test_data()
  trained <- iconic:::.auto_train_gan(idata, epochs = 5)
  idata2 <- iconic_data(X = idata$X, Y = idata$Y, M = idata$M,
                        W = idata$W, G = idata$G, Gm = idata$Gm,
                        covariates = idata$covariates,
                        trained_gan = trained)
  out <- capture.output(print(idata2))
  expect_true(any(grepl("GAN", out)))
})
