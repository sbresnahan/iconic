# Tests for infer_confounding(): data-calibrated confounding parameters.

# ── Helper: create a small iconic_data from run_single_iteration ──
.make_conf_data <- function(n = 100, n_features = 10, phi = 0.8, seed = 42) {
  dat <- run_single_iteration(
    n_synthetic_samples = n, n_features = n_features,
    mo_confounding = 0.8, phi = phi, seed = seed)
  iconic_data(Z = dat$Z, Y = t(dat$Y), M = dat$M, W = t(dat$W),
              G = dat$G[, 1], Gm = dat$Gm, covariates = dat$synthetic_data)
}

# ═══════════════════════════════════════════════════════════════
# Return structure
# ═══════════════════════════════════════════════════════════════

test_that("infer_confounding returns all parameter slots", {
  idata <- .make_conf_data()
  diag <- iconic_diagnose(idata)
  est <- iconic_estimate(idata, diagnosis = diag)
  conf <- infer_confounding(idata, diagnosis = diag, estimate = est)
  expect_s3_class(conf, "iconic_confounding")
  expect_true(all(c("conf_strength", "mo_confounding", "omega_1",
                     "omega_2", "k", "unavailable", "warnings") %in%
                   names(conf)))
})

test_that("infer_confounding each parameter has required fields", {
  idata <- .make_conf_data()
  diag <- iconic_diagnose(idata)
  est <- iconic_estimate(idata, diagnosis = diag)
  conf <- infer_confounding(idata, diagnosis = diag, estimate = est)
  for (param in c("conf_strength", "mo_confounding", "omega_1", "omega_2")) {
    p <- conf[[param]]
    expect_true(all(c("estimate", "method", "assumption",
                      "available", "warning") %in% names(p)),
                label = paste("param", param))
  }
  # k has ci too
  expect_true("ci" %in% names(conf$k))
})

# ═══════════════════════════════════════════════════════════════
# conf_strength inference (UNADJ-IV2SLS gap)
# ═══════════════════════════════════════════════════════════════

test_that("conf_strength is available with valid G", {
  idata <- .make_conf_data()
  diag <- iconic_diagnose(idata)
  est <- iconic_estimate(idata, diagnosis = diag)
  conf <- infer_confounding(idata, diagnosis = diag, estimate = est)
  expect_true(conf$conf_strength$available)
  expect_true(grepl("UNADJ-IV2SLS gap", conf$conf_strength$method))
  expect_true(conf$conf_strength$estimate > 0)
})

test_that("conf_strength unavailable without G", {
  dat <- run_single_iteration(n_synthetic_samples = 100, n_features = 10,
                              mo_confounding = 0.8, phi = 0, seed = 42)
  idata <- iconic_data(Z = dat$Z, Y = t(dat$Y), M = dat$M, W = t(dat$W),
                       covariates = dat$synthetic_data)
  conf <- infer_confounding(idata)
  expect_false(conf$conf_strength$available)
  expect_equal(conf$conf_strength$estimate, 0.8)  # default
})

# ═══════════════════════════════════════════════════════════════
# mo_confounding inference (IV2SLS-IV2SLS2 NIE gap)
# ═══════════════════════════════════════════════════════════════

test_that("mo_confounding is available with G + Gm", {
  idata <- .make_conf_data()
  diag <- iconic_diagnose(idata)
  est <- iconic_estimate(idata, diagnosis = diag)
  conf <- infer_confounding(idata, diagnosis = diag, estimate = est)
  expect_true(conf$mo_confounding$available)
  expect_equal(conf$mo_confounding$method, "IV2SLS-IV2SLS2 NIE gap")
})

test_that("mo_confounding unavailable without Gm", {
  dat <- run_single_iteration(n_synthetic_samples = 100, n_features = 10,
                              mo_confounding = 0.8, phi = 0, seed = 42)
  idata <- iconic_data(Z = dat$Z, Y = t(dat$Y), M = dat$M, W = t(dat$W),
                       G = dat$G[, 1], covariates = dat$synthetic_data)
  conf <- infer_confounding(idata)
  expect_false(conf$mo_confounding$available)
  expect_equal(conf$mo_confounding$estimate, 0.8)  # default
})

test_that("mo_confounding unavailable for non-mediation data", {
  dat <- run_single_iteration(n_synthetic_samples = 100, n_features = 10,
                              seed = 42)
  idata <- iconic_data(Z = dat$Z, Y = t(dat$Y), W = t(dat$W),
                       G = dat$G[, 1], covariates = dat$synthetic_data)
  conf <- infer_confounding(idata)
  expect_false(conf$mo_confounding$available)
})

# ═══════════════════════════════════════════════════════════════
# omega inference (W-Y residual R²)
# ═══════════════════════════════════════════════════════════════

test_that("omega_1 is available with W", {
  idata <- .make_conf_data()
  conf <- infer_confounding(idata)
  expect_true(conf$omega_1$available)
  expect_true(conf$omega_1$estimate >= 0 && conf$omega_1$estimate <= 1)
})

test_that("omega unavailable without W", {
  dat <- run_single_iteration(n_synthetic_samples = 100, n_features = 10,
                              mo_confounding = 0.8, phi = 0.8, seed = 42)
  idata <- iconic_data(Z = dat$Z, Y = t(dat$Y), M = dat$M,
                       G = dat$G[, 1], Gm = dat$Gm)
  conf <- infer_confounding(idata)
  expect_false(conf$omega_1$available)
  expect_equal(conf$omega_1$estimate, 0.7)  # default
})

test_that("omega has composite warning", {
  idata <- .make_conf_data()
  conf <- infer_confounding(idata)
  if (conf$omega_1$available)
    expect_true(grepl("composite", conf$omega_1$warning))
})

# ═══════════════════════════════════════════════════════════════
# k inference (parallel analysis)
# ═══════════════════════════════════════════════════════════════

test_that("k is available with >= 5 features", {
  idata <- .make_conf_data(n_features = 10)
  conf <- infer_confounding(idata)
  expect_true(conf$k$available)
  expect_equal(conf$k$method, "parallel analysis (Horn, 1965)")
  expect_true(conf$k$estimate >= 1)
  expect_length(conf$k$ci, 2)
})

test_that("k unavailable with < 5 features", {
  idata <- .make_conf_data(n_features = 3)
  conf <- infer_confounding(idata)
  expect_false(conf$k$available)
  expect_equal(conf$k$estimate, 1L)  # default
})

# ═══════════════════════════════════════════════════════════════
# rho_G1, rho_G2 always unavailable
# ═══════════════════════════════════════════════════════════════

test_that("rho_G1 and rho_G2 are always unavailable", {
  idata <- .make_conf_data()
  conf <- infer_confounding(idata)
  expect_true("rho_G1" %in% conf$unavailable)
  expect_true("rho_G2" %in% conf$unavailable)
})

# ═══════════════════════════════════════════════════════════════
# Bare data (most parameters unavailable)
# ═══════════════════════════════════════════════════════════════

test_that("bare data has most parameters unavailable", {
  bare <- iconic_data(Z = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50))
  conf <- infer_confounding(bare)
  expect_false(conf$conf_strength$available)
  expect_false(conf$mo_confounding$available)
  expect_false(conf$omega_1$available)
  expect_false(conf$k$available)  # < 5 features
})

# ═══════════════════════════════════════════════════════════════
# Print method
# ═══════════════════════════════════════════════════════════════

test_that("print.iconic_confounding produces output", {
  idata <- .make_conf_data()
  est <- iconic_estimate(idata)
  conf <- infer_confounding(idata, estimate = est)
  out <- capture.output(print(conf))
  expect_true(any(grepl("iconic_confounding", out)))
  expect_true(any(grepl("conf_strength", out)))
})

# ═══════════════════════════════════════════════════════════════
# Integration with iconic_sensitivity
# ═══════════════════════════════════════════════════════════════

test_that("iconic_sensitivity with confounding=inferred includes inferred_confounding", {
  idata <- .make_conf_data(n = 80, n_features = 5)
  diag <- iconic_diagnose(idata)
  sens <- iconic_sensitivity(idata, diagnosis = diag,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3,
                             confounding = "inferred")
  expect_true(!is.null(sens$inferred_confounding))
  expect_s3_class(sens$inferred_confounding, "iconic_confounding")
})

test_that("iconic_sensitivity with confounding=default has NULL inferred_confounding", {
  idata <- .make_conf_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3,
                             confounding = "default")
  expect_null(sens$inferred_confounding)
})
