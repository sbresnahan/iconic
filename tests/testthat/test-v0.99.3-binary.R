# tests/testthat/test-v0.99.3-binary.R
# Binary outcome support
#
# Test groups:
# 1. iconic_data() binary interface
# 2. Total-effect binary estimators (logor + riskdiff)
# 3. Mediation binary estimators (logor + riskdiff)
# 4. iconic_estimate() binary dispatch (total + mediation)
# 5. Binary simulation DGP (generate_toy_data + run_single_iteration)
# 6. Sensitivity / prospect binary threading
# 7. COCA returns NA for binary (structural incompatibility)
# 8. Backward compatibility with continuous outcomes

# ── Helper: generate binary data with known truth ──
# Mirrors .make_surv_data(): same exposure/instrument/mediator structure,
# with the linear predictor converted to a 0/1 outcome via a logistic
# model. true_NDE / true_NIE are on the conditional log-OR scale.
.make_bin_data <- function(n = 300, seed = 42, mo_confounding = 0,
                           phi = 0) {
  set.seed(seed)
  G <- rnorm(n)
  X <- 0.5 * G + rnorm(n, sd = 0.5)
  Gm <- if (phi > 0) rnorm(n) else NULL
  M <- 0.4 * X
  if (mo_confounding > 0) M <- M + mo_confounding * 0.5 * rnorm(n)
  if (phi > 0) M <- M + phi * Gm
  M <- M + rnorm(n, sd = 0.05)
  eta <- 0.2 * X + 0.15 * M
  eta_c <- eta - mean(eta)
  # intercept solved for ~50% prevalence
  f <- function(e0) mean(plogis(e0 + eta_c)) - 0.5
  e0 <- uniroot(f, c(-20, 20))$root
  p <- plogis(e0 + eta_c)
  list(
    y_bin = rbinom(n, 1, p),
    X = X, G = G, Gm = Gm, M = M,
    W = matrix(rnorm(n * 5), 5, n),
    true_NDE = 0.2, true_NIE = 0.15
  )
}

# ============================================================
# 1. iconic_data() binary interface
# ============================================================

test_that("iconic_data accepts outcome_type = binary", {
  d <- .make_bin_data(n = 100)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      G = d$G, W = d$W)
  expect_true(inherits(bdat, "iconic_data"))
  expect_equal(bdat$outcome_type, "binary")
  expect_equal(bdat$n, 100)
  expect_equal(bdat$Y_bin, d$y_bin)
  expect_null(bdat$Y)
  expect_equal(bdat$n_features, 1L)
})

test_that("iconic_data binary with mediator", {
  d <- .make_bin_data(n = 100, phi = 0.8)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      M = d$M, G = d$G, Gm = d$Gm, W = d$W)
  expect_true(bdat$is_mediation)
  expect_equal(bdat$outcome_type, "binary")
})

test_that("iconic_data binary validation: Y must be 0/1", {
  expect_error(
    iconic_data(X = rnorm(50), Y = c(rep(2, 25), rep(0, 25)),
                outcome_type = "binary"),
    "0/1"
  )
})

test_that("iconic_data binary validation: no NA in Y", {
  expect_error(
    iconic_data(X = rnorm(50), Y = c(rep(NA, 5), rbinom(45, 1, 0.5)),
                outcome_type = "binary"),
    "NA"
  )
})

test_that("iconic_data binary validation: both levels present", {
  expect_error(
    iconic_data(X = rnorm(50), Y = rep(1, 50), outcome_type = "binary"),
    "both"
  )
})

test_that("iconic_data binary: print method shows case count", {
  d <- .make_bin_data(n = 100)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      G = d$G, W = d$W)
  out <- capture.output(print(bdat))
  expect_true(any(grepl("binary", out, ignore.case = TRUE)))
  expect_true(any(grepl("case", out, ignore.case = TRUE)))
})

test_that("as_iconic_data passthrough works for binary", {
  d <- .make_bin_data(n = 50)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      G = d$G, W = d$W)
  bdat2 <- as_iconic_data(bdat)
  expect_equal(bdat2$outcome_type, "binary")
  expect_equal(bdat2$Y_bin, bdat$Y_bin)
})

# ============================================================
# 2. Total-effect binary estimators (logor + riskdiff)
# ============================================================

test_that("fit_unadj_bin returns correct structure (logor)", {
  d <- .make_bin_data(n = 200)
  r <- fit_unadj_bin(d$y_bin, d$X)
  expect_named(r, c("beta", "se", "pvalue"))
  expect_true(is.numeric(r$beta) || is.na(r$beta))
})

test_that("fit_unadj_bin returns correct structure (riskdiff)", {
  d <- .make_bin_data(n = 200)
  r <- fit_unadj_bin(d$y_bin, d$X, effect_scale = "riskdiff")
  expect_named(r, c("beta", "se", "pvalue"))
  expect_true(is.numeric(r$beta) || is.na(r$beta))
})

test_that("fit_direct_bin returns correct structure", {
  d <- .make_bin_data(n = 200)
  r <- fit_direct_bin(d$y_bin, d$X, d$G, t(d$W))
  expect_named(r, c("beta", "se", "pvalue"))
})

test_that("fit_iv2sls_bin returns correct structure", {
  d <- .make_bin_data(n = 300)
  r <- fit_iv2sls_bin(d$y_bin, d$X, d$G, t(d$W))
  expect_named(r, c("beta", "se", "pvalue"))
})

test_that("fit_pgc_bin returns correct structure", {
  d <- .make_bin_data(n = 300)
  r <- fit_pgc_bin(d$y_bin, d$X, d$G, t(d$W))
  expect_named(r, c("beta", "se", "pvalue"))
})

test_that("fit_coca_bin returns NA with reason", {
  d <- .make_bin_data(n = 200)
  r <- fit_coca_bin(d$y_bin, d$X, t(d$W)[, 1])
  expect_true(all(is.na(unlist(r))))
  expect_true(!is.null(attr(r, "reason")))
  expect_true(grepl("COCA", attr(r, "reason")))
})

test_that("fit_iv2sls_bin returns NA for weak instrument", {
  d <- .make_bin_data(n = 200)
  # Use noise as instrument (F < 10)
  r <- fit_iv2sls_bin(d$y_bin, d$X, rnorm(200), t(d$W), min_f = 10)
  expect_true(all(is.na(unlist(r))))
})

# ============================================================
# 3. Mediation binary estimators (logor + riskdiff)
# ============================================================

test_that("fit_unadj_mediation_bin returns correct structure", {
  d <- .make_bin_data(n = 200, phi = 0.8)
  r <- fit_unadj_mediation_bin(d$y_bin, d$X, d$M)
  expect_true(all(c("NDE", "NDE_se", "NDE_p", "NIE", "NIE_se", "NIE_p",
                    "alpha_M", "alpha_se", "beta_M", "beta_M_se") %in% names(r)))
})

test_that("fit_direct_mediation_bin returns correct structure", {
  d <- .make_bin_data(n = 200, phi = 0.8)
  r <- fit_direct_mediation_bin(d$y_bin, d$X, d$M, d$G, t(d$W))
  expect_true(all(c("NDE", "NIE", "alpha_M", "beta_M") %in% names(r)))
})

test_that("fit_iv2sls_mediation_bin returns correct structure", {
  d <- .make_bin_data(n = 300, phi = 0.8)
  r <- fit_iv2sls_mediation_bin(d$y_bin, d$X, d$M, d$G, t(d$W))
  expect_true(all(c("NDE", "NIE", "alpha_M", "beta_M") %in% names(r)))
  # stage-3 mediator column is the uninstrumented M; must not be NA
  expect_true(!is.na(r$NDE))
  expect_true(!is.na(r$NIE))
})

test_that("fit_iv2sls_mediation2_bin returns correct structure", {
  d <- .make_bin_data(n = 300, phi = 0.8)
  r <- fit_iv2sls_mediation2_bin(d$y_bin, d$X, d$M, d$G, d$Gm)
  expect_true(all(c("NDE", "NIE", "alpha_M", "beta_M") %in% names(r)))
})

test_that("fit_pgc_mediation2_bin returns correct structure", {
  d <- .make_bin_data(n = 300, phi = 0.8)
  r <- fit_pgc_mediation2_bin(d$y_bin, d$X, d$M, d$G,
                              t(d$W), t(d$W), gm = d$Gm)
  expect_true(all(c("NDE", "NIE", "alpha_M", "beta_M") %in% names(r)))
})

test_that("fit_pgc_mediation_bin returns correct structure", {
  d <- .make_bin_data(n = 300, phi = 0.8)
  r <- fit_pgc_mediation_bin(d$y_bin, d$X, d$M, d$G, t(d$W))
  expect_true(all(c("NDE", "NIE", "alpha_M", "beta_M") %in% names(r)))
  expect_true(is.numeric(r$NDE))
  expect_true(is.numeric(r$NIE))
})

test_that("fit_pgc_mediation_bin NIE = alpha_M * beta_M", {
  d <- .make_bin_data(n = 300, phi = 0.8)
  r <- fit_pgc_mediation_bin(d$y_bin, d$X, d$M, d$G, t(d$W))
  if (!is.na(r$NIE) && !is.na(r$alpha_M) && !is.na(r$beta_M)) {
    expect_equal(r$NIE, r$alpha_M * r$beta_M, tolerance = 1e-8)
  }
})

test_that("PGC and PGC2 produce distinct results in binary mediation", {
  # Use path-specific loadings so PGC (single-panel) and PGC2 (path-specific) differ
  dat <- iconic::generate_toy_data(
    n = 500, n_features = 1, beta_X = 0.25,
    alpha_M = 0.50, beta_M = 0.30, conf_str = 0.6,
    w_signal = 0.7, phi = 0.8, mo_confounding = 0.8, rho_G1 = 0.3,
    outcome_type = "binary", seed = 123)
  bdat <- iconic_data(X = dat$X, Y = dat$y_bin, outcome_type = "binary",
                      M = dat$M, G = dat$G[, 1], Gm = dat$Gm,
                      W = dat$W, W1 = dat$W1, W2 = dat$W2)
  est <- iconic_estimate(bdat, effect_scale = "logor")
  pgc <- est[est$method == "PGC", ]
  pgc2 <- est[est$method == "PGC2", ]
  # Both should produce non-NA estimates
  expect_true(!is.na(pgc$NDE))
  expect_true(!is.na(pgc2$NDE))
  # They should NOT be identical (different estimators)
  expect_true(pgc$NDE != pgc2$NDE)
  expect_true(pgc$NIE != pgc2$NIE)
})

test_that("fit_coca_mediation_bin returns NA with reason", {
  d <- .make_bin_data(n = 200, phi = 0.8)
  r <- fit_coca_mediation_bin(d$y_bin, d$X, d$M, t(d$W)[, 1])
  expect_true(all(is.na(unlist(r))))
  expect_true(!is.null(attr(r, "reason")))
})

test_that("NIE = alpha_M * beta_M for binary mediation (riskdiff)", {
  d <- .make_bin_data(n = 300, phi = 0.8)
  r <- fit_unadj_mediation_bin(d$y_bin, d$X, d$M,
                               effect_scale = "riskdiff")
  if (!is.na(r$NIE) && !is.na(r$alpha_M) && !is.na(r$beta_M)) {
    expect_equal(r$NIE, r$alpha_M * r$beta_M, tolerance = 1e-8)
  }
})

# ============================================================
# 4. iconic_estimate() binary dispatch
# ============================================================

test_that("iconic_estimate binary total-effect returns data.frame", {
  d <- .make_bin_data(n = 300)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      G = d$G, W = d$W)
  est <- iconic_estimate(bdat, effect_scale = "logor")
  expect_true(is.data.frame(est))
  expect_true(all(c("method", "beta", "se", "pvalue") %in% names(est)))
  expect_true("COCA" %in% est$method)
  expect_true(all(is.na(est$beta[est$method == "COCA"])))
})

test_that("iconic_estimate binary mediation returns data.frame", {
  d <- .make_bin_data(n = 300, phi = 0.8)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      M = d$M, G = d$G, Gm = d$Gm, W = d$W)
  est <- iconic_estimate(bdat, effect_scale = "logor")
  expect_true(is.data.frame(est))
  expect_true(all(c("method", "NDE", "NDE_se", "NIE", "NIE_se") %in% names(est)))
  # All 8 methods should be present
  expect_true(all(c("UNADJ", "DIRECT", "IV2SLS", "COCA", "PGC",
                    "IV2SLS2", "PGC2", "PGC2Gm") %in% est$method))
  # COCA should be NA
  expect_true(all(is.na(est$NDE[est$method == "COCA"])))
})

test_that("iconic_estimate binary riskdiff scale works", {
  d <- .make_bin_data(n = 300, phi = 0.8)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      M = d$M, G = d$G, Gm = d$Gm, W = d$W)
  est <- iconic_estimate(bdat, effect_scale = "riskdiff")
  expect_true(is.data.frame(est))
  expect_true(all(c("NDE", "NIE") %in% names(est)))
})

test_that("iconic_estimate effect_scale loghr remaps to logor for binary", {
  d <- .make_bin_data(n = 200)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      G = d$G, W = d$W)
  expect_message(iconic_estimate(bdat, effect_scale = "loghr"),
                 "logor")
})

# ============================================================
# 5. Binary simulation DGP
# ============================================================

test_that("generate_toy_data binary returns y_bin", {
  dat <- iconic:::generate_toy_data(n = 200, n_features = 5,
                                    outcome_type = "binary", seed = 42)
  expect_true("y_bin" %in% names(dat))
  expect_equal(dat$outcome_type, "binary")
  expect_length(dat$y_bin, 200)
  expect_true(all(dat$y_bin %in% c(0, 1)))
  # Prevalence should be reasonable (20-80%)
  prev <- mean(dat$y_bin)
  expect_true(prev > 0.2 && prev < 0.8)
})

test_that("generate_toy_data binary bin_prev controls prevalence", {
  dat_lo <- iconic:::generate_toy_data(n = 2000, outcome_type = "binary",
                                       bin_prev = 0.2, seed = 42)
  dat_hi <- iconic:::generate_toy_data(n = 2000, outcome_type = "binary",
                                       bin_prev = 0.8, seed = 42)
  expect_true(abs(mean(dat_lo$y_bin) - 0.2) < 0.05)
  expect_true(abs(mean(dat_hi$y_bin) - 0.8) < 0.05)
})

test_that("generate_toy_data binary dgp_params records outcome_type", {
  dat <- iconic:::generate_toy_data(n = 100, outcome_type = "binary", seed = 1)
  expect_equal(dat$dgp_params$outcome_type, "binary")
  expect_true("bin_prev" %in% names(dat$dgp_params))
})

test_that("run_single_iteration binary returns y_bin", {
  dat <- run_single_iteration(n_synthetic_samples = 200, n_features = 5,
                              outcome_type = "binary", seed = 42)
  expect_true("y_bin" %in% names(dat))
  expect_equal(dat$outcome_type, "binary")
  expect_length(dat$y_bin, 200)
  expect_true(all(dat$y_bin %in% c(0, 1)))
})

# ============================================================
# 6. Sensitivity / prospect binary threading
# ============================================================

test_that("gan_sensitivity binary returns summary", {
  sens <- gan_sensitivity(NULL, conf_grid = c(0.4, 0.8),
                          coverage_grid = c(0.5, 1), k_grid = 1,
                          n_iter = 2, n_samples = 200, n_features = 3,
                          outcome_type = "binary", base_seed = 100)
  expect_true(is.data.frame(sens$summary))
  expect_true("method" %in% names(sens$summary))
  expect_true("COCA" %in% sens$summary$method)
  # COCA should have NaN bias
  coca_bias <- sens$summary$bias[sens$summary$method == "COCA"]
  expect_true(all(is.nan(coca_bias)))
})

test_that("gan_mediation_sensitivity binary returns summary", {
  sens <- gan_mediation_sensitivity(NULL, conf_grid = c(0.4, 0.8),
                                    coverage_grid = c(0.5, 1), k_grid = 1,
                                    mo_confounding = 0.8, phi = 0.8,
                                    n_iter = 2, n_samples = 200,
                                    n_features = 3,
                                    outcome_type = "binary",
                                    base_seed = 200)
  expect_true(is.data.frame(sens$summary))
  expect_true(all(c("NDE_bias", "NIE_bias") %in% names(sens$summary)))
  # IV2SLS2 should be present (phi > 0)
  expect_true("IV2SLS2" %in% sens$summary$method)
})

test_that("iconic_sensitivity binary returns surface", {
  d <- .make_bin_data(n = 100, phi = 0.8)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      M = d$M, G = d$G, Gm = d$Gm, W = d$W)
  sens <- iconic_sensitivity(bdat, n_iter = 2, n_features = 3,
                             rho_G1_grid = c(0, 0.3),
                             rho_G2_grid = c(0, 0.3),
                             gan_epochs = 5, base_seed = 300)
  expect_true(is.data.frame(sens$surface))
  expect_true("tipped" %in% names(sens$surface))
  expect_true(is.data.frame(sens$tipping_points))
  # COCA tipping should have NA max_bias (not -Inf)
  coca_max <- sens$tipping_points$max_NDE_bias[
    sens$tipping_points$method == "COCA"]
  expect_true(is.na(coca_max))
})

test_that("iconic_prospect binary returns prospect object", {
  d <- .make_bin_data(n = 100, phi = 0)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      M = d$M)
  prospect <- iconic_prospect(bdat, n_iter = 2, n_features = 3,
                              gamma_G_grid = c(0.4, 0.8),
                              target_gamma_G = 0.6,
                              gan_epochs = 5, base_seed = 400)
  expect_true(inherits(prospect, "iconic_prospect"))
  expect_true(is.data.frame(prospect$strength_surface))
  expect_true(is.data.frame(prospect$prospective))
})

# ============================================================
# 7. COCA returns NA for binary (structural incompatibility)
# ============================================================

test_that("COCA total-effect binary: NA with reason", {
  d <- .make_bin_data(n = 200)
  r <- fit_coca_bin(d$y_bin, d$X, t(d$W)[, 1])
  expect_true(all(is.na(unlist(r))))
  expect_true(grepl("binary", attr(r, "reason"), ignore.case = TRUE))
})

test_that("COCA mediation binary: NA with reason", {
  d <- .make_bin_data(n = 200, phi = 0.8)
  r <- fit_coca_mediation_bin(d$y_bin, d$X, d$M, t(d$W)[, 1])
  expect_true(all(is.na(unlist(r))))
  expect_true(grepl("binary", attr(r, "reason"), ignore.case = TRUE))
})

test_that("COCA NA propagates through iconic_estimate total-effect", {
  d <- .make_bin_data(n = 200)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      G = d$G, W = d$W)
  est <- iconic_estimate(bdat)
  expect_true(all(is.na(est$beta[est$method == "COCA"])))
  expect_true(all(is.na(est$se[est$method == "COCA"])))
})

test_that("COCA NA propagates through iconic_estimate mediation", {
  d <- .make_bin_data(n = 200, phi = 0.8)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      M = d$M, G = d$G, Gm = d$Gm, W = d$W)
  est <- iconic_estimate(bdat)
  expect_true(all(is.na(est$NDE[est$method == "COCA"])))
  expect_true(all(is.na(est$NIE[est$method == "COCA"])))
})

# ============================================================
# 8. Backward compatibility with continuous outcomes
# ============================================================

test_that("iconic_estimate continuous total-effect still works", {
  dat <- iconic:::generate_toy_data(n = 200, seed = 42)
  cdat <- iconic_data(X = dat$X, Y = dat$Y, G = dat$G,
                      W = dat$W)
  est <- iconic_estimate(cdat)
  expect_true(is.data.frame(est))
  expect_true(all(c("method", "beta", "se", "pvalue") %in% names(est)))
  # COCA should NOT be NA for continuous (toy data has confounding structure)
  expect_true(!all(is.na(est$beta[est$method == "COCA"])))
})

test_that("iconic_estimate continuous mediation still works", {
  dat <- iconic:::generate_toy_data(n = 200, seed = 42, n_mediators = 1, phi = 0.8)
  cdat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M, G = dat$G, Gm = dat$Gm,
                      W = dat$W)
  est <- iconic_estimate(cdat)
  expect_true(is.data.frame(est))
  expect_true(all(c("NDE", "NIE") %in% names(est)))
})

test_that("infer_confounding errors informatively for binary", {
  d <- .make_bin_data(n = 100)
  bdat <- iconic_data(X = d$X, Y = d$y_bin, outcome_type = "binary",
                      G = d$G, W = d$W)
  expect_error(infer_confounding(bdat), "binary")
})
