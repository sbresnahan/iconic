# tests/testthat/test--survival.R
# Survival / time-to-event outcome support
#
# Test groups:
# 1. iconic_data() survival interface
# 2. Total-effect survival estimators (loghr + rmst)
# 3. Mediation survival estimators (loghr + rmst)
# 4. iconic_estimate() survival dispatch (total + mediation)
# 5. Survival simulation DGP (generate_toy_data + run_single_iteration)
# 6. Sensitivity / prospect survival threading
# 7. COCA returns NA for survival (structural incompatibility)
# 8. Backward compatibility with continuous outcomes

skip_if_not_installed("survival")

# ── Helper: generate survival data with known truth ──
.make_surv_data <- function(n = 300, seed = 42, mo_confounding = 0,
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
  h0 <- 0.1
  T_true <- -log(runif(n)) / (h0 * exp(eta - mean(eta)))
  C <- rexp(n, rate = h0 * (1 - 0.6) / 0.6)
  list(
    time = pmin(T_true, C),
    event = as.integer(T_true <= C),
    X = X, G = G, Gm = Gm, M = M,
    W = matrix(rnorm(n * 5), 5, n),
    true_NDE = 0.2, true_NIE = 0.15
  )
}

# ============================================================
# 1. iconic_data() survival interface
# ============================================================

test_that("iconic_data accepts outcome_type = survival", {
  d <- .make_surv_data(n = 100)
  sdat <- iconic_data(X = d$X, outcome_type = "survival",
                      surv_time = d$time, surv_event = d$event,
                      G = d$G, W = d$W)
  expect_true(inherits(sdat, "iconic_data"))
  expect_equal(sdat$outcome_type, "survival")
  expect_equal(sdat$n, 100)
  expect_equal(sdat$surv_time, d$time)
  expect_equal(sdat$surv_event, d$event)
  expect_null(sdat$Y)
})

test_that("iconic_data survival with mediator", {
  d <- .make_surv_data(n = 100, phi = 0.8)
  sdat <- iconic_data(X = d$X, outcome_type = "survival",
                      surv_time = d$time, surv_event = d$event,
                      M = d$M, G = d$G, Gm = d$Gm, W = d$W)
  expect_true(sdat$is_mediation)
  expect_equal(sdat$outcome_type, "survival")
})

test_that("iconic_data survival validation: positive time", {
  expect_error(
    iconic_data(X = rnorm(50), outcome_type = "survival",
                surv_time = c(rep(-1, 25), rexp(25)),
                surv_event = rbinom(50, 1, 0.5)),
    "positive"
  )
})

test_that("iconic_data survival validation: event 0/1", {
  expect_error(
    iconic_data(X = rnorm(50), outcome_type = "survival",
                surv_time = rexp(50), surv_event = c(rep(2, 25), rep(0, 25))),
    "0/1"
  )
})

test_that("iconic_data survival validation: no NA in time/event", {
  expect_error(
    iconic_data(X = rnorm(50), outcome_type = "survival",
                surv_time = c(rep(NA, 5), rexp(45)),
                surv_event = rbinom(50, 1, 0.5)),
    "NA"
  )
})

test_that("iconic_data survival: print method shows event count", {
  d <- .make_surv_data(n = 100)
  sdat <- iconic_data(X = d$X, outcome_type = "survival",
                      surv_time = d$time, surv_event = d$event,
                      G = d$G, W = d$W)
  out <- capture.output(print(sdat))
  expect_true(any(grepl("survival", out, ignore.case = TRUE)))
  expect_true(any(grepl("event", out, ignore.case = TRUE)))
})

test_that("as_iconic_data passthrough works for survival", {
  d <- .make_surv_data(n = 50)
  sdat <- iconic_data(X = d$X, outcome_type = "survival",
                      surv_time = d$time, surv_event = d$event,
                      G = d$G, W = d$W)
  sdat2 <- as_iconic_data(sdat)
  expect_equal(sdat2$outcome_type, "survival")
  expect_equal(sdat2$surv_time, sdat$surv_time)
})

# ============================================================
# 2. Total-effect survival estimators (loghr + rmst)
# ============================================================

test_that("fit_unadj_surv returns correct structure (loghr)", {
  d <- .make_surv_data(n = 200)
  r <- fit_unadj_surv(d$time, d$event, d$X)
  expect_named(r, c("beta", "se", "pvalue"))
  expect_true(is.numeric(r$beta) || is.na(r$beta))
})

test_that("fit_unadj_surv returns correct structure (rmst)", {
  d <- .make_surv_data(n = 200)
  r <- fit_unadj_surv(d$time, d$event, d$X, effect_scale = "rmst")
  expect_named(r, c("beta", "se", "pvalue"))
  expect_true(is.numeric(r$beta) || is.na(r$beta))
})

test_that("fit_direct_surv returns correct structure", {
  d <- .make_surv_data(n = 200)
  r <- fit_direct_surv(d$time, d$event, d$X, d$G, t(d$W))
  expect_named(r, c("beta", "se", "pvalue"))
})

test_that("fit_iv2sls_surv returns correct structure", {
  d <- .make_surv_data(n = 300)
  r <- fit_iv2sls_surv(d$time, d$event, d$X, d$G, t(d$W))
  expect_named(r, c("beta", "se", "pvalue"))
})

test_that("fit_pgc_surv returns correct structure", {
  d <- .make_surv_data(n = 300)
  r <- fit_pgc_surv(d$time, d$event, d$X, d$G, t(d$W))
  expect_named(r, c("beta", "se", "pvalue"))
})

test_that("fit_coca_surv returns NA with reason", {
  d <- .make_surv_data(n = 200)
  r <- fit_coca_surv(d$time, d$event, d$X, t(d$W)[, 1])
  expect_true(all(is.na(unlist(r))))
  expect_true(!is.null(attr(r, "reason")))
  expect_true(grepl("COCA", attr(r, "reason")))
})

test_that("fit_iv2sls_surv returns NA for weak instrument", {
  d <- .make_surv_data(n = 200)
  # Use noise as instrument (F < 10)
  r <- fit_iv2sls_surv(d$time, d$event, d$X, rnorm(200), t(d$W),
                       min_f = 10)
  expect_true(all(is.na(unlist(r))))
})

# ============================================================
# 3. Mediation survival estimators (loghr + rmst)
# ============================================================

test_that("fit_unadj_mediation_surv returns correct structure", {
  d <- .make_surv_data(n = 200, phi = 0.8)
  r <- fit_unadj_mediation_surv(d$time, d$event, d$X, d$M)
  expect_true(all(c("NDE", "NDE_se", "NDE_p", "NIE", "NIE_se", "NIE_p",
                    "alpha_M", "alpha_se", "beta_M", "beta_M_se") %in% names(r)))
})

test_that("fit_direct_mediation_surv returns correct structure", {
  d <- .make_surv_data(n = 200, phi = 0.8)
  r <- fit_direct_mediation_surv(d$time, d$event, d$X, d$M, d$G, t(d$W))
  expect_true(all(c("NDE", "NIE", "alpha_M", "beta_M") %in% names(r)))
})

test_that("fit_iv2sls_mediation_surv returns correct structure", {
  d <- .make_surv_data(n = 300, phi = 0.8)
  r <- fit_iv2sls_mediation_surv(d$time, d$event, d$X, d$M, d$G, t(d$W))
  expect_true(all(c("NDE", "NIE", "alpha_M", "beta_M") %in% names(r)))
})

test_that("fit_iv2sls_mediation2_surv returns correct structure", {
  d <- .make_surv_data(n = 300, phi = 0.8)
  r <- fit_iv2sls_mediation2_surv(d$time, d$event, d$X, d$M, d$G, d$Gm, t(d$W))
  expect_true(all(c("NDE", "NIE", "alpha_M", "beta_M") %in% names(r)))
})

test_that("fit_pgc_mediation2_surv returns correct structure", {
  d <- .make_surv_data(n = 300, phi = 0.8)
  r <- fit_pgc_mediation2_surv(d$time, d$event, d$X, d$M, d$G,
                               t(d$W), t(d$W), gm = d$Gm)
  expect_true(all(c("NDE", "NIE", "alpha_M", "beta_M") %in% names(r)))
})

test_that("fit_pgc_mediation_surv returns correct structure", {
  d <- .make_surv_data(n = 300, phi = 0.8)
  r <- fit_pgc_mediation_surv(d$time, d$event, d$X, d$M, d$G, t(d$W))
  expect_true(all(c("NDE", "NIE", "alpha_M", "beta_M") %in% names(r)))
  expect_true(is.numeric(r$NDE))
  expect_true(is.numeric(r$NIE))
})

test_that("fit_pgc_mediation_surv NIE = alpha_M * beta_M", {
  d <- .make_surv_data(n = 300, phi = 0.8)
  r <- fit_pgc_mediation_surv(d$time, d$event, d$X, d$M, d$G, t(d$W))
  if (!is.na(r$NIE) && !is.na(r$alpha_M) && !is.na(r$beta_M)) {
    expect_equal(r$NIE, r$alpha_M * r$beta_M, tolerance = 1e-8)
  }
})

test_that("PGC and PGC2 produce distinct results in survival mediation", {
  # Use path-specific loadings so PGC (single-panel) and PGC2 (path-specific) differ
  set.seed(123)
  dat <- iconic::generate_toy_data(
    n = 500, n_features = 1, beta_X = 0.25,
    alpha_M = 0.50, beta_M = 0.30, conf_str = 0.6,
    w_signal = 0.7, phi = 0.8, mo_confounding = 0.8, rho_G1 = 0.3,
    outcome_type = "survival",
    surv_event_frac = 0.6, seed = 123)
  sdat <- iconic_data(X = dat$X, outcome_type = "survival",
                      surv_time = dat$surv_time, surv_event = dat$surv_event,
                      M = dat$M, G = dat$G[, 1], Gm = dat$Gm,
                      W = dat$W, W1 = dat$W1, W2 = dat$W2)
  est <- iconic_estimate(sdat, effect_scale = "loghr")
  pgc <- est[est$method == "PGC", ]
  pgc2 <- est[est$method == "PGC2", ]
  # Both should produce non-NA estimates
  expect_true(!is.na(pgc$NDE))
  expect_true(!is.na(pgc2$NDE))
  # They should NOT be identical (different estimators)
  expect_true(pgc$NDE != pgc2$NDE)
  expect_true(pgc$NIE != pgc2$NIE)
})

test_that("fit_coca_mediation_surv returns NA with reason", {
  d <- .make_surv_data(n = 200, phi = 0.8)
  r <- fit_coca_mediation_surv(d$time, d$event, d$X, d$M, t(d$W)[, 1])
  expect_true(all(is.na(unlist(r))))
  expect_true(!is.null(attr(r, "reason")))
})

test_that("NIE = alpha_M * beta_M for survival mediation (rmst)", {
  d <- .make_surv_data(n = 300, phi = 0.8)
  r <- fit_unadj_mediation_surv(d$time, d$event, d$X, d$M,
                                effect_scale = "rmst")
  if (!is.na(r$NIE) && !is.na(r$alpha_M) && !is.na(r$beta_M)) {
    expect_equal(r$NIE, r$alpha_M * r$beta_M, tolerance = 1e-8)
  }
})

# ============================================================
# 4. iconic_estimate() survival dispatch
# ============================================================

test_that("iconic_estimate survival total-effect returns data.frame", {
  d <- .make_surv_data(n = 300)
  sdat <- iconic_data(X = d$X, outcome_type = "survival",
                      surv_time = d$time, surv_event = d$event,
                      G = d$G, W = d$W)
  est <- iconic_estimate(sdat, effect_scale = "loghr")
  expect_true(is.data.frame(est))
  expect_true(all(c("method", "beta", "se", "pvalue") %in% names(est)))
  expect_true("COCA" %in% est$method)
  expect_true(all(is.na(est$beta[est$method == "COCA"])))
})

test_that("iconic_estimate survival mediation returns data.frame", {
  d <- .make_surv_data(n = 300, phi = 0.8)
  sdat <- iconic_data(X = d$X, outcome_type = "survival",
                      surv_time = d$time, surv_event = d$event,
                      M = d$M, G = d$G, Gm = d$Gm, W = d$W)
  est <- iconic_estimate(sdat, effect_scale = "loghr")
  expect_true(is.data.frame(est))
  expect_true(all(c("method", "NDE", "NDE_se", "NIE", "NIE_se") %in% names(est)))
  # All 8 methods should be present
  expect_true(all(c("UNADJ", "DIRECT", "IV2SLS", "COCA", "PGC",
                    "IV2SLS2", "PGC2", "PGC2Gm") %in% est$method))
  # COCA should be NA
  expect_true(all(is.na(est$NDE[est$method == "COCA"])))
})

test_that("iconic_estimate survival rmst scale works", {
  d <- .make_surv_data(n = 300, phi = 0.8)
  sdat <- iconic_data(X = d$X, outcome_type = "survival",
                      surv_time = d$time, surv_event = d$event,
                      M = d$M, G = d$G, Gm = d$Gm, W = d$W)
  est <- iconic_estimate(sdat, effect_scale = "rmst")
  expect_true(is.data.frame(est))
  expect_true(all(c("NDE", "NIE") %in% names(est)))
})

test_that("iconic_estimate effect_scale rmst ignored for continuous", {
  cdat <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100), 1, 100),
                      G = rnorm(100), W = matrix(rnorm(100 * 5), 5, 100))
  expect_message(iconic_estimate(cdat, effect_scale = "rmst"),
                 "ignored")
})

# ============================================================
# 5. Survival simulation DGP
# ============================================================

test_that("generate_toy_data survival returns surv_time and surv_event", {
  dat <- iconic:::generate_toy_data(n = 200, n_features = 5,
                                     outcome_type = "survival", seed = 42)
  expect_true("surv_time" %in% names(dat))
  expect_true("surv_event" %in% names(dat))
  expect_equal(dat$outcome_type, "survival")
  expect_length(dat$surv_time, 200)
  expect_length(dat$surv_event, 200)
  expect_true(all(dat$surv_time >= 0))
  expect_true(all(dat$surv_event %in% c(0, 1)))
  # Event fraction should be reasonable (30-80%)
  frac <- mean(dat$surv_event)
  expect_true(frac > 0.3 && frac < 0.8)
})

test_that("generate_toy_data continuous backward compat", {
  dat <- iconic:::generate_toy_data(n = 100, n_features = 5, seed = 42)
  expect_false("surv_time" %in% names(dat))
  expect_equal(dat$outcome_type, "continuous")
  expect_equal(dim(dat$Y), c(100, 5))
})

test_that("run_single_iteration survival returns surv_time and surv_event", {
  dat <- run_single_iteration(n_synthetic_samples = 200, n_features = 5,
                              outcome_type = "survival", seed = 42)
  expect_true("surv_time" %in% names(dat))
  expect_true("surv_event" %in% names(dat))
  expect_equal(dat$outcome_type, "survival")
  expect_length(dat$surv_time, 200)
  expect_true(all(dat$surv_event %in% c(0, 1)))
})

test_that("run_single_iteration continuous backward compat", {
  dat <- run_single_iteration(n_synthetic_samples = 100, n_features = 5,
                              seed = 42)
  expect_false("surv_time" %in% names(dat))
  expect_equal(dat$outcome_type, "continuous")
})

test_that("generate_toy_data survival dgp_params records outcome_type", {
  dat <- iconic:::generate_toy_data(n = 100, outcome_type = "survival", seed = 1)
  expect_equal(dat$dgp_params$outcome_type, "survival")
  expect_true("surv_h0" %in% names(dat$dgp_params))
})

# ============================================================
# 6. Sensitivity / prospect survival threading
# ============================================================

test_that("gan_sensitivity survival returns summary", {
  sens <- gan_sensitivity(NULL, conf_grid = c(0.4, 0.8),
                          coverage_grid = c(0.5, 1), k_grid = 1,
                          n_iter = 2, n_samples = 200, n_features = 3,
                          outcome_type = "survival", base_seed = 100)
  expect_true(is.data.frame(sens$summary))
  expect_true("method" %in% names(sens$summary))
  expect_true("COCA" %in% sens$summary$method)
  # COCA should have NaN bias
  coca_bias <- sens$summary$bias[sens$summary$method == "COCA"]
  expect_true(all(is.nan(coca_bias)))
})

test_that("gan_mediation_sensitivity survival returns summary", {
  sens <- gan_mediation_sensitivity(NULL, conf_grid = c(0.4, 0.8),
                                     coverage_grid = c(0.5, 1), k_grid = 1,
                                     mo_confounding = 0.8, phi = 0.8,
                                     n_iter = 2, n_samples = 200,
                                     n_features = 3,
                                     outcome_type = "survival",
                                     base_seed = 200)
  expect_true(is.data.frame(sens$summary))
  expect_true(all(c("NDE_bias", "NIE_bias") %in% names(sens$summary)))
  # IV2SLS2 should be present (phi > 0)
  expect_true("IV2SLS2" %in% sens$summary$method)
})

test_that("iconic_sensitivity survival returns surface", {
  d <- .make_surv_data(n = 100, phi = 0.8)
  sdat <- iconic_data(X = d$X, outcome_type = "survival",
                      surv_time = d$time, surv_event = d$event,
                      M = d$M, G = d$G, Gm = d$Gm, W = d$W)
  sens <- iconic_sensitivity(sdat, n_iter = 2, n_features = 3,
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

test_that("iconic_prospect survival returns prospect object", {
  d <- .make_surv_data(n = 100, phi = 0)
  sdat <- iconic_data(X = d$X, outcome_type = "survival",
                      surv_time = d$time, surv_event = d$event,
                      M = d$M)
  prospect <- iconic_prospect(sdat, n_iter = 2, n_features = 3,
                              gamma_G_grid = c(0.4, 0.8),
                              target_gamma_G = 0.6,
                              gan_epochs = 5, base_seed = 400)
  expect_true(inherits(prospect, "iconic_prospect"))
  expect_true(is.data.frame(prospect$strength_surface))
  expect_true(is.data.frame(prospect$prospective))
})

# ============================================================
# 7. COCA returns NA for survival (structural incompatibility)
# ============================================================

test_that("COCA total-effect survival: NA with reason", {
  d <- .make_surv_data(n = 200)
  r <- fit_coca_surv(d$time, d$event, d$X, t(d$W)[, 1])
  expect_true(all(is.na(unlist(r))))
  expect_true(grepl("survival", attr(r, "reason"), ignore.case = TRUE))
})

test_that("COCA mediation survival: NA with reason", {
  d <- .make_surv_data(n = 200, phi = 0.8)
  r <- fit_coca_mediation_surv(d$time, d$event, d$X, d$M, t(d$W)[, 1])
  expect_true(all(is.na(unlist(r))))
  expect_true(grepl("survival", attr(r, "reason"), ignore.case = TRUE))
})

test_that("COCA NA propagates through iconic_estimate total-effect", {
  d <- .make_surv_data(n = 200)
  sdat <- iconic_data(X = d$X, outcome_type = "survival",
                      surv_time = d$time, surv_event = d$event,
                      G = d$G, W = d$W)
  est <- iconic_estimate(sdat)
  expect_true(all(is.na(est$beta[est$method == "COCA"])))
  expect_true(all(is.na(est$se[est$method == "COCA"])))
})

test_that("COCA NA propagates through iconic_estimate mediation", {
  d <- .make_surv_data(n = 200, phi = 0.8)
  sdat <- iconic_data(X = d$X, outcome_type = "survival",
                      surv_time = d$time, surv_event = d$event,
                      M = d$M, G = d$G, Gm = d$Gm, W = d$W)
  est <- iconic_estimate(sdat)
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

test_that("gan_sensitivity continuous backward compat", {
  sens <- gan_sensitivity(NULL, conf_grid = c(0.8),
                          coverage_grid = c(1), n_iter = 2,
                          n_samples = 100, n_features = 3,
                          base_seed = 500)
  expect_true(is.data.frame(sens$summary))
  # COCA should NOT be NaN for continuous
  coca_bias <- sens$summary$bias[sens$summary$method == "COCA"]
  expect_true(!all(is.nan(coca_bias)))
})
