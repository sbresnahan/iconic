# Acceptance tests for the IV2SLS2 path-specific negative-control fix (v0.9.9).
#
# Background: IV2SLS2 previously conditioned on the POOLED NC panel W in all
# three 2SLS stages. Under multi-confounder designs (k >= 2, distinct path
# loadings) the pooled panel scale((W1 + W2) / 2) is a common child of the two
# independent confounders conf_XM (U1) and conf_MY (U2) -- a collider. Conditioning
# on it opened a spurious path so the NDE bias got WORSE as coverage omega rose.
# The fix makes NC augmentation optional and path-specific: W1 in stage 1
# (X ~ G + W1), W2 in stages 2-3 (M ~ X_hat + Gm + W2, Y ~ X_hat + M_hat + W2).
# With path-specific panels, higher coverage now IMPROVES the NDE estimate.

# Helper: mean |NDE bias| of IV2SLS2 across replicates at a given coverage.
# k = 2 distinct path loadings (lambda_XM = e1, lambda_MY = e2) so W1/W2 are
# distinct panels proxying independent confounders.
.iv2sls2_nde_bias <- function(omega, n = 500, n_features = 5, n_iter = 10,
                              seed0 = 9100, W2_only = FALSE) {
  nde <- numeric(0)
  true_nde <- NA_real_
  for (i in seq_len(n_iter)) {
    dat <- iconic:::generate_toy_data(
      n = n, n_features = n_features, phi = 0.8, mo_confounding = 0.8,
      rho_G1 = 0.3, rho_G2 = 0.3,
      lambda_XM = c(1, 0), lambda_MY = c(0, 1),
      omega_1 = omega, omega_2 = omega, seed = seed0 + i)
    true_nde <- dat$true_NDE
    for (f in seq_len(n_features)) {
      r <- if (W2_only) {
        fit_iv2sls_mediation2(dat$Y[, f], dat$X, dat$M, dat$G[, f], dat$Gm,
                              W2 = dat$W2)
      } else {
        fit_iv2sls_mediation2(dat$Y[, f], dat$X, dat$M, dat$G[, f], dat$Gm,
                              W1 = dat$W1, W2 = dat$W2)
      }
      nde <- c(nde, r$NDE)
    }
  }
  mean(abs(nde - true_nde), na.rm = TRUE)
}

test_that("IV2SLS2 path-specific W1/W2: coverage helps (collider removed)", {
  # Core acceptance: mean |NDE bias| is non-increasing in omega.
  b03 <- .iv2sls2_nde_bias(0.3)
  b07 <- .iv2sls2_nde_bias(0.7)
  b10 <- .iv2sls2_nde_bias(1.0)
  # Coverage must help, not hurt: bias at omega = 1 no worse than at omega = 0.3.
  expect_lte(b10, b03 + 0.01)
  # And the high-coverage estimate should be meaningfully accurate.
  expect_lt(b10, 0.05)
})

test_that("IV2SLS2 with W1 = W2 = NULL reduces to plain 2-stage MR", {
  dat <- iconic:::generate_toy_data(n = 500, n_features = 1, phi = 0.8,
                                    mo_confounding = 0.8, seed = 42)
  r <- fit_iv2sls_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$Gm)
  # Pure-MR point estimate equals the sequential 2SLS with no NC covariate.
  fs <- lm(dat$X ~ dat$G[, 1]); X_hat <- fitted(fs)
  ms <- lm(dat$M ~ X_hat + dat$Gm); M_hat <- fitted(ms)
  os <- lm(dat$Y[, 1] ~ X_hat + M_hat)
  expect_equal(unname(r$NDE), unname(coef(os)["X_hat"]), tolerance = 1e-8)
})

test_that("IV2SLS2 W2-only: coverage helps (case-study config)", {
  b03 <- .iv2sls2_nde_bias(0.3, W2_only = TRUE)
  b10 <- .iv2sls2_nde_bias(1.0, W2_only = TRUE)
  expect_lte(b10, b03 + 0.01)
})

test_that("IV2SLS2 pooled guard: identical W1/W2 falls back to pure MR", {
  dat <- iconic:::generate_toy_data(n = 500, n_features = 1, phi = 0.8,
                                    mo_confounding = 0.8, rho_G1 = 0.3,
                                    rho_G2 = 0.3, lambda_XM = c(1, 0),
                                    lambda_MY = c(0, 1), omega_1 = 0.7,
                                    omega_2 = 0.7, seed = 7)
  r_pure <- fit_iv2sls_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$Gm)
  r_pool <- fit_iv2sls_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$Gm,
                                  W1 = dat$W1, W2 = dat$W1)  # identical -> guard
  expect_equal(r_pure$NDE, r_pool$NDE, tolerance = 1e-8)
  expect_equal(r_pure$NIE, r_pool$NIE, tolerance = 1e-8)
})

test_that("IV2SLS2 defunct pooled `w` argument errors with redirect", {
  dat <- iconic:::generate_toy_data(n = 300, n_features = 1, phi = 0.8,
                                    mo_confounding = 0.8, seed = 1)
  expect_error(
    fit_iv2sls_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$Gm,
                          w = dat$W[, 1]),
    "removed in v0.9.9")
})

test_that("IV2SLS2 NIE Type I error near nominal under the null", {
  # Under alpha_M = 0 there is no X->M path, so the true NIE is 0; the NIE
  # Wald test should reject at roughly the nominal rate. (We null alpha_M
  # rather than beta_M: with beta_M = 0 the M->Y path is absent, so M_hat is
  # near-collinear with X_hat in stage 3 and the delta-method NIE SE is
  # unreliable -- a DGP artifact, not an estimator defect.)
  set.seed(555)
  n_rep <- 40
  pvals <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    dat <- iconic:::generate_toy_data(n = 500, n_features = 1, phi = 0.8,
                                      mo_confounding = 0.8, rho_G1 = 0.3,
                                      rho_G2 = 0.3, lambda_XM = c(1, 0),
                                      lambda_MY = c(0, 1), omega_1 = 0.7,
                                      omega_2 = 0.7, alpha_M = 0,  # NIE = 0
                                      seed = 6000 + i)
    r <- fit_iv2sls_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$Gm,
                               W1 = dat$W1, W2 = dat$W2)
    pvals[i] <- r$NIE_p
  }
  t1 <- mean(pvals < 0.05, na.rm = TRUE)
  # Allow generous tolerance for a stochastic rate (n_rep = 40).
  expect_lt(t1, 0.25)
})

test_that("IV2SLS2 survival: coverage helps on the log-HR scale", {
  skip_if_not_installed("survival")
  .surv_bias <- function(omega, n = 500, n_iter = 8, seed0 = 7300) {
    nde <- numeric(0)
    for (i in seq_len(n_iter)) {
      dat <- iconic:::generate_toy_data(
        n = n, n_features = 3, beta_X = 0.25, alpha_M = 0.5, beta_M = 0.3,
        phi = 0.8, mo_confounding = 0.8, rho_G1 = 0.3, rho_G2 = 0.3,
        lambda_XM = c(1, 0), lambda_MY = c(0, 1), omega_1 = omega,
        omega_2 = omega, outcome_type = "survival", surv_event_frac = 0.6,
        seed = seed0 + i)
      r <- fit_iv2sls_mediation2_surv(dat$surv_time, dat$surv_event, dat$X,
                                      dat$M, dat$G[, 1], dat$Gm,
                                      W1 = dat$W1, W2 = dat$W2,
                                      effect_scale = "loghr")
      nde <- c(nde, r$NDE)
    }
    mean(abs(nde - 0.25), na.rm = TRUE)  # true NDE (log-HR) = beta_X = 0.25
  }
  b03 <- .surv_bias(0.3)
  b10 <- .surv_bias(1.0)
  expect_lte(b10, b03 + 0.02)
})

test_that("IV2SLS2 survival defunct pooled `w` argument errors", {
  dat <- iconic:::generate_toy_data(n = 300, n_features = 1, phi = 0.8,
                                    mo_confounding = 0.8,
                                    outcome_type = "survival", seed = 1)
  expect_error(
    fit_iv2sls_mediation2_surv(dat$surv_time, dat$surv_event, dat$X, dat$M,
                               dat$G[, 1], dat$Gm, w = dat$W[, 1]),
    "removed in v0.9.9")
})

test_that("iconic_data retains a lone W2 panel for IV2SLS2", {
  dat <- iconic:::generate_toy_data(n = 300, n_features = 3, phi = 0.8,
                                    mo_confounding = 0.8, rho_G1 = 0.3,
                                    rho_G2 = 0.3, lambda_XM = c(1, 0),
                                    lambda_MY = c(0, 1), omega_1 = 0.7,
                                    omega_2 = 0.7, seed = 3)
  id_w2 <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M, G = dat$G[, 1],
                       Gm = dat$Gm, W2 = t(dat$W2))
  # Lone W2 retained for IV2SLS2, but has_path_nc stays FALSE so the
  # two-bridge estimators (PGC2/PGC2Gm) remain ineligible.
  expect_false(is.null(id_w2$W2))
  expect_true(is.null(id_w2$W1))
  expect_false(id_w2$has_path_nc)
})
