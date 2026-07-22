# Ground truth regression test (v0.8.2)
#
# BLOCKING TEST: verifies that the imposed true_NDE and true_NIE are
# recoverable by an oracle estimator (OLS with the known confounder U)
# as n grows, at feat_cor = 0, 0.5, 0.8, and with GAN-learned correlations.
#
# If this test fails, every bias/RMSE/tipping-point measurement in the
# package is against a moving target and the feature must be reverted.
#
# Oracle estimator:
#   NDE = coef(lm(Y_f ~ Z + M + U1))[["Z"]]         (= beta_Z)
#   alpha_M = coef(lm(M ~ Z + U1))[["Z"]]            (= alpha_M)
#   beta_M = coef(lm(Y_f ~ Z + M + U1))[["M"]]      (= beta_M)
#   NIE = alpha_M * beta_M                            (= alpha_M * beta_M)
#
# The oracle is correct by construction under the linear DGP because
# conditioning on U (the true confounder) removes all confounding.
# Correlated noise is additive (MVN noise added to Y and W) and does
# not enter the structural equations for M or the coefficients beta_Z,
# alpha_M, beta_M.

# ── Helper: compute oracle bias for one dataset ──
.oracle_bias <- function(dat, n_features) {
  nde_bias <- numeric(n_features)
  nie_bias <- numeric(n_features)
  U1 <- dat$U1
  if (is.null(U1)) U1 <- dat$U[, 1]  # run_single_iteration stores U as matrix
  for (f in seq_len(n_features)) {
    df_y <- data.frame(y = dat$Y[, f], Z = dat$Z, M = dat$M, U1 = U1)
    fit_y <- lm(y ~ Z + M + U1, data = df_y)
    nde_bias[f] <- coef(fit_y)[["Z"]] - dat$true_NDE
    df_m <- data.frame(M = dat$M, Z = dat$Z, U1 = U1)
    fit_m <- lm(M ~ Z + U1, data = df_m)
    nie_bias[f] <- coef(fit_m)[["Z"]] * coef(fit_y)[["M"]] - dat$true_NIE
  }
  c(NDE_bias = mean(nde_bias), NIE_bias = mean(nie_bias))
}

# ═══════════════════════════════════════════════════════════════
# Test 1: Oracle recovers true_NDE/true_NIE at feat_cor = 0
# ═══════════════════════════════════════════════════════════════

test_that("oracle recovers true effects at feat_cor = 0, n = 10000", {
  dat <- generate_toy_data(n = 10000, n_features = 10, feat_cor = 0,
                           mo_confounding = 0.8, phi = 0.8, seed = 42)
  bias <- .oracle_bias(dat, 10)
  expect_lt(abs(bias["NDE_bias"]), 0.01)
  expect_lt(abs(bias["NIE_bias"]), 0.01)
})

# ═══════════════════════════════════════════════════════════════
# Test 2: Oracle recovers true_NDE/true_NIE at feat_cor = 0.5
# ═══════════════════════════════════════════════════════════════

test_that("oracle recovers true effects at feat_cor = 0.5, n = 10000", {
  dat <- generate_toy_data(n = 10000, n_features = 10, feat_cor = 0.5,
                           mo_confounding = 0.8, phi = 0.8, seed = 42)
  bias <- .oracle_bias(dat, 10)
  expect_lt(abs(bias["NDE_bias"]), 0.01)
  expect_lt(abs(bias["NIE_bias"]), 0.01)
})

# ═══════════════════════════════════════════════════════════════
# Test 3: Oracle recovers true_NDE/true_NIE at feat_cor = 0.8
# ═══════════════════════════════════════════════════════════════

test_that("oracle recovers true effects at feat_cor = 0.8, n = 10000", {
  dat <- generate_toy_data(n = 10000, n_features = 10, feat_cor = 0.8,
                           mo_confounding = 0.8, phi = 0.8, seed = 42)
  bias <- .oracle_bias(dat, 10)
  expect_lt(abs(bias["NDE_bias"]), 0.01)
  expect_lt(abs(bias["NIE_bias"]), 0.01)
})

# ═══════════════════════════════════════════════════════════════
# Test 4: Bias decreases as n grows (convergence to zero)
# ═══════════════════════════════════════════════════════════════

test_that("oracle bias decreases as n grows (feat_cor = 0.5)", {
  bias_500 <- .oracle_bias(
    generate_toy_data(n = 500, n_features = 10, feat_cor = 0.5,
                      mo_confounding = 0.8, phi = 0.8, seed = 42), 10)
  bias_10000 <- .oracle_bias(
    generate_toy_data(n = 10000, n_features = 10, feat_cor = 0.5,
                      mo_confounding = 0.8, phi = 0.8, seed = 42), 10)
  expect_lt(abs(bias_10000["NDE_bias"]), abs(bias_500["NDE_bias"]))
  expect_lt(abs(bias_10000["NIE_bias"]), abs(bias_500["NIE_bias"]))
})

# ═══════════════════════════════════════════════════════════════
# Test 5: GAN-learned correlations do not perturb ground truth
# ═══════════════════════════════════════════════════════════════

test_that("GAN-learned correlations do not perturb ground truth (n = 50000)", {
  skip_if_not_installed("MASS")
  real_data <- load_real_input_data(example = TRUE)
  trained_gan <- train_gan_on_real_data(real_data$gan_training_data,
                                         feature_correlations = real_data$feature_correlations,
                                         verbose = FALSE)
  dat <- run_single_iteration(trained_gan, n_synthetic_samples = 50000,
                              n_features = 30, n_confounders = 1,
                              mo_confounding = 0.8, phi = 0.8, feat_cor = 0,
                              seed = 42)
  bias <- .oracle_bias(dat, 30)
  expect_lt(abs(bias["NDE_bias"]), 0.01)
  expect_lt(abs(bias["NIE_bias"]), 0.01)
})

# ═══════════════════════════════════════════════════════════════
# Test 6: Feature correlations do not change oracle bias vs no correlations
# ═══════════════════════════════════════════════════════════════

test_that("feature correlations do not increase oracle bias vs feat_cor = 0", {
  bias_0 <- .oracle_bias(
    generate_toy_data(n = 10000, n_features = 10, feat_cor = 0,
                      mo_confounding = 0.8, phi = 0.8, seed = 42), 10)
  bias_8 <- .oracle_bias(
    generate_toy_data(n = 10000, n_features = 10, feat_cor = 0.8,
                      mo_confounding = 0.8, phi = 0.8, seed = 42), 10)
  # The difference in |bias| should be small (within Monte Carlo noise)
  expect_lt(abs(abs(bias_8["NDE_bias"]) - abs(bias_0["NDE_bias"])), 0.01)
  expect_lt(abs(abs(bias_8["NIE_bias"]) - abs(bias_0["NIE_bias"])), 0.01)
})
