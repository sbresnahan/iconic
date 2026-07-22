# Tests for ICONIC v0.9.3: composite null hypothesis test (JT-comp)
# Huang (2019) closed-form composite p-value for mediation NIE.
#
# These tests verify:
# 1. .np_cdf() — normal product distribution CDF
# 2. composite_p_value() — JT-comp formula
# 3. .apply_composite_pvalues() — two-pass post-processing
# 4. se_method = "composite" end-to-end via analyze_mediation_robust()
# 5. se_method = "composite" end-to-end via iconic_estimate()

# ── §1: Normal product distribution CDF ──

test_that(".np_cdf returns 1 at z = 0", {
  expect_equal(iconic:::.np_cdf(0), 1)
})

test_that(".np_cdf returns 1 for very small z", {
  expect_equal(iconic:::.np_cdf(1e-11), 1)
})

test_that(".np_cdf returns 0 for very large z", {
  # K0 underflows; the tail probability is negligibly small
  expect_equal(iconic:::.np_cdf(100), 0)
  expect_equal(iconic:::.np_cdf(50), 0)
})

test_that(".np_cdf is monotonically decreasing", {
  z_vals <- c(0.5, 1, 1.5, 2, 3, 5)
  cdf_vals <- iconic:::.np_cdf(z_vals)
  expect_true(all(diff(cdf_vals) < 0))
})

test_that(".np_cdf handles negative z (absolute value)", {
  expect_equal(iconic:::.np_cdf(-2), iconic:::.np_cdf(2))
})

test_that(".np_cdf matches Monte Carlo simulation", {
  # Under H0(1), if Z1, Z2 ~ N(0,1) independent, then
  # P(|Z1*Z2| >= |z|) = F(z).
  set.seed(42)
  n_sim <- 1e6
  z1 <- rnorm(n_sim)
  z2 <- rnorm(n_sim)
  prod <- z1 * z2
  z_test <- 1.96
  mc_prob <- mean(abs(prod) >= z_test)
  np_prob <- iconic:::.np_cdf(z_test)
  # Should match to ~2 decimal places with 1e6 samples
  expect_true(abs(mc_prob - np_prob) < 0.01)
})

test_that(".np_cdf is vectorized", {
  z_vec <- c(0, 1, 2, 5, 100)
  result <- iconic:::.np_cdf(z_vec)
  expect_length(result, 5)
  expect_true(all(result >= 0 & result <= 1))
})

# ── §2: Composite p-value formula ──

test_that("composite_p_value returns 1 under complete null (a=0, b=0)", {
  expect_equal(iconic::composite_p_value(0, 0), 1)
})

test_that("composite_p_value returns 1 under H0(2) (a!=0, b=0)", {
  # When b = 0, ab = 0, so F(0) + F(0) - F(0) = 1 + 1 - 1 = 1
  expect_equal(iconic::composite_p_value(5, 0), 1)
})

test_that("composite_p_value returns 1 under H0(3) (a=0, b!=0)", {
  expect_equal(iconic::composite_p_value(0, 5), 1)
})

test_that("composite_p_value returns small p under strong signal", {
  p <- iconic::composite_p_value(2, 2)
  expect_true(p < 0.05)
  expect_true(p > 0)
})

test_that("composite_p_value is bounded in [0, 1]", {
  # Test various combinations including extreme values
  a_vals <- c(-5, -2, -1, 0, 1, 2, 5, 10)
  b_vals <- c(-5, -2, -1, 0, 1, 2, 5, 10)
  for (a in a_vals) {
    for (b in b_vals) {
      p <- iconic::composite_p_value(a, b)
      expect_true(p >= 0 && p <= 1,
                  info = paste("a =", a, "b =", b, "p =", p))
    }
  }
})

test_that("composite_p_value is symmetric in sign of ab", {
  # p(a, b) should equal p(-a, b) = p(a, -b) = p(-a, -b)
  # because F depends on |ab|
  p1 <- iconic::composite_p_value(2, 3)
  p2 <- iconic::composite_p_value(-2, 3)
  p3 <- iconic::composite_p_value(2, -3)
  p4 <- iconic::composite_p_value(-2, -3)
  expect_equal(p1, p2)
  expect_equal(p1, p3)
  expect_equal(p1, p4)
})

test_that("composite_p_value with Var > 1 is less significant", {
  # Under H0(2)/H0(3), Var > 1, which should increase p-value
  p_var1 <- iconic::composite_p_value(2, 2, var_a = 1, var_b = 1)
  p_var15 <- iconic::composite_p_value(2, 2, var_a = 1.5, var_b = 1.5)
  expect_true(p_var15 >= p_var1)
})

test_that("composite_p_value handles extreme z-statistics", {
  # Very large z-statistics should give very small p-values
  p <- iconic::composite_p_value(10, 10)
  expect_true(p < 1e-5)
  # Should not error or produce NaN
  expect_true(is.finite(p))
})

# ── §3: .apply_composite_pvalues post-processing ──

test_that(".apply_composite_pvalues adds var_a and var_b columns", {
  set.seed(42)
  dat <- iconic::generate_toy_data(n = 200, n_features = 10, n_mediators = 1,
                                   seed = 42)
  res <- iconic:::analyze_mediation_robust(dat, se_method = "delta")
  # Add the z-statistic columns that .apply_composite_pvalues needs
  # (they are already present from the row() function in v0.9.3)
  res2 <- iconic:::.apply_composite_pvalues(res)
  expect_true("var_a" %in% names(res2))
  expect_true("var_b" %in% names(res2))
})

test_that(".apply_composite_pvalues replaces NIE_p but not NDE_p", {
  set.seed(42)
  dat <- iconic::generate_toy_data(n = 200, n_features = 10, n_mediators = 1,
                                   seed = 42)
  res_delta <- iconic:::analyze_mediation_robust(dat, se_method = "delta")
  res_comp  <- iconic:::.apply_composite_pvalues(res_delta)
  # NIE_p should change
  expect_true(!all(res_comp$NIE_p == res_delta$NIE_p, na.rm = TRUE))
  # NDE_p should NOT change
  expect_equal(res_delta$NDE_p, res_comp$NDE_p)
})

test_that(".apply_composite_pvalues preserves NIE_se", {
  set.seed(42)
  dat <- iconic::generate_toy_data(n = 200, n_features = 10, n_mediators = 1,
                                   seed = 42)
  res_delta <- iconic:::analyze_mediation_robust(dat, se_method = "delta")
  res_comp  <- iconic:::.apply_composite_pvalues(res_delta)
  expect_equal(res_delta$NIE_se, res_comp$NIE_se)
})

test_that(".apply_composite_pvalues clamps variance to [1, 1.5]", {
  set.seed(42)
  dat <- iconic::generate_toy_data(n = 200, n_features = 10, n_mediators = 1,
                                   seed = 42)
  res <- iconic:::.apply_composite_pvalues(
    iconic:::analyze_mediation_robust(dat, se_method = "delta"))
  # With a single mediator, Var(a) should be clamped to 1
  expect_true(all(res$var_a >= 1, na.rm = TRUE))
  expect_true(all(res$var_a <= 1.5, na.rm = TRUE))
  expect_true(all(res$var_b >= 1, na.rm = TRUE))
  expect_true(all(res$var_b <= 1.5, na.rm = TRUE))
})

test_that(".apply_composite_pvalues handles NULL and empty input", {
  expect_null(iconic:::.apply_composite_pvalues(NULL))
  empty <- data.frame(feature = integer(0), method = character(0),
                      NIE_p = numeric(0))
  result <- iconic:::.apply_composite_pvalues(empty)
  expect_equal(nrow(result), 0)
})

# ── §4: End-to-end via analyze_mediation_robust ──

test_that("analyze_mediation_robust accepts se_method = composite", {
  expect_true("se_method" %in% names(formals(iconic:::analyze_mediation_robust)))
  set.seed(42)
  dat <- iconic::generate_toy_data(n = 200, n_features = 5, n_mediators = 1,
                                   seed = 42)
  res <- iconic:::analyze_mediation_robust(dat, se_method = "composite")
  expect_s3_class(res, "data.frame")
  expect_true("NIE_p" %in% names(res))
  expect_true("var_a" %in% names(res))
  expect_true("var_b" %in% names(res))
})

test_that("composite NIE_p differs from delta NIE_p", {
  set.seed(42)
  dat <- iconic::generate_toy_data(n = 200, n_features = 10, n_mediators = 1,
                                   seed = 42)
  res_delta <- iconic:::analyze_mediation_robust(dat, se_method = "delta")
  res_comp  <- iconic:::analyze_mediation_robust(dat, se_method = "composite")
  # At least some NIE_p values should differ
  diffs <- res_delta$NIE_p != res_comp$NIE_p
  expect_true(sum(diffs, na.rm = TRUE) > 0)
})

test_that("composite preserves NDE_p and NIE_se", {
  set.seed(42)
  dat <- iconic::generate_toy_data(n = 200, n_features = 10, n_mediators = 1,
                                   seed = 42)
  res_delta <- iconic:::analyze_mediation_robust(dat, se_method = "delta")
  res_comp  <- iconic:::analyze_mediation_robust(dat, se_method = "composite")
  expect_equal(res_delta$NDE_p, res_comp$NDE_p)
  expect_equal(res_delta$NIE_se, res_comp$NIE_se)
})

test_that("composite NIE_significant flag is set", {
  set.seed(42)
  dat <- iconic::generate_toy_data(n = 200, n_features = 10, n_mediators = 1,
                                   seed = 42)
  res <- iconic:::analyze_mediation_robust(dat, se_method = "composite",
                                           alpha = 0.05)
  expect_true("NIE_significant" %in% names(res))
  expect_true(all(res$NIE_significant %in% c(TRUE, FALSE, NA)))
})

# ── §5: End-to-end via iconic_estimate ──

test_that("iconic_estimate accepts se_method = composite", {
  set.seed(42)
  n <- 200; p <- 10; q <- 10
  Z <- rnorm(n)
  Y <- matrix(rnorm(n * p), p, n)
  G <- rnorm(n)
  W <- matrix(rnorm(n * q), q, n)
  M <- matrix(rnorm(n), 1, n)
  idat <- iconic::iconic_data(Z = Z, Y = Y, G = G, W = W, M = M)
  est <- iconic::iconic_estimate(idat, se_method = "composite", n_cores = 1)
  expect_s3_class(est, "data.frame")
  expect_true("NIE_p" %in% names(est))
  expect_true("var_a" %in% names(est))
  expect_true("var_b" %in% names(est))
})

test_that("iconic_estimate composite preserves NDE_p", {
  set.seed(42)
  n <- 200; p <- 10; q <- 10
  Z <- rnorm(n)
  Y <- matrix(rnorm(n * p), p, n)
  G <- rnorm(n)
  W <- matrix(rnorm(n * q), q, n)
  M <- matrix(rnorm(n), 1, n)
  idat <- iconic::iconic_data(Z = Z, Y = Y, G = G, W = W, M = M)
  est_delta <- iconic::iconic_estimate(idat, se_method = "delta", n_cores = 1)
  est_comp  <- iconic::iconic_estimate(idat, se_method = "composite", n_cores = 1)
  expect_equal(est_delta$NDE_p, est_comp$NDE_p)
})

# ── §6: Type I error under point null ──

test_that("composite type I error is reasonable under point null", {
  # Under the true point null (alpha_M = 0, beta_M = 0, no confounding),
  # the composite test should have type I error close to nominal.
  # We use a small number of reps for speed; the check is that the
  # rate is not wildly inflated (e.g., < 0.15 at alpha = 0.05).
  set.seed(123)
  n_reps <- 20
  alpha <- 0.05
  sig_count <- 0
  n_total <- 0
  for (r in seq_len(n_reps)) {
    dat <- iconic::generate_toy_data(n = 200, mo_confounding = 0,
                                     seed = r + 5000, n_features = 20,
                                     n_mediators = 1, alpha_M = 0, beta_M = 0)
    res <- iconic:::analyze_mediation_robust(dat, se_method = "composite",
                                             alpha = alpha)
    sub <- res[res$method == "UNADJ", ]
    sig_count <- sig_count + sum(sub$NIE_significant, na.rm = TRUE)
    n_total <- n_total + sum(!is.na(sub$NIE_p))
  }
  rate <- sig_count / n_total
  # Should be below 0.15 (allowing for Monte Carlo noise and the known
  # slight inflation of JT-comp at small p-values)
  expect_true(rate < 0.15,
              info = paste("Type I error rate =", round(rate, 4)))
})

# ── §7: Real case study design (scalar Y, multiple mediators) ──
#
# ICONIC's real case studies (GDM -> placental isoforms -> birth weight;
# smoking -> tumor expression -> survival) have a scalar outcome Y and a
# panel of mediators M.  The composite test groups by method only, so
# Var(a) and Var(b) are estimated across mediators within each method.
# Each mediator has its own stage-1 regression M_m ~ Z, giving real
# variation in a = alpha_M / SE(alpha_M).

test_that(".apply_composite_pvalues groups by method, not method x mediator", {
  # With n_features = 1 and n_mediators = 20, each method has 20 rows
  # (one per mediator).  Grouping by method gives 20 tests per group,
  # enough to estimate Var(a) > 0.  Grouping by (method, mediator)
  # would give 1 test per group, falling back to Var = 1.
  set.seed(42)
  dat <- iconic::generate_toy_data(n = 200, n_features = 1, n_mediators = 20,
                                   alpha_M = 0.1, beta_M = 0.05, seed = 42)
  res <- iconic:::analyze_mediation_robust(dat, se_method = "delta")
  res_comp <- iconic:::.apply_composite_pvalues(res)
  # With 20 mediators, at least some methods should have var_a > 1
  # (i.e., not all clamped to the fallback value of 1)
  expect_true(any(res_comp$var_a > 1, na.rm = TRUE))
})

test_that("composite NIE_p differs from delta with multiple mediators", {
  set.seed(42)
  dat <- iconic::generate_toy_data(n = 200, n_features = 1, n_mediators = 20,
                                   alpha_M = 0.1, beta_M = 0.05, seed = 42)
  res_delta <- iconic:::analyze_mediation_robust(dat, se_method = "delta")
  res_comp  <- iconic:::analyze_mediation_robust(dat, se_method = "composite")
  diffs <- res_delta$NIE_p != res_comp$NIE_p
  expect_true(sum(diffs, na.rm = TRUE) > 0)
})

test_that("composite type I error under point null with multiple mediators", {
  # Real design: scalar Y, 20 mediators, point null
  set.seed(999)
  n_reps <- 20
  alpha <- 0.05
  sig_count <- 0
  n_total <- 0
  for (r in seq_len(n_reps)) {
    dat <- iconic::generate_toy_data(n = 200, mo_confounding = 0,
                                     seed = r + 9000, n_features = 1,
                                     n_mediators = 20, alpha_M = 0, beta_M = 0)
    res <- iconic:::analyze_mediation_robust(dat, se_method = "composite",
                                             alpha = alpha)
    sub <- res[res$method == "UNADJ", ]
    sig_count <- sig_count + sum(sub$NIE_significant, na.rm = TRUE)
    n_total <- n_total + sum(!is.na(sub$NIE_p))
  }
  rate <- sig_count / n_total
  expect_true(rate < 0.15,
              info = paste("Type I error rate (multi-mediator) =", round(rate, 4)))
})
