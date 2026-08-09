# Tests for the refactoring: format-agnostic estimation functions.
#
# iconic:::.estimate_total_feature() and iconic:::.estimate_mediation_feature() were
# extracted from iconic:::.analyze_feature() and iconic:::.analyze_mediation_feature()
# so that both the simulation pipeline and the real-data pipeline
# (iconic_estimate) share the same per-feature estimation code.
#
# These tests verify that the extracted functions produce identical
# results to the wrapper functions when given the same inputs.
#
# NOTE: Updated for the contract change. .analyze_feature()
# and .analyze_mediation_feature() now pass w = W_mat (the full W matrix,
# n x q) to .estimate_total_feature() / .estimate_mediation_feature(),
# not a single column. The tests below were updated to pass w = dat$W
# (full matrix) to match.

# ── Helper: create a dataset from run_single_iteration ──
.make_refac_dat <- function(n = 100, n_features = 5, seed = 42) {
  run_single_iteration(
    n_synthetic_samples = n, n_features = n_features,
    mo_confounding = 0.8, phi = 0.8,
    omega_1 = 0.7, omega_2 = 0.7,
    seed = seed)
}

# ═══════════════════════════════════════════════════════════════
# Total-effect: .estimate_total_feature vs .analyze_feature
# ═══════════════════════════════════════════════════════════════

test_that(".estimate_total_feature matches .analyze_feature for feature 1", {
  dat <- .make_refac_dat()
  W_avg <- rowMeans(dat$W)

  # Wrapper (simulation format)
  res_wrap <- iconic:::.analyze_feature(dat, f = 1, W_avg = W_avg)

  # Extracted (explicit vectors)
  res_ext <- iconic:::.estimate_total_feature(
    X = dat$X, y = dat$Y[, 1], g = dat$G[, 1], w = dat$W,
    W_mat = dat$W, W_avg = W_avg, covars = dat$synthetic_data,
    methods = iconic:::.methods_all, feature_idx = 1L)

  expect_equal(res_wrap, res_ext, tolerance = 1e-12)
})

test_that(".estimate_total_feature matches .analyze_feature for feature 3", {
  dat <- .make_refac_dat(n_features = 8)
  W_avg <- rowMeans(dat$W)

  res_wrap <- iconic:::.analyze_feature(dat, f = 3, W_avg = W_avg)
  res_ext <- iconic:::.estimate_total_feature(
    X = dat$X, y = dat$Y[, 3], g = dat$G[, 3], w = dat$W,
    W_mat = dat$W, W_avg = W_avg, covars = dat$synthetic_data,
    methods = iconic:::.methods_all, feature_idx = 3L)

  expect_equal(res_wrap, res_ext, tolerance = 1e-12)
})

test_that(".estimate_total_feature with multiple features (loop correctness)", {
  dat <- .make_refac_dat(n_features = 5)
  W_avg <- rowMeans(dat$W)

  # Run via wrapper for all features
  wrap_res <- do.call(rbind, lapply(1:5, function(f)
    iconic:::.analyze_feature(dat, f, W_avg)))

  # Run via extracted function for all features
  ext_res <- do.call(rbind, lapply(1:5, function(f)
    iconic:::.estimate_total_feature(
      X = dat$X, y = dat$Y[, f], g = dat$G[, f], w = dat$W,
      W_mat = dat$W, W_avg = W_avg, covars = dat$synthetic_data,
      methods = iconic:::.methods_all, feature_idx = f)))

  expect_equal(wrap_res, ext_res, tolerance = 1e-12)
})

# ═══════════════════════════════════════════════════════════════
# Mediation: .estimate_mediation_feature vs .analyze_mediation_feature
# ═══════════════════════════════════════════════════════════════

test_that(".estimate_mediation_feature matches .analyze_mediation_feature for feature 1", {
  dat <- .make_refac_dat()
  W_avg <- rowMeans(dat$W)

  res_wrap <- iconic:::.analyze_mediation_feature(dat, f = 1, W_avg = W_avg)
  res_ext <- iconic:::.estimate_mediation_feature(
    X = dat$X, y = dat$Y[, 1], M_vec = dat$M,
    g = dat$G[, 1], gm = dat$Gm,
    w = dat$W, W_mat = dat$W,
    W1_mat = dat$W1, W2_mat = dat$W2,
    W_avg = W_avg, covars = dat$synthetic_data,
    methods = NULL, feature_idx = 1L)

  expect_equal(res_wrap, res_ext, tolerance = 1e-12)
})

test_that(".estimate_mediation_feature matches .analyze_mediation_feature for feature 4", {
  dat <- .make_refac_dat(n_features = 6)
  W_avg <- rowMeans(dat$W)

  res_wrap <- iconic:::.analyze_mediation_feature(dat, f = 4, W_avg = W_avg)
  res_ext <- iconic:::.estimate_mediation_feature(
    X = dat$X, y = dat$Y[, 4], M_vec = dat$M,
    g = dat$G[, 4], gm = dat$Gm,
    w = dat$W, W_mat = dat$W,
    W1_mat = dat$W1, W2_mat = dat$W2,
    W_avg = W_avg, covars = dat$synthetic_data,
    methods = NULL, feature_idx = 4L)

  expect_equal(res_wrap, res_ext, tolerance = 1e-12)
})

test_that(".estimate_mediation_feature with multiple features (loop correctness)", {
  dat <- .make_refac_dat(n_features = 5)
  W_avg <- rowMeans(dat$W)

  wrap_res <- do.call(rbind, lapply(1:5, function(f)
    iconic:::.analyze_mediation_feature(dat, f, W_avg)))

  ext_res <- do.call(rbind, lapply(1:5, function(f)
    iconic:::.estimate_mediation_feature(
      X = dat$X, y = dat$Y[, f], M_vec = dat$M,
      g = dat$G[, f], gm = dat$Gm,
      w = dat$W, W_mat = dat$W,
      W1_mat = dat$W1, W2_mat = dat$W2,
      W_avg = W_avg, covars = dat$synthetic_data,
      methods = NULL, feature_idx = f)))

  expect_equal(wrap_res, ext_res, tolerance = 1e-12)
})

# ═══════════════════════════════════════════════════════════════
# Missing inputs: instrument-based methods return NA without G
# ═══════════════════════════════════════════════════════════════

test_that(".estimate_total_feature without G: instrument methods not present", {
  dat <- .make_refac_dat()
  W_avg <- rowMeans(dat$W)

  res <- iconic:::.estimate_total_feature(
    X = dat$X, y = dat$Y[, 1], g = NULL, w = dat$W,
    W_mat = dat$W, W_avg = W_avg, covars = dat$synthetic_data,
    methods = iconic:::.methods_all, feature_idx = 1L)

  methods_present <- unique(res$method)
  expect_true("UNADJ" %in% methods_present)
  expect_true("COCA" %in% methods_present)
  # DIRECT, IV2SLS, PGC require G
  expect_false("DIRECT" %in% methods_present)
  expect_false("IV2SLS" %in% methods_present)
  expect_false("PGC" %in% methods_present)
})

test_that(".estimate_mediation_feature without G: instrument methods not present", {
  dat <- .make_refac_dat()
  W_avg <- rowMeans(dat$W)

  res <- iconic:::.estimate_mediation_feature(
    X = dat$X, y = dat$Y[, 1], M_vec = dat$M,
    g = NULL, gm = dat$Gm,
    w = dat$W, W_mat = dat$W,
    W1_mat = dat$W1, W2_mat = dat$W2,
    W_avg = W_avg, covars = dat$synthetic_data,
    methods = NULL, feature_idx = 1L)

  methods_present <- unique(res$method)
  expect_true("UNADJ" %in% methods_present)
  expect_true("COCA" %in% methods_present)
  # DIRECT, IV2SLS, PGC, IV2SLS2, PGC2, PGC2Gm all require G
  expect_false("DIRECT" %in% methods_present)
  expect_false("IV2SLS" %in% methods_present)
  expect_false("PGC" %in% methods_present)
  expect_false("IV2SLS2" %in% methods_present)
  expect_false("PGC2" %in% methods_present)
  expect_false("PGC2Gm" %in% methods_present)
})

# ═══════════════════════════════════════════════════════════════
# Missing inputs: NC-based methods return NA without W
# ═══════════════════════════════════════════════════════════════

test_that(".estimate_total_feature without W: NC methods not present", {
  dat <- .make_refac_dat()

  res <- iconic:::.estimate_total_feature(
    X = dat$X, y = dat$Y[, 1], g = dat$G[, 1], w = NULL,
    W_mat = NULL, W_avg = NULL, covars = dat$synthetic_data,
    methods = iconic:::.methods_all, feature_idx = 1L)

  methods_present <- unique(res$method)
  expect_true("UNADJ" %in% methods_present)
  # COCA requires W, DIRECT requires W, PGC requires W_mat.
  # IV2SLS is identified by the instrument alone (W optional), so it runs.
  expect_false("COCA" %in% methods_present)
  expect_false("DIRECT" %in% methods_present)
  expect_true("IV2SLS" %in% methods_present)
  expect_false("PGC" %in% methods_present)
})

test_that(".estimate_mediation_feature without W: NC methods not present", {
  dat <- .make_refac_dat()

  res <- iconic:::.estimate_mediation_feature(
    X = dat$X, y = dat$Y[, 1], M_vec = dat$M,
    g = dat$G[, 1], gm = dat$Gm,
    w = NULL, W_mat = NULL,
    W1_mat = NULL, W2_mat = NULL,
    W_avg = NULL, covars = dat$synthetic_data,
    methods = NULL, feature_idx = 1L)

  methods_present <- unique(res$method)
  expect_true("UNADJ" %in% methods_present)
  # NC-based methods (COCA, DIRECT, PGC, PGC2, PGC2Gm) require W.
  # IV2SLS / IV2SLS2 are identified by the instrument(s) alone (W optional).
  expect_false("COCA" %in% methods_present)
  expect_false("DIRECT" %in% methods_present)
  expect_true("IV2SLS" %in% methods_present)
  expect_false("PGC" %in% methods_present)
  expect_true("IV2SLS2" %in% methods_present)
  expect_false("PGC2" %in% methods_present)
  expect_false("PGC2Gm" %in% methods_present)
})

# ═══════════════════════════════════════════════════════════════
# Consistency: iconic_estimate matches simulation pipeline
# ═══════════════════════════════════════════════════════════════

test_that("iconic_estimate total-effect matches analyze_methods_robust", {
  # Use total-effect only (no mediation)
  dat_te <- run_single_iteration(n_synthetic_samples = 100, n_features = 5,
                                 seed = 99)
  idata <- iconic_data(X = dat_te$X, Y = t(dat_te$Y), W = t(dat_te$W),
                       G = dat_te$G[, 1], covariates = dat_te$synthetic_data,
                       scale = FALSE)
  sim <- analyze_methods_robust(dat_te)
  real <- iconic_estimate(idata)
  for (m in unique(sim$method)) {
    s <- sim$beta[sim$method == m]
    r <- real$beta[real$method == m]
    expect_lt(max(abs(s - r)), 1e-10, label = paste("method", m))
  }
})

test_that("iconic_estimate mediation matches analyze_mediation_robust", {
  dat <- .make_refac_dat(seed = 99)
  # Pass W1 and W2 explicitly so PGC2/PGC2Gm use the path-specific NCs
  # from the DGP, not W1=W2=W (which happens when only W is supplied)
  idata <- iconic_data(X = dat$X, Y = t(dat$Y), M = dat$M,
                       W = t(dat$W), W1 = t(dat$W1), W2 = t(dat$W2),
                       G = dat$G[, 1], Gm = dat$Gm,
                       covariates = dat$synthetic_data, scale = FALSE)
  sim <- analyze_mediation_robust(dat)
  real <- iconic_estimate(idata)
  for (m in unique(sim$method)) {
    s <- sim$NDE[sim$method == m]
    r <- real$NDE[real$method == m]
    expect_lt(max(abs(s - r)), 1e-10, label = paste("NDE method", m))
    s <- sim$NIE[sim$method == m]
    r <- real$NIE[real$method == m]
    expect_lt(max(abs(s - r)), 1e-10, label = paste("NIE method", m))
  }
})

# ═══════════════════════════════════════════════════════════════
# Edge cases
# ═══════════════════════════════════════════════════════════════

test_that(".estimate_total_feature returns NULL with too few complete cases", {
  # Create data with many NAs
  X <- rnorm(50)
  y <- c(rnorm(15), rep(NA, 35))
  res <- iconic:::.estimate_total_feature(
    X = X, y = y, g = rnorm(50), w = rnorm(50),
    W_mat = matrix(rnorm(50 * 5), 50, 5), W_avg = rnorm(50),
    covars = NULL, methods = iconic:::.methods_all, feature_idx = 1L)
  expect_null(res)
})

test_that(".estimate_mediation_feature returns NULL with too few complete cases", {
  X <- rnorm(50)
  y <- c(rnorm(15), rep(NA, 35))
  res <- iconic:::.estimate_mediation_feature(
    X = X, y = y, M_vec = rnorm(50),
    g = rnorm(50), gm = rnorm(50),
    w = rnorm(50), W_mat = matrix(rnorm(50 * 5), 50, 5),
    W1_mat = matrix(rnorm(50 * 5), 50, 5),
    W2_mat = matrix(rnorm(50 * 5), 50, 5),
    W_avg = rnorm(50), covars = NULL,
    methods = NULL, feature_idx = 1L)
  expect_null(res)
})

test_that(".estimate_total_feature respects methods subset", {
  dat <- .make_refac_dat()
  W_avg <- rowMeans(dat$W)
  res <- iconic:::.estimate_total_feature(
    X = dat$X, y = dat$Y[, 1], g = dat$G[, 1], w = dat$W,
    W_mat = dat$W, W_avg = W_avg, covars = dat$synthetic_data,
    methods = c("UNADJ", "IV2SLS"), feature_idx = 1L)
  expect_true(all(res$method %in% c("UNADJ", "IV2SLS")))
  expect_equal(nrow(res), 2)
})
