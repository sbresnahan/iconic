test_that("generate_toy_data returns correct structure", {
  dat <- iconic:::generate_toy_data(n = 100, n_features = 5, seed = 1)
  expect_true(all(c("X", "G", "Y", "W", "U1", "M", "synthetic_data",
                      "true_total", "true_NDE", "true_NIE") %in% names(dat)))
  expect_length(dat$X, 100)
  expect_equal(dim(dat$Y), c(100, 5))
  expect_equal(dim(dat$W), c(100, 5))
  expect_equal(dat$true_total, 0.10 + 0.50 * 0.30)
  expect_equal(dat$true_NDE, 0.10)
  expect_equal(dat$true_NIE, 0.50 * 0.30)
})

test_that("fit_direct returns named list with numerics", {
  dat <- iconic:::generate_toy_data(n = 200, seed = 42)
  res <- fit_direct(dat$Y[, 1], dat$X, dat$G[, 1], dat$W[, 1])
  expect_named(res, c("beta", "se", "pvalue"))
  expect_true(is.numeric(res$beta))
  expect_true(res$pvalue >= 0 && res$pvalue <= 1)
})

test_that("fit_coca returns named list with numerics", {
  dat <- iconic:::generate_toy_data(n = 200, seed = 42)
  res <- fit_coca(dat$Y[, 1], dat$X, rowMeans(dat$W))
  expect_named(res, c("beta", "se", "pvalue"))
  # May return NA if unstable — that's OK, just check structure
  expect_true(is.numeric(res$beta) || is.na(res$beta))
})

test_that("fit_iv2sls returns named list with numerics", {
  dat <- iconic:::generate_toy_data(n = 300, seed = 42)
  res <- fit_iv2sls(dat$Y[, 1], dat$X, dat$G[, 1], dat$W[, 1])
  expect_named(res, c("beta", "se", "pvalue"))
  expect_true(is.numeric(res$beta) || is.na(res$beta))
  if (!is.na(res$pvalue)) {
    expect_true(res$pvalue >= 0 && res$pvalue <= 1)
  }
})

test_that("fit_pgc (matrix bridge) returns named list with numerics", {
  dat <- iconic:::generate_toy_data(n = 200, seed = 42)
  res <- fit_pgc(dat$Y[, 1], dat$X, dat$G[, 1], dat$W)
  expect_named(res, c("beta", "se", "pvalue"))
  expect_true(is.numeric(res$beta) || is.na(res$beta))
})

test_that("fit_pgc_scalar returns named list with numerics", {
  dat <- iconic:::generate_toy_data(n = 200, seed = 42)
  res <- fit_pgc_scalar(dat$Y[, 1], dat$X, dat$G[, 1], rowMeans(dat$W))
  expect_named(res, c("beta", "se", "pvalue"))
  expect_true(is.numeric(res$beta) || is.na(res$beta))
})

test_that("fit_pgc_scalar is approximately equal to IV2SLS", {
  # The scalar bridge is algebraically equivalent to IV/2SLS
  dat <- iconic:::generate_toy_data(n = 500, seed = 99)
  pgc_s <- fit_pgc_scalar(dat$Y[, 1], dat$X, dat$G[, 1], rowMeans(dat$W))
  iv <- fit_iv2sls(dat$Y[, 1], dat$X, dat$G[, 1], dat$W[, 1])
  if (!is.na(pgc_s$beta) && !is.na(iv$beta)) {
    expect_lt(abs(pgc_s$beta - iv$beta), 0.02)
  }
})

test_that("run_simulation returns correct structure with five methods", {
  res <- run_simulation(n_iter = 5, n_samples = 100, n_features = 3)
  expect_named(res, c("raw", "summary", "iter_bias", "true_total", "params"))
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC") %in% res$summary$method))
  expect_equal(res$true_total, 0.10 + 0.50 * 0.30)
})

test_that("run_null_sim Type I error structure is correct with five methods", {
  res <- run_null_sim(n_iter = 10, n_samples = 100, n_features = 3)
  expect_named(res, c("rates", "raw"))
  expect_true("type1_error" %in% names(res$rates))
  expect_true(all(res$rates$type1_error >= 0 & res$rates$type1_error <= 1))
  expect_equal(nrow(res$rates), 5)
})

test_that("IV2SLS is approximately unbiased under large N", {
  # Under large n with strong instrument, IV2SLS should recover true_total
  res <- run_simulation(n_iter = 30, n_samples = 1000, n_features = 5,
                        beta_X = 0.1, alpha_M = 0.5, beta_M = 0.3,
                        conf_str = 0.8)
  iv_bias <- res$summary$bias[res$summary$method == "IV2SLS"]
  expect_lt(abs(iv_bias), 0.05)
})
