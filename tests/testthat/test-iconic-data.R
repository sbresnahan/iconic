# Tests for the iconic_data() constructor and validation

test_that("iconic_data constructs from vectors", {
  d <- iconic_data(Z = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50))
  expect_s3_class(d, "iconic_data")
  expect_equal(d$n, 50)
  expect_equal(d$n_features, 3)
  expect_false(d$is_mediation)
  expect_false(d$has_instrument)
  expect_false(d$has_nc)
})

test_that("iconic_data auto-transposes samples x features Y", {
  Y <- matrix(rnorm(50 * 3), 50, 3) # samples x features
  d <- iconic_data(Z = rnorm(50), Y = Y)
  expect_equal(dim(d$Y), c(3, 50)) # stored as features x samples
})

test_that("iconic_data accepts features x samples Y", {
  Y <- matrix(rnorm(50 * 3), 3, 50) # features x samples
  d <- iconic_data(Z = rnorm(50), Y = Y)
  expect_equal(dim(d$Y), c(3, 50))
})

test_that("iconic_data handles mediation mode", {
  d <- iconic_data(Z = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50),
                   M = rnorm(50))
  expect_true(d$is_mediation)
  expect_equal(d$n_mediators, 1)
})

test_that("iconic_data handles matrix mediator", {
  M <- matrix(rnorm(50 * 2), 2, 50) # 2 mediators x 50 samples
  d <- iconic_data(Z = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50), M = M)
  expect_true(d$is_mediation)
  expect_equal(d$n_mediators, 2)
})

test_that("iconic_data sets W1=W2=W when only W supplied", {
  W <- matrix(rnorm(50 * 5), 5, 50)
  d <- iconic_data(Z = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50), W = W)
  expect_true(d$has_nc)
  expect_true(d$has_path_nc)
  expect_identical(d$W1, d$W)
  expect_identical(d$W2, d$W)
})

test_that("iconic_data accepts separate W1 and W2", {
  W1 <- matrix(rnorm(50 * 5), 5, 50)
  W2 <- matrix(rnorm(50 * 5), 5, 50)
  d <- iconic_data(Z = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50),
                   W1 = W1, W2 = W2)
  expect_true(d$has_path_nc)
  # W is derived from W1/W2 so single-panel estimators (COCA, PGC, IV2SLS)
  # and the instrument-strength check have a W matrix to work with.
  expect_true(d$has_nc)
  expect_false(is.null(d$W))
})

test_that("iconic_data rejects n < 20", {
  expect_error(iconic_data(Z = rnorm(10), Y = matrix(rnorm(10 * 3), 3, 10)),
               "At least 20 samples")
})

test_that("iconic_data rejects missing Z", {
  expect_error(iconic_data(Y = matrix(rnorm(50 * 3), 3, 50)), "Z.*required")
})

test_that("iconic_data rejects missing Y", {
  expect_error(iconic_data(Z = rnorm(50)), "Y.*required")
})

test_that("iconic_data encodes covariates", {
  cv <- data.frame(sex = rbinom(50, 1, 0.5), GA = rnorm(50))
  d <- iconic_data(Z = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50),
                   covariates = cv)
  expect_true(ncol(d$covariates) > 0)
})

test_that("iconic_data handles matrix Z (column means)", {
  Zmat <- matrix(rnorm(50 * 3), 3, 50)
  d <- iconic_data(Z = Zmat, Y = matrix(rnorm(50 * 3), 3, 50))
  expect_equal(d$n, 50)
  expect_length(d$Z, 50)
})

test_that("iconic_data handles matrix G (extracts first column)", {
  Gmat <- matrix(rnorm(50 * 5), 50, 5) # n x n_features, as from generate_toy_data
  d <- iconic_data(Z = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50), G = Gmat)
  expect_true(d$has_instrument)
  expect_length(d$G, 50)
  expect_equal(d$G, as.numeric(Gmat[, 1]))
})

test_that("print.iconic_data produces output", {
  d <- iconic_data(Z = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50),
                   G = rnorm(50), W = matrix(rnorm(50 * 5), 5, 50))
  out <- capture.output(print(d))
  expect_true(any(grepl("iconic_data", out)))
  expect_true(any(grepl("50", out)))
})

test_that("as_iconic_data bridges from load_real_input_data", {
  input <- load_real_input_data(example = TRUE)
  d <- as_iconic_data(input)
  expect_s3_class(d, "iconic_data")
  expect_true(d$n > 0)
})
