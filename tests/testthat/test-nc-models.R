test_that("built-in NC models return n x n_features matrices", {
  n <- 50; k <- 2; p <- 4
  U <- matrix(rnorm(n * k), n, k)
  cv <- data.frame(row.names = seq_len(n))

  Wp <- nc_proxy(U, cv, list(n_features = p, coverage = 0.7))
  Wc <- nc_cpg(U,   cv, list(n_features = p, coverage = 0.7))
  expect_equal(dim(Wp), c(n, p))
  expect_equal(dim(Wc), c(n, p))
  expect_true(all(is.finite(Wp)))
  expect_true(all(is.finite(Wc)))
})

test_that("nc_proxy coverage controls confounder signal captured", {
  set.seed(1)
  n <- 4000
  U <- matrix(rnorm(n), n, 1)
  hi <- nc_proxy(U, data.frame(row.names = seq_len(n)),
                 list(n_features = 1, coverage = 1, noise_sd = 0))
  lo <- nc_proxy(U, data.frame(row.names = seq_len(n)),
                 list(n_features = 1, coverage = 0, noise_sd = 0))
  expect_gt(abs(cor(hi[, 1], U[, 1])), abs(cor(lo[, 1], U[, 1])))
})

test_that("nc_proxy mode='shared' gives all columns the same signal", {
  set.seed(1)
  n <- 500; k <- 3; p <- 6
  U <- matrix(rnorm(n * k), n, k)
  W <- nc_proxy(U, data.frame(row.names = seq_len(n)),
                list(n_features = p, coverage = 1, noise_sd = 0, mode = "shared"))
  # All columns carry rowMeans(U) -> high inter-column correlation
  cors <- cor(W)
  expect_true(all(cors[lower.tri(cors)] > 0.8))
})

test_that("nc_proxy mode='distinct' gives columns different confounder signals", {
  set.seed(1)
  n <- 2000; k <- 3; p <- 6
  U <- matrix(rnorm(n * k), n, k)
  W <- nc_proxy(U, data.frame(row.names = seq_len(n)),
                list(n_features = p, coverage = 1, noise_sd = 0.1, mode = "distinct"))
  # Column f captures confounder ((f-1) %% k) + 1
  # So columns 1,4 capture U1; 2,5 capture U2; 3,6 capture U3
  expect_gt(abs(cor(W[, 1], U[, 1])), 0.5)
  expect_gt(abs(cor(W[, 2], U[, 2])), 0.5)
  expect_gt(abs(cor(W[, 3], U[, 3])), 0.5)
  # Cross-confounder correlations should be low
  expect_lt(abs(cor(W[, 1], U[, 2])), 0.2)
  expect_lt(abs(cor(W[, 2], U[, 3])), 0.2)
})

test_that("nc_proxy mode='distinct' with fewer columns than confounders is under-identified", {
  set.seed(1)
  n <- 1000; k <- 3; p <- 2
  U <- matrix(rnorm(n * k), n, k)
  W <- nc_proxy(U, data.frame(row.names = seq_len(n)),
                list(n_features = p, coverage = 1, noise_sd = 0.1, mode = "distinct"))
  # Only 2 columns -> captures confounders 1 and 2, misses confounder 3
  expect_gt(abs(cor(W[, 1], U[, 1])), 0.5)
  expect_gt(abs(cor(W[, 2], U[, 2])), 0.5)
  # No column captures confounder 3
  expect_lt(max(abs(cor(W[, 1], U[, 3])), abs(cor(W[, 2], U[, 3]))), 0.2)
})

test_that("registry and validator work; unknown model errors", {
  expect_true(all(c("proxy", "cpg") %in% list_nc_models()))
  expect_true(iconic:::.validate_nc_model("proxy"))
  expect_true(iconic:::.validate_nc_model(nc_cpg))
  expect_error(iconic:::.resolve_nc_model("nope"))
})

test_that("a user-supplied NC function honouring the contract is accepted", {
  my_nc <- function(U, covariates, params)
    matrix(rowMeans(U), nrow(U), params$n_features) + rnorm(nrow(U) * params$n_features)
  expect_true(iconic:::.validate_nc_model(my_nc))
})

test_that("simulate_single_genetic_instrument returns a standardised instrument", {
  gi <- simulate_single_genetic_instrument(500, seed = 1)
  expect_length(gi$G, 500)
  expect_lt(abs(mean(gi$G)), 0.05)
  expect_lt(abs(sd(gi$G) - 1), 0.05)
})
