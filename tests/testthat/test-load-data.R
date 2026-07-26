test_that("load_real_input_data returns the expected structure on example data", {
  inp <- load_real_input_data(example = TRUE)
  expect_named(inp, c("gan_training_data", "original_matrices", "covariates",
                      "feature_correlations", "feature_texture",
                      "sample_names", "feature_names", "n_samples", "n_features"))
  expect_true(all(c("exposure_level", "outcome_level", "mediator_level") %in%
                    names(inp$gan_training_data)))
  expect_equal(nrow(inp$gan_training_data), inp$n_samples)
  expect_true(all(vapply(inp$gan_training_data, is.numeric, logical(1))))
  # feature_texture should be present when M_matrix is supplied
  expect_s3_class(inp$feature_texture, "iconic_feature_texture")
  # M should NOT be in feature_correlations (replaced by copula)
  expect_null(inp$feature_correlations$M)
})

test_that("covariate encoding produces numeric columns and one-hot ethnicity", {
  inp <- load_real_input_data(example = TRUE)
  cv <- inp$covariates
  expect_true(all(vapply(cv, is.numeric, logical(1))))
  expect_true("sex" %in% names(cv))
  expect_true(any(grepl("^mother_ethnicity", names(cv)))) # one-hot dummies
  expect_false("sample_id" %in% names(cv))
})

test_that("covariate names colliding with reserved tokens are renamed", {
  n <- 40
  Zm <- matrix(rnorm(5 * n), 5, n)
  Ym <- matrix(rnorm(5 * n), 5, n)
  cov <- data.frame(Z = rnorm(n), g = rnorm(n), ok = rnorm(n)) # Z, g are reserved
  enc <- load_real_input_data(Z_matrix = Zm, Y_matrix = Ym, covariates_df = cov)$covariates
  expect_false(any(names(enc) %in% c("Z", "g")))
  expect_true(all(c("cov_Z", "cov_g", "ok") %in% names(enc)))
})

test_that("missing Y_matrix is an error", {
  expect_error(load_real_input_data(Z_matrix = matrix(1, 2, 2), example = FALSE),
               "required")
})
