test_that("train_gan_on_real_data errors without torch", {
  inp <- load_real_input_data(example = TRUE)
  if (!check_torch_setup()) {
    expect_error(train_gan_on_real_data(inp$gan_training_data, epochs = 5, verbose = FALSE),
                 "torch is required")
  } else {
    gan <- train_gan_on_real_data(inp$gan_training_data, epochs = 5, verbose = FALSE)
    expect_s3_class(gan, "iconic_gan")
    expect_equal(gan$model_type, "gan")
    tex <- sample_texture(gan, 20)
    expect_equal(nrow(tex), 20)
    expect_true(all(names(inp$gan_training_data) %in% names(tex)))
  }
})

test_that("run_single_iteration matches the generate_toy_data contract", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 300, n_features = 8, seed = 1)
  expect_true(all(c("Z", "G", "Y", "W", "U1", "M", "synthetic_data", "true_total") %in%
                    names(dat)))
  expect_equal(dim(dat$Y), c(300, 8))
  expect_equal(dim(dat$W), c(300, 8))
  expect_equal(dim(dat$G), c(300, 8))
  expect_length(dat$Z, 300)
  # single instrument replicated across columns
  expect_equal(dat$G[, 1], dat$G[, 8])
})

test_that("true_total is the closed-form causal effect and flows through the estimators", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 300, n_features = 5,
                              beta_Z = 0.1, alpha_M = 0.5, beta_M = 0.3, seed = 2)
  expect_equal(dat$true_total, 0.1 + 0.5 * 0.3)

  res <- analyze_methods_robust(dat)
  expect_true(all(c("feature", "method", "beta", "se", "pvalue", "significant") %in%
                    names(res)))
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC") %in% res$method))
})

test_that("effect_size = 0 gives a null (zero) true total", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 200, n_features = 4,
                              effect_size = 0, seed = 3)
  expect_equal(dat$true_total, 0)
})

test_that("column alignment holds when the NC model returns a different width", {
  narrow_nc <- function(U, covariates, params) # returns only 2 controls
    matrix(rowMeans(U), nrow(U), 2)
  dat <- run_single_iteration(NULL, n_synthetic_samples = 150, n_features = 6,
                              nc_model = narrow_nc, seed = 4)
  expect_equal(ncol(dat$W), 6) # recycled to match Y
  expect_silent(iconic:::run_methods(dat, 6))
})

test_that("IV2SLS is approximately unbiased on synthetic data", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 1500, n_features = 8,
                              conf_strength = 0.8, seed = 5)
  res <- iconic:::run_methods(dat, 8)
  iv <- mean(res$beta[res$method == "IV2SLS"], na.rm = TRUE)
  expect_lt(abs(iv - dat$true_total), 0.06)
})

test_that("torch GAN trains when available (skipped otherwise)", {
  skip_on_cran()
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  inp <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(inp$gan_training_data,
                                feature_correlations = inp$feature_correlations,
                                feature_texture = inp$feature_texture,
                                epochs = 5, verbose = FALSE)
  expect_equal(gan$model_type, "gan")
  expect_equal(nrow(sample_texture(gan, 10)), 10)
  # feature_texture should be stored
  expect_s3_class(gan$feature_texture, "iconic_feature_texture")
})

# ---- Binary-column handling tests ------------------------------------------

test_that("binary columns are detected and stored in the iconic_gan object", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  inp <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(inp$gan_training_data, epochs = 5, verbose = FALSE)
  expect_true("binary_cols" %in% names(gan))
  expect_true("onehot_groups" %in% names(gan))
  # sex is 0/1 in the example data, so it should be flagged as binary
  expect_true("sex" %in% gan$binary_cols)
  # one-hot ethnicity dummies should also be flagged as binary
  eth_cols <- grep("^mother_ethnicity", gan$binary_cols, value = TRUE)
  expect_true(length(eth_cols) >= 1)
})

test_that("sample_texture returns 0/1 values for binary columns", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  inp <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(inp$gan_training_data, epochs = 5, verbose = FALSE)
  tex <- sample_texture(gan, 200)
  for (cl in gan$binary_cols) {
    vals <- unique(tex[[cl]])
    expect_true(all(vals %in% c(0, 1)),
                info = paste("Column", cl, "has values:", paste(vals, collapse = ", ")))
  }
})

test_that("one-hot dummy groups are mutually exclusive in synthetic draws", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  inp <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(inp$gan_training_data, epochs = 5, verbose = FALSE)
  tex <- sample_texture(gan, 500)
  for (grp in gan$onehot_groups) {
    present <- intersect(grp, names(tex))
    if (length(present) >= 2) {
      rowsums <- rowSums(tex[, present, drop = FALSE])
      # Each row should have at most one active dummy (0 or 1, never 2+)
      expect_true(all(rowsums <= 1),
                  info = paste("Group", paste(grp, collapse = "/"),
                               "has rows with multiple active dummies"))
    }
  }
})

test_that("plot_gan_diagnostics handles binary columns without error", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  inp <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(inp$gan_training_data, epochs = 5, verbose = FALSE)
  # Should produce a plot with both bar charts (binary) and density overlays (continuous)
  png(tempfile(fileext = ".png"))
  expect_silent(plot_gan_diagnostics(gan, inp$gan_training_data))
  dev.off()
})

test_that("continuous columns are not flagged as binary", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  inp <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(inp$gan_training_data, epochs = 5, verbose = FALSE)
  # exposure_level, outcome_level, and mediator_level are continuous
  expect_false("exposure_level" %in% gan$binary_cols)
  expect_false("outcome_level" %in% gan$binary_cols)
  expect_false("mediator_level" %in% gan$binary_cols)
  # GA is z-scored continuous
  expect_false("GA" %in% gan$binary_cols)
})

test_that("backward compatibility: sample_texture works with gan lacking binary_cols", {
  skip_on_cran()
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  # Simulate an old iconic_gan object without the new fields
  inp <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(inp$gan_training_data, epochs = 5, verbose = FALSE)
  gan_old <- gan
  gan_old$binary_cols <- NULL
  gan_old$onehot_groups <- NULL
  # Should still work (no rounding, but no error)
  tex <- expect_silent(sample_texture(gan_old, 20))
  expect_equal(nrow(tex), 20)
})
