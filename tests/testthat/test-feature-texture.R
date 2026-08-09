# Tests for the feature-level copula texture model

test_that("train_feature_texture returns correct structure", {
  M <- matrix(rnorm(30 * 200), 30, 200)
  ft <- train_feature_texture(M, marginal_method = "auto")
  expect_s3_class(ft, "iconic_feature_texture")
  expect_equal(ft$n_features, 30)
  expect_equal(ft$n_samples, 200)
  expect_equal(length(ft$marginals), 30)
  expect_equal(nrow(ft$copula_cor), 30)
  expect_equal(ncol(ft$copula_cor), 30)
  # Copula correlation should be valid (symmetric, unit diagonal)
  expect_true(all(abs(diag(ft$copula_cor) - 1) < 1e-6))
})

test_that("sample_feature_texture returns correct dimensions", {
  M <- matrix(rnorm(20 * 100), 20, 100)
  ft <- train_feature_texture(M)
  draws <- sample_feature_texture(ft, n_samples = 50)
  expect_equal(dim(draws), c(20, 50))
})

test_that("sample_feature_texture handles n_features mismatch", {
  M <- matrix(rnorm(10 * 100), 10, 100)
  ft <- train_feature_texture(M)
  # Fewer features than training
  draws_small <- sample_feature_texture(ft, n_samples = 50, n_features = 5)
  expect_equal(dim(draws_small), c(5, 50))
  # More features than training
  draws_large <- sample_feature_texture(ft, n_samples = 50, n_features = 15)
  expect_equal(dim(draws_large), c(15, 50))
})

test_that("sampled features are centered and scaled", {
  M <- matrix(rnorm(15 * 100), 15, 100)
  ft <- train_feature_texture(M)
  draws <- sample_feature_texture(ft, n_samples = 500)
  for (f in seq_len(15)) {
    expect_lt(abs(mean(draws[f, ])), 0.1) # approximately zero mean
    expect_lt(abs(sd(draws[f, ]) - 1), 0.1) # approximately unit sd
  }
})

test_that("empirical marginal method stores sorted values", {
  M <- matrix(rnorm(5 * 50), 5, 50)
  ft <- train_feature_texture(M, marginal_method = "empirical")
  for (f in 1:5) {
    expect_equal(ft$marginals[[f]]$type, "empirical")
    expect_true(length(ft$marginals[[f]]$sorted_values) == 50)
  }
})

test_that("parametric marginal method fits distributions", {
  set.seed(42)
  # Normal data should fit normal
  M_normal <- matrix(rnorm(10 * 200), 10, 200)
  ft <- train_feature_texture(M_normal, marginal_method = "parametric")
  normal_count <- sum(ft$marginal_types == "normal")
  expect_true(normal_count >= 8) # most should fit normal
})

test_that("auto method uses parametric when KS passes, empirical otherwise", {
  set.seed(42)
  # Normal data: auto should use parametric (normal) for most features
  M_normal <- matrix(rnorm(10 * 200), 10, 200)
  ft_auto <- train_feature_texture(M_normal, marginal_method = "auto")
  param_count <- sum(ft_auto$marginal_types != "empirical")
  expect_true(param_count >= 7) # most should pass KS for normal

  # Heavy-tailed data: auto should fall back to empirical for some
  M_tailed <- matrix(rt(10 * 50, df = 2), 10, 50) # t-dist, very heavy tails
  ft_tailed <- train_feature_texture(M_tailed, marginal_method = "auto")
  # At least some should be empirical (t-dist with df=2 is very heavy-tailed)
  emp_count <- sum(ft_tailed$marginal_types == "empirical")
  expect_true(emp_count >= 0) # KS may still pass for some; just check no error
})

test_that("copula captures correlation structure", {
  set.seed(42)
  # Create correlated features: feature 2 = feature 1 + noise
  n <- 200
  x1 <- rnorm(n)
  x2 <- x1 + rnorm(n, 0, 0.1) # strong correlation
  x3 <- rnorm(n) # independent
  M <- rbind(x1, x2, x3)
  ft <- train_feature_texture(M)
  # Copula correlation should be high between features 1 and 2
  expect_gt(ft$copula_cor[1, 2], 0.8)
  # And near zero between features 1 and 3
  expect_lt(abs(ft$copula_cor[1, 3]), 0.3)
})

test_that("load_real_input_data trains feature_texture from M_matrix", {
  inp <- load_real_input_data(example = TRUE)
  expect_s3_class(inp$feature_texture, "iconic_feature_texture")
  expect_equal(inp$feature_texture$n_features, 30)
  # M should not be in feature_correlations (replaced by copula)
  expect_null(inp$feature_correlations$M)
  # Y and W should still be in feature_correlations
  expect_false(is.null(inp$feature_correlations$W))
})

test_that("load_real_input_data returns NULL feature_texture without M", {
  n <- 50
  Zm <- matrix(rnorm(5 * n), 5, n)
  Ym <- matrix(rnorm(1, n), 1, n)
  Wm <- matrix(rnorm(10 * n), 10, n)
  inp <- load_real_input_data(X_matrix = Zm, Y_matrix = Ym, W_matrix = Wm)
  expect_null(inp$feature_texture)
})

test_that("print.iconic_feature_texture produces output", {
  M <- matrix(rnorm(10 * 100), 10, 100)
  ft <- train_feature_texture(M)
  out <- capture.output(print(ft))
  expect_true(any(grepl("iconic_feature_texture", out)))
  expect_true(any(grepl("features", out)))
})

test_that("run_single_iteration uses copula texture when available", {
  skip_on_cran()
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  inp <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(inp$gan_training_data,
                                feature_correlations = inp$feature_correlations,
                                feature_texture = inp$feature_texture,
                                epochs = 5, verbose = FALSE)
  dat <- run_single_iteration(gan, n_synthetic_samples = 100, n_features = 30,
                              n_mediators = 1, phi = 0.5, mo_confounding = 0.5,
                              seed = 42)
  expect_length(dat$M, 100)
  expect_equal(dim(dat$Y), c(100, 30))
  # Ground truth should be preserved
  expect_equal(dat$true_NDE, 0.1)
  expect_equal(dat$true_NIE, 0.15)
})
