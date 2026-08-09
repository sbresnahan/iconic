# v0.9.8 revision tests (plan §4)

# ── 1. Rename ──
test_that("iconic_data(X = ...) works and Z = ... errors", {
  dat <- generate_toy_data(n = 100, n_features = 3, seed = 1)
  idat <- iconic_data(X = dat$X, Y = dat$Y, covariates = dat$synthetic_data)
  expect_s3_class(idat, "iconic_data")
  expect_error(iconic_data(Z = dat$X, Y = dat$Y), "renamed to `X`")
})

test_that("generate_toy_data(beta_X = ...) sets exposure effect; beta_Z errors", {
  dat <- generate_toy_data(n = 100, n_features = 3, beta_X = 0.25, seed = 2)
  expect_equal(dat$dgp_params$beta_X, 0.25)
  expect_error(generate_toy_data(n = 100, n_features = 3, beta_Z = 0.25),
               "renamed to `beta_X`")
})

# ── 2. Completeness wiring ──
test_that("iconic_diagnose flags under-identification when valid NCs < k", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 400, n_features = 10,
                              n_confounders = 2, mo_confounding = 0.8, phi = 0.8,
                              rho_G1 = 0, rho_G2 = 0, seed = 3)
  # supply only 1 NC (first column of the n x q W) so valid count < k = 2
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G[, 1], Gm = dat$Gm, W = dat$W[, 1, drop = FALSE],
                      covariates = dat$synthetic_data)
  diag <- iconic_diagnose(idat, k = 2)
  expect_equal(diag$k, 2)
  expect_false(diag$k_inferred)
  expect_equal(diag$completeness$completeness, "under-identified")
})

test_that("k falls back to 1 when <5 outcome features", {
  dat <- generate_toy_data(n = 200, n_features = 3, mo_confounding = 0.8,
                           phi = 0.8, seed = 4)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G, Gm = dat$Gm, W = dat$W,
                      covariates = dat$synthetic_data)
  diag <- iconic_diagnose(idat)
  expect_equal(diag$k, 1)
})

test_that("capture screen runs and returns a verdict", {
  dat <- generate_toy_data(n = 300, n_features = 8, mo_confounding = 0.8,
                           phi = 0.8, rho_G1 = 0, rho_G2 = 0, seed = 5)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G, Gm = dat$Gm, W = dat$W,
                      covariates = dat$synthetic_data)
  diag <- iconic_diagnose(idat)
  expect_true(!is.null(diag$completeness$capture))
  expect_true(diag$completeness$capture$capture_verdict %in%
                c("strong", "weak", "negligible"))
})

# ── 3. Omega sweep ──
test_that("iconic_sensitivity omega sweep returns omega columns", {
  dat <- generate_toy_data(n = 200, n_features = 5, mo_confounding = 0.8,
                           phi = 0.8, seed = 6)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G, Gm = dat$Gm, W = dat$W,
                      covariates = dat$synthetic_data)
  fake_gan <- list(feature_correlations = NULL, feature_texture = NULL)
  class(fake_gan) <- "iconic_gan"
  sens <- iconic_sensitivity(idat, trained_gan = fake_gan,
    rho_G1_grid = c(0, 0.3), rho_G2_grid = c(0, 0.3),
    omega_1 = c(0.3, 1.0), omega_2 = 0.7,
    n_iter = 3, n_cores = 1)
  expect_true("omega_1" %in% names(sens$surface))
  expect_true("omega_2" %in% names(sens$surface))
  expect_true(sens$omega_swept)
  # 2 rho1 x 2 rho2 x 2 omega1 x 1 omega2 = 8 cells x 8 methods = 64 rows
  expect_equal(nrow(sens$surface), 8 * 8)
})

test_that("iconic_prospect Phase 1 surface carries omega", {
  dat <- generate_toy_data(n = 200, n_features = 5, mo_confounding = 0.8,
                           phi = 0.8, seed = 7)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G, Gm = dat$Gm, W = dat$W,
                      covariates = dat$synthetic_data)
  fake_gan <- list(feature_correlations = NULL, feature_texture = NULL)
  class(fake_gan) <- "iconic_gan"
  prosp <- iconic_prospect(idat, trained_gan = fake_gan,
    gamma_G_grid = c(0.4, 0.8),
    omega_1 = c(0.3, 1.0), omega_2 = 0.7,
    n_iter = 3)
  expect_true("omega_1" %in% names(prosp$strength_surface))
  expect_true(prosp$omega_swept)
})

# ── 4. Tier removal ──
test_that("iconic_recommend has no tier field and per-estimand robustness", {
  dat <- generate_toy_data(n = 200, n_features = 5, mo_confounding = 0.8,
                           phi = 0.8, seed = 8)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G, Gm = dat$Gm, W = dat$W,
                      covariates = dat$synthetic_data)
  diag <- iconic_diagnose(idat)
  est <- iconic_estimate(idat, diagnosis = diag)
  fake_gan <- list(feature_correlations = NULL, feature_texture = NULL)
  class(fake_gan) <- "iconic_gan"
  sens <- iconic_sensitivity(idat, trained_gan = fake_gan,
    rho_G1_grid = c(0, 0.3), rho_G2_grid = c(0, 0.3),
    n_iter = 3, n_cores = 1)
  rec <- iconic_recommend(idat, diagnosis = diag, estimate = est, sensitivity = sens)
  expect_null(rec$ranking$tier)
  expect_null(rec$recommended_tier)
  expect_true("robustness_NDE" %in% names(rec$ranking))
  expect_true("robustness_NIE" %in% names(rec$ranking))
  expect_true(!is.na(rec$recommended))
})

# ── 5. TE output ──
test_that("mediation estimate returns TE, TE_se, TE_p", {
  dat <- generate_toy_data(n = 200, n_features = 3, mo_confounding = 0.8,
                           phi = 0.8, seed = 9)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G, Gm = dat$Gm, W = dat$W,
                      covariates = dat$synthetic_data)
  est <- iconic_estimate(idat)
  expect_true(all(c("TE", "TE_se", "TE_p") %in% names(est)))
  # TE = NDE + NIE within tolerance
  ok <- !is.na(est$TE) & !is.na(est$NDE) & !is.na(est$NIE)
  expect_true(all(abs(est$TE[ok] - (est$NDE[ok] + est$NIE[ok])) < 1e-8))
})

# ── 6. Scaling ──
test_that("iconic_data scales all inputs and records parameters", {
  set.seed(10)
  X_raw <- rnorm(100, mean = 5, sd = 3)
  Y_raw <- matrix(rnorm(100 * 3, mean = 2, sd = 4), 3, 100)
  idat <- iconic_data(X = X_raw, Y = Y_raw)
  expect_equal(mean(idat$X), 0, tolerance = 1e-8)
  expect_equal(sd(idat$X), 1, tolerance = 1e-8)
  expect_equal(mean(idat$Y[1, ]), 0, tolerance = 1e-8)
  expect_true(!is.null(idat$scaling$X))
  expect_true(!is.null(idat$scaling$Y))
  # back-transformation
  bt <- idat$X * idat$scaling$X$scale + idat$scaling$X$center
  expect_true(all(abs(bt - X_raw) < 1e-8))
})

# ── 7. Print / summary ──
test_that("assigned iconic_diagnose result has explicit summary()", {
  dat <- generate_toy_data(n = 200, n_features = 5, mo_confounding = 0.8,
                           phi = 0.8, seed = 11)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G, Gm = dat$Gm, W = dat$W,
                      covariates = dat$synthetic_data)
  expect_message(diag <- iconic_diagnose(idat), "summary")
  out <- capture.output(summary(diag))
  expect_true(any(grepl("iconic_diagnosis", out)))
})

# ── 8. separate_U removal ──
test_that("separate_U errors informatively", {
  expect_error(generate_toy_data(n = 100, n_features = 3, separate_U = TRUE),
               "separate_U")
  expect_error(run_single_iteration(NULL, n_synthetic_samples = 100,
                                    n_features = 3, separate_U = TRUE),
               "separate_U")
})

test_that("default shared loadings reproduce old separate_U = FALSE", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 300, n_features = 5,
                              n_confounders = 2, mo_confounding = 0.8, phi = 0.8,
                              rho_G1 = 0.1, rho_G2 = 0, seed = 12)
  # shared loadings: conf_XM and conf_MY are the same composite
  expect_equal(cor(dat$conf_XM, dat$conf_MY), 1.0, tolerance = 1e-8)
})

test_that("path-specific loadings reproduce old separate_U = TRUE", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 300, n_features = 5,
                              n_confounders = 2, mo_confounding = 0.8, phi = 0.8,
                              lambda_XM = c(1, 0), lambda_MY = c(0, 1),
                              rho_G1 = 0.1, rho_G2 = 0, seed = 13)
  # e1/e2 loadings: conf_XM and conf_MY are distinct confounder directions
  expect_true(abs(cor(dat$conf_XM, dat$conf_MY)) < 0.3)
})

test_that("overlapping-but-distinct loadings work (new capability)", {
  dat <- run_single_iteration(NULL, n_synthetic_samples = 300, n_features = 5,
                              n_confounders = 3, mo_confounding = 0.8, phi = 0.8,
                              lambda_XM = c(1, 1, 0), lambda_MY = c(0, 1, 1),
                              rho_G1 = 0.1, rho_G2 = 0, seed = 14)
  # partial overlap: correlation between 0 and 1
  r <- cor(dat$conf_XM, dat$conf_MY)
  expect_true(r > 0.1 && r < 0.9)
})

test_that("plot_nc_configuration_comparison is removed", {
  expect_false(exists("plot_nc_configuration_comparison"))
})

test_that("scenario_manifest no longer records separate_U", {
  dat <- generate_toy_data(n = 100, n_features = 3, seed = 15)
  m <- scenario_manifest(dat)
  expect_false("separate_U" %in% m$modifiable_parameters$parameter)
  expect_false("separate_U" %in% m$fixed_parameters$parameter)
})

test_that("IV estimators are eligible and run without W (instrument-only)", {
  # IV2SLS / IV2SLS2 are identified by the instrument(s) alone; W is an
  # optional proximal augmentation. With G + Gm and no W, both must be
  # eligible and must produce finite estimates.
  dat <- generate_toy_data(n = 300, n_features = 5, mo_confounding = 0.8,
                           phi = 0.8, rho_G1 = 0, rho_G2 = 0, seed = 3)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G[, 1], Gm = dat$Gm,
                      covariates = dat$synthetic_data)
  expect_false(idat$has_nc)
  diag <- iconic_diagnose(idat, n_cores = 1)
  elig <- diag$eligibility
  expect_true(elig$eligible[elig$estimator == "IV2SLS"])
  expect_true(elig$eligible[elig$estimator == "IV2SLS2"])
  # proximal bridge estimators still require W
  expect_false(elig$eligible[elig$estimator == "PGC"])
  expect_false(elig$eligible[elig$estimator == "PGC2"])
  expect_false(elig$eligible[elig$estimator == "PGC2Gm"])
  # and the estimators actually run (finite NDE)
  est <- iconic_estimate(idat, diagnosis = diag, n_cores = 1)
  s <- if (!is.null(est$summary)) est$summary else est
  expect_true(all(c("UNADJ", "IV2SLS", "IV2SLS2") %in% s$method))
  expect_true(all(is.finite(s$NDE[s$method %in% c("IV2SLS", "IV2SLS2")])))
})

test_that("IV2SLS with G only (no Gm, no W) is eligible; IV2SLS2 is not", {
  dat <- generate_toy_data(n = 300, n_features = 5, mo_confounding = 0.8,
                           phi = 0.8, rho_G1 = 0, rho_G2 = 0, seed = 7)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G[, 1],
                      covariates = dat$synthetic_data)
  diag <- iconic_diagnose(idat, n_cores = 1)
  elig <- diag$eligibility
  expect_true(elig$eligible[elig$estimator == "IV2SLS"])
  expect_false(elig$eligible[elig$estimator == "IV2SLS2"])
})

test_that("infer_confounding subsets a wide mediator panel to max_infer_tasks", {
  # 120 mediators, single vector outcome: the estimate auto-run should use
  # only 50 mediators, and the subset should be recorded.
  dat <- generate_toy_data(n = 200, n_features = 1, n_mediators = 120,
                           mo_confounding = 0.8, phi = 0.8, rho_G1 = 0,
                           rho_G2 = 0, seed = 5)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G[, 1], Gm = dat$Gm, W = dat$W,
                      covariates = dat$synthetic_data)
  diag <- iconic_diagnose(idat, n_cores = 1)
  conf <- infer_confounding(idat, diagnosis = diag, estimate = NULL,
                            n_cores = 1, max_infer_tasks = 50)
  is <- conf$inference_subset
  expect_true(is$subsetted)
  expect_equal(is$n_mediators_used, 50)
  expect_equal(is$n_mediators_total, 120)
  # gap-based conf_strength should be inferred (strong instrument) and the
  # method string should disclose the random subset
  expect_true(conf$conf_strength$available)
  expect_match(conf$conf_strength$method, "random subset")
})

test_that("infer_confounding subsets both mediators and features when wide", {
  dat <- generate_toy_data(n = 150, n_features = 80, n_mediators = 90,
                           mo_confounding = 0.8, phi = 0.8, rho_G1 = 0,
                           rho_G2 = 0, seed = 6)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G[, 1], Gm = dat$Gm, W = dat$W,
                      covariates = dat$synthetic_data)
  sub <- iconic:::.subset_data_for_inference(idat, 50)
  is <- attr(sub, "inference_subset")
  expect_true(is$subsetted)
  expect_equal(is$n_mediators_used, 50)
  expect_equal(is$n_features_used, 50)
  expect_equal(nrow(sub$M), 50)
  expect_equal(nrow(sub$Y), 50)
})

test_that("infer_confounding uses the full grid when the panel is small", {
  dat <- generate_toy_data(n = 150, n_features = 5, n_mediators = 10,
                           mo_confounding = 0.8, phi = 0.8, rho_G1 = 0,
                           rho_G2 = 0, seed = 7)
  idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G[, 1], Gm = dat$Gm, W = dat$W,
                      covariates = dat$synthetic_data)
  sub <- iconic:::.subset_data_for_inference(idat, 50)
  is <- attr(sub, "inference_subset")
  expect_false(is$subsetted)
  expect_equal(is$n_mediators_used, 10)
  expect_equal(is$n_features_used, 5)
})

test_that("robustness scoring is not poisoned by an all-NA estimator (survival COCA)", {
  # COCA is structurally incompatible with survival outcomes, so its
  # sensitivity-surface bias column is all NA -> max_bias = -Inf. The
  # combined-criterion normalization must average over finite values only,
  # so one such estimator does not turn every robustness score into NaN/-Inf.
  mk <- function(m, nde_mu, nie_mu) {
    data.frame(method = m, rho_G1 = 0, rho_G2 = 0, omega_1 = 0.7, omega_2 = 0.7,
               NDE_bias = rnorm(3, nde_mu, 0.02), NIE_bias = rnorm(3, nie_mu, 0.02),
               NDE_coverage = runif(3, 0.9, 0.96), NIE_coverage = runif(3, 0.9, 0.96),
               stringsAsFactors = FALSE)
  }
  set.seed(1)
  surface <- rbind(
    mk("UNADJ", 0.26, 0.42), mk("DIRECT", 0.33, 0.30),
    data.frame(method = "COCA", rho_G1 = 0, rho_G2 = 0, omega_1 = 0.7, omega_2 = 0.7,
               NDE_bias = NA_real_, NIE_bias = NA_real_,
               NDE_coverage = NA_real_, NIE_coverage = NA_real_,
               stringsAsFactors = FALSE),
    mk("IV2SLS", 0.19, 0.44), mk("PGC", 0.21, 0.24), mk("IV2SLS2", 0.20, 0.11),
    mk("PGC2", 0.31, 0.26), mk("PGC2Gm", 0.18, 0.11))
  rob <- iconic:::.extract_robustness(list(surface = surface), criterion = "combined")
  # every estimator except COCA must have a finite NDE robustness score
  finite_scores <- rob$score_NDE[rob$method != "COCA"]
  expect_true(all(is.finite(finite_scores)))
  # COCA (all-NA bias) must rank worst
  expect_equal(rob$score_NDE[rob$method == "COCA"], -Inf)
  # and the best (lowest max bias) estimator must score highest
  best <- rob$method[which.max(rob$score_NDE)]
  expect_equal(best, rob$method[which.min(ifelse(is.finite(rob$max_bias_NDE),
                                                 rob$max_bias_NDE, Inf))])
})
