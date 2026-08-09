# Tests for iconic_sensitivity(): degradation surface calibrated to user data.
#
# The sensitivity surface sweeps a 2D grid of instrument-independence
# violations (rho_G1 x rho_G2) and reports how each estimator's NDE/NIE
# bias degrades.

# ── Helper: create a small iconic_data from run_single_iteration ──
.make_sens_data <- function(n = 100, n_features = 5, phi = 0.8, seed = 42) {
  dat <- run_single_iteration(
    n_synthetic_samples = n, n_features = n_features,
    mo_confounding = 0.8, phi = phi, seed = seed)
  iconic_data(X = dat$X, Y = t(dat$Y), M = dat$M, W = t(dat$W),
              G = dat$G[, 1], Gm = dat$Gm, covariates = dat$synthetic_data)
}

# ═══════════════════════════════════════════════════════════════
# Surface structure
# ═══════════════════════════════════════════════════════════════

test_that("iconic_sensitivity mediation surface includes IV2SLS2 and PGC2Gm", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0, 0.3),
                             rho_G2_grid = c(0, 0.3),
                             n_iter = 3, n_features = 3)
  methods_present <- unique(sens$surface$method)
  expect_true("IV2SLS2" %in% methods_present)
  expect_true("PGC2Gm" %in% methods_present)
})

test_that("iconic_sensitivity surface includes all 8 methods at origin", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0),
                             rho_G2_grid = c(0),
                             n_iter = 5, n_features = 3)
  methods_present <- unique(sens$surface$method)
  # All 8 estimators should be present when G, Gm, W are all supplied
  for (m in c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
              "IV2SLS2", "PGC2", "PGC2Gm")) {
    expect_true(m %in% methods_present, label = paste("method", m, "present"))
  }
})

test_that("iconic_sensitivity rejects non-mediation data", {
  idata <- .make_sens_data()
  bare <- iconic_data(X = idata$X, Y = idata$Y, G = idata$G, W = idata$W)
  expect_error(iconic_sensitivity(bare), "mediation data")
})

# ═══════════════════════════════════════════════════════════════
# Custom grids
# ═══════════════════════════════════════════════════════════════

test_that("iconic_sensitivity custom grid dimensions match", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  r1 <- c(0, 0.1, 0.2, 0.3, 0.5)
  r2 <- c(0, 0.2, 0.5)
  # Pin omega to a single value so the row count isolates the rho grid.
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = r1, rho_G2_grid = r2,
                             omega_1 = 0.7, omega_2 = 0.7,
                             n_iter = 2, n_features = 3)
  expect_equal(length(unique(sens$surface$rho_G1)), length(r1))
  expect_equal(length(unique(sens$surface$rho_G2)), length(r2))
  # Total cells = length(r1) * length(r2) * n_methods
  n_methods <- length(unique(sens$surface$method))
  expect_equal(nrow(sens$surface), length(r1) * length(r2) * n_methods)
})

test_that("iconic_sensitivity single-cell grid works", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0.3),
                             rho_G2_grid = c(0.3),
                             n_iter = 3, n_features = 3)
  expect_equal(length(unique(sens$surface$rho_G1)), 1)
  expect_equal(length(unique(sens$surface$rho_G2)), 1)
})

# ═══════════════════════════════════════════════════════════════
# Tipping points
# ═══════════════════════════════════════════════════════════════

test_that("iconic_sensitivity tipping_points has expected structure", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0, 0.3, 0.5),
                             rho_G2_grid = c(0, 0.3, 0.5),
                             n_iter = 5, n_features = 3)
  expect_true(is.data.frame(sens$tipping_points))
  expect_true("method" %in% names(sens$tipping_points))
  expect_true("max_NDE_bias" %in% names(sens$tipping_points))
  expect_true("max_NIE_bias" %in% names(sens$tipping_points))
  expect_true("tip_rho_G2_NDE" %in% names(sens$tipping_points))
  expect_true("tip_rho_G1_NDE" %in% names(sens$tipping_points))
  expect_true(nrow(sens$tipping_points) > 0)
})

test_that("iconic_sensitivity tipping points are within grid range", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  r1 <- c(0, 0.1, 0.2, 0.3, 0.5)
  r2 <- c(0, 0.1, 0.2, 0.3, 0.5)
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = r1, rho_G2_grid = r2,
                             n_iter = 5, n_features = 3)
  tp <- sens$tipping_points
  # Non-NA tipping points should be within the grid
  nde_tips <- tp$tip_rho_G2_NDE[!is.na(tp$tip_rho_G2_NDE)]
  if (length(nde_tips) > 0)
    expect_true(all(nde_tips %in% r2))
  g1_tips <- tp$tip_rho_G1_NDE[!is.na(tp$tip_rho_G1_NDE)]
  if (length(g1_tips) > 0)
    expect_true(all(g1_tips %in% r1))
})

# ═══════════════════════════════════════════════════════════════
# Surface content
# ═══════════════════════════════════════════════════════════════

test_that("iconic_sensitivity surface has NDE_bias and tipped columns", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0, 0.3),
                             rho_G2_grid = c(0, 0.3),
                             n_iter = 3, n_features = 3)
  expect_true("NDE_bias" %in% names(sens$surface))
  expect_true("NIE_bias" %in% names(sens$surface))
  expect_true("tipped" %in% names(sens$surface))
  expect_true("tipped_NDE" %in% names(sens$surface))
  expect_true("tipped_NIE" %in% names(sens$surface))
})

test_that("iconic_sensitivity PGC2Gm robust at origin", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0),
                             rho_G2_grid = c(0),
                             n_iter = 10, n_features = 5)
  pgc2gm <- sens$surface[sens$surface$method == "PGC2Gm" &
                         sens$surface$rho_G1 == 0 &
                         sens$surface$rho_G2 == 0, ]
  expect_lt(abs(pgc2gm$NDE_bias), 0.1)
})

test_that("iconic_sensitivity IV2SLS2 robust at origin", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0),
                             rho_G2_grid = c(0),
                             n_iter = 10, n_features = 5)
  iv2sls2 <- sens$surface[sens$surface$method == "IV2SLS2" &
                          sens$surface$rho_G1 == 0 &
                          sens$surface$rho_G2 == 0, ]
  expect_lt(abs(iv2sls2$NDE_bias), 0.1)
})

# ═══════════════════════════════════════════════════════════════
# Consistency with gan_mediation_sensitivity
# ═══════════════════════════════════════════════════════════════

test_that("iconic_sensitivity bias values in same range as gan_mediation_sensitivity", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data(n = 200, n_features = 5)
  # Use the same parameters for both
  common_args <- list(
    n_iter = 5, n_samples = 200, n_features = 5,
    mo_confounding = 0.8, phi = 0.8,
    omega_1 = 0.7, omega_2 = 0.7
  )
  sens <- do.call(iconic_sensitivity, c(list(data = idata,
                                             rho_G1_grid = c(0, 0.3),
                                             rho_G2_grid = c(0, 0.3)),
                                        common_args))
  # gan_mediation_sensitivity at rho_G1=0.3, rho_G2=0.3
  gms <- gan_mediation_sensitivity(
    trained_gan = NULL,
    conf_grid = 0.8, coverage_grid = 0.7, k_grid = 1,
    mo_confounding = 0.8, phi = 0.8,
    rho_G1 = 0.3, rho_G2 = 0.3,
    omega_1 = 0.7, omega_2 = 0.7,
    n_iter = 5, n_samples = 200, n_features = 5)

  # Compare NDE bias for IV2SLS2 and PGC2Gm at the same violation level
  sens_iv <- sens$surface[sens$surface$method == "IV2SLS2" &
                          sens$surface$rho_G1 == 0.3 &
                          sens$surface$rho_G2 == 0.3, ]$NDE_bias
  gms_iv <- gms$summary[gms$summary$method == "IV2SLS2", ]$NDE_bias

  sens_pg <- sens$surface[sens$surface$method == "PGC2Gm" &
                          sens$surface$rho_G1 == 0.3 &
                          sens$surface$rho_G2 == 0.3, ]$NDE_bias
  gms_pg <- gms$summary[gms$summary$method == "PGC2Gm", ]$NDE_bias

  # Both should be in a similar range (within 0.1 of each other)
  # They won't be identical due to different seeds, but should be close
  expect_lt(abs(sens_iv - gms_iv), 0.15)
  expect_lt(abs(sens_pg - gms_pg), 0.15)
})

# ═══════════════════════════════════════════════════════════════
# Print method
# ═══════════════════════════════════════════════════════════════

test_that("print.iconic_sensitivity produces output", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3)
  out <- capture.output(print(sens))
  expect_true(any(grepl("iconic_sensitivity", out)))
  expect_true(any(grepl("Grid", out)))
  expect_true(any(grepl("bias", out, ignore.case = TRUE)))
})

# ═══════════════════════════════════════════════════════════════
# Calibration
# ═══════════════════════════════════════════════════════════════

test_that("iconic_sensitivity calibrates n_samples to data", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data(n = 150)
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3)
  expect_equal(sens$n_samples, 150)
})

test_that("iconic_sensitivity respects explicit n_samples override", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data(n = 100)
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3,
                             n_samples = 300)
  expect_equal(sens$n_samples, 300)
})

test_that("iconic_sensitivity respects bias_threshold override", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0, 0.3),
                             rho_G2_grid = c(0, 0.3),
                             n_iter = 3, n_features = 3,
                             bias_threshold = 0.05)
  expect_equal(sens$bias_threshold, 0.05)
  # Tipped flags should use the 0.05 threshold
  expect_true(all(abs(sens$surface$NDE_bias[sens$surface$tipped_NDE]) > 0.05))
})

test_that("iconic_sensitivity returns S3 class", {
  skip_if_not_installed("torch")
  skip_if_not(check_torch_setup())
  idata <- .make_sens_data()
  sens <- iconic_sensitivity(idata,
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3)
  expect_s3_class(sens, "iconic_sensitivity")
})
