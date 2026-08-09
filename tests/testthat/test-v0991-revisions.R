# Tests for the v0.9.9.1 revisions:
#  1. Lone-panel (W2-only / W1-only) data now derives the pooled W so
#     DIRECT / COCA / PGC become eligible and the NC validity screens run.
#     Applies to BOTH continuous and survival outcomes (the W/W1/W2
#     derivation in iconic_data() is outcome-type-agnostic).
#  2. recycle_lone_panel = TRUE opts in to using the lone panel as BOTH
#     path-specific bridges, making PGC2 / PGC2Gm eligible (with a warning).
#  3. iconic_sensitivity(confounding = "inferred") uses infer_confounding()'s
#     documented random subset (max_infer_tasks) and accepts a precomputed
#     iconic_confounding object.
#  4. User-supplied omega_1 / omega_2 sweep vectors take precedence over
#     inferred scalar omegas (the inferred values only fill in defaults).

# ── Helper: lone-W2 mediation data (GUSTO case-study wiring), continuous ──
.make_lone_w2 <- function(n = 200, n_features = 3, seed = 3) {
  dat <- iconic:::generate_toy_data(
    n = n, n_features = n_features, phi = 0.8, mo_confounding = 0.8,
    rho_G1 = 0.3, rho_G2 = 0.3, lambda_XM = c(1, 0), lambda_MY = c(0, 1),
    omega_1 = 0.7, omega_2 = 0.7, seed = seed)
  list(dat = dat,
       idata = iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                           G = dat$G[, 1], Gm = dat$Gm, W2 = t(dat$W2)))
}

# ── Helper: lone-W2 survival mediation data ──
.make_lone_w2_surv <- function(n = 300, n_features = 3, seed = 11,
                               recycle = FALSE) {
  dat <- iconic:::generate_toy_data(
    n = n, n_features = n_features, phi = 0.8, mo_confounding = 0.8,
    rho_G1 = 0.3, rho_G2 = 0.3, lambda_XM = c(1, 0), lambda_MY = c(0, 1),
    omega_1 = 0.7, omega_2 = 0.7, outcome_type = "survival", seed = seed)
  idata <- iconic_data(X = dat$X, outcome_type = "survival",
                       surv_time = dat$surv_time, surv_event = dat$surv_event,
                       M = dat$M, G = dat$G[, 1], Gm = dat$Gm,
                       W2 = t(dat$W2), recycle_lone_panel = recycle)
  list(dat = dat, idata = idata)
}

# A mock iconic_gan that makes sample_texture() fail, so the sweep falls
# back to default texture without requiring torch.
.mock_gan <- function() structure(list(model_type = "mock"), class = "iconic_gan")

# ═══════════════════════════════════════════════════════════════
# 1. Lone-panel W derivation -> single-panel estimator eligibility
# ═══════════════════════════════════════════════════════════════

test_that("lone W2 derives pooled W and enables DIRECT/COCA/PGC (continuous)", {
  id <- .make_lone_w2()$idata
  expect_true(id$has_nc)
  expect_false(is.null(id$W))
  expect_equal(id$W, id$W2)
  expect_false(id$has_path_nc)          # two-bridge estimators stay gated
  expect_false(isTRUE(id$recycled_lone_panel))

  d <- iconic_diagnose(id, min_f = 5)
  el <- d$eligibility
  elig <- el$estimator[el$eligible]
  # 6/8: single-panel + IV estimators eligible; PGC2/PGC2Gm still gated.
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC", "IV2SLS2")
                  %in% elig))
  expect_false("PGC2" %in% elig)
  expect_false("PGC2Gm" %in% elig)
})

test_that("lone W1 derives pooled W symmetrically (continuous)", {
  dat <- iconic:::generate_toy_data(
    n = 200, n_features = 3, phi = 0.8, mo_confounding = 0.8,
    rho_G1 = 0.3, rho_G2 = 0.3, lambda_XM = c(1, 0), lambda_MY = c(0, 1),
    omega_1 = 0.7, omega_2 = 0.7, seed = 3)
  id <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                    G = dat$G[, 1], Gm = dat$Gm, W1 = t(dat$W1))
  # Wiring: the lone W1 panel derives the pooled W (has_nc) and stays
  # single-panel (has_path_nc FALSE), exactly as for lone W2.
  expect_true(id$has_nc)
  expect_equal(id$W, id$W1)
  expect_false(id$has_path_nc)
  d <- iconic_diagnose(id, min_f = 5)
  elig <- d$eligibility$estimator[d$eligibility$eligible]
  # DIRECT and COCA are gated only on G + W / W, so they become eligible.
  # (PGC additionally requires the completeness screen, which a W1 panel
  # proxying the X->M confounder may fail because its features correlate
  # with X by construction -- a data-dependent gate, not a wiring issue.)
  expect_true(all(c("DIRECT", "COCA") %in% elig))
  expect_false("PGC2" %in% elig)
  expect_false("PGC2Gm" %in% elig)
})

# ═══════════════════════════════════════════════════════════════
# 1b. Survival outcomes get the same lone-panel fix
# ═══════════════════════════════════════════════════════════════

test_that("lone W2 derives pooled W and enables DIRECT/COCA/PGC (survival)", {
  id <- .make_lone_w2_surv()$idata
  expect_true(id$has_nc)
  expect_false(is.null(id$W))
  expect_equal(id$W, id$W2)
  expect_false(id$has_path_nc)

  d <- iconic_diagnose(id, min_f = 5)
  elig <- d$eligibility$estimator[d$eligibility$eligible]
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC", "IV2SLS2")
                  %in% elig))
  expect_false("PGC2" %in% elig)
  expect_false("PGC2Gm" %in% elig)
})

test_that("survival estimation runs newly-eligible single-panel estimators", {
  id <- .make_lone_w2_surv()$idata
  est <- iconic_estimate(id, min_f = 5, effect_scale = "loghr")
  methods <- sort(unique(est$method))
  expect_true(all(c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC", "IV2SLS2")
                  %in% methods))
  expect_false("PGC2" %in% methods)
})

test_that("survival recycle enables PGC2/PGC2Gm estimation", {
  expect_warning(
    id <- .make_lone_w2_surv(recycle = TRUE)$idata,
    "recycle_lone_panel")
  expect_true(id$has_path_nc)
  est <- iconic_estimate(id, min_f = 5, effect_scale = "loghr")
  methods <- sort(unique(est$method))
  expect_true(all(c("PGC2", "PGC2Gm") %in% methods))
})

# ═══════════════════════════════════════════════════════════════
# 2. Opt-in recycle -> PGC2 / PGC2Gm eligible (continuous)
# ═══════════════════════════════════════════════════════════════

test_that("recycle_lone_panel = TRUE makes PGC2/PGC2Gm eligible with warning", {
  dat <- .make_lone_w2()$dat
  expect_warning(
    id <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                      G = dat$G[, 1], Gm = dat$Gm, W2 = t(dat$W2),
                      recycle_lone_panel = TRUE),
    "recycle_lone_panel")
  expect_true(id$has_path_nc)
  expect_true(id$recycled_lone_panel)
  expect_equal(id$W1, id$W2)

  d <- iconic_diagnose(id, min_f = 5)
  el <- d$eligibility
  elig <- el$estimator[el$eligible]
  expect_true(all(c("PGC2", "PGC2Gm") %in% elig))
  expect_equal(sum(el$eligible), 8L)
  # Reason strings disclose the shared recycled panel.
  expect_match(el$reason[el$estimator == "PGC2"], "shared recycled panel")
  expect_match(el$reason[el$estimator == "PGC2Gm"], "shared recycled panel")
})

test_that("recycle_lone_panel = FALSE (default) keeps PGC2 gated", {
  id <- .make_lone_w2()$idata
  expect_false(id$has_path_nc)
  d <- iconic_diagnose(id, min_f = 5)
  expect_false("PGC2" %in% d$eligibility$estimator[d$eligibility$eligible])
})

# ═══════════════════════════════════════════════════════════════
# 3. No regression on the pooled-W and W1+W2 paths
# ═══════════════════════════════════════════════════════════════

test_that("pooled W still sets W1=W2=W (backward compatible)", {
  W <- matrix(rnorm(50 * 5), 5, 50)
  d <- iconic_data(X = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50), W = W)
  expect_true(d$has_nc)
  expect_true(d$has_path_nc)
  expect_false(isTRUE(d$recycled_lone_panel))
})

test_that("separate W1 and W2 still derive combined W", {
  W1 <- matrix(rnorm(50 * 5), 5, 50)
  W2 <- matrix(rnorm(50 * 5), 5, 50)
  d <- iconic_data(X = rnorm(50), Y = matrix(rnorm(50 * 3), 3, 50),
                   W1 = W1, W2 = W2)
  expect_true(d$has_path_nc)
  expect_true(d$has_nc)
  expect_false(isTRUE(d$recycled_lone_panel))
})

# ═══════════════════════════════════════════════════════════════
# 4. iconic_sensitivity: subset path + precomputed confounding + omega precedence
# ═══════════════════════════════════════════════════════════════

# Build a >50-mediator panel so the inference subset is exercised.
.make_big_med_data <- function(n = 120, n_mediators = 60, seed = 42) {
  dat <- iconic:::generate_toy_data(n = n, n_features = 5, phi = 0.8,
                                    mo_confounding = 0.8, seed = seed)
  M <- matrix(rep(dat$M, n_mediators), nrow = n_mediators, byrow = TRUE) +
    matrix(rnorm(n_mediators * n, 0, 0.05), n_mediators, n)
  Gm <- matrix(rep(dat$Gm, n_mediators), nrow = n_mediators, byrow = TRUE)
  iconic_data(X = dat$X, Y = t(dat$Y), M = M, W = t(dat$W),
              G = dat$G[, 1], Gm = Gm, covariates = dat$synthetic_data)
}

test_that("iconic_sensitivity(confounding='inferred') uses the random subset", {
  idata <- .make_big_med_data()
  diag <- iconic_diagnose(idata, min_f = 5)
  sens <- iconic_sensitivity(idata, diagnosis = diag, confounding = "inferred",
                             trained_gan = .mock_gan(),
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3)
  ic <- sens$inferred_confounding
  expect_s3_class(ic, "iconic_confounding")
  # Gap-based calibration ran on a random subset, not the full 60 mediators.
  expect_true(isTRUE(ic$inference_subset$subsetted))
  expect_equal(ic$inference_subset$n_mediators_used, 50L)
  expect_equal(ic$inference_subset$n_mediators_total, 60L)
})

test_that("iconic_sensitivity accepts a precomputed iconic_confounding object", {
  idata <- .make_big_med_data()
  diag <- iconic_diagnose(idata, min_f = 5)
  conf <- infer_confounding(idata, diagnosis = diag, n_cores = 1)
  sens <- iconic_sensitivity(idata, diagnosis = diag, confounding = conf,
                             trained_gan = .mock_gan(),
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3)
  expect_s3_class(sens$inferred_confounding, "iconic_confounding")
})

test_that("user-supplied omega sweep takes precedence over inferred omegas", {
  idata <- .make_big_med_data()
  diag <- iconic_diagnose(idata, min_f = 5)
  conf <- infer_confounding(idata, diagnosis = diag, n_cores = 1)
  sens <- iconic_sensitivity(idata, diagnosis = diag, confounding = conf,
                             trained_gan = .mock_gan(),
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3,
                             omega_1 = c(0.3, 0.7, 1.0),
                             omega_2 = c(0.3, 0.7, 1.0))
  # The 3-level diagonal sweep must be preserved, not collapsed to a scalar.
  expect_equal(length(unique(sens$surface$omega_1)), 3L)
})

test_that("default omega is filled by the inferred scalar", {
  idata <- .make_big_med_data()
  diag <- iconic_diagnose(idata, min_f = 5)
  conf <- infer_confounding(idata, diagnosis = diag, n_cores = 1)
  sens <- iconic_sensitivity(idata, diagnosis = diag, confounding = conf,
                             trained_gan = .mock_gan(),
                             rho_G1_grid = c(0), rho_G2_grid = c(0),
                             n_iter = 2, n_features = 3)
  # omega left at default -> inferred scalar -> single coverage level.
  expect_equal(length(unique(sens$surface$omega_1)), 1L)
})
