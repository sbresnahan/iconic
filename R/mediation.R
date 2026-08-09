# ============================================================
# Mediation estimators: natural direct (NDE) and natural
# indirect (NIE) effects under unmeasured confounding.
#
# Each estimator returns list(NDE, NDE_se, NDE_p, NIE, NIE_se, NIE_p).
# NDE = natural direct effect (beta_X, X -> Y not through M)
# NIE = natural indirect effect (alpha_M * beta_M, X -> M -> Y)
#
# The data-generating process may include mediator-outcome (M-O)
# confounding via a shared unmeasured confounder U1 that affects
# both M and Y. Under M-O confounding, natural effects are not
# point-identified by a single genetic instrument (Rudolph et al.,
# 2024); the estimators below are approximations whose bias is the
# central quantity the simulation benchmarks.
#
# fit_pgc_mediation() now uses a MATRIX bridge (regresses
# X_resid on the full W matrix). The original scalar-bridge version
# is retained as fit_pgc_scalar_mediation().
#
# fit_iv2sls_mediation2() implements a 2-stage MR mediation
# estimator that uses TWO instruments -- G for X and Gm for M --
# making NDE and NIE point-identified under M-O confounding when both
# instruments are valid and strong. This resolves the identification
# failure that limits the single-instrument estimators above.
#
# References:
# Baron & Kenny (1986); Robins & Greenland (1992);
# VanderWeele (2015) -- mediation foundations
# Tchetgen Tchetgen (2014) -- COCA
# Miao, Geng & Tchetgen Tchetgen (2018) -- proximal inference
# Rudolph et al. (2024) -- natural effects with a single IV
# Loh et al. (2024) -- M-O confounding distorts NDE and NIE
# ============================================================


# ── Internal helpers ──

#' Delta-method SE for a product alpha * beta (internal)
#'
#' @param alpha Point estimate of alpha.
#' @param alpha_se SE of alpha.
#' @param beta Point estimate of beta.
#' @param beta_se SE of beta.
#' @param cov_ab Estimated covariance of alpha and beta. Default 0.
#' @return Scalar SE of the product alpha * beta.
#' @keywords internal
delta_se_product <- function(alpha, alpha_se, beta, beta_se, cov_ab = 0) {
  grad <- c(beta, alpha)
  V <- matrix(c(alpha_se^2, cov_ab, cov_ab, beta_se^2), 2, 2)
  sqrt(as.numeric(t(grad) %*% V %*% grad))
}

#' Bootstrap SE for a mediation estimator
#'
#' Resamples the data `n_boot` times and returns the SD of the bootstrap
#' NDE/NIE distribution. This is an opt-in alternative to the delta-method
#' SE for users who want a non-parametric SE (
#' bootstrapped reasonably"). Slower than the delta method.
#'
#' The estimator is supplied as a closure `estimator_fn(idx)` that captures
#' all model-specific data (outcome, exposure, mediator, instruments,
#' negative controls, covariates) and subsets each by `idx` internally.
#' This guarantees every resampled draw uses a synchronised bootstrap
#' sample across all variables — critical for the instrumented estimators
#' (IV2SLS, PGC, PGC2Gm) whose instruments and NC panels must be resampled
#' in lockstep with y, X, and M.
#'
#' @param estimator_fn A closure `function(idx)` returning a fit list with
#' `NDE` and `NIE` (and optionally `NDE_se`, `NIE_se`).
#' @param n Sample size (length of the resampling index).
#' @param n_boot Number of bootstrap resamples. Default 500.
#' @return A list with NDE_boot_se, NIE_boot_se, NDE_boot_dist, NIE_boot_dist.
#' @keywords internal
bootstrap_mediation_se <- function(estimator_fn, n, n_boot = 500) {
  nde_dist <- numeric(n_boot)
  nie_dist <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, replace = TRUE)
    fit <- tryCatch(estimator_fn(idx), error = function(e) NULL)
    if (!is.null(fit) && !is.na(fit$NDE)) {
      nde_dist[b] <- fit$NDE
      nie_dist[b] <- fit$NIE
    } else {
      nde_dist[b] <- NA_real_
      nie_dist[b] <- NA_real_
    }
  }
  list(
    NDE_boot_se = sd(nde_dist, na.rm = TRUE),
    NIE_boot_se = sd(nie_dist, na.rm = TRUE),
    NDE_boot_dist = nde_dist,
    NIE_boot_dist = nie_dist
  )
}

# Reuse .covar_str and .bind_covars from estimators.R


# ── Composite null hypothesis test ──
#
# The NIE = alpha_M * beta_M is a product of two coefficients. The null
# H0: alpha_M * beta_M = 0 is a COMPOSITE null with three cases:
# H0(1): alpha_M = 0 AND beta_M = 0
# H0(2): alpha_M != 0, beta_M = 0
# H0(3): alpha_M = 0, beta_M != 0
# The Wald/Sobel test assumes the product is approximately normal, which
# only holds under H0(2)/H0(3). Under H0(1) (sparse signals) the product
# follows a normal product distribution (density f(z) = K0(z)/pi), making
# the Sobel test conservative.
#
# Huang (2019, Annals of Applied Statistics) proposed a closed-form
# composite p-value that accounts for all three null cases without
# estimating their proportions:
#
# p_comp = F(ab / sqrt(Var(a))) + F(ab / sqrt(Var(b))) - F(ab)
#
# where a = alpha_hat / SE(alpha_hat), b = beta_hat / SE(beta_hat) are
# the standardized z-statistics, Var(a)/Var(b) are their variances across
# the collection of tests, and F(z) is the CDF of the normal product
# distribution:
#
# F(z) = 2 * integral_{|z|}^{Inf} K0(x) / pi dx
#
# K0 is the modified Bessel function of the second kind (besselK(x, nu=0)).
# The formula is derived via a Taylor-series approximation (Theorem 3.3);
# the error term delta_N -> 0 rapidly when signals are sparse and Var(a),
# Var(b) < 1.5 (approximately n < 2000).
#
# References:
# Huang, Y.-T. (2019). Genome-wide analyses of sparse mediation effects
# under composite null hypotheses. Annals of Applied Statistics, 13(1),
# 60-84.
# Du, J. et al. (2023). Methods for large-scale single mediator
# hypothesis testing. Genetic Epidemiology, 47(2), 167-184.

#' Normal product distribution CDF
#'
#' Computes F(z) = 2 * \eqn{integral_{|z|}^{Inf}} K0(x)/pi dx, the two-sided
#' tail probability of the standard normal product distribution.
#' Under H0(1), if Z1, Z2 ~ N(0,1) independent, then P(|Z1*Z2| >= |z|)
#' = F(z). K0 is the modified Bessel function (besselK(x, nu=0)).
#'
#' @param z Numeric scalar or vector.
#' @return Numeric scalar or vector of CDF values in \eqn{[0, 1]}.
#' @keywords internal
.np_cdf <- function(z) {
  z <- abs(z)
  vapply(z, function(zz) {
    if (zz < 1e-10) return(1)
    # For large z, K0(x) underflows to 0 well before the upper limit,
    # causing integrate() to hit non-finite values. The tail probability
    # is negligibly small (F(10) ~ 1e-5, F(20) ~ 1e-9), so return 0.
    if (zz > 50) return(0)
    2 * integrate(function(x) besselK(x, nu = 0) / pi,
                  lower = zz, upper = Inf,
                  rel.tol = 1e-10, subdivisions = 500L)$value
  }, numeric(1))
}

#' Composite null hypothesis test p-value
#'
#' Computes the Huang (2019) JT-comp composite p-value for testing
#' H0: alpha * beta = 0 against H1: alpha * beta != 0.
#'
#' The composite null decomposes into three cases (alpha=beta=0,
#' alpha!=0/beta=0, alpha=0/beta!=0). The Wald/Sobel test is
#' conservative under the first case because the product of two
#' independent normals follows a normal product distribution, not a
#' normal. This function computes a closed-form p-value that accounts
#' for all three cases without estimating their proportions.
#'
#' @param a Standardized z-statistic for alpha (alpha_hat / SE).
#' @param b Standardized z-statistic for beta (beta_hat / SE).
#' @param var_a Variance of the z-statistic for alpha across the
#' collection of tests. Under the point null this is ~1; under
#' H0(2)/H0(3) it is >1. Default 1 (conservative).
#' @param var_b Variance of the z-statistic for beta across the
#' collection of tests. Default 1.
#' @return Numeric scalar in \eqn{[0, 1]}.
#' @export
#'
#' @references Huang, Y.-T. (2019). Genome-wide analyses of sparse
#' mediation effects under composite null hypotheses. \emph{Annals of
#' Applied Statistics}, 13(1), 60-84.
#'
#' @examples
#' composite_p_value(0, 0) # null: p = 1
#' composite_p_value(2, 2) # strong signal: small p
#' composite_p_value(2, 0) # H0(2): p = 1 (no mediation)
composite_p_value <- function(a, b, var_a = 1, var_b = 1) {
  ab <- a * b
  p <- .np_cdf(ab / sqrt(var_a)) +
       .np_cdf(ab / sqrt(var_b)) -
       .np_cdf(ab)
  pmin(pmax(p, 0), 1)
}

#' Apply composite null p-values to a mediation results data frame (internal)
#'
#' Two-pass post-processing for \code{se_method = "composite"}.
#'
#' The JT-comp test (Huang 2019) requires Var(a) and Var(b), the variances
#' of the standardized z-statistics \code{a = alpha_M / alpha_se} and
#' \code{b = beta_M / beta_M_se} across the collection of tests sharing
#' the same estimator. This function:
#'
#' 1. Groups rows by \code{method} (and \code{mediator} when present).
#' 2. For each group, computes the z-statistics from the delta-method
#' estimates already stored in the data frame.
#' 3. Estimates Var(a) / Var(b) as the sample variance of the z-statistics.
#' When fewer than 5 tests are available in a group, falls back to
#' Var = 1 (the point-null value, conservative).
#' 4. Replaces \code{NIE_p} with the JT-comp composite p-value.
#'
#' NDE_p is left unchanged (the NDE is a single coefficient, not a product,
#' so the composite null does not apply). NIE_se is left as the delta-method
#' SE (used for CI construction); only the p-value is replaced.
#'
#' @param res Data frame from \code{.estimate_mediation_feature()} /
#' \code{.estimate_mediation_driver()} / \code{run_mediation_methods()}.
#' Must contain columns: \code{method}, \code{NIE}, \code{NIE_p},
#' \code{alpha_M}, \code{alpha_se}, \code{beta_M}, \code{beta_M_se}.
#' @return The same data frame with \code{NIE_p} replaced by composite
#' p-values. Also adds columns \code{var_a}, \code{var_b} for
#' transparency.
#' @keywords internal
.apply_composite_pvalues <- function(res) {
  if (is.null(res) || nrow(res) == 0) return(res)

  # Group by method only. Var(a) and Var(b) are estimated across all
  # tests sharing the same estimator — i.e., across mediators (and
  # outcome features, when present) within each method. This is the
  # correct grouping for ICONIC's design: in the real case studies Y is
  # a scalar outcome and M is a panel of mediators, so the variation in
  # a = alpha_M / SE(alpha_M) comes from the different stage-1
  # regressions M_m ~ X for each mediator m. Grouping by (method,
  # mediator) instead would give one test per group when n_features = 1,
  # falling back to Var = 1 and making the composite test a no-op.
  groups <- factor(res$method, levels = unique(res$method))

  res$var_a <- NA_real_
  res$var_b <- NA_real_

  for (g in levels(groups)) {
    idx <- which(groups == g)
    if (!length(idx)) next

    a <- res$alpha_M[idx] / res$alpha_se[idx]
    b <- res$beta_M[idx] / res$beta_M_se[idx]

    # Keep only finite z-statistics
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 5) {
      # Too few tests to estimate variance; use point-null Var = 1
      var_a <- 1
      var_b <- 1
    } else {
      var_a <- stats::var(a[ok])
      var_b <- stats::var(b[ok])
    }

    # Variance clamping:
    #
    # 1. Lower bound (>= 1): Huang (2019) assumes independent tests, so
    # under the point null Var(a) = Var(b) = 1. In ICONIC, the tests
    # share the same samples and are positively correlated, which
    # deflates the sample variance below 1. Using Var < 1 makes the
    # test anti-conservative (inflated type I error). Clamping to
    # >= 1 restores the point-null calibration.
    #
    # 2. Upper bound (<= 1.5): Huang (2019) recommends the JT-comp
    # approximation only when Var(a), Var(b) < 1.5 (approximately
    # n < 2000 with sparse signals). Beyond this the Taylor-series
    # error term grows.
    if (!is.finite(var_a) || var_a < 1) var_a <- 1
    if (!is.finite(var_b) || var_b < 1) var_b <- 1
    if (var_a > 1.5) var_a <- 1.5
    if (var_b > 1.5) var_b <- 1.5

    res$var_a[idx] <- var_a
    res$var_b[idx] <- var_b

    # Compute composite p-value for each row in the group
    for (i in idx) {
      if (!is.na(res$alpha_M[i]) && !is.na(res$alpha_se[i]) &&
          res$alpha_se[i] != 0 && !is.na(res$beta_M[i]) &&
          !is.na(res$beta_M_se[i]) && res$beta_M_se[i] != 0) {
        a_i <- res$alpha_M[i] / res$alpha_se[i]
        b_i <- res$beta_M[i] / res$beta_M_se[i]
        res$NIE_p[i] <- composite_p_value(a_i, b_i, var_a, var_b)
      }
    }
  }

  res
}


# ── 1. UNADJ mediation ──

#' UNADJ mediation estimator: naive Baron-Kenny style
#'
#' Stage 1: \code{M ~ X} (estimate alpha_M).
#' Stage 2: \code{Y ~ X + M} (estimate NDE = beta_X, beta_M).
#' NIE = alpha_M * beta_M.
#'
#' Does not adjust for unmeasured confounding. Provided as a bias
#' reference floor, analogous to UNADJ in the total-effect setting.
#'
#' @param y Numeric outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param covars Optional data frame of additional covariates (n rows).
#'
#' @return Named list: \code{NDE}, \code{NDE_se}, \code{NDE_p},
#' \code{NIE}, \code{NIE_se}, \code{NIE_p}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, mo_confounding = 0.8, seed = 1)
#' fit_unadj_mediation(dat$Y[, 1], dat$X, dat$M)
#' }
fit_unadj_mediation <- function(y, X, M, covars = NULL) {
  NA_res <- list(NDE = NA_real_, NDE_se = NA_real_, NDE_p = NA_real_,
                 NIE = NA_real_, NIE_se = NA_real_, NIE_p = NA_real_,
                 alpha_M = NA_real_, alpha_se = NA_real_, beta_M = NA_real_, beta_M_se = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  # Stage 1: M ~ X
  d1 <- .bind_covars(data.frame(M = M, X = X), covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ X", cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"X" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["X"])
  alpha_se <- as.numeric(s1["X", 2])

  # Stage 2: Y ~ X + M
  d2 <- .bind_covars(data.frame(y = y, X = X, M = M), covars)
  fit2 <- tryCatch(lm(as.formula(paste0("y ~ X + M", cs)), data = d2),
                   error = function(e) NULL)
  if (is.null(fit2)) return(NA_res)
  s2 <- summary(fit2)$coefficients
  if (!"X" %in% rownames(s2) || !"M" %in% rownames(s2)) return(NA_res)
  beta_X <- as.numeric(coef(fit2)["X"])
  beta_X_se <- as.numeric(s2["X", 2])
  beta_M <- as.numeric(coef(fit2)["M"])
  beta_M_se <- as.numeric(s2["M", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_X, NDE_se = beta_X_se, NDE_p = as.numeric(s2["X", 4]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 2. DIRECT mediation ──

#' DIRECT mediation estimator: OLS with instrument and NC as covariates
#'
#' Adjusts for the genetic instrument G and negative-control W in both
#' the mediator and outcome regressions. Like \code{\link{fit_direct}},
#' this is a naive adjustment that does not correct for unmeasured
#' confounding via a ratio or IV approach.
#'
#' @param y Numeric outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param g Numeric instrument vector (length n).
#' @param w Numeric negative-control vector (length n).
#' @param covars Optional data frame of additional covariates (n rows).
#'
#' @return Named list: \code{NDE}, \code{NDE_se}, \code{NDE_p},
#' \code{NIE}, \code{NIE_se}, \code{NIE_p}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, mo_confounding = 0.8, seed = 1)
#' fit_direct_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$W[, 1])
#' }
fit_direct_mediation <- function(y, X, M, g, w, covars = NULL) {
  NA_res <- list(NDE = NA_real_, NDE_se = NA_real_, NDE_p = NA_real_,
                 NIE = NA_real_, NIE_se = NA_real_, NIE_p = NA_real_,
                 alpha_M = NA_real_, alpha_se = NA_real_, beta_M = NA_real_, beta_M_se = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  # Stage 1: M ~ X + G + W (full W panel)
  d1 <- data.frame(M = M, X = X, g = g)
  d1 <- cbind(d1, we$df)
  d1 <- .bind_covars(d1, covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ X + g + ", we$frag, cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"X" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["X"])
  alpha_se <- as.numeric(s1["X", 2])

  # Stage 2: Y ~ X + M + G + W (full W panel)
  d2 <- data.frame(y = y, X = X, M = M, g = g)
  d2 <- cbind(d2, we$df)
  d2 <- .bind_covars(d2, covars)
  fit2 <- tryCatch(lm(as.formula(paste0("y ~ X + M + g + ", we$frag, cs)), data = d2),
                   error = function(e) NULL)
  if (is.null(fit2)) return(NA_res)
  s2 <- summary(fit2)$coefficients
  if (!"X" %in% rownames(s2) || !"M" %in% rownames(s2)) return(NA_res)
  beta_X <- as.numeric(coef(fit2)["X"])
  beta_X_se <- as.numeric(s2["X", 2])
  beta_M <- as.numeric(coef(fit2)["M"])
  beta_M_se <- as.numeric(s2["M", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_X, NDE_se = beta_X_se, NDE_p = as.numeric(s2["X", 4]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 3. COCA mediation ──

#' COCA mediation estimator: negative-control calibration of both stages
#'
#' Uses the negative-control outcome W to calibrate both the mediator
#' and outcome regressions, extending \code{\link{fit_coca}} to the
#' mediation setting.
#'
#' Stage 1: \code{W ~ M + X} -> calibrated alpha_M = -beta_X / beta_M.
#' Stage 2: \code{W ~ Y + X + M} -> calibrated NDE = -beta_X / beta_Y,
#' calibrated beta_M = -beta_M / beta_Y.
#' NIE = alpha_M * beta_M (both calibrated).
#'
#' @param y Numeric primary outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param w Numeric negative-control outcome vector (length n).
#' Recommended: pass \code{rowMeans(W_matrix)} for stability.
#' @param covars Optional data frame of additional covariates (n rows).
#' @param ratio_cap Maximum absolute value of any ratio estimate before
#' flagging as unstable. Default 10.
#'
#' @return Named list: \code{NDE}, \code{NDE_se}, \code{NDE_p},
#' \code{NIE}, \code{NIE_se}, \code{NIE_p}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, mo_confounding = 0.8, seed = 1)
#' fit_coca_mediation(dat$Y[, 1], dat$X, dat$M, rowMeans(dat$W))
#' }
fit_coca_mediation <- function(y, X, M, w, covars = NULL, ratio_cap = 10) {
  NA_res <- list(NDE = NA_real_, NDE_se = NA_real_, NDE_p = NA_real_,
                 NIE = NA_real_, NIE_se = NA_real_, NIE_p = NA_real_,
                 alpha_M = NA_real_, alpha_se = NA_real_, beta_M = NA_real_, beta_M_se = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  # Stage 1: W ~ M + X -> alpha_M calibrated as -beta_X / beta_M
  d1 <- .bind_covars(data.frame(w = w, M = M, X = X), covars)
  fit1 <- tryCatch(lm(as.formula(paste0("w ~ M + X", cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  b1 <- coef(fit1)
  bZ1 <- b1["X"]; bM1 <- b1["M"]
  if (is.na(bM1) || abs(bM1) < 1e-8) return(NA_res)
  alpha <- -bZ1 / bM1
  if (abs(alpha) > ratio_cap) return(NA_res)
  V1 <- vcov(fit1)[c("X", "M"), c("X", "M")]
  grad1 <- c(-1 / bM1, bZ1 / bM1^2)
  alpha_se <- sqrt(as.numeric(t(grad1) %*% V1 %*% grad1))

  # Stage 2: W ~ Y + X + M -> NDE and beta_M calibrated
  d2 <- .bind_covars(data.frame(w = w, y = y, X = X, M = M), covars)
  fit2 <- tryCatch(lm(as.formula(paste0("w ~ y + X + M", cs)), data = d2),
                   error = function(e) NULL)
  if (is.null(fit2)) return(NA_res)
  b2 <- coef(fit2)
  bY2 <- b2["y"]; bZ2 <- b2["X"]; bM2 <- b2["M"]
  if (is.na(bY2) || abs(bY2) < 1e-8) return(NA_res)
  NDE <- -bZ2 / bY2
  beta_M_cal <- -bM2 / bY2
  if (abs(NDE) > ratio_cap) return(NA_res)

  V2 <- vcov(fit2)[c("X", "y"), c("X", "y")]
  grad2 <- c(-1 / bY2, bZ2 / bY2^2)
  NDE_se <- sqrt(as.numeric(t(grad2) %*% V2 %*% grad2))

  V2m <- vcov(fit2)[c("M", "y"), c("M", "y")]
  grad_bm <- c(-1 / bY2, bM2 / bY2^2)
  beta_M_se <- sqrt(as.numeric(t(grad_bm) %*% V2m %*% grad_bm))
  NIE <- alpha * beta_M_cal
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M_cal, beta_M_se)

  list(
    NDE = NDE, NDE_se = NDE_se, NDE_p = 2 * pnorm(-abs(NDE / NDE_se)),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M_cal, beta_M_se = beta_M_se
  )
}


# ── 4. IV2SLS mediation (single instrument) ──

#' IV2SLS mediation estimator: instrumented exposure in both stages
#'
#' Uses the genetic instrument G to purge U1 from X, then estimates
#' the mediator and outcome regressions with the cleaned exposure.
#'
#' Strategy:
#' \enumerate{
#' \item \code{X ~ G + W} -> X_hat (purge U1 from X).
#' \item \code{M ~ X_hat} -> alpha_M (clean effect of X on M).
#' \item \code{Y ~ X_hat + M + W} -> NDE = beta_X, beta_M (OLS).
#' }
#'
#' The IV cleans X of U1 confounding, but M remains endogenous via
#' U1 -> M. With a single instrument, natural effects are not fully
#' identified (Rudolph et al., 2024); NDE and NIE are approximations
#' whose bias from M-O confounding is the key finding the simulation
#' demonstrates.
#'
#' @param y Numeric outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param g Numeric instrument vector (length n).
#' @param w Numeric negative-control vector (length n).
#' @param covars Optional data frame of additional covariates (n rows).
#' @param min_f Minimum acceptable partial F-statistic for the excluded
#' instrument. Default 10.
#'
#' @return Named list: \code{NDE}, \code{NDE_se}, \code{NDE_p},
#' \code{NIE}, \code{NIE_se}, \code{NIE_p}.
#' @export
#'
#' @references
#' Rudolph, K. E., et al. (2024). Natural direct and indirect effects
#' with an instrumental variable. *Biometrics*.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 300, mo_confounding = 0.8, seed = 1)
#' fit_iv2sls_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$W[, 1])
#' }
fit_iv2sls_mediation <- function(y, X, M, g, w, covars = NULL, min_f = 10) {
  NA_res <- list(NDE = NA_real_, NDE_se = NA_real_, NDE_p = NA_real_,
                 NIE = NA_real_, NIE_se = NA_real_, NIE_p = NA_real_,
                 alpha_M = NA_real_, alpha_se = NA_real_, beta_M = NA_real_, beta_M_se = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  # Weak instrument check: X ~ G (+ W when present)
  d_fs <- data.frame(X = X, g = g)
  d_fs <- .bind_covars(d_fs, we$df)
  d_fs <- .bind_covars(d_fs, covars)
  fs <- tryCatch(lm(as.formula(paste0("X ~ g", .plus_frag(we$frag), cs)), data = d_fs),
                 error = function(e) NULL)
  if (is.null(fs)) return(NA_res)
  Fst <- .partial_F(fs, "g")
  if (is.na(Fst) || Fst < min_f) return(NA_res)

  # Stage 1: X ~ G (+ W) -> X_hat (purge U1 from X)
  X_hat <- fitted(fs)

  # Stage 2 (mediator): M ~ X_hat -> alpha_M
  d1 <- .bind_covars(data.frame(M = M, X_hat = X_hat), covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ X_hat", cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"X_hat" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["X_hat"])
  alpha_se <- as.numeric(s1["X_hat", 2])

  # Stage 3 (outcome): Y ~ X_hat + M (+ W when present; X cleaned, M observed)
  d2 <- data.frame(y = y, X_hat = X_hat, M = M)
  d2 <- .bind_covars(d2, we$df)
  d2 <- .bind_covars(d2, covars)
  fit2 <- tryCatch(lm(as.formula(paste0("y ~ X_hat + M", .plus_frag(we$frag), cs)), data = d2),
                   error = function(e) NULL)
  if (is.null(fit2)) return(NA_res)
  s2 <- summary(fit2)$coefficients
  if (!"X_hat" %in% rownames(s2) || !"M" %in% rownames(s2)) return(NA_res)
  beta_X <- as.numeric(coef(fit2)["X_hat"])
  beta_X_se <- as.numeric(s2["X_hat", 2])
  beta_M <- as.numeric(coef(fit2)["M"])
  beta_M_se <- as.numeric(s2["M", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_X, NDE_se = beta_X_se, NDE_p = as.numeric(s2["X_hat", 4]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 5. IV2SLS2 mediation (two instruments: 2-stage MR) ──

#' IV2SLS2 mediation estimator: 2-stage MR with instruments for both X and M
#'
#' Uses TWO genetic instruments: G for the exposure X and Gm for the
#' mediator M. This is the key extension: by instrumenting both
#' endogenous variables, NDE and NIE become \strong{point-identified} even
#' under mediator-outcome (M-O) confounding, resolving the identification
#' failure that limits the single-instrument estimators.
#'
#' The motivating example is placental eQTLs: fetal-genotype-derived eQTLs
#' instrument placental isoform expression (M), while a PFAS-metabolism PRS
#' instruments the exposure (X). The mediator set must be restricted to
#' isoforms for which eQTLs have been identified.
#'
#' Strategy (sequential 2SLS, three OLS stages):
#' \enumerate{
#' \item \code{X ~ G + W + covars} -> X_hat (purge U1 from X).
#' Weak-IV check: partial F for G >= \code{min_f}.
#' \item \code{M ~ X_hat + Gm + W + covars} -> M_hat, alpha_M = coef on X_hat.
#' Weak-IV check: partial F for Gm >= \code{min_f}.
#' \item \code{Y ~ X_hat + M_hat + W + covars} -> NDE = beta_X_hat, beta_M = coef on M_hat.
#' }
#'
#' NIE = alpha_M * beta_M (delta-method SE).
#'
#' Including X_hat in the M first-stage (stage 2) is essential: M
#' structurally depends on X, so the X -> M path must be captured for
#' alpha_M and the NIE to be correctly estimated. Gm provides the
#' exogenous variation that identifies the M -> Y effect net of U1
#' confounding.
#'
#' The sequential 2SLS implementation matches the pattern of
#' \code{\link{fit_iv2sls_mediation}} (three separate OLS stages).
#' Standard errors in the outcome stage do not account for the
#' two-stage estimation of M_hat, consistent with the existing
#' estimators; the simulation benchmarks bias, not SE accuracy. A unit
#' test cross-validates this estimator's NDE and beta_M against
#' \code{AER::ivreg} on the just-identified system
#' (\code{Y ~ X + M + W | G + Gm + W}).
#'
#' @param y Numeric outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param g Numeric instrument for X (length n).
#' @param gm Numeric instrument for M (length n).
#' @param w Numeric negative-control vector (length n).
#' @param covars Optional data frame of additional covariates (n rows).
#' @param min_f Minimum acceptable partial F-statistic for each excluded
#' instrument. Default 10.
#'
#' @return Named list: \code{NDE}, \code{NDE_se}, \code{NDE_p},
#' \code{NIE}, \code{NIE_se}, \code{NIE_p}. Returns all-\code{NA} if
#' either first-stage partial F is below \code{min_f}.
#' @export
#'
#' @references
#' Rudolph, K. E., et al. (2024). Natural direct and indirect effects
#' with an instrumental variable. *Biometrics*.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 500, mo_confounding = 0.8,
#' phi = 0.8, seed = 1)
#' fit_iv2sls_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$Gm,
#' dat$W[, 1])
#' }
fit_iv2sls_mediation2 <- function(y, X, M, g, gm, w, covars = NULL, min_f = 10) {
  NA_res <- list(NDE = NA_real_, NDE_se = NA_real_, NDE_p = NA_real_,
                 NIE = NA_real_, NIE_se = NA_real_, NIE_p = NA_real_,
                 alpha_M = NA_real_, alpha_se = NA_real_, beta_M = NA_real_, beta_M_se = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  # ── Stage 1: X ~ G (+ W) + covars -> X_hat (purge U1 from X) ──
  d_fs <- data.frame(X = X, g = g)
  d_fs <- .bind_covars(d_fs, we$df)
  d_fs <- .bind_covars(d_fs, covars)
  fs <- tryCatch(lm(as.formula(paste0("X ~ g", .plus_frag(we$frag), cs)), data = d_fs),
                 error = function(e) NULL)
  if (is.null(fs)) return(NA_res)
  Fst_g <- .partial_F(fs, "g")
  if (is.na(Fst_g) || Fst_g < min_f) return(NA_res)
  X_hat <- fitted(fs)

  # ── Stage 2: M ~ X_hat + Gm (+ W) + covars -> M_hat, alpha_M ──
  # X_hat captures the X -> M path; Gm provides exogenous variation
  # that identifies M net of U1 confounding.
  d_ms <- data.frame(M = M, X_hat = X_hat, gm = gm)
  d_ms <- .bind_covars(d_ms, we$df)
  d_ms <- .bind_covars(d_ms, covars)
  ms <- tryCatch(lm(as.formula(paste0("M ~ X_hat + gm", .plus_frag(we$frag), cs)), data = d_ms),
                 error = function(e) NULL)
  if (is.null(ms)) return(NA_res)
  Fst_gm <- .partial_F(ms, "gm")
  if (is.na(Fst_gm) || Fst_gm < min_f) return(NA_res)

  s_ms <- summary(ms)$coefficients
  if (!"X_hat" %in% rownames(s_ms)) return(NA_res)
  alpha <- as.numeric(coef(ms)["X_hat"])
  alpha_se <- as.numeric(s_ms["X_hat", 2])
  M_hat <- fitted(ms)

  # ── Stage 3: Y ~ X_hat + M_hat (+ W) + covars -> NDE, beta_M ──
  d_os <- data.frame(y = y, X_hat = X_hat, M_hat = M_hat)
  d_os <- .bind_covars(d_os, we$df)
  d_os <- .bind_covars(d_os, covars)
  os <- tryCatch(lm(as.formula(paste0("y ~ X_hat + M_hat", .plus_frag(we$frag), cs)), data = d_os),
                 error = function(e) NULL)
  if (is.null(os)) return(NA_res)
  s_os <- summary(os)$coefficients
  if (!"X_hat" %in% rownames(s_os) || !"M_hat" %in% rownames(s_os)) return(NA_res)
  beta_X <- as.numeric(coef(os)["X_hat"])
  beta_X_se <- as.numeric(s_os["X_hat", 2])
  beta_M <- as.numeric(coef(os)["M_hat"])
  beta_M_se <- as.numeric(s_os["M_hat", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_X, NDE_se = beta_X_se, NDE_p = as.numeric(s_os["X_hat", 4]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 6a. PGC mediation (matrix bridge) -- DEFAULT ──

#' PGC mediation estimator: bridge-function-adjusted natural effects (matrix bridge)
#'
#' Extends \code{\link{fit_pgc}} (matrix bridge) to the mediation setting
#' by constructing a confounding proxy \eqn{\hat W} from the full W matrix
#' and including it in both the mediator and outcome regressions.
#'
#' Steps:
#' \enumerate{
#' \item Residualise X on G -> X_resid (U-driven residual).
#' \item Bridge X_resid on the FULL W matrix -> W_hat (proxy for U).
#' \item \code{M ~ X + W_hat} -> alpha_M (adjusted for confounding proxy).
#' \item \code{Y ~ X + M + W_hat} -> NDE = beta_X, beta_M (adjusted).
#' }
#'
#' The matrix bridge requires \code{ncol(W) >= k} (proximal completeness)
#' for the bridge to span the confounder subspace.
#'
#' @param y Numeric outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param g Numeric instrument vector (length n).
#' @param W Numeric negative-control matrix (n x q) or vector
#' (length n). If a matrix, the bridge uses all q columns.
#' @param covars Optional data frame of additional covariates (n rows).
#'
#' @return Named list: \code{NDE}, \code{NDE_se}, \code{NDE_p},
#' \code{NIE}, \code{NIE_se}, \code{NIE_p}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, mo_confounding = 0.8, seed = 1)
#' fit_pgc_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$W)
#' }
fit_pgc_mediation <- function(y, X, M, g, W, covars = NULL) {
  NA_res <- list(NDE = NA_real_, NDE_se = NA_real_, NDE_p = NA_real_,
                 NIE = NA_real_, NIE_se = NA_real_, NIE_p = NA_real_,
                 alpha_M = NA_real_, alpha_se = NA_real_, beta_M = NA_real_, beta_M_se = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  # Ensure W is a matrix
  if (!is.matrix(W)) W <- as.matrix(W)

  # Step 1: residualise X on G -> U-driven residual
  d_r <- .bind_covars(data.frame(Xc = X, g = g), covars)
  fit_resid <- tryCatch(lm(as.formula(paste0("Xc ~ g", cs)), data = d_r),
                        error = function(e) NULL)
  if (is.null(fit_resid)) return(NA_res)
  X_resid <- residuals(fit_resid)

  # Step 2: bridge X_resid on the FULL W matrix -> W_hat
  d_b <- data.frame(X_resid = X_resid)
  d_b <- cbind(d_b, as.data.frame(W))
  if (!is.null(covars)) d_b <- cbind(d_b, covars)
  w_names <- paste0("W", seq_len(ncol(W)))
  names(d_b)[2:(ncol(W) + 1)] <- w_names
  fml_b <- as.formula(paste0("X_resid ~ ",
                             paste(w_names, collapse = " + "), cs))
  fit_b <- tryCatch(lm(fml_b, data = d_b), error = function(e) NULL)
  if (is.null(fit_b)) return(NA_res)
  W_hat <- fitted(fit_b)

  # Stage 1: M ~ X + W_hat
  d1 <- .bind_covars(data.frame(M = M, X = X, W_hat = W_hat), covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ X + W_hat", cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"X" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["X"])
  alpha_se <- as.numeric(s1["X", 2])

  # Stage 2: Y ~ X + M + W_hat
  d2 <- .bind_covars(data.frame(y = y, X = X, M = M, W_hat = W_hat), covars)
  fit2 <- tryCatch(lm(as.formula(paste0("y ~ X + M + W_hat", cs)), data = d2),
                   error = function(e) NULL)
  if (is.null(fit2)) return(NA_res)
  s2 <- summary(fit2)$coefficients
  if (!"X" %in% rownames(s2) || !"M" %in% rownames(s2)) return(NA_res)
  beta_X <- as.numeric(coef(fit2)["X"])
  beta_X_se <- as.numeric(s2["X", 2])
  beta_M <- as.numeric(coef(fit2)["M"])
  beta_M_se <- as.numeric(s2["M", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_X, NDE_se = beta_X_se, NDE_p = as.numeric(s2["X", 4]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 7. PGC-2 mediation (two-stage proximal, path-specific bridges) ──

#' PGC-2 mediation estimator: two-stage proximal mediation with path-specific bridges
#'
#' Extends the proximal-inference approach to the two-linked-DAG mediation
#' setting. Uses \strong{path-specific} negative controls — W1
#' capturing the X->M confounder conf_XM and W2 capturing the M->Y confounder
#' conf_MY — to purge confounding at both stages via bridge functions.
#'
#' When \code{gm} is \code{NULL} (no mediator instrument), Stage 2 residualises
#' M on the cleaned exposure X_hat to isolate the conf_MY-driven component, then
#' bridges on W2. This is pure negative-control identification at both stages
#' — no mediator instrument is required.
#'
#' When \code{gm} is supplied (mediator instrument present but its exogeneity
#' may be in doubt), Stage 2 uses Gm to help isolate conf_MY's effect on M before
#' bridging on W2. The bridge W_hat_M does the confounding removal, so the
#' estimator is robust to Gm being correlated with conf_MY — the residual bias
#' from Gm-U correlation is smaller than IV2SLS2's.
#'
#' Strategy (three stages):
#' \enumerate{
#' \item \strong{Bridge for X} (purge conf_XM from X):
#' \code{X_resid = residuals(X ~ G1 + C)};
#' \code{W_hat_X = bridge(X_resid ~ W1)} (fitted values, proxy for conf_XM);
#' \code{X_hat = fitted(X ~ G1 + W_hat_X + C)}.
#' Weak-IV check: partial F for G1 >= \code{min_f}.
#' \item \strong{Bridge for M} (purge conf_MY from M):
#' \code{M_resid = residuals(M ~ X_hat + C)} (\code{gm = NULL}), or
#' \code{M_resid = residuals(M ~ Gm + C)} (\code{gm} supplied);
#' \code{W_hat_M = bridge(M_resid ~ W2)} (fitted values, proxy for conf_MY);
#' \code{M_hat = fitted(M ~ X_hat + W_hat_M + C)} (\code{gm = NULL}), or
#' \code{M_hat = fitted(M ~ X_hat + Gm + W_hat_M + C)} (\code{gm} supplied);
#' \code{alpha_M = coefficient on X_hat}.
#' \item \strong{Outcome}:
#' \code{Y ~ X_hat + M_hat + W_hat_X + W_hat_M + C};
#' \code{NDE = coefficient on X_hat}, \code{beta_M = coefficient on M_hat}.
#' }
#'
#' NIE = alpha_M * beta_M (delta-method SE).
#'
#' The key advantage over \code{\link{fit_iv2sls_mediation2}}: PGC-2 does not
#' require instrument exogeneity. When the mediator instrument Gm is
#' correlated with the confounder conf_MY (rho_G2 > 0), IV2SLS2 is biased but
#' PGC-2's bridge absorbs conf_MY regardless of the instrument violation. The
#' tipping-point simulation maps where PGC-2 bias crosses below IV2SLS2 bias.
#'
#' @param y Numeric outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param g Numeric instrument for X (length n).
#' @param W1 Numeric negative-control matrix (n x q) or vector for the
#' X->M path (captures conf_XM). If a matrix, the bridge uses
#' all q columns.
#' @param W2 Numeric negative-control matrix (n x q) or vector for the
#' M->Y path (captures conf_MY). If a matrix, the bridge uses
#' all q columns.
#' @param gm Optional numeric mediator instrument vector (length n).
#' When \code{NULL} (default), Stage 2 uses pure NC
#' identification. When supplied, Gm helps isolate conf_MY
#' before bridging — robust to Gm-U correlation.
#' @param covars Optional data frame of additional covariates (n rows).
#' @param min_f Minimum acceptable partial F-statistic for the excluded
#' instrument G1. Default 10.
#'
#' @return Named list: \code{NDE}, \code{NDE_se}, \code{NDE_p},
#' \code{NIE}, \code{NIE_se}, \code{NIE_p}. Returns all-\code{NA} if
#' the first-stage partial F for G1 is below \code{min_f}.
#' @export
#'
#' @references
#' Miao, W., Geng, Z., & Tchetgen Tchetgen, E. (2018). Identifying causal
#' effects with proxy variables of an unmeasured confounder. *Biometrika*,
#' 105(4), 987-993.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 500, mo_confounding = 0.8,
#' rho_G2 = 0.3, lambda_XM = c(1, 0), lambda_MY = c(0, 1),
#' omega_1 = 0.7, omega_2 = 0.7, seed = 1)
#' # Without mediator instrument (pure NC identification)
#' fit_pgc_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G1, dat$W1, dat$W2)
#' # With (possibly imperfect) mediator instrument
#' fit_pgc_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G1, dat$W1, dat$W2,
#' gm = dat$Gm)
#' }
fit_pgc_mediation2 <- function(y, X, M, g, W1, W2, gm = NULL,
                               covars = NULL, min_f = 10) {
  NA_res <- list(NDE = NA_real_, NDE_se = NA_real_, NDE_p = NA_real_,
                 NIE = NA_real_, NIE_se = NA_real_, NIE_p = NA_real_,
                 alpha_M = NA_real_, alpha_se = NA_real_, beta_M = NA_real_, beta_M_se = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  # Ensure W1, W2 are matrices
  if (!is.matrix(W1)) W1 <- as.matrix(W1)
  if (!is.matrix(W2)) W2 <- as.matrix(W2)

  # ── Stage 1: Bridge for X (purge conf_XM from X) ──
  # Weak-IV check: partial F for G1 in X ~ G1 + W1 + covars
  d_fs <- .bind_covars(data.frame(X = X, g = g), covars)
  # Include W1 in the first stage so the partial F is conditional on W1
  d_fs <- cbind(d_fs, as.data.frame(W1))
  w1_names <- paste0("W1_", seq_len(ncol(W1)))
  names(d_fs)[(ncol(d_fs) - ncol(W1) + 1):ncol(d_fs)] <- w1_names
  fs_fml <- as.formula(paste0("X ~ g + ",
                              paste(w1_names, collapse = " + "), cs))
  fs <- tryCatch(lm(fs_fml, data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_res)
  Fst_g <- .partial_F(fs, "g")
  if (is.na(Fst_g) || Fst_g < min_f) return(NA_res)

  # Residualise X on G1 -> conf_XM-driven residual
  d_r <- .bind_covars(data.frame(Xc = X, g = g), covars)
  fit_resid <- tryCatch(lm(as.formula(paste0("Xc ~ g", cs)), data = d_r),
                        error = function(e) NULL)
  if (is.null(fit_resid)) return(NA_res)
  X_resid <- residuals(fit_resid)

  # Bridge X_resid on W1 -> W_hat_X (proxy for conf_XM)
  d_b1 <- data.frame(X_resid = X_resid)
  d_b1 <- cbind(d_b1, as.data.frame(W1))
  if (!is.null(covars)) d_b1 <- cbind(d_b1, covars)
  w1b_names <- paste0("W1b_", seq_len(ncol(W1)))
  names(d_b1)[2:(ncol(W1) + 1)] <- w1b_names
  fml_b1 <- as.formula(paste0("X_resid ~ ",
                              paste(w1b_names, collapse = " + "), cs))
  fit_b1 <- tryCatch(lm(fml_b1, data = d_b1), error = function(e) NULL)
  if (is.null(fit_b1)) return(NA_res)
  W_hat_X <- fitted(fit_b1)

  # X_hat = fitted(X ~ G1 + W_hat_X + C)
  d_zh <- .bind_covars(data.frame(X = X, g = g, W_hat_X = W_hat_X), covars)
  fit_zh <- tryCatch(lm(as.formula(paste0("X ~ g + W_hat_X", cs)), data = d_zh),
                     error = function(e) NULL)
  if (is.null(fit_zh)) return(NA_res)
  X_hat <- fitted(fit_zh)

  # ── Stage 2: Bridge for M (purge conf_MY from M) ──
  # Branch on gm: residualise on X_hat (pure NC) or on Gm (NC-augmented)
  if (is.null(gm)) {
    # Pure NC: M_resid = residuals(M ~ X_hat + C)
    # Isolates conf_MY's effect on M + noise (X path removed)
    d_mr <- .bind_covars(data.frame(M = M, X_hat = X_hat), covars)
    fit_mr <- tryCatch(lm(as.formula(paste0("M ~ X_hat", cs)), data = d_mr),
                       error = function(e) NULL)
    if (is.null(fit_mr)) return(NA_res)
    M_resid <- residuals(fit_mr)
  } else {
    # NC-augmented: M_resid = residuals(M ~ Gm + C)
    # Gm helps isolate conf_MY; bridge W_hat_M does confounding removal
    d_mr <- .bind_covars(data.frame(M = M, gm = gm), covars)
    fit_mr <- tryCatch(lm(as.formula(paste0("M ~ gm", cs)), data = d_mr),
                       error = function(e) NULL)
    if (is.null(fit_mr)) return(NA_res)
    M_resid <- residuals(fit_mr)
  }

  # Bridge M_resid on W2 -> W_hat_M (proxy for conf_MY)
  d_b2 <- data.frame(M_resid = M_resid)
  d_b2 <- cbind(d_b2, as.data.frame(W2))
  if (!is.null(covars)) d_b2 <- cbind(d_b2, covars)
  w2b_names <- paste0("W2b_", seq_len(ncol(W2)))
  names(d_b2)[2:(ncol(W2) + 1)] <- w2b_names
  fml_b2 <- as.formula(paste0("M_resid ~ ",
                              paste(w2b_names, collapse = " + "), cs))
  fit_b2 <- tryCatch(lm(fml_b2, data = d_b2), error = function(e) NULL)
  if (is.null(fit_b2)) return(NA_res)
  W_hat_M <- fitted(fit_b2)

  # M_hat and alpha_M
  if (is.null(gm)) {
    # M_hat = fitted(M ~ X_hat + W_hat_M + C)
    d_mh <- .bind_covars(data.frame(M = M, X_hat = X_hat, W_hat_M = W_hat_M), covars)
    fit_mh <- tryCatch(lm(as.formula(paste0("M ~ X_hat + W_hat_M", cs)),
                          data = d_mh), error = function(e) NULL)
  } else {
    # M_hat = fitted(M ~ X_hat + Gm + W_hat_M + C)
    d_mh <- .bind_covars(data.frame(M = M, X_hat = X_hat, gm = gm,
                                    W_hat_M = W_hat_M), covars)
    fit_mh <- tryCatch(lm(as.formula(paste0("M ~ X_hat + gm + W_hat_M", cs)),
                          data = d_mh), error = function(e) NULL)
  }
  if (is.null(fit_mh)) return(NA_res)
  s_mh <- summary(fit_mh)$coefficients
  if (!"X_hat" %in% rownames(s_mh)) return(NA_res)
  alpha <- as.numeric(coef(fit_mh)["X_hat"])
  alpha_se <- as.numeric(s_mh["X_hat", 2])
  M_hat <- fitted(fit_mh)

  # ── Stage 3: Outcome ──
  # Y ~ X_hat + M_hat + W_hat_X + W_hat_M + C
  d_os <- .bind_covars(data.frame(y = y, X_hat = X_hat, M_hat = M_hat,
                                  W_hat_X = W_hat_X, W_hat_M = W_hat_M), covars)
  os <- tryCatch(lm(as.formula(paste0("y ~ X_hat + M_hat + W_hat_X + W_hat_M", cs)),
                    data = d_os), error = function(e) NULL)
  if (is.null(os)) return(NA_res)
  s_os <- summary(os)$coefficients
  if (!"X_hat" %in% rownames(s_os) || !"M_hat" %in% rownames(s_os)) return(NA_res)
  beta_X <- as.numeric(coef(os)["X_hat"])
  beta_X_se <- as.numeric(s_os["X_hat", 2])
  beta_M <- as.numeric(coef(os)["M_hat"])
  beta_M_se <- as.numeric(s_os["M_hat", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_X, NDE_se = beta_X_se, NDE_p = as.numeric(s_os["X_hat", 4]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 6b. PGC mediation (scalar bridge) -- LEGACY / FALLBACK ──

#' PGC mediation estimator: bridge-function-adjusted natural effects (scalar bridge)
#'
#' The original ICONIC PGC mediation implementation, which summarises the
#' negative-control panel as a scalar (\code{rowMeans(W)}) before bridging.
#' Like \code{\link{fit_pgc_scalar}}, this version is algebraically
#' equivalent to the IV2SLS mediation estimator when the instrument is
#' valid. Use \code{\link{fit_pgc_mediation}} (matrix bridge) when the
#' completeness condition is of interest.
#'
#' @param y Numeric outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param g Numeric instrument vector (length n).
#' @param w Numeric negative-control vector (length n).
#' Pass \code{rowMeans(W_matrix)} for stability.
#' @param covars Optional data frame of additional covariates (n rows).
#'
#' @return Named list: \code{NDE}, \code{NDE_se}, \code{NDE_p},
#' \code{NIE}, \code{NIE_se}, \code{NIE_p}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, mo_confounding = 0.8, seed = 1)
#' fit_pgc_scalar_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], rowMeans(dat$W))
#' }
fit_pgc_scalar_mediation <- function(y, X, M, g, w, covars = NULL) {
  NA_res <- list(NDE = NA_real_, NDE_se = NA_real_, NDE_p = NA_real_,
                 NIE = NA_real_, NIE_se = NA_real_, NIE_p = NA_real_,
                 alpha_M = NA_real_, alpha_se = NA_real_, beta_M = NA_real_, beta_M_se = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  # Step 1: residualise X on G -> U-driven residual
  d_r <- .bind_covars(data.frame(Xc = X, g = g), covars)
  fit_resid <- tryCatch(lm(as.formula(paste0("Xc ~ g", cs)), data = d_r),
                        error = function(e) NULL)
  if (is.null(fit_resid)) return(NA_res)
  X_resid <- residuals(fit_resid)

  # Step 2: bridge W (scalar) on X_resid -> W_hat
  d_b <- .bind_covars(data.frame(w = w, X_resid = X_resid), covars)
  fit_b <- tryCatch(lm(as.formula(paste0("w ~ X_resid", cs)), data = d_b),
                    error = function(e) NULL)
  if (is.null(fit_b)) return(NA_res)
  W_hat <- fitted(fit_b)

  # Stage 1: M ~ X + W_hat
  d1 <- .bind_covars(data.frame(M = M, X = X, W_hat = W_hat), covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ X + W_hat", cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"X" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["X"])
  alpha_se <- as.numeric(s1["X", 2])

  # Stage 2: Y ~ X + M + W_hat
  d2 <- .bind_covars(data.frame(y = y, X = X, M = M, W_hat = W_hat), covars)
  fit2 <- tryCatch(lm(as.formula(paste0("y ~ X + M + W_hat", cs)), data = d2),
                   error = function(e) NULL)
  if (is.null(fit2)) return(NA_res)
  s2 <- summary(fit2)$coefficients
  if (!"X" %in% rownames(s2) || !"M" %in% rownames(s2)) return(NA_res)
  beta_X <- as.numeric(coef(fit2)["X"])
  beta_X_se <- as.numeric(s2["X", 2])
  beta_M <- as.numeric(coef(fit2)["M"])
  beta_M_se <- as.numeric(s2["M", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_X, NDE_se = beta_X_se, NDE_p = as.numeric(s2["X", 4]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── Driver: apply all mediation methods across features ──

# The base set of methods always run. IV2SLS2 is added when a mediator
# instrument (Gm) is available in the dataset. PGC2 and PGC2Gm are added
# when path-specific negative controls (W1, W2) are available.
.mediation_methods_all <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC")
.mediation_methods_with_gm <- c(.mediation_methods_all, "IV2SLS2")
.mediation_methods_with_v05 <- c(.mediation_methods_all, "PGC2")
.mediation_methods_with_v05_gm <- c(.mediation_methods_with_gm, "PGC2", "PGC2Gm")

#' Estimate mediation effects for a single feature with specified methods (internal)
#'
#' Format-agnostic per-feature mediation estimation: takes explicit vectors
#' and covariates (no dependency on the `dat` list format). Called by both
#' the simulation driver (.analyze_mediation_feature) and the real-data
#' driver (iconic_estimate).
#'
#' @param X Numeric exposure vector (length n).
#' @param y Numeric outcome vector (length n).
#' @param M_vec Numeric mediator vector (length n).
#' @param g Numeric instrument vector (length n), or NULL.
#' @param gm Numeric mediator instrument vector (length n), or NULL.
#' @param w Numeric NC vector (length n) or matrix (n x q), or NULL
#' (for DIRECT, IV2SLS, IV2SLS2 — full panel as covariates).
#' @param W_mat Numeric NC matrix (n x q), or NULL (for PGC matrix bridge).
#' @param W1_mat Numeric path-specific NC matrix (n x q) for X->M, or NULL.
#' @param W2_mat Numeric path-specific NC matrix (n x q) for M->Y, or NULL.
#' @param W_avg Numeric vector (length n): row means of NC panel for COCA.
#' If NULL but W_mat present, computed inline.
#' @param covars Optional data frame of covariates (n rows).
#' @param methods Character vector of methods to run. Default: all applicable.
#' @param feature_idx Integer or character: feature identifier. Default 1L.
#' @param se_method "delta" (default), "bootstrap", or "composite".
#' When "bootstrap", NDE_se/NIE_se are replaced by the SD of `n_boot`
#' nonparametric bootstrap resamples. Point estimates (NDE/NIE)
#' and p-values are unchanged; only the SE column is swapped.
#' When "composite", NIE_p is replaced by the Huang (2019) JT-comp
#' composite null p-value (post-processed across features by the driver).
#' @param n_boot Number of bootstrap resamples when `se_method="bootstrap"`.
#' @return Data frame: `feature`, `method`, `NDE`, `NDE_se`, `NDE_p`,
#' `NIE`, `NIE_se`, `NIE_p`. Returns NULL if < 20 complete cases.
#' When `se_method="composite"`, also includes `alpha_M`, `alpha_se`,
#' `beta_M`, `beta_M_se` for downstream composite p-value computation.
#' @keywords internal
.estimate_mediation_feature <- function(X, y, M_vec, g = NULL, gm = NULL,
                                        w = NULL, W_mat = NULL,
                                        W1_mat = NULL, W2_mat = NULL,
                                        W_avg = NULL, covars = NULL,
                                        methods = NULL,
                                        feature_idx = 1L, min_f = 10,
                                        se_method = c("delta", "bootstrap",
                                                      "composite"),
                                        n_boot = 500) {
  se_method <- match.arg(se_method)
  na <- list(NDE = NA_real_, NDE_se = NA_real_, NDE_p = NA_real_,
             NIE = NA_real_, NIE_se = NA_real_, NIE_p = NA_real_,
             alpha_M = NA_real_, alpha_se = NA_real_,
             beta_M = NA_real_, beta_M_se = NA_real_)
  row <- function(m, r) {
    # Total effect = NDE + NIE on the estimand scale. SE from the delta
    # method on the sum, treating NDE and NIE as independent (they are
    # estimated from the same outcome model but the covariance is not
    # returned by the fit_* functions; the independence approximation is
    # conservative for the sum).
    TE <- if (!is.na(r$NDE) && !is.na(r$NIE)) r$NDE + r$NIE else NA_real_
    TE_se <- if (!is.na(r$NDE_se) && !is.na(r$NIE_se))
      sqrt(r$NDE_se^2 + r$NIE_se^2) else NA_real_
    TE_p <- if (!is.na(TE) && !is.na(TE_se) && TE_se > 0)
      2 * pnorm(-abs(TE / TE_se)) else NA_real_
    data.frame(
      feature = feature_idx, method = m,
      NDE = as.numeric(r$NDE), NDE_se = as.numeric(r$NDE_se),
      NDE_p = as.numeric(r$NDE_p),
      NIE = as.numeric(r$NIE), NIE_se = as.numeric(r$NIE_se),
      NIE_p = as.numeric(r$NIE_p),
      TE = TE, TE_se = TE_se, TE_p = TE_p,
      alpha_M = as.numeric(r$alpha_M), alpha_se = as.numeric(r$alpha_se),
      beta_M = as.numeric(r$beta_M), beta_M_se = as.numeric(r$beta_M_se),
      stringsAsFactors = FALSE)
  }

  # Default methods: all that can run. The IV estimators (IV2SLS, IV2SLS2)
  # are identified by the instrument(s) alone; W is an optional proximal
  # augmentation, so they run with or without W.
  if (is.null(methods)) {
    methods <- "UNADJ"
    if (!is.null(g)) methods <- c(methods, "IV2SLS")
    if (!is.null(g) && !is.null(w)) methods <- c(methods, "DIRECT")
    if (!is.null(w)) methods <- c(methods, "COCA")
    if (!is.null(g) && !is.null(W_mat)) methods <- c(methods, "PGC")
    if (!is.null(g) && !is.null(gm)) methods <- c(methods, "IV2SLS2")
    if (!is.null(g) && !is.null(W1_mat) && !is.null(W2_mat)) methods <- c(methods, "PGC2")
    if (!is.null(g) && !is.null(W1_mat) && !is.null(W2_mat) && !is.null(gm))
      methods <- c(methods, "PGC2Gm")
  }

  # Determine which methods can actually run
  can_run <- character(0)
  if ("UNADJ" %in% methods) can_run <- c(can_run, "UNADJ")
  if ("DIRECT" %in% methods && !is.null(g) && !is.null(w)) can_run <- c(can_run, "DIRECT")
  if ("COCA" %in% methods && !is.null(w)) can_run <- c(can_run, "COCA")
  if ("IV2SLS" %in% methods && !is.null(g)) can_run <- c(can_run, "IV2SLS")
  if ("PGC" %in% methods && !is.null(g) && !is.null(W_mat)) can_run <- c(can_run, "PGC")
  if ("IV2SLS2" %in% methods && !is.null(g) && !is.null(gm))
    can_run <- c(can_run, "IV2SLS2")
  if ("PGC2" %in% methods && !is.null(g) && !is.null(W1_mat) && !is.null(W2_mat))
    can_run <- c(can_run, "PGC2")
  if ("PGC2Gm" %in% methods && !is.null(g) && !is.null(W1_mat) && !is.null(W2_mat) && !is.null(gm))
    can_run <- c(can_run, "PGC2Gm")

  if (!length(can_run)) return(NULL)

  # Complete cases
  needed <- cbind(y, X, M_vec)
  if (!is.null(g)) needed <- cbind(needed, g)
  if (!is.null(w)) needed <- cbind(needed, as.matrix(w))
  if (!is.null(gm)) needed <- cbind(needed, gm)
  if (!is.null(covars) && ncol(covars) > 0) needed <- cbind(needed, covars)
  ok <- stats::complete.cases(needed)
  if (sum(ok) < 20) return(NULL)

  X_f <- X[ok]
  y_f <- y[ok]
  M_f <- M_vec[ok]
  g_f <- if (!is.null(g)) g[ok] else NULL
  gm_f <- if (!is.null(gm)) gm[ok] else NULL
  w_f <- if (!is.null(w)) as.matrix(w)[ok, , drop = FALSE] else NULL
  cv_f <- if (!is.null(covars)) covars[ok, , drop = FALSE] else NULL
  W_mat_f <- if (!is.null(W_mat)) W_mat[ok, , drop = FALSE] else NULL
  W1_f <- if (!is.null(W1_mat)) W1_mat[ok, , drop = FALSE] else NULL
  W2_f <- if (!is.null(W2_mat)) W2_mat[ok, , drop = FALSE] else NULL
  Wa_f <- if (!is.null(W_avg)) W_avg[ok] else if (!is.null(W_mat_f)) rowMeans(W_mat_f) else NULL

  rows <- list()

  # when se_method == "bootstrap", run the estimator once for the
  # point estimate / p-value, then resample n_boot times and replace the
  # delta-method SE with the bootstrap SD. The closure `boot_fn(idx)`
  # subsets every model-specific input by `idx` so instruments, NCs, and
  # covariates are resampled in lockstep with y, X, M.
  .with_se <- function(method, fit_call, boot_call) {
    r <- tryCatch(fit_call, error = function(e) na)
    if (se_method == "bootstrap" && !is.na(r$NDE)) {
      bs <- tryCatch(bootstrap_mediation_se(boot_call, n = sum(ok),
                                            n_boot = n_boot),
                     error = function(e) NULL)
      if (!is.null(bs) && is.finite(bs$NDE_boot_se)) {
        r$NDE_se <- bs$NDE_boot_se
        r$NIE_se <- bs$NIE_boot_se
        # Recompute p-values as Wald z-tests using the bootstrap SE.
        # The delta-method p-values (t-test for NDE from summary.lm(),
        # z-test for NIE) were computed inside each fit_* function with
        # the delta SE; they are stale once the bootstrap SE is swapped in.
        r$NDE_p <- 2 * pnorm(-abs(r$NDE / r$NDE_se))
        r$NIE_p <- 2 * pnorm(-abs(r$NIE / r$NIE_se))
      }
    }
    row(method, r)
  }

  if ("UNADJ" %in% can_run) {
    rows[["UNADJ"]] <- .with_se("UNADJ",
      fit_unadj_mediation(y_f, X_f, M_f, cv_f),
      function(idx) fit_unadj_mediation(y_f[idx], X_f[idx], M_f[idx],
                                        cv_f[idx, , drop = FALSE]))
  }
  if ("DIRECT" %in% can_run) {
    rows[["DIRECT"]] <- .with_se("DIRECT",
      fit_direct_mediation(y_f, X_f, M_f, g_f, w_f, cv_f),
      function(idx) fit_direct_mediation(y_f[idx], X_f[idx], M_f[idx],
                                         g_f[idx], w_f[idx, , drop = FALSE],
                                         cv_f[idx, , drop = FALSE]))
  }
  if ("COCA" %in% can_run) {
    rows[["COCA"]] <- .with_se("COCA",
      fit_coca_mediation(y_f, X_f, M_f, Wa_f, cv_f),
      function(idx) fit_coca_mediation(y_f[idx], X_f[idx], M_f[idx],
                                       Wa_f[idx], cv_f[idx, , drop = FALSE]))
  }
  if ("IV2SLS" %in% can_run) {
    rows[["IV2SLS"]] <- .with_se("IV2SLS",
      fit_iv2sls_mediation(y_f, X_f, M_f, g_f, w_f, cv_f, min_f = min_f),
      function(idx) fit_iv2sls_mediation(y_f[idx], X_f[idx], M_f[idx],
                                         g_f[idx], w_f[idx, , drop = FALSE],
                                         cv_f[idx, , drop = FALSE],
                                         min_f = min_f))
  }
  if ("PGC" %in% can_run) {
    rows[["PGC"]] <- .with_se("PGC",
      fit_pgc_mediation(y_f, X_f, M_f, g_f, W_mat_f, cv_f),
      function(idx) fit_pgc_mediation(y_f[idx], X_f[idx], M_f[idx],
                                      g_f[idx],
                                      W_mat_f[idx, , drop = FALSE],
                                      cv_f[idx, , drop = FALSE]))
  }
  if ("IV2SLS2" %in% can_run) {
    rows[["IV2SLS2"]] <- .with_se("IV2SLS2",
      fit_iv2sls_mediation2(y_f, X_f, M_f, g_f, gm_f, w_f, cv_f, min_f = min_f),
      function(idx) fit_iv2sls_mediation2(y_f[idx], X_f[idx], M_f[idx],
                                          g_f[idx], gm_f[idx],
                                          w_f[idx, , drop = FALSE],
                                          cv_f[idx, , drop = FALSE],
                                          min_f = min_f))
  }
  if ("PGC2" %in% can_run) {
    rows[["PGC2"]] <- .with_se("PGC2",
      fit_pgc_mediation2(y_f, X_f, M_f, g_f, W1_f, W2_f,
                         gm = NULL, covars = cv_f, min_f = min_f),
      function(idx) fit_pgc_mediation2(y_f[idx], X_f[idx], M_f[idx],
                                       g_f[idx], W1_f[idx, , drop = FALSE],
                                       W2_f[idx, , drop = FALSE],
                                       gm = NULL,
                                       covars = cv_f[idx, , drop = FALSE],
                                       min_f = min_f))
  }
  if ("PGC2Gm" %in% can_run) {
    rows[["PGC2Gm"]] <- .with_se("PGC2Gm",
      fit_pgc_mediation2(y_f, X_f, M_f, g_f, W1_f, W2_f,
                         gm = gm_f, covars = cv_f, min_f = min_f),
      function(idx) fit_pgc_mediation2(y_f[idx], X_f[idx], M_f[idx],
                                       g_f[idx], W1_f[idx, , drop = FALSE],
                                       W2_f[idx, , drop = FALSE],
                                       gm = gm_f[idx],
                                       covars = cv_f[idx, , drop = FALSE],
                                       min_f = min_f))
  }

  do.call(rbind, rows)
}

#' Apply all mediation estimators to a single feature (internal)
#'
#' Thin wrapper around .estimate_mediation_feature() that extracts vectors
#' from the `dat` list format used by the simulation pipeline.
#'
#' When \code{dat$Gm} is present (i.e. a mediator instrument was supplied),
#' the 2-stage MR estimator \code{\link{fit_iv2sls_mediation2}} is also
#' run, producing a sixth row (method = "IV2SLS2"). When \code{dat$W1}
#' and \code{dat$W2} are present, the two-stage proximal
#' estimator \code{\link{fit_pgc_mediation2}} is run as "PGC2" (without
#' Gm) and "PGC2Gm" (with Gm, when Gm is also present).
#'
#' @param dat Dataset list (from generate_toy_data / run_single_iteration).
#' @param f Feature (column) index.
#' @param W_avg Row means of the full negative-control panel (for COCA).
#' @param W_valid Optional: validity-screened W matrix for matrix-bridge PGC.
#' @param se_method "delta" (default) or "bootstrap".
#' @param n_boot Number of bootstrap resamples when `se_method="bootstrap"`.
#' @return Data frame of five to eight rows (one per method) or NULL.
#' When \code{dat$M} is a matrix (n_mediators > 1), results are returned
#' for each mediator separately with a \code{mediator} column.
#' @keywords internal
.analyze_mediation_feature <- function(dat, f, W_avg, W_valid = NULL,
                                       se_method = "delta", n_boot = 500) {
  X <- dat$X
  cv <- dat$synthetic_data
  M <- dat$M
  y <- dat$Y[, f]
  g <- if (!is.null(dat$G)) dat$G[, f] else NULL
  gm <- dat$Gm # NULL when no mediator instrument
  W1 <- dat$W1 # NULL when not active
  W2 <- dat$W2 # NULL when not active

  W_mat <- if (!is.null(W_valid)) W_valid else dat$W
  # pass the full W panel (n x q) to all estimators.
  w <- W_mat # full matrix (NULL when no W)

  # when n_mediators > 1, M is a matrix (n_mediators x n).
  # Estimate mediation for each mediator separately (one M_m at a time).
  if (is.matrix(M)) {
    nm <- nrow(M)
    rows <- list()
    for (m in seq_len(nm)) {
      M_vec <- M[m, ]
      gm_m <- if (is.matrix(gm)) gm[m, ] else gm
      res <- .estimate_mediation_feature(
        X = X, y = y, M_vec = M_vec, g = g, gm = gm_m,
        w = w, W_mat = W_mat, W1_mat = W1, W2_mat = W2,
        W_avg = W_avg, covars = cv,
        methods = NULL, feature_idx = f,
        se_method = se_method, n_boot = n_boot
      )
      if (!is.null(res)) {
        res$mediator <- m
        rows[[m]] <- res
      }
    }
    do.call(rbind, rows)
  } else {
    .estimate_mediation_feature(
      X = X, y = y, M_vec = M, g = g, gm = gm,
      w = w, W_mat = W_mat, W1_mat = W1, W2_mat = W2,
      W_avg = W_avg, covars = cv,
      methods = NULL, feature_idx = f,
      se_method = se_method, n_boot = n_boot
    )
  }
}

#' Apply all mediation estimators across features (internal)
#'
#' @param dat List returned by \code{generate_toy_data()} / \code{run_single_iteration()}.
#' @param n_features Number of outcome columns to process.
#' @param W_valid Optional: validity-screened W matrix for matrix-bridge PGC.
#' @param n_cores Number of parallel workers. Default 1 (sequential).
#' @param se_method "delta" (default), "bootstrap", or
#' "composite".
#' @param n_boot Number of bootstrap resamples when `se_method="bootstrap"`.
#' @return Data frame with columns: feature, method, NDE, NDE_se, NDE_p, NIE, NIE_se, NIE_p.
#' When `se_method="composite"`, also includes alpha_M, alpha_se, beta_M,
#' beta_M_se, var_a, var_b.
#' @keywords internal
run_mediation_methods <- function(dat, n_features = ncol(dat$Y), W_valid = NULL,
                                  n_cores = 1, se_method = "delta", n_boot = 500) {
  W_avg <- rowMeans(dat$W)
  results <- .parallel_lapply(seq_len(n_features),
                    function(f) .analyze_mediation_feature(dat, f, W_avg, W_valid,
                                                           se_method = se_method,
                                                           n_boot = n_boot),
                    n_cores = n_cores)
  out <- do.call(rbind, Filter(Negate(is.null), results))
  if (se_method == "composite" && !is.null(out))
    out <- .apply_composite_pvalues(out)
  out
}


#' Run all mediation estimators on one synthetic dataset
#'
#' Reference-named entry point: applies UNADJ, DIRECT, COCA, IV2SLS, and PGC
#' (matrix bridge) mediation estimators to each tested outcome feature and
#' returns tidy per-feature results with significance flags for both NDE
#' and NIE. When the dataset includes a mediator instrument (\code{Gm}),
#' the 2-stage MR estimator \code{\link{fit_iv2sls_mediation2}} (method
#' "IV2SLS2") is also run. When the dataset includes path-specific
#' negative controls (\code{W1}, \code{W2};), the two-stage
#' proximal estimator \code{\link{fit_pgc_mediation2}} is run as "PGC2"
#' (without Gm) and "PGC2Gm" (with Gm, when Gm is also present).
#'
#' Negative controls are always summarised over the full W panel for COCA,
#' so restricting \code{test_features} does not change its control summary.
#' The matrix-bridge PGC uses the full W panel.
#'
#' A scalar-bridge variant ([fit_pgc_scalar_mediation()]) is exported for
#' standalone use but is not included in the default pipeline.
#'
#' @param iteration_data Dataset list from \code{\link{run_single_iteration}()}
#' (or \code{generate_toy_data()}). When \code{Gm} is present, the
#' 2-stage MR estimator is included.
#' @param test_features Optional integer indices of outcome features to test.
#' Default \code{NULL} (all features).
#' @param alpha Significance threshold for the significance flags. Default 0.05.
#' @param n_cores Number of parallel workers. Default 1 (sequential).
#' Uses \code{parallel::mclapply} on Unix and a PSOCK cluster on Windows.
#' @param se_method "delta" (default), "bootstrap", or
#' "composite". See
#' \code{\link{iconic_estimate}()} for details.
#' @param n_boot Number of bootstrap resamples when
#' \code{se_method = "bootstrap"}. Default 500.
#'
#' @return Data frame: \code{feature}, \code{method}, \code{NDE}, \code{NDE_se},
#' \code{NDE_p}, \code{NIE}, \code{NIE_se}, \code{NIE_p},
#' \code{NDE_significant}, \code{NIE_significant}.
#' @export
#'
#' @examples
#' \donttest{
#' dat <- iconic:::generate_toy_data(n = 300, mo_confounding = 0.8, seed = 1)
#' analyze_mediation_robust(dat)
#' }
analyze_mediation_robust <- function(iteration_data, test_features = NULL,
                                     alpha = 0.05, n_cores = 1,
                                     se_method = c("delta", "bootstrap",
                                                   "composite"),
                                     n_boot = 500) {
  se_method <- match.arg(se_method)
  feats <- if (is.null(test_features)) seq_len(ncol(iteration_data$Y)) else test_features
  W_avg <- rowMeans(iteration_data$W)
  res <- .parallel_lapply(feats,
                    function(f) .analyze_mediation_feature(iteration_data, f, W_avg,
                                                           se_method = se_method,
                                                           n_boot = n_boot),
                    n_cores = n_cores, progress = "analyze_mediation_robust")
  out <- do.call(rbind, Filter(Negate(is.null), res))
  if (is.null(out)) return(out)
  if (se_method == "composite")
    out <- .apply_composite_pvalues(out)
  out$NDE_significant <- out$NDE_p < alpha
  out$NIE_significant <- out$NIE_p < alpha
  out
}


#' Summarise mediation simulation results across features (internal)
#'
#' Iterates over the methods present in \code{combined} (rather than a
#' hardcoded list), so it handles both the 5-method (no Gm) and 6-method
#' (with Gm) cases.
#'
#' @param combined Data frame from \code{run_mediation_methods()}.
#' @param true_NDE Scalar true natural direct effect.
#' @param true_NIE Scalar true natural indirect effect.
#' @return Data frame with one row per method: NDE/NIE mean, bias, sd, rmse,
#' Type I error rates, and counts.
#' @keywords internal
summarise_mediation_results <- function(combined, true_NDE, true_NIE) {
  # Use the methods actually present in the data, preserving canonical order
  methods_present <- unique(combined$method)
  ordered <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC", "IV2SLS2",
               "PGC2", "PGC2Gm")
  methods <- intersect(ordered, methods_present)

  rows <- lapply(methods, function(m) {
    sub <- combined[combined$method == m, ]
    # CI coverage, %bias, and mean SE.
    # Wald CI coverage: fraction of replicates where true value falls in
    # estimate +/- 1.96 * SE. Requires NDE_se / NIE_se columns.
    nde_cov <- if ("NDE_se" %in% names(sub)) {
      ci_lower <- sub$NDE - 1.96 * sub$NDE_se
      ci_upper <- sub$NDE + 1.96 * sub$NDE_se
      mean(true_NDE >= ci_lower & true_NDE <= ci_upper, na.rm = TRUE)
    } else NA_real_
    nie_cov <- if ("NIE_se" %in% names(sub)) {
      ci_lower <- sub$NIE - 1.96 * sub$NIE_se
      ci_upper <- sub$NIE + 1.96 * sub$NIE_se
      mean(true_NIE >= ci_lower & true_NIE <= ci_upper, na.rm = TRUE)
    } else NA_real_
    nde_pct_bias <- if (abs(true_NDE) > 1e-8)
      (mean(sub$NDE, na.rm = TRUE) - true_NDE) / abs(true_NDE) else NA_real_
    nie_pct_bias <- if (abs(true_NIE) > 1e-8)
      (mean(sub$NIE, na.rm = TRUE) - true_NIE) / abs(true_NIE) else NA_real_

    data.frame(
      method = m,
      NDE_mean = mean(sub$NDE, na.rm = TRUE),
      NDE_bias = mean(sub$NDE, na.rm = TRUE) - true_NDE,
      NDE_pct_bias = nde_pct_bias,
      NDE_sd = sd(sub$NDE, na.rm = TRUE),
      NDE_rmse = sqrt(mean((sub$NDE - true_NDE)^2, na.rm = TRUE)),
      NDE_mean_se = if ("NDE_se" %in% names(sub)) mean(sub$NDE_se, na.rm = TRUE) else NA_real_,
      NDE_coverage = nde_cov,
      NIE_mean = mean(sub$NIE, na.rm = TRUE),
      NIE_bias = mean(sub$NIE, na.rm = TRUE) - true_NIE,
      NIE_pct_bias = nie_pct_bias,
      NIE_sd = sd(sub$NIE, na.rm = TRUE),
      NIE_rmse = sqrt(mean((sub$NIE - true_NIE)^2, na.rm = TRUE)),
      NIE_mean_se = if ("NIE_se" %in% names(sub)) mean(sub$NIE_se, na.rm = TRUE) else NA_real_,
      NIE_coverage = nie_cov,
      NIE_type1 = mean(sub$NIE_p < 0.05, na.rm = TRUE),
      NDE_type1 = mean(sub$NDE_p < 0.05, na.rm = TRUE),
      n_NDE = sum(!is.na(sub$NDE)),
      n_NIE = sum(!is.na(sub$NIE)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
