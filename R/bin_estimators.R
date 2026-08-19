# ============================================================
# Binary (0/1) total-effect estimators
#
# Mirrors R/estimators.R but replaces the OLS outcome stage with a
# logistic regression (stats::glm, family = binomial) when
# effect_scale = "logor" (the default), or with an OLS regression on
# the 0/1 outcome — a linear probability model — when
# effect_scale = "riskdiff", a collapsible scale on which the
# NDE/NIE product decomposition is exact. This mirrors the survival
# package's Cox log-HR / RMST duality.
#
# First-stage regressions (X ~ G + W, X_resid ~ W) remain OLS: X and
# the NC panel are continuous, so the 2SPS (two-stage predictor
# substitution) structure is identical to the continuous estimators.
# Only the outcome stage changes.
#
# COCA is unsupported for binary outcomes: it regresses W on Y
# (W ~ y + X) and recovers the effect as the ratio -beta_X / beta_Y,
# an identification argument that assumes a linear structural outcome
# model. With a binary outcome the structural model is nonlinear
# (logistic), so the linear COCA ratio recovers neither the causal
# log-OR nor the risk difference. fit_coca_bin() therefore returns NA
# with an informative attribute.
# ============================================================


# ── Internal helpers ──

# Fit the binary outcome stage: logistic regression (logor scale) or
# OLS on the 0/1 outcome (riskdiff scale; linear probability model).
# Returns NULL on fit failure (callers convert to an NA result).
.fit_bin_model <- function(fml, data, effect_scale) {
  if (effect_scale == "logor") {
    tryCatch(stats::glm(fml, data = data, family = stats::binomial()),
             error = function(e) NULL)
  } else {
    tryCatch(stats::lm(fml, data = data), error = function(e) NULL)
  }
}

# Coefficient extraction is scale-agnostic: .extract_surv_coef() (in
# R/surv_estimators.R) reads the Estimate / Std. Error / p-value columns
# from any model with a coefficient table, selecting Pr(>|z|) for glm
# fits and Pr(>|t|) for lm fits. The binary estimators reuse it
# unchanged.


# ── 1. UNADJ (binary) ──

#' UNADJ binary estimator: unadjusted logistic / linear-probability regression
#'
#' Regresses the binary outcome on the exposure X (plus covariates)
#' with no instrument or negative-control adjustment. Bias reference.
#'
#' When \code{effect_scale = "logor"} (default), fits a logistic
#' regression and returns the conditional log-odds ratio for X. When
#' \code{effect_scale = "riskdiff"}, fits a linear probability model
#' (OLS on the 0/1 outcome), returning a risk difference — a
#' collapsible, linear alternative to the non-collapsible odds ratio.
#'
#' @param y Numeric 0/1 outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param covars Optional data frame of covariates (n rows).
#' @param effect_scale Character: \code{"logor"} (logistic) or
#' \code{"riskdiff"} (linear probability model). Default \code{"logor"}.
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- generate_toy_data(n = 200, outcome_type = "binary", seed = 1)
#' fit_unadj_bin(dat$y_bin, dat$X)
#' fit_unadj_bin(dat$y_bin, dat$X, effect_scale = "riskdiff")
fit_unadj_bin <- function(y, X, covars = NULL,
                          effect_scale = c("logor", "riskdiff")) {
  effect_scale <- match.arg(effect_scale)
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  d <- .bind_covars(data.frame(y = y, X = X), covars)
  fml <- as.formula(paste0("y ~ X", cs))
  fit <- .fit_bin_model(fml, d, effect_scale)
  if (is.null(fit)) return(NA_result)

  ex <- .extract_surv_coef(fit, "X")
  if (is.null(ex)) return(NA_result)
  list(beta = ex$b, se = ex$se, pvalue = ex$p)
}


# ── 2. DIRECT (binary) ──

#' DIRECT binary estimator: logistic / LPM with instrument and NC covariates
#'
#' Regresses the binary outcome on X plus the genetic instrument G, the
#' negative-control panel W, and covariates. Naive adjustment; does not
#' correct for unmeasured confounding via a ratio or IV approach.
#'
#' @param y Numeric 0/1 outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param g Numeric instrument vector (length n).
#' @param w Numeric NC vector (length n) or matrix (n x q).
#' @param covars Optional data frame of covariates (n rows).
#' @param effect_scale Character: \code{"logor"} or \code{"riskdiff"}.
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- generate_toy_data(n = 200, outcome_type = "binary", seed = 1)
#' fit_direct_bin(dat$y_bin, dat$X, dat$G[, 1], dat$W[, 1])
fit_direct_bin <- function(y, X, g, w, covars = NULL,
                           effect_scale = c("logor", "riskdiff")) {
  effect_scale <- match.arg(effect_scale)
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  d <- data.frame(y = y, X = X, g = g)
  d <- cbind(d, we$df)
  d <- .bind_covars(d, covars)
  fml <- as.formula(paste0("y ~ X + g + ", we$frag, cs))
  fit <- .fit_bin_model(fml, d, effect_scale)
  if (is.null(fit)) return(NA_result)

  ex <- .extract_surv_coef(fit, "X")
  if (is.null(ex)) return(NA_result)
  list(beta = ex$b, se = ex$se, pvalue = ex$p)
}


# ── 3. IV2SLS (binary, 2SPS) ──

#' IV2SLS binary estimator: two-stage predictor substitution with logistic / LPM
#'
#' Two-stage predictor substitution (2SPS): the first stage regresses X on
#' the instrument G (plus W and covariates) via OLS, producing fitted
#' \eqn{\hat X}; the second stage regresses the binary outcome on
#' \eqn{\hat X} (plus W and covariates) via logistic regression (log-OR)
#' or a linear probability model (risk difference).
#'
#' A weak-instrument check (partial F for the excluded instrument G,
#' Stock & Yogo 2005) is applied to the OLS first stage. If the partial
#' F is below \code{min_f}, the function returns NA.
#'
#' @param y Numeric 0/1 outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param g Numeric instrument vector (length n).
#' @param w Numeric NC vector (length n) or matrix (n x q).
#' @param covars Optional data frame of covariates (n rows).
#' @param min_f Minimum partial F for the excluded instrument. Default 10.
#' @param effect_scale Character: \code{"logor"} or \code{"riskdiff"}.
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- generate_toy_data(n = 200, outcome_type = "binary", seed = 1)
#' fit_iv2sls_bin(dat$y_bin, dat$X, dat$G[, 1], dat$W[, 1])
fit_iv2sls_bin <- function(y, X, g, w, covars = NULL, min_f = 10,
                           effect_scale = c("logor", "riskdiff")) {
  effect_scale <- match.arg(effect_scale)
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  # ── First stage (OLS): X ~ g (+ W) + covars -> X_hat ──
  d_fs <- data.frame(X = X, g = g)
  d_fs <- .bind_covars(d_fs, we$df)
  d_fs <- .bind_covars(d_fs, covars)
  fml_fs <- as.formula(paste0("X ~ g", .plus_frag(we$frag), cs))
  fs <- tryCatch(lm(fml_fs, data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_result)
  Fst <- .partial_F(fs, "g")
  if (is.na(Fst) || Fst < min_f) return(NA_result)
  X_hat <- fitted(fs)

  # ── Second stage (logistic / LPM): y ~ X_hat (+ W) + covars ──
  d_os <- data.frame(y = y, X_hat = X_hat)
  d_os <- .bind_covars(d_os, we$df)
  d_os <- .bind_covars(d_os, covars)
  fml_os <- as.formula(paste0("y ~ X_hat", .plus_frag(we$frag), cs))
  fit <- .fit_bin_model(fml_os, d_os, effect_scale)
  if (is.null(fit)) return(NA_result)

  ex <- .extract_surv_coef(fit, "X_hat")
  if (is.null(ex)) return(NA_result)
  list(beta = ex$b, se = ex$se, pvalue = ex$p)
}


# ── 4. PGC (binary, matrix bridge) ──

#' PGC binary estimator: proxy G-component correction with logistic / LPM
#'
#' Three-step bridge-function estimator with a binary outcome stage:
#' \enumerate{
#' \item Residualise X on G -> X_resid (OLS).
#' \item Bridge X_resid on the FULL W matrix -> W_hat (OLS).
#' \item Regress the binary outcome on X + W_hat via logistic
#' regression (log-OR) or a linear probability model (risk difference).
#' }
#'
#' @param y Numeric 0/1 outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param g Numeric instrument vector (length n).
#' @param W Numeric NC matrix (n x q) or vector (length n).
#' @param covars Optional data frame of covariates (n rows).
#' @param effect_scale Character: \code{"logor"} or \code{"riskdiff"}.
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- generate_toy_data(n = 200, outcome_type = "binary", seed = 1)
#' fit_pgc_bin(dat$y_bin, dat$X, dat$G[, 1], dat$W)
fit_pgc_bin <- function(y, X, g, W, covars = NULL,
                        effect_scale = c("logor", "riskdiff")) {
  effect_scale <- match.arg(effect_scale)
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  if (!is.matrix(W)) W <- as.matrix(W)

  # Step 1: residualise X on G -> X_resid (OLS)
  d_r <- .bind_covars(data.frame(Xc = X, g = g), covars)
  fit_resid <- tryCatch(
    lm(as.formula(paste0("Xc ~ g", cs)), data = d_r),
    error = function(e) NULL
  )
  if (is.null(fit_resid)) return(NA_result)
  X_resid <- residuals(fit_resid)

  # Step 2: bridge X_resid on the FULL W matrix -> W_hat (OLS)
  d_b <- data.frame(X_resid = X_resid)
  d_b <- cbind(d_b, as.data.frame(W))
  if (!is.null(covars)) d_b <- cbind(d_b, covars)
  w_names <- paste0("W", seq_len(ncol(W)))
  names(d_b)[2:(ncol(W) + 1)] <- w_names
  fml_b <- as.formula(paste0("X_resid ~ ", paste(w_names, collapse = " + "), cs))
  fit_b <- tryCatch(lm(fml_b, data = d_b), error = function(e) NULL)
  if (is.null(fit_b)) return(NA_result)
  W_hat <- fitted(fit_b)

  # Step 3: binary outcome stage with W_hat as confounder proxy
  d_f <- .bind_covars(data.frame(y = y, X = X, W_hat = W_hat), covars)
  fml_f <- as.formula(paste0("y ~ X + W_hat", cs))
  fit <- .fit_bin_model(fml_f, d_f, effect_scale)
  if (is.null(fit)) return(NA_result)

  ex <- .extract_surv_coef(fit, "X")
  if (is.null(ex)) return(NA_result)
  list(beta = ex$b, se = ex$se, pvalue = ex$p)
}


# ── 5. COCA (binary) — unsupported ──

#' COCA binary estimator: NOT supported (returns NA)
#'
#' COCA (Correlated Outcome Control Approach) regresses the negative
#' control W on the outcome Y (\code{W ~ y + X}) and recovers the causal
#' effect as a ratio \eqn{-\hat\beta_X / \hat\beta_Y}. That
#' identification argument assumes a linear structural outcome model:
#' the W-Y regression coefficient must be proportional to the
#' confounder-Outcome association on the same scale as the causal
#' effect. With a binary outcome the structural model is nonlinear
#' (logistic), so the linear COCA ratio recovers neither the causal
#' log-odds ratio nor the risk difference. COCA is therefore
#' unsupported for binary outcomes and always returns
#' \code{list(beta=NA, se=NA, pvalue=NA)}.
#'
#' @param y Numeric 0/1 outcome vector (length n).
#' @param X Numeric exposure vector (length n).
#' @param w Numeric NC vector (length n).
#' @param covars Optional data frame of covariates (n rows).
#' @param ... Ignored (accepted for signature compatibility).
#'
#' @return \code{list(beta = NA, se = NA, pvalue = NA)} with an
#' informative \code{"reason"} attribute.
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- generate_toy_data(n = 200, outcome_type = "binary", seed = 1)
#' fit_coca_bin(dat$y_bin, dat$X, dat$W[, 1])
#' # $beta [1] NA (COCA unsupported for binary outcomes)
fit_coca_bin <- function(y, X, w, covars = NULL, ...) {
  out <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  attr(out, "reason") <- paste(
    "COCA regresses W on Y (W ~ y + X) and recovers the effect as a",
    "ratio, an identification argument that assumes a linear structural",
    "outcome model. With a binary (nonlinear, logistic) outcome the",
    "ratio recovers neither the log-OR nor the risk difference. COCA is",
    "unsupported for binary outcomes."
  )
  out
}
