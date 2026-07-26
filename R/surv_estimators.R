# ============================================================
# Survival (time-to-event) total-effect estimators — v0.9.4
#
# Mirrors R/estimators.R but replaces the OLS outcome stage with a
# Cox proportional-hazards model (survival::coxph + Surv) or, when
# effect_scale = "rmst", an OLS regression on leave-one-out RMST
# pseudo-observations (Graw et al. 2009), which restores a linear,
# collapsible scale.
#
# First-stage regressions (Z ~ G + W, Z_resid ~ W) remain OLS: Z and
# the NC panel are continuous, so the 2SPS (two-stage predictor
# substitution) structure is identical to the continuous estimators.
# Only the outcome stage changes.
#
# COCA is structurally incompatible with survival: it regresses W on Y
# (W ~ y + Z), placing the outcome on the RHS — impossible with a Surv
# object. fit_coca_surv() therefore returns NA with an informative
# attribute, and is documented as unsupported for survival outcomes.
# ============================================================


# ── Internal helpers ──

# RMST pseudo-observations (leave-one-out jackknife of the Kaplan-Meier
# estimate of restricted mean survival time up to tau).
#
# For subject i, the pseudo-observation is
#   psi_i = n * RMST_full - (n - 1) * RMST_{-i}
# where RMST is the area under the Kaplan-Meier curve up to tau.
# Pseudo-observations are approximately unbiased for E[min(T, tau)] and
# can be regressed via OLS / GEE to obtain collapsible, linear effect
# estimates (Andersen et al. 2003; Graw et al. 2009).
#
# @param time   numeric, observed/censored follow-up times (length n).
# @param event  numeric 0/1, event indicator (1 = observed event).
# @param tau    numeric, restriction time horizon. Default: 90th
#               percentile of `time` (avoids extrapolation beyond the
#               bulk of observed follow-up).
# @return numeric vector of length n (the pseudo-observations).
.rmst_pseudo <- function(time, event, tau = NULL) {
  time  <- as.numeric(time)
  event <- as.numeric(event)
  n <- length(time)
  if (is.null(tau)) tau <- as.numeric(quantile(time, 0.90))
  if (tau >= max(time)) tau <- max(time) * 0.999

  # Full-sample RMST (area under KM curve up to tau)
  km_full <- survival::survfit(survival::Surv(time, event) ~ 1)
  rmst_full <- .km_rmst(km_full, tau)

  # Leave-one-out jackknife
  psi <- numeric(n)
  for (i in seq_len(n)) {
    km_i <- survival::survfit(survival::Surv(time[-i], event[-i]) ~ 1)
    rmst_i <- .km_rmst(km_i, tau)
    psi[i] <- n * rmst_full - (n - 1) * rmst_i
  }
  psi
}

# Area under a Kaplan-Meier curve up to tau (restricted mean survival time).
# @param km   a survfit object (KM estimate).
# @param tau  restriction time.
.km_rmst <- function(km, tau) {
  tt <- c(0, km$time)
  ss <- c(1, km$surv)
  # Keep only the part of the curve up to tau
  idx <- tt <= tau
  tt_t <- c(tt[idx], tau)
  ss_t <- c(ss[idx], ss[max(which(idx))])
  # Trapezoidal integration of S(t) from 0 to tau
  sum(diff(tt_t) * (ss_t[-length(ss_t)] + ss_t[-1]) / 2)
}

# Build the outcome-stage response for a survival estimator.
# Returns either a survival::Surv object (loghr scale) or a numeric
# vector of RMST pseudo-observations (rmst scale).
# @param time, event  survival inputs.
# @param effect_scale "loghr" or "rmst".
# @param tau           RMST horizon (used only when rmst).
.make_surv_response <- function(time, event, effect_scale, tau = NULL) {
  if (effect_scale == "rmst") {
    .rmst_pseudo(time, event, tau)
  } else {
    survival::Surv(time, event)
  }
}

# Extract (beta, se, pvalue) for `term` from a fitted model.
# Works for both coxph (loghr) and lm (rmst pseudo-observations).
.extract_surv_coef <- function(fit, term) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) return(NULL)
  list(
    b  = as.numeric(coef(fit)[term]),
    se = as.numeric(sm[term, 2]),
    p  = as.numeric(sm[term, if ("Pr(>|z|)" %in% colnames(sm)) "Pr(>|z|)" else "Pr(>|t|)"])
  )
}


# ── 1. UNADJ (survival) ──

#' UNADJ survival estimator: unadjusted Cox / RMST regression
#'
#' Regresses the survival outcome on the exposure Z (plus covariates)
#' with no instrument or negative-control adjustment.  Bias reference.
#'
#' When \code{effect_scale = "loghr"} (default), fits a Cox
#' proportional-hazards model and returns the log-hazard ratio for Z.
#' When \code{effect_scale = "rmst"}, regresses leave-one-out RMST
#' pseudo-observations (Graw et al. 2009) on Z via OLS, returning an
#' effect on the restricted-mean-survival-time (time) scale — a
#' collapsible, linear alternative to the non-collapsible hazard ratio.
#'
#' @param time         Numeric follow-up time vector (length n).
#' @param event        Numeric 0/1 event indicator (length n).
#' @param Z            Numeric exposure vector (length n).
#' @param covars       Optional data frame of covariates (n rows).
#' @param effect_scale Character: \code{"loghr"} (Cox) or \code{"rmst"}
#'   (pseudo-observation OLS). Default \code{"loghr"}.
#' @param tau          RMST restriction time horizon. Default \code{NULL}
#'   (90th percentile of follow-up). Used only when
#'   \code{effect_scale = "rmst"}.
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
#' fit_unadj_surv(dat$surv_time, dat$surv_event, dat$Z)
#' fit_unadj_surv(dat$surv_time, dat$surv_event, dat$Z, effect_scale = "rmst")
#' }
fit_unadj_surv <- function(time, event, Z, covars = NULL,
                           effect_scale = c("loghr", "rmst"), tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  resp <- .make_surv_response(time, event, effect_scale, tau)
  d <- .bind_covars(data.frame(y = resp, Z = Z), covars)
  fml <- as.formula(paste0("y ~ Z", cs))

  fit <- if (effect_scale == "loghr") {
    tryCatch(survival::coxph(fml, data = d), error = function(e) NULL)
  } else {
    tryCatch(lm(fml, data = d), error = function(e) NULL)
  }
  if (is.null(fit)) return(NA_result)

  ex <- .extract_surv_coef(fit, "Z")
  if (is.null(ex)) return(NA_result)
  list(beta = ex$b, se = ex$se, pvalue = ex$p)
}


# ── 2. DIRECT (survival) ──

#' DIRECT survival estimator: Cox / RMST with instrument and NC covariates
#'
#' Regresses the survival outcome on Z plus the genetic instrument G, the
#' negative-control panel W, and covariates.  Naive adjustment; does not
#' correct for unmeasured confounding via a ratio or IV approach.
#'
#' @param time         Numeric follow-up time vector (length n).
#' @param event        Numeric 0/1 event indicator (length n).
#' @param Z            Numeric exposure vector (length n).
#' @param g            Numeric instrument vector (length n).
#' @param w            Numeric NC vector (length n) or matrix (n x q).
#' @param covars       Optional data frame of covariates (n rows).
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau          RMST horizon (rmst only).
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
#' fit_direct_surv(dat$surv_time, dat$surv_event, dat$Z, dat$G[, 1], dat$W[, 1])
#' }
fit_direct_surv <- function(time, event, Z, g, w, covars = NULL,
                            effect_scale = c("loghr", "rmst"), tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  resp <- .make_surv_response(time, event, effect_scale, tau)
  d <- data.frame(y = resp, Z = Z, g = g)
  d <- cbind(d, we$df)
  d <- .bind_covars(d, covars)
  fml <- as.formula(paste0("y ~ Z + g + ", we$frag, cs))

  fit <- if (effect_scale == "loghr") {
    tryCatch(survival::coxph(fml, data = d), error = function(e) NULL)
  } else {
    tryCatch(lm(fml, data = d), error = function(e) NULL)
  }
  if (is.null(fit)) return(NA_result)

  ex <- .extract_surv_coef(fit, "Z")
  if (is.null(ex)) return(NA_result)
  list(beta = ex$b, se = ex$se, pvalue = ex$p)
}


# ── 3. IV2SLS (survival, 2SPS) ──

#' IV2SLS survival estimator: two-stage predictor substitution with Cox / RMST
#'
#' Two-stage predictor substitution (2SPS): the first stage regresses Z on
#' the instrument G (plus W and covariates) via OLS, producing fitted
#' \eqn{\hat Z}; the second stage regresses the survival outcome on
#' \eqn{\hat Z} (plus W and covariates) via Cox (log-HR) or RMST
#' pseudo-observation OLS.
#'
#' A weak-instrument check (partial F for the excluded instrument G,
#' Stock & Yogo 2005) is applied to the OLS first stage.  If the partial
#' F is below \code{min_f}, the function returns NA.
#'
#' @param time         Numeric follow-up time vector (length n).
#' @param event        Numeric 0/1 event indicator (length n).
#' @param Z            Numeric exposure vector (length n).
#' @param g            Numeric instrument vector (length n).
#' @param w            Numeric NC vector (length n) or matrix (n x q).
#' @param covars       Optional data frame of covariates (n rows).
#' @param min_f        Minimum partial F for the excluded instrument. Default 10.
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau          RMST horizon (rmst only).
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
#' fit_iv2sls_surv(dat$surv_time, dat$surv_event, dat$Z, dat$G[, 1], dat$W[, 1])
#' }
fit_iv2sls_surv <- function(time, event, Z, g, w, covars = NULL, min_f = 10,
                            effect_scale = c("loghr", "rmst"), tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  # ── First stage (OLS): Z ~ g + W + covars -> Z_hat ──
  d_fs <- data.frame(Z = Z, g = g)
  d_fs <- cbind(d_fs, we$df)
  d_fs <- .bind_covars(d_fs, covars)
  fml_fs <- as.formula(paste0("Z ~ g + ", we$frag, cs))
  fs <- tryCatch(lm(fml_fs, data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_result)
  Fst <- .partial_F(fs, "g")
  if (is.na(Fst) || Fst < min_f) return(NA_result)
  Z_hat <- fitted(fs)

  # ── Second stage (Cox / RMST): Surv(t,e) ~ Z_hat + W + covars ──
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d_os <- data.frame(y = resp, Z_hat = Z_hat)
  d_os <- cbind(d_os, we$df)
  d_os <- .bind_covars(d_os, covars)
  fml_os <- as.formula(paste0("y ~ Z_hat + ", we$frag, cs))

  fit <- if (effect_scale == "loghr") {
    tryCatch(survival::coxph(fml_os, data = d_os), error = function(e) NULL)
  } else {
    tryCatch(lm(fml_os, data = d_os), error = function(e) NULL)
  }
  if (is.null(fit)) return(NA_result)

  ex <- .extract_surv_coef(fit, "Z_hat")
  if (is.null(ex)) return(NA_result)
  list(beta = ex$b, se = ex$se, pvalue = ex$p)
}


# ── 4. PGC (survival, matrix bridge) ──

#' PGC survival estimator: proxy G-component correction with Cox / RMST
#'
#' Three-step bridge-function estimator with a survival outcome stage:
#' \enumerate{
#'   \item Residualise Z on G -> Z_resid (OLS).
#'   \item Bridge Z_resid on the FULL W matrix -> W_hat (OLS).
#'   \item Regress the survival outcome on Z + W_hat via Cox (log-HR) or
#'         RMST pseudo-observation OLS.
#' }
#'
#' @param time         Numeric follow-up time vector (length n).
#' @param event        Numeric 0/1 event indicator (length n).
#' @param Z            Numeric exposure vector (length n).
#' @param g            Numeric instrument vector (length n).
#' @param W            Numeric NC matrix (n x q) or vector (length n).
#' @param covars       Optional data frame of covariates (n rows).
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau          RMST horizon (rmst only).
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
#' fit_pgc_surv(dat$surv_time, dat$surv_event, dat$Z, dat$G[, 1], dat$W)
#' }
fit_pgc_surv <- function(time, event, Z, g, W, covars = NULL,
                         effect_scale = c("loghr", "rmst"), tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  if (!is.matrix(W)) W <- as.matrix(W)

  # Step 1: residualise Z on G -> Z_resid (OLS)
  d_r <- .bind_covars(data.frame(Zc = Z, g = g), covars)
  fit_resid <- tryCatch(
    lm(as.formula(paste0("Zc ~ g", cs)), data = d_r),
    error = function(e) NULL
  )
  if (is.null(fit_resid)) return(NA_result)
  Z_resid <- residuals(fit_resid)

  # Step 2: bridge Z_resid on the FULL W matrix -> W_hat (OLS)
  d_b <- data.frame(Z_resid = Z_resid)
  d_b <- cbind(d_b, as.data.frame(W))
  if (!is.null(covars)) d_b <- cbind(d_b, covars)
  w_names <- paste0("W", seq_len(ncol(W)))
  names(d_b)[2:(ncol(W) + 1)] <- w_names
  fml_b <- as.formula(paste0("Z_resid ~ ", paste(w_names, collapse = " + "), cs))
  fit_b <- tryCatch(lm(fml_b, data = d_b), error = function(e) NULL)
  if (is.null(fit_b)) return(NA_result)
  W_hat <- fitted(fit_b)

  # Step 3: survival outcome stage with W_hat as confounder proxy
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d_f <- .bind_covars(data.frame(y = resp, Z = Z, W_hat = W_hat), covars)
  fml_f <- as.formula(paste0("y ~ Z + W_hat", cs))

  fit <- if (effect_scale == "loghr") {
    tryCatch(survival::coxph(fml_f, data = d_f), error = function(e) NULL)
  } else {
    tryCatch(lm(fml_f, data = d_f), error = function(e) NULL)
  }
  if (is.null(fit)) return(NA_result)

  ex <- .extract_surv_coef(fit, "Z")
  if (is.null(ex)) return(NA_result)
  list(beta = ex$b, se = ex$se, pvalue = ex$p)
}


# ── 5. COCA (survival) — structurally unsupported ──

#' COCA survival estimator: NOT supported (returns NA)
#'
#' COCA (Correlated Outcome Control Approach) regresses the negative
#' control W on the outcome Y (\code{W ~ y + Z}) and recovers the causal
#' effect as a ratio \eqn{-\hat\beta_Z / \hat\beta_Y}.  This places the
#' outcome on the right-hand side of the regression, which is
#' structurally impossible when the outcome is a \code{survival::Surv}
#' object (time-to-event).  COCA is therefore unsupported for survival
#' outcomes and always returns \code{list(beta=NA, se=NA, pvalue=NA)}.
#'
#' @param time   Numeric follow-up time vector (length n).
#' @param event  Numeric 0/1 event indicator (length n).
#' @param Z      Numeric exposure vector (length n).
#' @param w      Numeric NC vector (length n).
#' @param covars Optional data frame of covariates (n rows).
#' @param ...    Ignored (accepted for signature compatibility).
#'
#' @return \code{list(beta = NA, se = NA, pvalue = NA)} with an
#'   informative \code{"reason"} attribute.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
#' fit_coca_surv(dat$surv_time, dat$surv_event, dat$Z, dat$W[, 1])
#' # $beta [1] NA  (COCA unsupported for survival outcomes)
#' }
fit_coca_surv <- function(time, event, Z, w, covars = NULL, ...) {
  out <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  attr(out, "reason") <- paste(
    "COCA regresses W on Y (W ~ y + Z), placing the outcome on the RHS,",
    "which is impossible with a Surv object. COCA is unsupported for",
    "survival outcomes."
  )
  out
}
