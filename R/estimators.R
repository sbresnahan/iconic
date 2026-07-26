# ============================================================
# Five exported causal estimators (plus UNADJ reference).
# Each takes a tidy data frame / vectors and returns a named
# list(beta, se, pvalue).
#
# fit_pgc() now uses a MATRIX bridge (regresses Z_resid
# on the full W matrix), making the proximal completeness
# condition (dim(W_valid) >= k) binding. The original scalar-
# bridge version is retained as fit_pgc_scalar(), which is
# algebraically equivalent to IV/2SLS when the instrument is
# valid.
# ============================================================


#Internal helpers

.covar_str <- function(covar_names) {
  if (length(covar_names) == 0) return("")
  paste0(" + ", paste(covar_names, collapse = " + "))
}

# cbind() a covariate frame onto `d` only when it is non-NULL. Guards against
# R (>= 4.4) treating `cbind(data.frame(...), NULL)` as a zero-row argument.
.bind_covars <- function(d, covars) {
  if (is.null(covars)) d else cbind(d, covars)
}

# Expand a negative-control input (vector or matrix) into a data frame
# with named columns W1, W2, ..., Wq, and return the formula fragment
# "W1 + W2 + ... + Wq" for inclusion in a regression formula. When w is
# a single vector, it is treated as one column (W1), preserving backward
# compatibility for callers that pass a scalar NC.
#
# @param w Numeric vector (length n) or matrix (n x q) of NC features.
# @return List with elements:
# $df — data.frame with q columns named W1..Wq
# $frag — character string "W1 + W2 + ... + Wq"
.expand_w <- function(w) {
  if (is.null(w)) return(NULL)
  if (!is.matrix(w)) w <- as.matrix(w)
  q <- ncol(w)
  w_names <- paste0("W", seq_len(q))
  df <- as.data.frame(w)
  names(df) <- w_names
  list(df = df, frag = paste(w_names, collapse = " + "))
}

.extract_coef <- function(fit, term) {
  sm <- summary(fit)$coefficients
  list(b = coef(fit)[term],
       se = sm[term, 2],
       p = sm[term, 4])
}

# Partial F-statistic for an excluded instrument (internal)
#
# Computes the partial F for `term` in a first-stage regression, i.e. the
# squared t-statistic for that coefficient. This is the correct weak-
# instrument diagnostic (Stock & Yogo 2005) when the first stage includes
# included instruments / covariates (e.g. W) that inflate the overall model
# F regardless of the excluded instrument's relevance.
#
# @param fs An lm object from the first-stage regression.
# @param term Name of the excluded instrument coefficient.
# @return Scalar partial F (t^2), or NA if the term is absent.
.partial_F <- function(fs, term) {
  sm <- summary(fs)$coefficients
  if (!term %in% rownames(sm)) return(NA_real_)
  as.numeric(sm[term, "t value"]^2)
}


#1. DIRECT

#' DIRECT estimator: OLS with instrument and negative-control as covariates
#'
#' Regresses Y on Z plus the genetic instrument G, the negative-control W,
#' and any additional covariates. This is a "naive" adjustment that uses
#' whatever observables are available but does NOT correct for unmeasured
#' confounding via a ratio or IV approach.
#'
#' @param y Numeric outcome vector (length n).
#' @param Z Numeric exposure vector (length n), assumed pre-scaled.
#' @param g Numeric instrument vector (length n).
#' @param w Numeric negative-control vector (length n) or matrix
#' (n x q). When a matrix, all q columns are included as
#' separate covariates.
#' @param covars Optional data frame of additional covariates (n rows).
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, seed = 1)
#' fit_direct(dat$Y[, 1], dat$Z, dat$G[, 1], dat$W[, 1])
#' }


fit_direct <- function(y, Z, g, w, covars = NULL) {
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  we <- .expand_w(w)
  d <- data.frame(y = y, Z = Z, g = g)
  d <- cbind(d, we$df)
  d <- .bind_covars(d, covars)
  fml <- as.formula(paste0("y ~ Z + g + ", we$frag, .covar_str(cnames)))

  fit <- tryCatch(lm(fml, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NA_result)

  sm <- summary(fit)$coefficients
  if (!"Z" %in% rownames(sm)) return(NA_result)

  list(
    beta = as.numeric(coef(fit)["Z"]),
    se = as.numeric(sm["Z", 2]),
    pvalue = as.numeric(sm["Z", 4])
  )
}


#2. COCA

#' COCA estimator: Negative-Control Outcome Correction via ratio
#'
#' Implements the Correlated Outcome Control Approach (COCA). Fits
#' \code{w ~ y + Z + covars} and recovers the causal effect as
#' \eqn{\hat\beta = -\hat\beta_Z / \hat\beta_Y}. Standard errors are
#' obtained via the delta method.
#'
#' The negative-control W should be an outcome that shares the same
#' unmeasured confounders as Y but has no direct causal path from Z.
#'
#' @param y Numeric primary outcome vector (length n).
#' @param Z Numeric exposure vector (length n).
#' @param w Numeric negative-control outcome vector (length n).
#' Recommended: pass \code{rowMeans(W_matrix)} for stability.
#' @param covars Optional data frame of additional covariates (n rows).
#' @param ratio_cap Maximum absolute value of the ratio estimate before
#' flagging as unstable and returning \code{NA}. Default 10.
#' @param se_cap Maximum SE before flagging as unstable. Default 5.
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' Returns \code{list(beta=NA, se=NA, pvalue=NA)} if estimation is
#' unstable (near-zero \eqn{\hat\beta_Y} or extreme ratio).
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, seed = 1)
#' fit_coca(dat$Y[, 1], dat$Z, rowMeans(dat$W))
#' }


fit_coca <- function(y, Z, w, covars = NULL, ratio_cap = 10, se_cap = 5) {
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  d <- .bind_covars(data.frame(w = w, y = y, Z = Z), covars)
  fml <- as.formula(paste0("w ~ y + Z", .covar_str(cnames)))

  fit <- tryCatch(lm(fml, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NA_result)

  b <- coef(fit)
  bZ <- b["Z"]
  bY <- b["y"]

  if (is.na(bY) || abs(bY) < 1e-8) return(NA_result)

  bhat <- -bZ / bY
  V <- vcov(fit)[c("Z", "y"), c("Z", "y")]
  grad <- c(-1 / bY, bZ / bY^2)
  se_h <- sqrt(as.numeric(t(grad) %*% V %*% grad))

  if (abs(bhat) > ratio_cap || is.nan(se_h) || se_h > se_cap) return(NA_result)

  list(
    beta = bhat,
    se = se_h,
    pvalue = 2 * (1 - pnorm(abs(bhat / se_h)))
  )
}


#3. IV2SLS
#' IV2SLS estimator: Two-Stage Least Squares with genetic instrument
#'
#' Uses the genetic instrument G to instrument for the exposure Z,
#' controlling for the negative-control W and any additional covariates.
#' Requires \pkg{AER}.
#'
#' A weak-instrument check is applied using the partial F-statistic for the
#' excluded instrument G (testing G conditional on W and covariates in the
#' first stage), following Stock & Yogo (2005). If the partial F is below
#' \code{min_f}, the function returns \code{NA}.
#'
#' @param y Numeric outcome vector (length n).
#' @param Z Numeric exposure vector (length n).
#' @param g Numeric instrument vector (length n).
#' @param w Numeric negative-control vector (length n).
#' @param covars Optional data frame of additional covariates (n rows).
#' @param min_f Minimum acceptable partial F-statistic for the excluded
#' instrument. Default 10.
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, seed = 1)
#' fit_iv2sls(dat$Y[, 1], dat$Z, dat$G[, 1], dat$W[, 1])
#' }


fit_iv2sls <- function(y, Z, g, w, covars = NULL, min_f = 10) {
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  # First-stage weak-instrument check: partial F for the excluded instrument G
  # (not the overall model F, which is inflated by W sharing confounders with Z)
  d_fs <- data.frame(Z = Z, g = g)
  d_fs <- cbind(d_fs, we$df)
  d_fs <- .bind_covars(d_fs, covars)
  fml_fs <- as.formula(paste0("Z ~ g + ", we$frag, cs))
  fs <- tryCatch(lm(fml_fs, data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_result)
  Fst <- .partial_F(fs, "g")
  if (is.na(Fst) || Fst < min_f) return(NA_result)

  # 2SLS via AER::ivreg
  d_iv <- data.frame(y = y, Z = Z, G_inst = g)
  d_iv <- cbind(d_iv, we$df)
  d_iv <- .bind_covars(d_iv, covars)
  fml_2sls <- as.formula(
    paste0("y ~ Z + ", we$frag, cs, " | G_inst + ", we$frag, cs)
  )
  fit <- tryCatch(
    AER::ivreg(fml_2sls, data = d_iv),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_result)
  if (!"Z" %in% names(coef(fit))) return(NA_result)

  sm <- summary(fit)$coefficients
  if (!"Z" %in% rownames(sm)) return(NA_result)

  list(
    beta = as.numeric(coef(fit)["Z"]),
    se = as.numeric(sm["Z", 2]),
    pvalue = as.numeric(sm["Z", 4])
  )
}


#4a. PGC (matrix bridge) -- DEFAULT

#' PGC estimator: Proxy G-Component Correction (matrix bridge)
#'
#' A three-step bridge-function estimator:
#' \enumerate{
#' \item Residualises Z on G to isolate the U-driven component Z_resid.
#' \item Regresses Z_resid on the FULL W matrix to construct W_hat,
#' a proxy for unmeasured confounding. This step requires
#' \code{ncol(W) >= k} (the proximal completeness condition):
#' if W has fewer valid columns than confounders, the bridge
#' cannot span the confounder subspace and the estimator is
#' under-identified.
#' \item Fits Y ~ Z + W_hat to absorb confounding bias.
#' }
#'
#' Unlike the scalar version (\code{\link{fit_pgc_scalar}}), which
#' collapses W to \code{rowMeans(W)} and is algebraically equivalent
#' to IV/2SLS, the matrix bridge preserves the dimensional structure
#' of W and is the estimator for which the completeness condition is
#' binding.
#'
#' @param y Numeric outcome vector (length n).
#' @param Z Numeric exposure vector (length n).
#' @param g Numeric instrument vector (length n).
#' @param W Numeric negative-control matrix (n x q) or vector
#' (length n). If a matrix, the bridge uses all q
#' columns. Pass only validity-screened columns for
#' the completeness condition to be meaningful.
#' @param covars Optional data frame of additional covariates (n rows).
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, seed = 1)
#' fit_pgc(dat$Y[, 1], dat$Z, dat$G[, 1], dat$W)
#' }
fit_pgc <- function(y, Z, g, W, covars = NULL) {
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  # Ensure W is a matrix
  if (!is.matrix(W)) W <- as.matrix(W)

  # Step 1: residualise Z on G -> U-driven residual
  d_r <- .bind_covars(data.frame(Zc = Z, g = g), covars)
  fit_resid <- tryCatch(
    lm(as.formula(paste0("Zc ~ g", cs)), data = d_r),
    error = function(e) NULL
  )
  if (is.null(fit_resid)) return(NA_result)
  Z_resid <- residuals(fit_resid)

  # Step 2: bridge Z_resid on the FULL W matrix
  # (not rowMeans -- this is what makes completeness binding)
  d_b <- data.frame(Z_resid = Z_resid)
  d_b <- cbind(d_b, as.data.frame(W))
  if (!is.null(covars)) d_b <- cbind(d_b, covars)
  w_names <- paste0("W", seq_len(ncol(W)))
  names(d_b)[2:(ncol(W) + 1)] <- w_names
  fml_b <- as.formula(paste0("Z_resid ~ ",
                             paste(w_names, collapse = " + "), cs))
  fit_b <- tryCatch(lm(fml_b, data = d_b), error = function(e) NULL)
  if (is.null(fit_b)) return(NA_result)
  W_hat <- fitted(fit_b)

  # Step 3: final outcome regression with W_hat as confounder proxy
  d_f <- .bind_covars(data.frame(y = y, Z = Z, W_hat = W_hat), covars)
  fit_f <- tryCatch(
    lm(as.formula(paste0("y ~ Z + W_hat", cs)), data = d_f),
    error = function(e) NULL
  )
  if (is.null(fit_f)) return(NA_result)
  if (!"Z" %in% names(coef(fit_f))) return(NA_result)

  sm <- summary(fit_f)$coefficients
  if (!"Z" %in% rownames(sm)) return(NA_result)

  list(
    beta = as.numeric(coef(fit_f)["Z"]),
    se = as.numeric(sm["Z", 2]),
    pvalue = as.numeric(sm["Z", 4])
  )
}


#4b. PGC (scalar bridge) -- LEGACY / FALLBACK

#' PGC estimator: Proxy G-Component Correction (scalar bridge)
#'
#' The original ICONIC PGC implementation, which summarises the
#' negative-control panel as a scalar (\code{rowMeans(W)}) before
#' bridging. This version is numerically stable and works in small
#' samples, but the scalar bridge produces \eqn{\hat W} proportional
#' to the G-residualised exposure by construction, making the
#' estimator algebraically equivalent to IV/2SLS when the instrument
#' is valid. The proximal completeness condition
#' (\code{dim(W_valid) >= k}) is NOT binding for this estimator.
#'
#' Use \code{\link{fit_pgc}} (matrix bridge) when the completeness
#' condition is of interest. Use this function as a stable fallback
#' when \code{n} is small relative to the number of valid controls,
#' or as an IV-equivalent benchmark.
#'
#' @param y Numeric outcome vector (length n).
#' @param Z Numeric exposure vector (length n).
#' @param g Numeric instrument vector (length n).
#' @param w Numeric negative-control vector (length n).
#' Pass \code{rowMeans(W_matrix)} for stability.
#' @param covars Optional data frame of additional covariates (n rows).
#'
#' @return Named list: \code{beta}, \code{se}, \code{pvalue}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, seed = 1)
#' fit_pgc_scalar(dat$Y[, 1], dat$Z, dat$G[, 1], rowMeans(dat$W))
#' }
fit_pgc_scalar <- function(y, Z, g, w, covars = NULL) {
  NA_result <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  # Step 1: residualise Z on G -> U-driven residual
  d_r <- .bind_covars(data.frame(Zc = Z, g = g), covars)
  fit_resid <- tryCatch(
    lm(as.formula(paste0("Zc ~ g", cs)), data = d_r),
    error = function(e) NULL
  )
  if (is.null(fit_resid)) return(NA_result)
  Z_resid <- residuals(fit_resid)

  # Step 2: bridge W (scalar) on the U-component
  d_b <- .bind_covars(data.frame(w = w, Z_resid = Z_resid), covars)
  fit_b <- tryCatch(
    lm(as.formula(paste0("w ~ Z_resid", cs)), data = d_b),
    error = function(e) NULL
  )
  if (is.null(fit_b)) return(NA_result)
  W_hat <- fitted(fit_b)

  # Step 3: final outcome regression with W_hat as confounder proxy
  d_f <- .bind_covars(data.frame(y = y, Z = Z, W_hat = W_hat), covars)
  fit_f <- tryCatch(
    lm(as.formula(paste0("y ~ Z + W_hat", cs)), data = d_f),
    error = function(e) NULL
  )
  if (is.null(fit_f)) return(NA_result)
  if (!"Z" %in% names(coef(fit_f))) return(NA_result)

  sm <- summary(fit_f)$coefficients
  if (!"Z" %in% rownames(sm)) return(NA_result)

  list(
    beta = as.numeric(coef(fit_f)["Z"]),
    se = as.numeric(sm["Z", 2]),
    pvalue = as.numeric(sm["Z", 4])
  )
}
