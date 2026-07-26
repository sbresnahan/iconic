# ============================================================
# Survival (time-to-event) mediation estimators — v0.9.4
#
# Mirrors the continuous mediation estimators in R/mediation.R using
# two-stage predictor substitution (2SPS): the mediator first stages
# (M ~ Z / M ~ Z_hat + Gm) and the confounder bridges (Z_resid ~ W,
# M_resid ~ W2) remain OLS, because M and the NC panels are continuous.
# Only the outcome stage (stage 3) switches to a Cox proportional-
# hazards model (log-HR scale) or an OLS regression on RMST pseudo-
# observations (collapsible time scale).
#
# NIE = alpha_M * beta_M (product of the stage-2 alpha and the stage-3
# coefficient on M_hat).  On the log-HR scale this product is an
# approximation because the hazard ratio is non-collapsible; on the RMST
# scale the decomposition is exact.  See the man-pages and the v0.9.4
# methods notes for the non-collapsibility caveat.
#
# COCA mediation is structurally unsupported for survival (it regresses
# W on Y in both stages) and returns NA.
# ============================================================


# Shared NA-result template for mediation survival estimators.
.surv_med_NA <- function() {
  list(NDE = NA_real_, NDE_se = NA_real_, NDE_p = NA_real_,
       NIE = NA_real_, NIE_se = NA_real_, NIE_p = NA_real_,
       alpha_M = NA_real_, alpha_se = NA_real_,
       beta_M = NA_real_, beta_M_se = NA_real_)
}

# Fit the outcome stage (stage 3) on the survival response and extract
# NDE (coef on Z_hat) and beta_M (coef on M_hat).
# @param resp     Surv object (loghr) or numeric pseudo-observations (rmst).
# @param d_os     data frame with columns y, Z_hat, M_hat, [W_hat_*], [covars].
# @param fml_os   outcome-stage formula (character).
# @param scale    "loghr" or "rmst".
# @return list(NDE, NDE_se, NDE_p, beta_M, beta_M_se) or NULL on failure.
.fit_surv_outcome_stage <- function(resp, d_os, fml_os, scale,
                                     nde_name = "Z_hat", med_name = "M_hat") {
  d_os$y <- resp
  fit <- if (scale == "loghr") {
    tryCatch(survival::coxph(fml_os, data = d_os), error = function(e) NULL)
  } else {
    tryCatch(lm(fml_os, data = d_os), error = function(e) NULL)
  }
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)$coefficients
  pcol <- if ("Pr(>|z|)" %in% colnames(sm)) "Pr(>|z|)" else "Pr(>|t|)"
  if (!nde_name %in% rownames(sm) || !med_name %in% rownames(sm)) return(NULL)
  list(
    NDE       = as.numeric(coef(fit)[nde_name]),
    NDE_se    = as.numeric(sm[nde_name, 2]),
    NDE_p     = as.numeric(sm[nde_name, pcol]),
    beta_M    = as.numeric(coef(fit)[med_name]),
    beta_M_se = as.numeric(sm[med_name, 2])
  )
}


# ── 1. UNADJ mediation (survival) ──

#' UNADJ survival mediation estimator: naive Baron-Kenny with Cox / RMST
#'
#' Stage 1 (OLS): \code{M ~ Z} -> alpha_M.  Stage 2 (Cox / RMST):
#' survival outcome ~ Z + M -> NDE (coef on Z), beta_M (coef on M).
#' \code{NIE = alpha_M * beta_M}.  No confounding adjustment; bias
#' reference.
#'
#' @param time         Numeric follow-up time vector (length n).
#' @param event        Numeric 0/1 event indicator (length n).
#' @param Z            Numeric exposure vector (length n).
#' @param M            Numeric mediator vector (length n).
#' @param covars       Optional data frame of covariates (n rows).
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau          RMST horizon (rmst only).
#'
#' @return Named list: \code{NDE}, \code{NDE_se}, \code{NDE_p},
#'   \code{NIE}, \code{NIE_se}, \code{NIE_p}, \code{alpha_M},
#'   \code{alpha_se}, \code{beta_M}, \code{beta_M_se}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival",
#'                                   mo_confounding = 0.8, seed = 1)
#' fit_unadj_mediation_surv(dat$surv_time, dat$surv_event, dat$Z, dat$M)
#' }
fit_unadj_mediation_surv <- function(time, event, Z, M, covars = NULL,
                                     effect_scale = c("loghr", "rmst"),
                                     tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  # Stage 1 (OLS): M ~ Z -> alpha_M
  d1 <- .bind_covars(data.frame(M = M, Z = Z), covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ Z", cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"Z" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["Z"])
  alpha_se <- as.numeric(s1["Z", 2])

  # Stage 2 (Cox / RMST): outcome ~ Z + M
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d2 <- .bind_covars(data.frame(Z = Z, M = M), covars)
  d2$y <- resp
  fml2 <- as.formula(paste0("y ~ Z + M", cs))
  fit2 <- if (effect_scale == "loghr") {
    tryCatch(survival::coxph(fml2, data = d2), error = function(e) NULL)
  } else {
    tryCatch(lm(fml2, data = d2), error = function(e) NULL)
  }
  if (is.null(fit2)) return(NA_res)
  sm <- summary(fit2)$coefficients
  pcol <- if ("Pr(>|z|)" %in% colnames(sm)) "Pr(>|z|)" else "Pr(>|t|)"
  if (!"Z" %in% rownames(sm) || !"M" %in% rownames(sm)) return(NA_res)
  beta_Z <- as.numeric(coef(fit2)["Z"])
  beta_Z_se <- as.numeric(sm["Z", 2])
  beta_M <- as.numeric(coef(fit2)["M"])
  beta_M_se <- as.numeric(sm["M", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_Z, NDE_se = beta_Z_se, NDE_p = as.numeric(sm["Z", pcol]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 2. DIRECT mediation (survival) ──

#' DIRECT survival mediation estimator: Cox / RMST with G and W covariates
#'
#' Adjusts for the instrument G and negative-control W in both the
#' mediator (OLS) and outcome (Cox / RMST) stages.  Naive adjustment.
#'
#' @param time         Numeric follow-up time vector (length n).
#' @param event        Numeric 0/1 event indicator (length n).
#' @param Z            Numeric exposure vector (length n).
#' @param M            Numeric mediator vector (length n).
#' @param g            Numeric instrument vector (length n).
#' @param w            Numeric NC vector (length n) or matrix (n x q).
#' @param covars       Optional data frame of covariates (n rows).
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau          RMST horizon (rmst only).
#'
#' @return Named list (same fields as \code{\link{fit_unadj_mediation_surv}}).
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival",
#'                                   mo_confounding = 0.8, seed = 1)
#' fit_direct_mediation_surv(dat$surv_time, dat$surv_event, dat$Z, dat$M,
#'                           dat$G[, 1], dat$W[, 1])
#' }
fit_direct_mediation_surv <- function(time, event, Z, M, g, w, covars = NULL,
                                      effect_scale = c("loghr", "rmst"),
                                      tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  # Stage 1 (OLS): M ~ Z + g + W
  d1 <- data.frame(M = M, Z = Z, g = g)
  d1 <- cbind(d1, we$df)
  d1 <- .bind_covars(d1, covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ Z + g + ", we$frag, cs)),
                      data = d1), error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"Z" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["Z"])
  alpha_se <- as.numeric(s1["Z", 2])

  # Stage 2 (Cox / RMST): outcome ~ Z + M + g + W
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d2 <- data.frame(Z = Z, M = M, g = g)
  d2 <- cbind(d2, we$df)
  d2 <- .bind_covars(d2, covars)
  d2$y <- resp
  fml2 <- as.formula(paste0("y ~ Z + M + g + ", we$frag, cs))
  fit2 <- if (effect_scale == "loghr") {
    tryCatch(survival::coxph(fml2, data = d2), error = function(e) NULL)
  } else {
    tryCatch(lm(fml2, data = d2), error = function(e) NULL)
  }
  if (is.null(fit2)) return(NA_res)
  sm <- summary(fit2)$coefficients
  pcol <- if ("Pr(>|z|)" %in% colnames(sm)) "Pr(>|z|)" else "Pr(>|t|)"
  if (!"Z" %in% rownames(sm) || !"M" %in% rownames(sm)) return(NA_res)
  beta_Z <- as.numeric(coef(fit2)["Z"])
  beta_Z_se <- as.numeric(sm["Z", 2])
  beta_M <- as.numeric(coef(fit2)["M"])
  beta_M_se <- as.numeric(sm["M", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_Z, NDE_se = beta_Z_se, NDE_p = as.numeric(sm["Z", pcol]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 3. IV2SLS mediation (survival, single instrument, 2SPS) ──

#' IV2SLS survival mediation estimator: single-instrument 2SPS with Cox / RMST
#'
#' Two-stage predictor substitution with a single instrument (G for Z, no
#' mediator instrument):
#' \enumerate{
#'   \item OLS: \code{Z ~ g + W + covars} -> Z_hat (purge U1 from Z).
#'   \item OLS: \code{M ~ Z_hat + covars} -> alpha_M.
#'   \item Cox / RMST: \code{Surv(t,e) ~ Z_hat + M + W + covars} ->
#'         NDE (coef on Z_hat), beta_M (coef on M).
#' }
#' \code{NIE = alpha_M * beta_M}.  Weak-instrument gate (partial F for G)
#' applies to the OLS first stage.  Unlike
#' \code{\link{fit_iv2sls_mediation2_surv}}, this estimator does not
#' instrument the mediator and is NOT point-identified under M-O
#' confounding — it is included for parity with the continuous
#' \code{\link{fit_iv2sls_mediation}}.
#'
#' @param time         Numeric follow-up time vector (length n).
#' @param event        Numeric 0/1 event indicator (length n).
#' @param Z            Numeric exposure vector (length n).
#' @param M            Numeric mediator vector (length n).
#' @param g            Numeric instrument for Z (length n).
#' @param w            Numeric NC vector (length n) or matrix (n x q).
#' @param covars       Optional data frame of covariates (n rows).
#' @param min_f        Minimum partial F for the excluded instrument. Default 10.
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau          RMST horizon (rmst only).
#'
#' @return Named list (same fields as \code{\link{fit_unadj_mediation_surv}}).
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
#' fit_iv2sls_mediation_surv(dat$surv_time, dat$surv_event, dat$Z, dat$M,
#'                           dat$G[, 1], dat$W[, 1])
#' }
fit_iv2sls_mediation_surv <- function(time, event, Z, M, g, w,
                                      covars = NULL, min_f = 10,
                                      effect_scale = c("loghr", "rmst"),
                                      tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  # Stage 1 (OLS): Z ~ g + W + covars -> Z_hat
  d_fs <- data.frame(Z = Z, g = g)
  d_fs <- cbind(d_fs, we$df)
  d_fs <- .bind_covars(d_fs, covars)
  fs <- tryCatch(lm(as.formula(paste0("Z ~ g + ", we$frag, cs)),
                    data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_res)
  Fst <- .partial_F(fs, "g")
  if (is.na(Fst) || Fst < min_f) return(NA_res)
  Z_hat <- fitted(fs)

  # Stage 2 (OLS): M ~ Z_hat -> alpha_M
  d1 <- .bind_covars(data.frame(M = M, Z_hat = Z_hat), covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ Z_hat", cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"Z_hat" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["Z_hat"])
  alpha_se <- as.numeric(s1["Z_hat", 2])

  # Stage 3 (Cox / RMST): outcome ~ Z_hat + M + W + covars
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d_os <- data.frame(Z_hat = Z_hat, M = M)
  d_os <- cbind(d_os, we$df)
  d_os <- .bind_covars(d_os, covars)
  d_os$y <- resp
  fml_os <- as.formula(paste0("y ~ Z_hat + M + ", we$frag, cs))
  fit_os <- if (effect_scale == "loghr") {
    tryCatch(survival::coxph(fml_os, data = d_os), error = function(e) NULL)
  } else {
    tryCatch(lm(fml_os, data = d_os), error = function(e) NULL)
  }
  if (is.null(fit_os)) return(NA_res)
  sm <- summary(fit_os)$coefficients
  pcol <- if ("Pr(>|z|)" %in% colnames(sm)) "Pr(>|z|)" else "Pr(>|t|)"
  if (!"Z_hat" %in% rownames(sm) || !"M" %in% rownames(sm)) return(NA_res)
  beta_Z <- as.numeric(coef(fit_os)["Z_hat"])
  beta_Z_se <- as.numeric(sm["Z_hat", 2])
  beta_M <- as.numeric(coef(fit_os)["M"])
  beta_M_se <- as.numeric(sm["M", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_Z, NDE_se = beta_Z_se, NDE_p = as.numeric(sm["Z_hat", pcol]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 4. IV2SLS2 mediation (survival, 2SPS) ──

#' IV2SLS2 survival mediation estimator: 2-stage MR with Cox / RMST outcome
#'
#' Two-stage Mendelian-randomization mediation with a survival outcome
#' stage (2SPS):
#' \enumerate{
#'   \item OLS: \code{Z ~ g + W + covars} -> Z_hat (purge U1 from Z).
#'   \item OLS: \code{M ~ Z_hat + gm + W + covars} -> M_hat, alpha_M.
#'   \item Cox / RMST: \code{Surv(t,e) ~ Z_hat + M_hat + W + covars} ->
#'         NDE (coef on Z_hat), beta_M (coef on M_hat).
#' }
#' \code{NIE = alpha_M * beta_M}.  Weak-instrument gates (partial F for
#' G and Gm) apply to the OLS first stages.
#'
#' @param time         Numeric follow-up time vector (length n).
#' @param event        Numeric 0/1 event indicator (length n).
#' @param Z            Numeric exposure vector (length n).
#' @param M            Numeric mediator vector (length n).
#' @param g            Numeric instrument for Z (length n).
#' @param gm           Numeric instrument for M (length n).
#' @param w            Numeric NC vector (length n) or matrix (n x q).
#' @param covars       Optional data frame of covariates (n rows).
#' @param min_f        Minimum partial F for each excluded instrument. Default 10.
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau          RMST horizon (rmst only).
#'
#' @return Named list (same fields as \code{\link{fit_unadj_mediation_surv}}).
#'   Returns all-NA if either first-stage partial F is below \code{min_f}.
#' @export
#'
#' @references
#' Rudolph, K. E., et al. (2024). Natural direct and indirect effects
#' with an instrumental variable. *Biometrics*.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 500, outcome_type = "survival",
#'                                   mo_confounding = 0.8, phi = 0.8, seed = 1)
#' fit_iv2sls_mediation2_surv(dat$surv_time, dat$surv_event, dat$Z, dat$M,
#'                            dat$G[, 1], dat$Gm, dat$W[, 1])
#' }
fit_iv2sls_mediation2_surv <- function(time, event, Z, M, g, gm, w,
                                       covars = NULL, min_f = 10,
                                       effect_scale = c("loghr", "rmst"),
                                       tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  # Stage 1 (OLS): Z ~ g + W + covars -> Z_hat
  d_fs <- data.frame(Z = Z, g = g)
  d_fs <- cbind(d_fs, we$df)
  d_fs <- .bind_covars(d_fs, covars)
  fs <- tryCatch(lm(as.formula(paste0("Z ~ g + ", we$frag, cs)),
                    data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_res)
  Fst_g <- .partial_F(fs, "g")
  if (is.na(Fst_g) || Fst_g < min_f) return(NA_res)
  Z_hat <- fitted(fs)

  # Stage 2 (OLS): M ~ Z_hat + gm + W + covars -> M_hat, alpha_M
  d_ms <- data.frame(M = M, Z_hat = Z_hat, gm = gm)
  d_ms <- cbind(d_ms, we$df)
  d_ms <- .bind_covars(d_ms, covars)
  ms <- tryCatch(lm(as.formula(paste0("M ~ Z_hat + gm + ", we$frag, cs)),
                    data = d_ms), error = function(e) NULL)
  if (is.null(ms)) return(NA_res)
  Fst_gm <- .partial_F(ms, "gm")
  if (is.na(Fst_gm) || Fst_gm < min_f) return(NA_res)
  s_ms <- summary(ms)$coefficients
  if (!"Z_hat" %in% rownames(s_ms)) return(NA_res)
  alpha <- as.numeric(coef(ms)["Z_hat"])
  alpha_se <- as.numeric(s_ms["Z_hat", 2])
  M_hat <- fitted(ms)

  # Stage 3 (Cox / RMST): outcome ~ Z_hat + M_hat + W + covars
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d_os <- data.frame(Z_hat = Z_hat, M_hat = M_hat)
  d_os <- cbind(d_os, we$df)
  d_os <- .bind_covars(d_os, covars)
  fml_os <- as.formula(paste0("y ~ Z_hat + M_hat + ", we$frag, cs))
  os <- .fit_surv_outcome_stage(resp, d_os, fml_os, effect_scale)
  if (is.null(os)) return(NA_res)

  NIE <- alpha * os$beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, os$beta_M, os$beta_M_se)

  list(
    NDE = os$NDE, NDE_se = os$NDE_se, NDE_p = os$NDE_p,
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se,
    beta_M = os$beta_M, beta_M_se = os$beta_M_se
  )
}


# ── 4. PGC mediation (survival, single-panel matrix bridge) ──

#' PGC survival mediation estimator: single-panel bridge with Cox / RMST
#'
#' Bridge-function-adjusted natural direct and indirect effects with a
#' survival outcome stage, using a single negative-control panel W for
#' both the Z->M and M->Y confounding paths.  This is the survival
#' analogue of \code{\link{fit_pgc_mediation}} (continuous outcome):
#' \enumerate{
#'   \item OLS: residualise Z on G -> Z_resid.
#'   \item OLS: bridge Z_resid on the FULL W matrix -> W_hat (proxy for U).
#'   \item OLS: \code{M ~ Z + W_hat + covars} -> alpha_M.
#'   \item Cox / RMST: \code{Surv(t,e) ~ Z + M + W_hat + covars} ->
#'         NDE (coef on Z), beta_M (coef on M).
#' }
#' \code{NIE = alpha_M * beta_M}.  Unlike
#' \code{\link{fit_pgc_mediation2_surv}}, which uses path-specific W1/W2
#' bridges, this estimator uses a single combined W panel and is
#' appropriate when separate U_XM / U_MY confounders are not assumed.
#'
#' @param time         Numeric follow-up time vector (length n).
#' @param event        Numeric 0/1 event indicator (length n).
#' @param Z            Numeric exposure vector (length n).
#' @param M            Numeric mediator vector (length n).
#' @param g            Numeric instrument for Z (length n).
#' @param W            Numeric NC matrix (n x q) or vector (length n).
#' @param covars       Optional data frame of covariates (n rows).
#' @param min_f        Minimum partial F for G. Default 10.
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau          RMST horizon (rmst only).
#'
#' @return Named list (same fields as \code{\link{fit_unadj_mediation_surv}}).
#'   Returns all-NA if the first-stage partial F for G is below \code{min_f}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
#' fit_pgc_mediation_surv(dat$surv_time, dat$surv_event, dat$Z, dat$M,
#'                        dat$G[, 1], dat$W)
#' }
fit_pgc_mediation_surv <- function(time, event, Z, M, g, W, covars = NULL,
                                   min_f = 10,
                                   effect_scale = c("loghr", "rmst"),
                                   tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  if (!is.matrix(W)) W <- as.matrix(W)

  # Weak-instrument check: partial F for G in Z ~ G + W + covars
  d_fs <- .bind_covars(data.frame(Z = Z, g = g), covars)
  d_fs <- cbind(d_fs, as.data.frame(W))
  w_fs_names <- paste0("Wfs_", seq_len(ncol(W)))
  names(d_fs)[(ncol(d_fs) - ncol(W) + 1):ncol(d_fs)] <- w_fs_names
  fs_fml <- as.formula(paste0("Z ~ g + ",
                              paste(w_fs_names, collapse = " + "), cs))
  fs <- tryCatch(lm(fs_fml, data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_res)
  Fst_g <- .partial_F(fs, "g")
  if (is.na(Fst_g) || Fst_g < min_f) return(NA_res)

  # Step 1: residualise Z on G -> Z_resid (OLS)
  d_r <- .bind_covars(data.frame(Zc = Z, g = g), covars)
  fit_resid <- tryCatch(lm(as.formula(paste0("Zc ~ g", cs)), data = d_r),
                        error = function(e) NULL)
  if (is.null(fit_resid)) return(NA_res)
  Z_resid <- residuals(fit_resid)

  # Step 2: bridge Z_resid on the FULL W matrix -> W_hat (OLS)
  d_b <- data.frame(Z_resid = Z_resid)
  d_b <- cbind(d_b, as.data.frame(W))
  if (!is.null(covars)) d_b <- cbind(d_b, covars)
  w_names <- paste0("W", seq_len(ncol(W)))
  names(d_b)[2:(ncol(W) + 1)] <- w_names
  fml_b <- as.formula(paste0("Z_resid ~ ", paste(w_names, collapse = " + "), cs))
  fit_b <- tryCatch(lm(fml_b, data = d_b), error = function(e) NULL)
  if (is.null(fit_b)) return(NA_res)
  W_hat <- fitted(fit_b)

  # Stage 1 (OLS): M ~ Z + W_hat + covars -> alpha_M
  d1 <- .bind_covars(data.frame(M = M, Z = Z, W_hat = W_hat), covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ Z + W_hat", cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"Z" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["Z"])
  alpha_se <- as.numeric(s1["Z", 2])

  # Stage 2 (Cox / RMST): outcome ~ Z + M + W_hat + covars
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d_os <- .bind_covars(data.frame(Z = Z, M = M, W_hat = W_hat), covars)
  fml_os <- as.formula(paste0("y ~ Z + M + W_hat", cs))
  os <- .fit_surv_outcome_stage(resp, d_os, fml_os, effect_scale,
                                nde_name = "Z", med_name = "M")
  if (is.null(os)) return(NA_res)

  # Numerical-stability guard (see fit_pgc_mediation2_surv for rationale).
  if (!is.finite(os$NDE) || abs(os$NDE) > 10 ||
      !is.finite(os$beta_M) || abs(os$beta_M) > 10)
    return(NA_res)

  NIE <- alpha * os$beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, os$beta_M, os$beta_M_se)

  list(
    NDE = os$NDE, NDE_se = os$NDE_se, NDE_p = os$NDE_p,
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se,
    beta_M = os$beta_M, beta_M_se = os$beta_M_se
  )
}


# ── 5. PGC2 / PGC2Gm mediation (survival, path-specific bridges) ──

#' PGC2 / PGC2Gm survival mediation estimator: path-specific bridges with Cox / RMST
#'
#' Two-stage proximal mediation with path-specific negative controls and
#' a survival outcome stage.  Stages 1-2 (Z bridge on W1, M bridge on W2,
#' Z_hat / M_hat construction) are identical to \code{\link{fit_pgc_mediation2}}
#' and remain OLS.  Only stage 3 switches to Cox (log-HR) or RMST
#' pseudo-observation OLS:
#' \code{Surv(t,e) ~ Z_hat + M_hat + W_hat_Z + W_hat_M + covars}.
#'
#' When \code{gm = NULL} (PGC2), stage 2 uses pure NC identification.
#' When \code{gm} is supplied (PGC2Gm), the mediator instrument helps
#' isolate U_MY before bridging W2.
#'
#' @param time         Numeric follow-up time vector (length n).
#' @param event        Numeric 0/1 event indicator (length n).
#' @param Z            Numeric exposure vector (length n).
#' @param M            Numeric mediator vector (length n).
#' @param g            Numeric instrument for Z (length n).
#' @param W1           Numeric NC matrix (n x q) or vector for the Z->M path.
#' @param W2           Numeric NC matrix (n x q) or vector for the M->Y path.
#' @param gm           Optional numeric mediator instrument (length n).
#' @param covars       Optional data frame of covariates (n rows).
#' @param min_f        Minimum partial F for G1. Default 10.
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau          RMST horizon (rmst only).
#'
#' @return Named list (same fields as \code{\link{fit_unadj_mediation_surv}}).
#'   Returns all-NA if the first-stage partial F for G1 is below \code{min_f}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 500, outcome_type = "survival",
#'                                   mo_confounding = 0.8, rho_G2 = 0.3,
#'                                   separate_U = TRUE, seed = 1)
#' fit_pgc_mediation2_surv(dat$surv_time, dat$surv_event, dat$Z, dat$M,
#'                         dat$G[, 1], dat$W1, dat$W2, gm = dat$Gm)
#' }
fit_pgc_mediation2_surv <- function(time, event, Z, M, g, W1, W2, gm = NULL,
                                    covars = NULL, min_f = 10,
                                    effect_scale = c("loghr", "rmst"),
                                    tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  if (!is.matrix(W1)) W1 <- as.matrix(W1)
  if (!is.matrix(W2)) W2 <- as.matrix(W2)

  # ── Stage 1: Bridge for Z (purge U_XM) — OLS, identical to continuous ──
  d_fs <- .bind_covars(data.frame(Z = Z, g = g), covars)
  d_fs <- cbind(d_fs, as.data.frame(W1))
  w1_names <- paste0("W1_", seq_len(ncol(W1)))
  names(d_fs)[(ncol(d_fs) - ncol(W1) + 1):ncol(d_fs)] <- w1_names
  fs_fml <- as.formula(paste0("Z ~ g + ",
                              paste(w1_names, collapse = " + "), cs))
  fs <- tryCatch(lm(fs_fml, data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_res)
  Fst_g <- .partial_F(fs, "g")
  if (is.na(Fst_g) || Fst_g < min_f) return(NA_res)

  d_r <- .bind_covars(data.frame(Zc = Z, g = g), covars)
  fit_resid <- tryCatch(lm(as.formula(paste0("Zc ~ g", cs)), data = d_r),
                        error = function(e) NULL)
  if (is.null(fit_resid)) return(NA_res)
  Z_resid <- residuals(fit_resid)

  d_b1 <- data.frame(Z_resid = Z_resid)
  d_b1 <- cbind(d_b1, as.data.frame(W1))
  if (!is.null(covars)) d_b1 <- cbind(d_b1, covars)
  w1b_names <- paste0("W1b_", seq_len(ncol(W1)))
  names(d_b1)[2:(ncol(W1) + 1)] <- w1b_names
  fml_b1 <- as.formula(paste0("Z_resid ~ ",
                              paste(w1b_names, collapse = " + "), cs))
  fit_b1 <- tryCatch(lm(fml_b1, data = d_b1), error = function(e) NULL)
  if (is.null(fit_b1)) return(NA_res)
  W_hat_Z <- fitted(fit_b1)

  d_zh <- .bind_covars(data.frame(Z = Z, g = g, W_hat_Z = W_hat_Z), covars)
  fit_zh <- tryCatch(lm(as.formula(paste0("Z ~ g + W_hat_Z", cs)),
                        data = d_zh), error = function(e) NULL)
  if (is.null(fit_zh)) return(NA_res)
  Z_hat <- fitted(fit_zh)

  # ── Stage 2: Bridge for M (purge U_MY) — OLS, identical to continuous ──
  if (is.null(gm)) {
    d_mr <- .bind_covars(data.frame(M = M, Z_hat = Z_hat), covars)
    fit_mr <- tryCatch(lm(as.formula(paste0("M ~ Z_hat", cs)),
                          data = d_mr), error = function(e) NULL)
  } else {
    d_mr <- .bind_covars(data.frame(M = M, gm = gm), covars)
    fit_mr <- tryCatch(lm(as.formula(paste0("M ~ gm", cs)),
                          data = d_mr), error = function(e) NULL)
  }
  if (is.null(fit_mr)) return(NA_res)
  M_resid <- residuals(fit_mr)

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

  if (is.null(gm)) {
    d_mh <- .bind_covars(data.frame(M = M, Z_hat = Z_hat,
                                    W_hat_M = W_hat_M), covars)
    fit_mh <- tryCatch(lm(as.formula(paste0("M ~ Z_hat + W_hat_M", cs)),
                          data = d_mh), error = function(e) NULL)
  } else {
    d_mh <- .bind_covars(data.frame(M = M, Z_hat = Z_hat, gm = gm,
                                    W_hat_M = W_hat_M), covars)
    fit_mh <- tryCatch(lm(as.formula(paste0("M ~ Z_hat + gm + W_hat_M", cs)),
                          data = d_mh), error = function(e) NULL)
  }
  if (is.null(fit_mh)) return(NA_res)
  s_mh <- summary(fit_mh)$coefficients
  if (!"Z_hat" %in% rownames(s_mh)) return(NA_res)
  alpha <- as.numeric(coef(fit_mh)["Z_hat"])
  alpha_se <- as.numeric(s_mh["Z_hat", 2])
  M_hat <- fitted(fit_mh)

  # ── Stage 3 (Cox / RMST): outcome ~ Z_hat + M_hat + W_hat_Z + W_hat_M ──
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d_os <- data.frame(Z_hat = Z_hat, M_hat = M_hat,
                     W_hat_Z = W_hat_Z, W_hat_M = W_hat_M)
  d_os <- .bind_covars(d_os, covars)
  fml_os <- as.formula(paste0("y ~ Z_hat + M_hat + W_hat_Z + W_hat_M", cs))
  os <- .fit_surv_outcome_stage(resp, d_os, fml_os, effect_scale)
  if (is.null(os)) return(NA_res)

  # Numerical-stability guard: the path-specific bridge can produce
  # extreme generated regressors that cause the Cox partial likelihood
  # to return wildly unstable coefficients (RMSE > 4 in benchmarks).
  # Reject estimates whose absolute value exceeds a threshold that is
  # implausible on either scale (log-HR > 10 ≈ HR 22,000; RMST shift
  # > 10 time units is equally extreme for standardised data).
  if (!is.finite(os$NDE) || abs(os$NDE) > 10 ||
      !is.finite(os$beta_M) || abs(os$beta_M) > 10)
    return(NA_res)

  NIE <- alpha * os$beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, os$beta_M, os$beta_M_se)

  list(
    NDE = os$NDE, NDE_se = os$NDE_se, NDE_p = os$NDE_p,
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se,
    beta_M = os$beta_M, beta_M_se = os$beta_M_se
  )
}


# ── 5. COCA mediation (survival) — structurally unsupported ──

#' COCA survival mediation estimator: NOT supported (returns NA)
#'
#' COCA mediation regresses the negative control W on the outcome Y in
#' both stages (\code{W ~ M + Z}, \code{W ~ y + Z + M}), placing the
#' outcome on the right-hand side.  This is structurally impossible when
#' the outcome is a \code{survival::Surv} object.  COCA mediation is
#' therefore unsupported for survival outcomes and always returns NA.
#'
#' @param time   Numeric follow-up time vector (length n).
#' @param event  Numeric 0/1 event indicator (length n).
#' @param Z      Numeric exposure vector (length n).
#' @param M      Numeric mediator vector (length n).
#' @param w      Numeric NC vector (length n).
#' @param covars Optional data frame of covariates (n rows).
#' @param ...    Ignored (signature compatibility).
#'
#' @return All-NA list with an informative \code{"reason"} attribute.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
#' fit_coca_mediation_surv(dat$surv_time, dat$surv_event, dat$Z, dat$M, dat$W[, 1])
#' # all NA (COCA mediation unsupported for survival outcomes)
#' }
fit_coca_mediation_surv <- function(time, event, Z, M, w, covars = NULL, ...) {
  out <- .surv_med_NA()
  attr(out, "reason") <- paste(
    "COCA mediation regresses W on Y (W ~ y + Z + M), placing the outcome",
    "on the RHS, which is impossible with a Surv object. COCA mediation",
    "is unsupported for survival outcomes."
  )
  out
}
