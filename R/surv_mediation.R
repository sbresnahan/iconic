# ============================================================
# Survival (time-to-event) mediation estimators
#
# Mirrors the continuous mediation estimators in R/mediation.R using
# two-stage predictor substitution (2SPS): the mediator first stages
# (M ~ X / M ~ X_hat + Gm) and the confounder bridges (X_resid ~ W,
# M_resid ~ W2) remain OLS, because M and the NC panels are continuous.
# Only the outcome stage (stage 3) switches to a Cox proportional-
# hazards model (log-HR scale) or an OLS regression on RMST pseudo-
# observations (collapsible time scale).
#
# NIE = alpha_M * beta_M (product of the stage-2 alpha and the stage-3
# coefficient on M_hat). On the log-HR scale this product is an
# approximation because the hazard ratio is non-collapsible; on the RMST
# scale the decomposition is exact. See the man-pages and the
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
# NDE (coef on X_hat) and beta_M (coef on M_hat).
# @param resp Surv object (loghr) or numeric pseudo-observations (rmst).
# @param d_os data frame with columns y, X_hat, M_hat, [W_hat_*], [covars].
# @param fml_os outcome-stage formula (character).
# @param scale "loghr" or "rmst".
# @return list(NDE, NDE_se, NDE_p, beta_M, beta_M_se) or NULL on failure.
.fit_surv_outcome_stage <- function(resp, d_os, fml_os, scale,
                                     nde_name = "X_hat", med_name = "M_hat") {
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
    NDE = as.numeric(coef(fit)[nde_name]),
    NDE_se = as.numeric(sm[nde_name, 2]),
    NDE_p = as.numeric(sm[nde_name, pcol]),
    beta_M = as.numeric(coef(fit)[med_name]),
    beta_M_se = as.numeric(sm[med_name, 2])
  )
}


# ── 1. UNADJ mediation (survival) ──

#' UNADJ survival mediation estimator: naive Baron-Kenny with Cox / RMST
#'
#' Stage 1 (OLS): \code{M ~ X} -> alpha_M. Stage 2 (Cox / RMST):
#' survival outcome ~ X + M -> NDE (coef on X), beta_M (coef on M).
#' \code{NIE = alpha_M * beta_M}. No confounding adjustment; bias
#' reference.
#'
#' @param time Numeric follow-up time vector (length n).
#' @param event Numeric 0/1 event indicator (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param covars Optional data frame of covariates (n rows).
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau RMST horizon (rmst only).
#'
#' @return Named list: \code{NDE}, \code{NDE_se}, \code{NDE_p},
#' \code{NIE}, \code{NIE_se}, \code{NIE_p}, \code{alpha_M},
#' \code{alpha_se}, \code{beta_M}, \code{beta_M_se}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival",
#' mo_confounding = 0.8, seed = 1)
#' fit_unadj_mediation_surv(dat$surv_time, dat$surv_event, dat$X, dat$M)
#' }
fit_unadj_mediation_surv <- function(time, event, X, M, covars = NULL,
                                     effect_scale = c("loghr", "rmst"),
                                     tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  # Stage 1 (OLS): M ~ X -> alpha_M
  d1 <- .bind_covars(data.frame(M = M, X = X), covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ X", cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"X" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["X"])
  alpha_se <- as.numeric(s1["X", 2])

  # Stage 2 (Cox / RMST): outcome ~ X + M
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d2 <- .bind_covars(data.frame(X = X, M = M), covars)
  d2$y <- resp
  fml2 <- as.formula(paste0("y ~ X + M", cs))
  fit2 <- if (effect_scale == "loghr") {
    tryCatch(survival::coxph(fml2, data = d2), error = function(e) NULL)
  } else {
    tryCatch(lm(fml2, data = d2), error = function(e) NULL)
  }
  if (is.null(fit2)) return(NA_res)
  sm <- summary(fit2)$coefficients
  pcol <- if ("Pr(>|z|)" %in% colnames(sm)) "Pr(>|z|)" else "Pr(>|t|)"
  if (!"X" %in% rownames(sm) || !"M" %in% rownames(sm)) return(NA_res)
  beta_X <- as.numeric(coef(fit2)["X"])
  beta_X_se <- as.numeric(sm["X", 2])
  beta_M <- as.numeric(coef(fit2)["M"])
  beta_M_se <- as.numeric(sm["M", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_X, NDE_se = beta_X_se, NDE_p = as.numeric(sm["X", pcol]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 2. DIRECT mediation (survival) ──

#' DIRECT survival mediation estimator: Cox / RMST with G and W covariates
#'
#' Adjusts for the instrument G and negative-control W in both the
#' mediator (OLS) and outcome (Cox / RMST) stages. Naive adjustment.
#'
#' @param time Numeric follow-up time vector (length n).
#' @param event Numeric 0/1 event indicator (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param g Numeric instrument vector (length n).
#' @param w Numeric NC vector (length n) or matrix (n x q).
#' @param covars Optional data frame of covariates (n rows).
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau RMST horizon (rmst only).
#'
#' @return Named list (same fields as \code{\link{fit_unadj_mediation_surv}}).
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival",
#' mo_confounding = 0.8, seed = 1)
#' fit_direct_mediation_surv(dat$surv_time, dat$surv_event, dat$X, dat$M,
#' dat$G[, 1], dat$W[, 1])
#' }
fit_direct_mediation_surv <- function(time, event, X, M, g, w, covars = NULL,
                                      effect_scale = c("loghr", "rmst"),
                                      tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  # Stage 1 (OLS): M ~ X + g + W
  d1 <- data.frame(M = M, X = X, g = g)
  d1 <- cbind(d1, we$df)
  d1 <- .bind_covars(d1, covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ X + g + ", we$frag, cs)),
                      data = d1), error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"X" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["X"])
  alpha_se <- as.numeric(s1["X", 2])

  # Stage 2 (Cox / RMST): outcome ~ X + M + g + W
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d2 <- data.frame(X = X, M = M, g = g)
  d2 <- cbind(d2, we$df)
  d2 <- .bind_covars(d2, covars)
  d2$y <- resp
  fml2 <- as.formula(paste0("y ~ X + M + g + ", we$frag, cs))
  fit2 <- if (effect_scale == "loghr") {
    tryCatch(survival::coxph(fml2, data = d2), error = function(e) NULL)
  } else {
    tryCatch(lm(fml2, data = d2), error = function(e) NULL)
  }
  if (is.null(fit2)) return(NA_res)
  sm <- summary(fit2)$coefficients
  pcol <- if ("Pr(>|z|)" %in% colnames(sm)) "Pr(>|z|)" else "Pr(>|t|)"
  if (!"X" %in% rownames(sm) || !"M" %in% rownames(sm)) return(NA_res)
  beta_X <- as.numeric(coef(fit2)["X"])
  beta_X_se <- as.numeric(sm["X", 2])
  beta_M <- as.numeric(coef(fit2)["M"])
  beta_M_se <- as.numeric(sm["M", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_X, NDE_se = beta_X_se, NDE_p = as.numeric(sm["X", pcol]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 3. IV2SLS mediation (survival, single instrument, 2SPS) ──

#' IV2SLS survival mediation estimator: single-instrument 2SPS with Cox / RMST
#'
#' Two-stage predictor substitution with a single instrument (G for X, no
#' mediator instrument):
#' \enumerate{
#' \item OLS: \code{X ~ g + W + covars} -> X_hat (purge U1 from X).
#' \item OLS: \code{M ~ X_hat + covars} -> alpha_M.
#' \item Cox / RMST: \code{Surv(t,e) ~ X_hat + M + W + covars} ->
#' NDE (coef on X_hat), beta_M (coef on M).
#' }
#' \code{NIE = alpha_M * beta_M}. Weak-instrument gate (partial F for G)
#' applies to the OLS first stage. Unlike
#' \code{\link{fit_iv2sls_mediation2_surv}}, this estimator does not
#' instrument the mediator and is NOT point-identified under M-O
#' confounding — it is included for parity with the continuous
#' \code{\link{fit_iv2sls_mediation}}.
#'
#' @param time Numeric follow-up time vector (length n).
#' @param event Numeric 0/1 event indicator (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param g Numeric instrument for X (length n).
#' @param w Numeric NC vector (length n) or matrix (n x q).
#' @param covars Optional data frame of covariates (n rows).
#' @param min_f Minimum partial F for the excluded instrument. Default 10.
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau RMST horizon (rmst only).
#'
#' @return Named list (same fields as \code{\link{fit_unadj_mediation_surv}}).
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
#' fit_iv2sls_mediation_surv(dat$surv_time, dat$surv_event, dat$X, dat$M,
#' dat$G[, 1], dat$W[, 1])
#' }
fit_iv2sls_mediation_surv <- function(time, event, X, M, g, w,
                                      covars = NULL, min_f = 10,
                                      effect_scale = c("loghr", "rmst"),
                                      tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we <- .expand_w(w)

  # Stage 1 (OLS): X ~ g (+ W) + covars -> X_hat
  d_fs <- data.frame(X = X, g = g)
  d_fs <- .bind_covars(d_fs, we$df)
  d_fs <- .bind_covars(d_fs, covars)
  fs <- tryCatch(lm(as.formula(paste0("X ~ g", .plus_frag(we$frag), cs)),
                    data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_res)
  Fst <- .partial_F(fs, "g")
  if (is.na(Fst) || Fst < min_f) return(NA_res)
  X_hat <- fitted(fs)

  # Stage 2 (OLS): M ~ X_hat -> alpha_M
  d1 <- .bind_covars(data.frame(M = M, X_hat = X_hat), covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ X_hat", cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"X_hat" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["X_hat"])
  alpha_se <- as.numeric(s1["X_hat", 2])

  # Stage 3 (Cox / RMST): outcome ~ X_hat + M (+ W) + covars
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d_os <- data.frame(X_hat = X_hat, M = M)
  d_os <- .bind_covars(d_os, we$df)
  d_os <- .bind_covars(d_os, covars)
  d_os$y <- resp
  fml_os <- as.formula(paste0("y ~ X_hat + M", .plus_frag(we$frag), cs))
  fit_os <- if (effect_scale == "loghr") {
    tryCatch(survival::coxph(fml_os, data = d_os), error = function(e) NULL)
  } else {
    tryCatch(lm(fml_os, data = d_os), error = function(e) NULL)
  }
  if (is.null(fit_os)) return(NA_res)
  sm <- summary(fit_os)$coefficients
  pcol <- if ("Pr(>|z|)" %in% colnames(sm)) "Pr(>|z|)" else "Pr(>|t|)"
  if (!"X_hat" %in% rownames(sm) || !"M" %in% rownames(sm)) return(NA_res)
  beta_X <- as.numeric(coef(fit_os)["X_hat"])
  beta_X_se <- as.numeric(sm["X_hat", 2])
  beta_M <- as.numeric(coef(fit_os)["M"])
  beta_M_se <- as.numeric(sm["M", 2])

  NIE <- alpha * beta_M
  NIE_se <- delta_se_product(alpha, alpha_se, beta_M, beta_M_se)

  list(
    NDE = beta_X, NDE_se = beta_X_se, NDE_p = as.numeric(sm["X_hat", pcol]),
    NIE = NIE, NIE_se = NIE_se, NIE_p = 2 * pnorm(-abs(NIE / NIE_se)),
    alpha_M = alpha, alpha_se = alpha_se, beta_M = beta_M, beta_M_se = beta_M_se
  )
}


# ── 4. IV2SLS2 mediation (survival, 2SPS) ──

#' IV2SLS2 survival mediation estimator: 2-stage MR with Cox / RMST outcome
#'
#' Two-stage Mendelian-randomization mediation with a survival outcome
#' stage (2SPS), with optional \strong{path-specific} negative-control (NC)
#' augmentation:
#' \enumerate{
#' \item OLS: \code{X ~ g (+ W1) + covars} -> X_hat (purge U1 from X).
#' \item OLS: \code{M ~ X_hat + gm (+ W2) + covars} -> M_hat, alpha_M.
#' \item Cox / RMST: \code{Surv(t,e) ~ X_hat + M_hat (+ W2) + covars} ->
#' NDE (coef on X_hat), beta_M (coef on M_hat).
#' }
#' \code{NIE = alpha_M * beta_M}. Weak-instrument gates (partial F for
#' G and Gm) apply to the OLS first stages.
#'
#' \strong{Path-specific NC augmentation (optional).} \code{W1} proxies the
#' exposure-mediator confounder (U1, X->M path) and is added to stage 1
#' only; \code{W2} proxies the mediator-outcome confounder (U2, M->Y path)
#' and is added to stages 2 and 3. Either panel may be omitted; with both
#' \code{NULL} the estimator reduces to plain two-instrument 2-stage MR.
#' Conditioning on a \emph{pooled} panel in all three stages is a collider
#' under multi-confounder designs and is not supported: identical
#' \code{W1}/\code{W2} are treated as absent (pure MR). When \code{W1} and
#' \code{W2} are distinct-noise proxies of the \emph{same} latent composite
#' (single-confounder design), their column spaces are near-collinear and
#' the estimator likewise falls back to plain two-instrument 2-stage MR.
#'
#' @param time Numeric follow-up time vector (length n).
#' @param event Numeric 0/1 event indicator (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param g Numeric instrument for X (length n).
#' @param gm Numeric instrument for M (length n).
#' @param covars Optional data frame of covariates (n rows).
#' @param min_f Minimum partial F for each excluded instrument. Default 10.
#' @param W1 Optional NC panel (vector length n or matrix n x q) proxying
#' the exposure-mediator confounder (X->M path); added to stage 1 only.
#' Default \code{NULL}.
#' @param W2 Optional NC panel (vector length n or matrix n x q) proxying
#' the mediator-outcome confounder (M->Y path); added to stages 2 and 3.
#' Default \code{NULL}.
#' @param w Defunct. The pooled single-panel argument was removed in
#' v0.9.9 (collider under multi-confounder designs). Use \code{W1} and/or
#' \code{W2} instead.
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau RMST horizon (rmst only).
#'
#' @return Named list (same fields as \code{\link{fit_unadj_mediation_surv}}).
#' Returns all-NA if either first-stage partial F is below \code{min_f}.
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
#' mo_confounding = 0.8, phi = 0.8, lambda_XM = c(1, 0),
#' lambda_MY = c(0, 1), omega_1 = 0.7, omega_2 = 0.7, seed = 1)
#' fit_iv2sls_mediation2_surv(dat$surv_time, dat$surv_event, dat$X, dat$M,
#' dat$G[, 1], dat$Gm, W1 = dat$W1, W2 = dat$W2)
#' }
fit_iv2sls_mediation2_surv <- function(time, event, X, M, g, gm,
                                       covars = NULL, min_f = 10,
                                       W1 = NULL, W2 = NULL, w = NULL,
                                       effect_scale = c("loghr", "rmst"),
                                       tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()

  # Defunct-argument trap: the pooled single-panel `w` was removed in v0.9.9.
  if (!is.null(w))
    stop("argument `w` was removed in v0.9.9. IV2SLS2 now takes optional ",
         "path-specific negative-control panels `W1` (stage 1, X->M path) ",
         "and `W2` (stages 2-3, M->Y path); a pooled panel conditioned on in ",
         "all three stages is a collider under multi-confounder designs. ",
         "Use `W1 = ...` and/or `W2 = ...`, or omit both for plain 2-stage MR.",
         call. = FALSE)

  # Pooled guard: identical W1 and W2 panels are a pooled panel in disguise
  # (a common child of the two confounders). Treat both as absent -> pure MR.
  if (!is.null(W1) && !is.null(W2) && identical(W1, W2)) {
    W1 <- NULL
    W2 <- NULL
  }

  # Shared-composite guard (k=1 fallback): distinct-noise panels of the SAME
  # latent composite (single-confounder design) are near-collinear; only one
  # confounder exists, so path-specific augmentation is undefined. Fall back to
  # plain 2-stage MR (drop both panels), as for an identical pooled panel.
  if (!is.null(W1) && !is.null(W2) && .w_panels_collinear(W1, W2)) {
    W1 <- NULL
    W2 <- NULL
  }

  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)
  we1 <- .expand_w(W1)   # stage 1 panel (X->M path)
  we2 <- .expand_w(W2)   # stages 2-3 panel (M->Y path)

  # Stage 1 (OLS): X ~ g (+ W1) + covars -> X_hat
  d_fs <- data.frame(X = X, g = g)
  d_fs <- .bind_covars(d_fs, we1$df)
  d_fs <- .bind_covars(d_fs, covars)
  fs <- tryCatch(lm(as.formula(paste0("X ~ g", .plus_frag(we1$frag), cs)),
                    data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_res)
  Fst_g <- .partial_F(fs, "g")
  if (is.na(Fst_g) || Fst_g < min_f) return(NA_res)
  X_hat <- fitted(fs)

  # Stage 2 (OLS): M ~ X_hat + gm (+ W2) + covars -> M_hat, alpha_M
  d_ms <- data.frame(M = M, X_hat = X_hat, gm = gm)
  d_ms <- .bind_covars(d_ms, we2$df)
  d_ms <- .bind_covars(d_ms, covars)
  ms <- tryCatch(lm(as.formula(paste0("M ~ X_hat + gm", .plus_frag(we2$frag), cs)),
                    data = d_ms), error = function(e) NULL)
  if (is.null(ms)) return(NA_res)
  Fst_gm <- .partial_F(ms, "gm")
  if (is.na(Fst_gm) || Fst_gm < min_f) return(NA_res)
  s_ms <- summary(ms)$coefficients
  if (!"X_hat" %in% rownames(s_ms)) return(NA_res)
  alpha <- as.numeric(coef(ms)["X_hat"])
  alpha_se <- as.numeric(s_ms["X_hat", 2])
  M_hat <- fitted(ms)

  # Stage 3 (Cox / RMST): outcome ~ X_hat + M_hat (+ W2) + covars
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d_os <- data.frame(X_hat = X_hat, M_hat = M_hat)
  d_os <- .bind_covars(d_os, we2$df)
  d_os <- .bind_covars(d_os, covars)
  fml_os <- as.formula(paste0("y ~ X_hat + M_hat", .plus_frag(we2$frag), cs))
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
#' both the X->M and M->Y confounding paths. This is the survival
#' analogue of \code{\link{fit_pgc_mediation}} (continuous outcome):
#' \enumerate{
#' \item OLS: residualise X on G -> X_resid.
#' \item OLS: bridge X_resid on the FULL W matrix -> W_hat (proxy for U).
#' \item OLS: \code{M ~ X + W_hat + covars} -> alpha_M.
#' \item Cox / RMST: \code{Surv(t,e) ~ X + M + W_hat + covars} ->
#' NDE (coef on X), beta_M (coef on M).
#' }
#' \code{NIE = alpha_M * beta_M}. Unlike
#' \code{\link{fit_pgc_mediation2_surv}}, which uses path-specific W1/W2
#' bridges, this estimator uses a single combined W panel and is
#' appropriate when separate conf_XM / conf_MY confounders are not assumed.
#'
#' @param time Numeric follow-up time vector (length n).
#' @param event Numeric 0/1 event indicator (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param g Numeric instrument for X (length n).
#' @param W Numeric NC matrix (n x q) or vector (length n).
#' @param covars Optional data frame of covariates (n rows).
#' @param min_f Minimum partial F for G. Default 10.
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau RMST horizon (rmst only).
#'
#' @return Named list (same fields as \code{\link{fit_unadj_mediation_surv}}).
#' Returns all-NA if the first-stage partial F for G is below \code{min_f}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
#' fit_pgc_mediation_surv(dat$surv_time, dat$surv_event, dat$X, dat$M,
#' dat$G[, 1], dat$W)
#' }
fit_pgc_mediation_surv <- function(time, event, X, M, g, W, covars = NULL,
                                   min_f = 10,
                                   effect_scale = c("loghr", "rmst"),
                                   tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  if (!is.matrix(W)) W <- as.matrix(W)

  # Weak-instrument check: partial F for G in X ~ G + W + covars
  d_fs <- .bind_covars(data.frame(X = X, g = g), covars)
  d_fs <- cbind(d_fs, as.data.frame(W))
  w_fs_names <- paste0("Wfs_", seq_len(ncol(W)))
  names(d_fs)[(ncol(d_fs) - ncol(W) + 1):ncol(d_fs)] <- w_fs_names
  fs_fml <- as.formula(paste0("X ~ g + ",
                              paste(w_fs_names, collapse = " + "), cs))
  fs <- tryCatch(lm(fs_fml, data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_res)
  Fst_g <- .partial_F(fs, "g")
  if (is.na(Fst_g) || Fst_g < min_f) return(NA_res)

  # Step 1: residualise X on G -> X_resid (OLS)
  d_r <- .bind_covars(data.frame(Xc = X, g = g), covars)
  fit_resid <- tryCatch(lm(as.formula(paste0("Xc ~ g", cs)), data = d_r),
                        error = function(e) NULL)
  if (is.null(fit_resid)) return(NA_res)
  X_resid <- residuals(fit_resid)

  # Step 2: bridge X_resid on the FULL W matrix -> W_hat (OLS)
  d_b <- data.frame(X_resid = X_resid)
  d_b <- cbind(d_b, as.data.frame(W))
  if (!is.null(covars)) d_b <- cbind(d_b, covars)
  w_names <- paste0("W", seq_len(ncol(W)))
  names(d_b)[2:(ncol(W) + 1)] <- w_names
  fml_b <- as.formula(paste0("X_resid ~ ", paste(w_names, collapse = " + "), cs))
  fit_b <- tryCatch(lm(fml_b, data = d_b), error = function(e) NULL)
  if (is.null(fit_b)) return(NA_res)
  W_hat <- fitted(fit_b)

  # Stage 1 (OLS): M ~ X + W_hat + covars -> alpha_M
  d1 <- .bind_covars(data.frame(M = M, X = X, W_hat = W_hat), covars)
  fit1 <- tryCatch(lm(as.formula(paste0("M ~ X + W_hat", cs)), data = d1),
                   error = function(e) NULL)
  if (is.null(fit1)) return(NA_res)
  s1 <- summary(fit1)$coefficients
  if (!"X" %in% rownames(s1)) return(NA_res)
  alpha <- as.numeric(coef(fit1)["X"])
  alpha_se <- as.numeric(s1["X", 2])

  # Stage 2 (Cox / RMST): outcome ~ X + M + W_hat + covars
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d_os <- .bind_covars(data.frame(X = X, M = M, W_hat = W_hat), covars)
  fml_os <- as.formula(paste0("y ~ X + M + W_hat", cs))
  os <- .fit_surv_outcome_stage(resp, d_os, fml_os, effect_scale,
                                nde_name = "X", med_name = "M")
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
#' a survival outcome stage. Stages 1-2 (X bridge on W1, M bridge on W2,
#' X_hat / M_hat construction) are identical to \code{\link{fit_pgc_mediation2}}
#' and remain OLS. Only stage 3 switches to Cox (log-HR) or RMST
#' pseudo-observation OLS:
#' \code{Surv(t,e) ~ X_hat + M_hat + W_hat_X + W_hat_M + covars}.
#'
#' When \code{gm = NULL} (PGC2), stage 2 uses pure NC identification.
#' When \code{gm} is supplied (PGC2Gm), the mediator instrument helps
#' isolate conf_MY before bridging W2.
#'
#' @param time Numeric follow-up time vector (length n).
#' @param event Numeric 0/1 event indicator (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param g Numeric instrument for X (length n).
#' @param W1 Numeric NC matrix (n x q) or vector for the X->M path.
#' @param W2 Numeric NC matrix (n x q) or vector for the M->Y path.
#' @param gm Optional numeric mediator instrument (length n).
#' @param covars Optional data frame of covariates (n rows).
#' @param min_f Minimum partial F for G1. Default 10.
#' @param effect_scale Character: \code{"loghr"} or \code{"rmst"}.
#' @param tau RMST horizon (rmst only).
#'
#' @return Named list (same fields as \code{\link{fit_unadj_mediation_surv}}).
#' Returns all-NA if the first-stage partial F for G1 is below \code{min_f}.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 500, outcome_type = "survival",
#' mo_confounding = 0.8, rho_G2 = 0.3,
#' lambda_XM = c(1, 0), lambda_MY = c(0, 1), seed = 1)
#' fit_pgc_mediation2_surv(dat$surv_time, dat$surv_event, dat$X, dat$M,
#' dat$G[, 1], dat$W1, dat$W2, gm = dat$Gm)
#' }
fit_pgc_mediation2_surv <- function(time, event, X, M, g, W1, W2, gm = NULL,
                                    covars = NULL, min_f = 10,
                                    effect_scale = c("loghr", "rmst"),
                                    tau = NULL) {
  effect_scale <- match.arg(effect_scale)
  NA_res <- .surv_med_NA()
  cnames <- if (!is.null(covars)) names(covars) else character(0)
  cs <- .covar_str(cnames)

  if (!is.matrix(W1)) W1 <- as.matrix(W1)
  if (!is.matrix(W2)) W2 <- as.matrix(W2)

  # ── Stage 1: Bridge for X (purge conf_XM) — OLS, identical to continuous ──
  d_fs <- .bind_covars(data.frame(X = X, g = g), covars)
  d_fs <- cbind(d_fs, as.data.frame(W1))
  w1_names <- paste0("W1_", seq_len(ncol(W1)))
  names(d_fs)[(ncol(d_fs) - ncol(W1) + 1):ncol(d_fs)] <- w1_names
  fs_fml <- as.formula(paste0("X ~ g + ",
                              paste(w1_names, collapse = " + "), cs))
  fs <- tryCatch(lm(fs_fml, data = d_fs), error = function(e) NULL)
  if (is.null(fs)) return(NA_res)
  Fst_g <- .partial_F(fs, "g")
  if (is.na(Fst_g) || Fst_g < min_f) return(NA_res)

  d_r <- .bind_covars(data.frame(Xc = X, g = g), covars)
  fit_resid <- tryCatch(lm(as.formula(paste0("Xc ~ g", cs)), data = d_r),
                        error = function(e) NULL)
  if (is.null(fit_resid)) return(NA_res)
  X_resid <- residuals(fit_resid)

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

  d_zh <- .bind_covars(data.frame(X = X, g = g, W_hat_X = W_hat_X), covars)
  fit_zh <- tryCatch(lm(as.formula(paste0("X ~ g + W_hat_X", cs)),
                        data = d_zh), error = function(e) NULL)
  if (is.null(fit_zh)) return(NA_res)
  X_hat <- fitted(fit_zh)

  # ── Stage 2: Bridge for M (purge conf_MY) — OLS, identical to continuous ──
  if (is.null(gm)) {
    d_mr <- .bind_covars(data.frame(M = M, X_hat = X_hat), covars)
    fit_mr <- tryCatch(lm(as.formula(paste0("M ~ X_hat", cs)),
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
    d_mh <- .bind_covars(data.frame(M = M, X_hat = X_hat,
                                    W_hat_M = W_hat_M), covars)
    fit_mh <- tryCatch(lm(as.formula(paste0("M ~ X_hat + W_hat_M", cs)),
                          data = d_mh), error = function(e) NULL)
  } else {
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

  # ── Stage 3 (Cox / RMST): outcome ~ X_hat + M_hat + W_hat_X + W_hat_M ──
  resp <- .make_surv_response(time, event, effect_scale, tau)
  d_os <- data.frame(X_hat = X_hat, M_hat = M_hat,
                     W_hat_X = W_hat_X, W_hat_M = W_hat_M)
  d_os <- .bind_covars(d_os, covars)
  fml_os <- as.formula(paste0("y ~ X_hat + M_hat + W_hat_X + W_hat_M", cs))
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
#' both stages (\code{W ~ M + X}, \code{W ~ y + X + M}), placing the
#' outcome on the right-hand side. This is structurally impossible when
#' the outcome is a \code{survival::Surv} object. COCA mediation is
#' therefore unsupported for survival outcomes and always returns NA.
#'
#' @param time Numeric follow-up time vector (length n).
#' @param event Numeric 0/1 event indicator (length n).
#' @param X Numeric exposure vector (length n).
#' @param M Numeric mediator vector (length n).
#' @param w Numeric NC vector (length n).
#' @param covars Optional data frame of covariates (n rows).
#' @param ... Ignored (signature compatibility).
#'
#' @return All-NA list with an informative \code{"reason"} attribute.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- iconic:::generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
#' fit_coca_mediation_surv(dat$surv_time, dat$surv_event, dat$X, dat$M, dat$W[, 1])
#' # all NA (COCA mediation unsupported for survival outcomes)
#' }
fit_coca_mediation_surv <- function(time, event, X, M, w, covars = NULL, ...) {
  out <- .surv_med_NA()
  attr(out, "reason") <- paste(
    "COCA mediation regresses W on Y (W ~ y + X + M), placing the outcome",
    "on the RHS, which is impossible with a Surv object. COCA mediation",
    "is unsupported for survival outcomes."
  )
  out
}
