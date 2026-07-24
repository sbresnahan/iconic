# ============================================================
# iconic_estimate: Real-data estimation layer for the v0.6.0
# model selection workflow.
#
# Fits all eligible estimators on the user's real data (supplied
# as an iconic_data object) and returns per-feature (and
# per-mediator) point estimates, SEs, and p-values.
#
# This is the real-data analogue of analyze_methods_robust() /
# analyze_mediation_robust(), which operate on the internal `dat`
# list from run_single_iteration().  Both share the same
# format-agnostic per-feature estimation functions
# (.estimate_total_feature / .estimate_mediation_feature).
# ============================================================

#' Fit all eligible estimators on real data
#'
#' Applies all eligible causal estimators to the user's real data
#' (supplied as an \code{\link{iconic_data}} object) and returns tidy
#' per-feature (and per-mediator) point estimates, standard errors, and
#' p-values.
#'
#' If a \code{diagnosis} from \code{\link{iconic_diagnose}()} is supplied,
#' only eligible estimators are run by default.  If \code{methods} is
#' supplied, only those methods are run (user override).  If neither is
#' supplied, all estimators whose required inputs are present are run.
#'
#' @param data      An \code{iconic_data} object from \code{\link{iconic_data}()}.
#' @param methods   Optional character vector of estimator names to run
#'   (e.g. \code{c("IV2SLS", "PGC")}).  Default \code{NULL}: run all
#'   eligible.
#' @param diagnosis Optional \code{iconic_diagnosis} object from
#'   \code{\link{iconic_diagnose}()}.  When supplied, only eligible
#'   estimators are run (unless \code{methods} or \code{run_all} overrides).
#' @param alpha     Significance threshold for significance flags. Default 0.05.
#' @param n_cores   Number of parallel workers for per-feature (and
#'   per-mediator) estimation.  Default 1 (sequential).  Uses
#'   \code{parallel::mclapply} on Unix and a PSOCK cluster on Windows.
#' @param min_f     Minimum partial F for the per-transcript weak-instrument
#'   gate inside IV2SLS, IV2SLS2, and PGC2Gm.  Default \code{NULL}: inherits
#'   from \code{diagnosis$min_f} when a diagnosis is supplied, otherwise 10.
#'   Pass an explicit value to override.
#' @param run_all   Logical.  Default \code{FALSE}.  When \code{TRUE},
#'   overrides diagnosis-based eligibility and runs every method whose
#'   required data exists (via \code{.auto_eligible_methods()}).  The
#'   per-transcript \code{min_f} gate still applies.  This is the "force
#'   run" escape hatch for exploratory analysis.
#' @param se_method Character (v0.9.2, JYH #867): \code{"delta"} (default),
#'   \code{"bootstrap"}, or \code{"composite"} (v0.9.3). When
#'   \code{"bootstrap"}, mediation NDE_se/NIE_se are replaced by the SD of
#'   \code{n_boot} nonparametric bootstrap resamples. When \code{"composite"},
#'   mediation NIE_p is replaced by the Huang (2019) JT-comp composite null
#'   p-value, which accounts for the three-case structure of H0: alpha*beta=0
#'   and provides higher power than the Sobel/Wald test when signals are
#'   sparse. Only applies in mediation mode. NDE_p and NIE_se are unchanged.
#' @param n_boot    Integer (v0.9.2): number of bootstrap resamples when
#'   \code{se_method = "bootstrap"}. Default 500.
#' @param effect_scale Character (v0.9.4): \code{"loghr"} (default) or
#'   \code{"rmst"}.  Only used when \code{data$outcome_type = "survival"}.
#'   \code{"loghr"} fits Cox proportional-hazards models and reports
#'   log-hazard ratios.  \code{"rmst"} regresses leave-one-out RMST
#'   pseudo-observations (Graw et al. 2009) via OLS, reporting effects on
#'   the restricted-mean-survival-time (time) scale — a collapsible
#'   alternative where the NDE/NIE product decomposition is exact.  Ignored
#'   (with a message) when \code{outcome_type = "continuous"}.
#' @param tau       Numeric (v0.9.4): RMST restriction time horizon.  Default
#'   \code{NULL} (90th percentile of follow-up).  Used only when
#'   \code{effect_scale = "rmst"}.
#'
#' @return A data frame.  For total-effect mode: \code{feature}, \code{method},
#'   \code{beta}, \code{se}, \code{pvalue}, \code{significant}.  For mediation
#'   mode: \code{feature}, \code{mediator}, \code{method}, \code{NDE},
#'   \code{NDE_se}, \code{NDE_p}, \code{NIE}, \code{NIE_se}, \code{NIE_p},
#'   \code{NDE_significant}, \code{NIE_significant}.  When
#'   \code{outcome_type = "survival"}, estimates are on the log-HR scale
#'   (\code{effect_scale = "loghr"}) or the RMST/time scale
#'   (\code{effect_scale = "rmst"}).
#' @export
#'
#' @examples
#' \dontrun{
#' data <- iconic_data(Z = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
#'                     G = rnorm(100), W = matrix(rnorm(100 * 10), 10, 100))
#' diag <- iconic_diagnose(data)
#' est <- iconic_estimate(data, diagnosis = diag, n_cores = 4)
#' head(est)
#'
#' # Survival outcome (v0.9.4)
#' sdat <- iconic_data(Z = rnorm(100), outcome_type = "survival",
#'                     surv_time = rexp(100), surv_event = rbinom(100, 1, 0.6),
#'                     G = rnorm(100), W = matrix(rnorm(100 * 10), 10, 100))
#' est <- iconic_estimate(sdat, effect_scale = "loghr")
#' est_rmst <- iconic_estimate(sdat, effect_scale = "rmst")
#' }
iconic_estimate <- function(data, methods = NULL, diagnosis = NULL,
                            alpha = 0.05, n_cores = 1, min_f = NULL,
                            run_all = FALSE,
                            se_method = c("delta", "bootstrap", "composite"),
                            n_boot = 500,
                            effect_scale = c("loghr", "rmst"),
                            tau = NULL) {
  se_method <- match.arg(se_method)
  effect_scale <- match.arg(effect_scale)
  if (!inherits(data, "iconic_data"))
    stop("data must be an iconic_data object from iconic_data().")

  # effect_scale only applies to survival outcomes
  if (data$outcome_type != "survival" && effect_scale == "rmst")
    message("effect_scale = 'rmst' is ignored for continuous outcomes.")

  # Inherit min_f from diagnosis when not explicitly set; fall back to 10.
  if (is.null(min_f)) {
    min_f <- if (!is.null(diagnosis) && !is.null(diagnosis$min_f))
      diagnosis$min_f else 10
  }

  # Determine which methods to run
  # run_all overrides diagnosis eligibility (runs every method whose data
  # exists) but the per-transcript min_f gate still applies.
  if (run_all) {
    methods_to_run <- .auto_eligible_methods(data)
  } else if (!is.null(methods)) {
    methods_to_run <- methods
  } else if (!is.null(diagnosis)) {
    elig <- diagnosis$eligibility
    methods_to_run <- elig$estimator[elig$eligible]
  } else {
    # Auto-eligibility: run all whose inputs are present
    methods_to_run <- .auto_eligible_methods(data)
  }

  if (!length(methods_to_run))
    stop("No eligible estimators for this data. ",
         "Supply instruments (G) or negative controls (W) ",
         "or specify methods explicitly.")

  # v0.9.4: branch on outcome_type.  Survival outcomes use dedicated
  # survival drivers that call the fit_*_surv / fit_*_mediation_surv
  # estimators.  Continuous outcomes use the existing drivers unchanged.
  if (data$outcome_type == "survival") {
    if (data$is_mediation) {
      .estimate_mediation_surv_driver(data, methods_to_run, alpha, n_cores,
                                      min_f, se_method, n_boot,
                                      effect_scale, tau)
    } else {
      .estimate_total_surv_driver(data, methods_to_run, alpha, n_cores,
                                  min_f, effect_scale, tau)
    }
  } else if (data$is_mediation) {
    .estimate_mediation_driver(data, methods_to_run, alpha, n_cores, min_f,
                               se_method = se_method, n_boot = n_boot)
  } else {
    .estimate_total_driver(data, methods_to_run, alpha, n_cores, min_f)
  }
}


#' Determine eligible methods from data alone (internal)
#'
#' Lightweight auto-eligibility when no diagnosis is supplied.
#' @param data iconic_data object.
#' @return Character vector of method names.
#' @keywords internal
.auto_eligible_methods <- function(data) {
  methods <- "UNADJ"
  if (data$has_instrument && data$has_nc)
    methods <- c(methods, "DIRECT", "IV2SLS")
  if (data$has_nc)
    methods <- c(methods, "COCA")
  if (data$has_instrument && data$has_nc)
    methods <- c(methods, "PGC")
  if (data$is_mediation && data$has_instrument &&
      data$has_mediator_instrument && data$has_nc)
    methods <- c(methods, "IV2SLS2")
  if (data$is_mediation && data$has_instrument && data$has_path_nc)
    methods <- c(methods, "PGC2")
  if (data$is_mediation && data$has_instrument &&
      data$has_path_nc && data$has_mediator_instrument)
    methods <- c(methods, "PGC2Gm")
  methods
}


#' Total-effect estimation driver for iconic_data (internal)
#'
#' Loops over outcome features, calling .estimate_total_feature().
#' @keywords internal
.estimate_total_driver <- function(data, methods, alpha, n_cores = 1,
                                   min_f = 10) {
  n <- data$n
  nf <- data$n_features
  Z <- data$Z
  cv <- data$covariates
  G <- data$G
  W <- data$W
  # iconic_data stores W as q x n (features x samples); estimation functions
  # expect n x q (samples x NC features), matching the simulation pipeline.
  W_mat <- if (!is.null(W)) t(W) else NULL
  W_avg <- if (!is.null(W)) colMeans(W) else NULL  # per-sample averages

  results <- .parallel_lapply(seq_len(nf), function(f) {
    y <- data$Y[f, ]  # features x samples -> row f is feature f
    # v0.8.4: pass the full W panel (n x q) to all estimators.
    # DIRECT, IV2SLS, and IV2SLS2 now use all NC features as covariates.
    w <- W_mat  # full n x q matrix (NULL when no W)
    .estimate_total_feature(
      Z = Z, y = y, g = G, w = w, W_mat = W_mat,
      W_avg = W_avg, covars = cv,
      methods = methods, feature_idx = f, min_f = min_f
    )
  }, n_cores = n_cores, progress = "Estimating features")

  # mclapply swallows errors as try-error objects; surface the first one.
  err_idx <- which(vapply(results, inherits, logical(1), "try-error"))
  if (length(err_idx)) {
    stop("Estimation failed for task ", err_idx[1], ". Underlying error:\n",
         conditionMessage(attr(results[[err_idx[1]]], "condition")))
  }

  out <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(out)) return(out)

  out$significant <- out$pvalue < alpha
  # Map feature indices to names
  out$feature <- data$feature_names[out$feature]
  out
}


#' Mediation estimation driver for iconic_data (internal)
#'
#' Loops over mediators x outcome features, calling
#' .estimate_mediation_feature().
#' @keywords internal
.estimate_mediation_driver <- function(data, methods, alpha, n_cores = 1,
                                       min_f = 10, se_method = "delta",
                                       n_boot = 500) {
  n <- data$n
  nf <- data$n_features
  nm <- data$n_mediators
  Z <- data$Z
  cv <- data$covariates
  G <- data$G
  W <- data$W
  # Transpose from q x n to n x q (see .estimate_total_driver)
  W_mat <- if (!is.null(W)) t(W) else NULL
  W_avg <- if (!is.null(W)) colMeans(W) else NULL
  W1_mat <- if (!is.null(data$W1)) t(data$W1) else NULL
  W2_mat <- if (!is.null(data$W2)) t(data$W2) else NULL

  # Flatten the mediator x feature grid into a single task list so it
  # can be parallelised in one .parallel_lapply() call.
  pairs <- expand.grid(m = seq_len(nm), f = seq_len(nf),
                        KEEP.OUT.ATTRS = FALSE)

  results <- .parallel_lapply(seq_len(nrow(pairs)), function(idx) {
    m <- pairs$m[idx]
    f <- pairs$f[idx]
    M_vec <- data$M[m, ]
    gm <- if (!is.null(data$Gm)) {
      if (is.matrix(data$Gm) && nrow(data$Gm) == nm) data$Gm[m, ]
      else if (is.matrix(data$Gm) && ncol(data$Gm) == 1) rep(data$Gm[1, ], n)
      else if (is.matrix(data$Gm)) data$Gm[min(m, nrow(data$Gm)), ]  # clamp
      else data$Gm  # scalar vector
    } else NULL

    y <- data$Y[f, ]
    # v0.8.4: pass the full W panel (n x q) to all estimators.
    w <- W_mat  # full n x q matrix (NULL when no W)

    res <- .estimate_mediation_feature(
      Z = Z, y = y, M_vec = M_vec, g = G, gm = gm,
      w = w, W_mat = W_mat, W1_mat = W1_mat, W2_mat = W2_mat,
      W_avg = W_avg, covars = cv,
      methods = methods, feature_idx = f, min_f = min_f,
      se_method = se_method, n_boot = n_boot
    )
    if (!is.null(res)) {
      res$mediator <- m
      res
    } else {
      NULL
    }
  }, n_cores = n_cores, progress = "Estimating mediation effects")

  # mclapply swallows errors as try-error objects; surface the first one
  # instead of producing a confusing downstream error.
  err_idx <- which(vapply(results, inherits, logical(1), "try-error"))
  if (length(err_idx)) {
    stop("Estimation failed for task ", err_idx[1], ". Underlying error:\n",
         conditionMessage(attr(results[[err_idx[1]]], "condition")))
  }

  results <- Filter(Negate(is.null), results)

  out <- do.call(rbind, results)
  if (is.null(out)) return(out)

  # v0.9.3: when se_method == "composite", replace NIE_p with the
  # Huang (2019) JT-comp composite null p-value.  This is a two-pass
  # operation: the per-feature fits above already collected the
  # z-statistics (alpha_M/alpha_se, beta_M/beta_M_se); here we estimate
  # Var(a)/Var(b) per method (and per mediator) across features and
  # compute the composite p-value.
  if (se_method == "composite")
    out <- .apply_composite_pvalues(out)

  out$NDE_significant <- out$NDE_p < alpha
  out$NIE_significant <- out$NIE_p < alpha

  # Map indices to names
  out$feature <- data$feature_names[out$feature]
  out$mediator <- data$mediator_names[out$mediator]

  # Reorder columns
  front <- c("feature", "mediator", "method")
  out <- out[, c(front, setdiff(names(out), front))]
  out
}


# ============================================================
# v0.9.4: Survival outcome estimation drivers
#
# These mirror .estimate_total_driver / .estimate_mediation_driver but
# call the fit_*_surv / fit_*_mediation_surv estimators and pass
# surv_time / surv_event / effect_scale / tau.  The per-feature loop
# runs once (n_features = 1 for survival); the per-mediator loop runs
# over all mediators.
# ============================================================

#' Total-effect survival estimation driver (internal, v0.9.4)
#'
#' Runs all eligible survival estimators on a single time-to-event
#' outcome.  Returns a data frame with the same columns as
#' \code{.estimate_total_driver}.
#' @keywords internal
.estimate_total_surv_driver <- function(data, methods, alpha, n_cores = 1,
                                        min_f = 10, effect_scale = "loghr",
                                        tau = NULL) {
  time <- data$surv_time
  event <- data$surv_event
  Z <- data$Z
  cv <- data$covariates
  G <- data$G
  W <- data$W
  W_mat <- if (!is.null(W)) t(W) else NULL  # n x q

  run_one <- function(method) {
    r <- switch(method,
      UNADJ  = fit_unadj_surv(time, event, Z, cv, effect_scale, tau),
      DIRECT = fit_direct_surv(time, event, Z, G, W_mat, cv, effect_scale, tau),
      COCA   = fit_coca_surv(time, event, Z, if (!is.null(W_mat)) W_mat[, 1] else NULL, cv),
      IV2SLS = fit_iv2sls_surv(time, event, Z, G, W_mat, cv, min_f, effect_scale, tau),
      PGC    = fit_pgc_surv(time, event, Z, G, W_mat, cv, effect_scale, tau),
      NULL
    )
    if (is.null(r)) return(NULL)
    data.frame(feature = 1L, method = method,
               beta = r$beta, se = r$se, pvalue = r$pvalue,
               stringsAsFactors = FALSE)
  }

  results <- lapply(methods, run_one)
  out <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(out)) return(out)

  out$significant <- out$pvalue < alpha
  out$feature <- data$feature_names[out$feature]
  out
}


#' Mediation survival estimation driver (internal, v0.9.4)
#'
#' Loops over mediators, calling the survival mediation estimators for
#' each.  Returns a data frame with the same columns as
#' \code{.estimate_mediation_driver}.
#' @keywords internal
.estimate_mediation_surv_driver <- function(data, methods, alpha,
                                            n_cores = 1, min_f = 10,
                                            se_method = "delta",
                                            n_boot = 500,
                                            effect_scale = "loghr",
                                            tau = NULL) {
  time <- data$surv_time
  event <- data$surv_event
  n <- data$n
  nm <- data$n_mediators
  Z <- data$Z
  cv <- data$covariates
  G <- data$G
  W <- data$W
  W_mat <- if (!is.null(W)) t(W) else NULL
  W1_mat <- if (!is.null(data$W1)) t(data$W1) else NULL
  W2_mat <- if (!is.null(data$W2)) t(data$W2) else NULL

  run_one <- function(m, method) {
    M_vec <- data$M[m, ]
    gm <- if (!is.null(data$Gm)) {
      if (is.matrix(data$Gm) && nrow(data$Gm) == nm) data$Gm[m, ]
      else if (is.matrix(data$Gm) && ncol(data$Gm) == 1) rep(data$Gm[1, ], n)
      else if (is.matrix(data$Gm)) data$Gm[min(m, nrow(data$Gm)), ]
      else data$Gm
    } else NULL

    w <- W_mat
    w1 <- W1_mat
    w2 <- W2_mat

    r <- switch(as.character(method),
      UNADJ   = fit_unadj_mediation_surv(time, event, Z, M_vec, cv,
                                         effect_scale, tau),
      DIRECT  = fit_direct_mediation_surv(time, event, Z, M_vec, G, w, cv,
                                          effect_scale, tau),
      COCA    = fit_coca_mediation_surv(time, event, Z, M_vec,
                                        if (!is.null(w)) w[, 1] else NULL, cv),
      IV2SLS  = fit_iv2sls_mediation_surv(time, event, Z, M_vec, G, w,
                                          cv, min_f, effect_scale, tau),
      IV2SLS2 = fit_iv2sls_mediation2_surv(time, event, Z, M_vec, G, gm, w,
                                           cv, min_f, effect_scale, tau),
      PGC     = fit_pgc_mediation_surv(time, event, Z, M_vec, G, w,
                                       cv, min_f, effect_scale, tau),
      PGC2    = fit_pgc_mediation2_surv(time, event, Z, M_vec, G, w1, w2,
                                        gm = NULL, cv, min_f,
                                        effect_scale, tau),
      PGC2Gm  = fit_pgc_mediation2_surv(time, event, Z, M_vec, G, w1, w2,
                                        gm = gm, cv, min_f,
                                        effect_scale, tau),
      NULL
    )
    if (is.null(r)) return(NULL)
    data.frame(
      feature = 1L, mediator = m, method = method,
      NDE = r$NDE, NDE_se = r$NDE_se, NDE_p = r$NDE_p,
      NIE = r$NIE, NIE_se = r$NIE_se, NIE_p = r$NIE_p,
      alpha_M = r$alpha_M, alpha_se = r$alpha_se,
      beta_M = r$beta_M, beta_M_se = r$beta_M_se,
      stringsAsFactors = FALSE
    )
  }

  # Flatten mediator x method grid
  grid <- expand.grid(m = seq_len(nm), method = methods,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  results <- .parallel_lapply(seq_len(nrow(grid)), function(i) {
    run_one(grid$m[i], grid$method[i])
  }, n_cores = n_cores, progress = "Estimating survival mediation effects")

  err_idx <- which(vapply(results, inherits, logical(1), "try-error"))
  if (length(err_idx)) {
    stop("Estimation failed for task ", err_idx[1], ". Underlying error:\n",
         conditionMessage(attr(results[[err_idx[1]]], "condition")))
  }

  results <- Filter(Negate(is.null), results)
  out <- do.call(rbind, results)
  if (is.null(out)) return(out)

  # Composite p-values (same two-pass logic as continuous driver)
  if (se_method == "composite")
    out <- .apply_composite_pvalues(out)

  out$NDE_significant <- out$NDE_p < alpha
  out$NIE_significant <- out$NIE_p < alpha

  out$feature <- data$feature_names[out$feature]
  out$mediator <- data$mediator_names[out$mediator]

  front <- c("feature", "mediator", "method")
  out <- out[, c(front, setdiff(names(out), front))]
  out
}
