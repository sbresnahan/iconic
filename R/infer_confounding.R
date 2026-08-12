# ============================================================
# infer_confounding: Data-calibrated confounding parameters.
#
# Estimates the held-fixed confounding parameters (delta, mo_confounding,
# omega, k) from the user's data, so that iconic_sensitivity() and
# iconic_prospect() can run on data-calibrated rather than default values.
#
# Inference methods:
# conf_strength (delta) -- UNADJ--IV2SLS gap (requires valid G)
# mo_confounding (delta_mo) -- IV2SLS--IV2SLS2 NIE gap (requires G + Gm)
# omega (omega) -- sqrt(R^2) of W on residualized Y (requires W)
# k (k) -- parallel analysis on residualized Y (requires >=5 features)
#
# Parameters NOT inferred:
# rho_G1, rho_G2 -- unestimable (U unobserved); always remain sweep variables.
#
# Each inferred value carries its method, the assumption required, and
# any warnings about reliability. Parameters that cannot be inferred
# (missing instruments, NCs, or insufficient features) are marked
# unavailable and fall back to defaults.
#
# The circularity caveat: inference uses estimator validity (e.g., IV2SLS
# is unbiased) to calibrate a benchmark whose purpose is to test estimator
# validity. This is documented in the return object and should be
# acknowledged in any analysis using confounding = "inferred".
# ============================================================

#' Infer confounding parameters from the user's data
#'
#' Estimates the held-fixed confounding parameters that
#' \code{\link{iconic_sensitivity}()} and \code{\link{iconic_prospect}()}
#' use when \code{confounding = "inferred"}.
#'
#' @section Inference methods:
#' \itemize{
#' \item \strong{conf_strength} (delta): the gap between the unadjusted
#' OLS estimate and the IV2SLS estimate, averaged across features.
#' Requires a valid exposure instrument G (F >= 10). Assumes the
#' exclusion restriction holds -- using IV2SLS validity to calibrate
#' a benchmark that tests IV2SLS validity is circular, so the
#' estimate should be interpreted as a best-case calibration.
#' \item \strong{mo_confounding} (delta_mo): the gap between the
#' IV2SLS NIE and the IV2SLS2 NIE, averaged across features.
#' Requires valid G + Gm (both F >= 10). Same circularity caveat.
#' \item \strong{omega} (omega_1, omega_2): the square root of the
#' R-squared from regressing each negative-control feature on the
#' outcome residualized on X + C, averaged across features.
#' Conflates NC coverage with confounder strength -- reported as a
#' composite, not pure coverage.
#' \item \strong{k}: the number of latent confounders, estimated via
#' parallel analysis (Horn, 1965) on the correlation matrix of
#' outcomes residualized on X + C. Requires at least 5 outcome
#' features. Returns a point estimate and a bootstrap confidence
#' interval.
#' }
#'
#' @param data An \code{iconic_data} object.
#' @param diagnosis Optional \code{iconic_diagnosis} from
#' \code{\link{iconic_diagnose}()}. Used to check instrument strength
#' before using estimator gaps.
#' @param estimate Optional estimate data frame from
#' \code{\link{iconic_estimate}()}. If \code{NULL}, estimates are
#' computed internally.
#' @param n_cores Number of parallel workers for omega inference and
#' k permutation. Default 1 (sequential). Uses
#' \code{parallel::mclapply} on Unix and a PSOCK cluster on Windows.
#'
#' @section Defaults (when inference is unavailable):
#' When a parameter cannot be inferred from the data, the following
#' defaults are used (from the simulation calibration):
#' \tabular{lll}{
#' \strong{Parameter} \tab \strong{Default} \tab \strong{Source} \cr
#' \code{conf_strength} (delta) \tab 0.8 \tab Simulation calibration \cr
#' \code{mo_confounding} (delta_mo) \tab 0.8 \tab Simulation calibration \cr
#' \code{omega_1, omega_2} \tab 0.7 \tab NC coverage (simulation calibration) \cr
#' \code{k} \tab 1 \tab Single-confounder assumption (typical mediation) \cr
#' \code{phi} \tab 0.8 \tab Strong mediator instrument assumption \cr
#' \code{lambda_XM}, \code{lambda_MY} \tab shared \tab Per-path confounder loadings \cr
#' }
#' These defaults are reported with a warning so the user knows the
#' value was not inferred from their data.
#'
#' @return An \code{iconic_confounding} S3 object: a named list with
#' \code{$conf_strength}, \code{$mo_confounding}, \code{$omega_1},
#' \code{$omega_2}, \code{$k}, \code{$unavailable} (character vector
#' of parameters that could not be inferred), and \code{$warnings}
#' (character vector of accumulated warnings). Each parameter slot
#' is a list with \code{estimate}, \code{method}, \code{assumption},
#' \code{available} (logical), and \code{warning} (character or NULL).
#' @export
#'
#' @examples
#' data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100*10), 10, 100),
#' M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
#' W = matrix(rnorm(100*10), 10, 100))
#' diag <- iconic_diagnose(data)
#' est <- iconic_estimate(data, diagnosis = diag)
#' conf <- infer_confounding(data, diagnosis = diag, estimate = est,
#'                           max_infer_tasks = 5)
#' print(conf)
#'
#' @param max_infer_tasks Cap on the number of mediators and of outcome
#' features used when \code{estimate = NULL} and estimates must be computed
#' internally. The confounding-strength gaps are averages across the
#' mediator x feature grid, so a random subset of at most
#' \code{max_infer_tasks} mediators and \code{max_infer_tasks} features
#' gives an unbiased Monte Carlo estimate at much lower cost. Default 50.
#' Panels smaller than the cap are used in full.
infer_confounding <- function(data, diagnosis = NULL, estimate = NULL,
                              n_cores = 1, max_infer_tasks = 50) {
  if (!inherits(data, "iconic_data"))
    stop("data must be an iconic_data object from iconic_data().")

  warnings <- character(0)
  unavailable <- c("rho_G1", "rho_G2") # always unestimable

  # -- Compute estimates if not supplied --
  # conf_strength (delta) and mo_confounding (delta_mo) are simple averages
  # of per-cell estimator gaps across the mediator x feature grid, so a
  # random subset of mediators/features gives an unbiased Monte Carlo
  # estimate of the full-panel value at a fraction of the cost. Subset the
  # grid to at most max_infer_tasks mediators and max_infer_tasks features
  # before running iconic_estimate(); small panels are used in full.
  inference_subset <- NULL
  if (is.null(estimate)) {
    sub_data <- .subset_data_for_inference(data, max_infer_tasks)
    inference_subset <- attr(sub_data, "inference_subset")
    estimate <- tryCatch(
      iconic_estimate(sub_data, diagnosis = diagnosis),
      error = function(e) NULL)
    if (is.null(estimate)) {
      warnings <- c(warnings,
        "Could not compute estimates for confounding inference; ",
        "all parameters will fall back to defaults.")
    }
  }

  # -- Instrument strength from diagnosis --
  # F_G and F_Gm can be vectors (one F-stat per mediator) when the
  # instrument is a panel. Use the median as a scalar summary for the
  # gap-based inference checks below.
  F_G <- NA_real_
  F_Gm <- NA_real_
  if (!is.null(diagnosis)) {
    F_G <- diagnosis$instrument_strength$F_G
    F_Gm <- diagnosis$instrument_strength$F_Gm
    if (length(F_G) > 1)
      F_G <- diagnosis$instrument_strength$F_G_median
    if (length(F_Gm) > 1)
      F_Gm <- diagnosis$instrument_strength$F_Gm_median
  }

  # === conf_strength (delta): UNADJ--IV2SLS gap ===
  conf_strength <- .infer_conf_strength(estimate, data, F_G, warnings)
  warnings <- conf_strength$warnings

  # === mo_confounding (delta_mo): IV2SLS--IV2SLS2 NIE gap ===
  mo_confounding <- .infer_mo_confounding(estimate, data, F_G, F_Gm, warnings)
  warnings <- mo_confounding$warnings

  # === omega_1, omega_2: W--Y residual R^2 ===
  omega_1 <- .infer_omega(data, path = 1, warnings, n_cores = n_cores)
  warnings <- omega_1$warnings
  omega_2 <- .infer_omega(data, path = 2, warnings, n_cores = n_cores)
  warnings <- omega_2$warnings

  # === k: parallel analysis on residualized outcomes ===
  k_inf <- .infer_k(data, warnings, n_cores = n_cores)
  warnings <- k_inf$warnings

  # Track unavailable parameters
  if (!conf_strength$available) unavailable <- c(unavailable, "conf_strength")
  if (!mo_confounding$available) unavailable <- c(unavailable, "mo_confounding")
  if (!omega_1$available) unavailable <- c(unavailable, "omega_1")
  if (!omega_2$available) unavailable <- c(unavailable, "omega_2")
  if (!k_inf$available) unavailable <- c(unavailable, "k")

  # Annotate the gap-based parameters when they were estimated on a random
  # subset of the mediator x feature grid (the estimate = NULL auto-run).
  if (!is.null(inference_subset) && isTRUE(inference_subset$subsetted)) {
    note <- sprintf("random subset (%d/%d mediators, %d/%d features)",
                    inference_subset$n_mediators_used, inference_subset$n_mediators_total,
                    inference_subset$n_features_used, inference_subset$n_features_total)
    if (!is.null(conf_strength$method))
      conf_strength$method <- paste0(conf_strength$method, "; ", note)
    if (!is.null(mo_confounding$method))
      mo_confounding$method <- paste0(mo_confounding$method, "; ", note)
  }

  obj <- list(
    conf_strength = conf_strength,
    mo_confounding = mo_confounding,
    omega_1 = omega_1,
    omega_2 = omega_2,
    k = k_inf,
    unavailable = unavailable,
    warnings = warnings,
    inference_subset = inference_subset
  )
  class(obj) <- c("iconic_confounding", "list")
  obj
}


#' Infer conf_strength (delta) from UNADJ--IV2SLS gap (internal)
#'
#' The gap between the unadjusted OLS estimate and the IV2SLS estimate
#' reflects the total confounding bias. Averaged across features and
#' scaled by the standard deviation of X.
#'
#' @keywords internal
#' @noRd
.infer_conf_strength <- function(estimate, data, F_G, warnings) {
  default <- list(
    estimate = 0.8, method = "default (no inference)",
    assumption = "none", available = FALSE, warning = NULL)

  if (is.null(estimate)) {
    warnings <- c(warnings,
      "conf_strength: no estimates available, using default 0.8.")
    default$warning <- "no estimates available"
    return(c(default, list(warnings = warnings)))
  }

  if (!data$has_instrument) {
    warnings <- c(warnings,
      "conf_strength: no exposure instrument (G), using default 0.8.")
    default$warning <- "no instrument available"
    return(c(default, list(warnings = warnings)))
  }

  if (is.na(F_G) || F_G < 10) {
    warnings <- c(warnings,
      sprintf("conf_strength: weak instrument (F_G=%.1f < 10), inference unreliable, using default 0.8.",
              ifelse(is.na(F_G), NA, F_G)))
    default$warning <- "weak instrument"
    return(c(default, list(warnings = warnings)))
  }

  # For total-effect data: use beta column (UNADJ vs IV2SLS)
  # For mediation data: use NDE column (UNADJ vs IV2SLS NDE)
  if ("beta" %in% names(estimate)) {
    unadj <- estimate$beta[estimate$method == "UNADJ"]
    iv2sls <- estimate$beta[estimate$method == "IV2SLS"]
    effect_label <- "total effect"
  } else if ("NDE" %in% names(estimate)) {
    unadj <- estimate$NDE[estimate$method == "UNADJ"]
    iv2sls <- estimate$NDE[estimate$method == "IV2SLS"]
    effect_label <- "NDE"
  } else {
    warnings <- c(warnings,
      "conf_strength: no beta or NDE column in estimates, using default 0.8.")
    default$warning <- "no effect estimates"
    return(c(default, list(warnings = warnings)))
  }

  if (length(unadj) == 0 || length(iv2sls) == 0 ||
      all(is.na(unadj)) || all(is.na(iv2sls))) {
    warnings <- c(warnings,
      "conf_strength: UNADJ or IV2SLS estimates missing, using default 0.8.")
    default$warning <- "estimates missing"
    return(c(default, list(warnings = warnings)))
  }

  # Gap = |UNADJ - IV2SLS|, averaged across features, scaled by sd(X)
  gap <- mean(abs(unadj - iv2sls), na.rm = TRUE)
  sd_X <- sd(data$X, na.rm = TRUE)
  if (sd_X == 0 || is.na(sd_X)) sd_X <- 1
  delta_est <- min(gap / sd_X, 1.5) # cap at 1.5 for stability

  warning_txt <- NULL
  if (delta_est < 0.1) {
    warning_txt <- paste0("inferred delta is very low; may indicate weak ",
                          "confounding or a strong instrument that already ",
                          "removes most bias")
  }

  list(
    estimate = as.numeric(delta_est),
    method = paste0("UNADJ-IV2SLS gap (", effect_label, ")"),
    assumption = "valid G (F>=10, exclusion restriction)",
    available = TRUE,
    warning = warning_txt,
    warnings = warnings
  )
}


#' Infer mo_confounding (delta_mo) from IV2SLS--IV2SLS2 NIE gap (internal)
#'
#' The gap between the IV2SLS NIE and the IV2SLS2 NIE reflects
#' mediator-outcome confounding. Averaged across features.
#'
#' @keywords internal
#' @noRd
.infer_mo_confounding <- function(estimate, data, F_G, F_Gm, warnings) {
  default <- list(
    estimate = 0.8, method = "default (no inference)",
    assumption = "none", available = FALSE, warning = NULL)

  if (is.null(estimate) || !("NIE" %in% names(estimate))) {
    warnings <- c(warnings,
      "mo_confounding: no mediation estimates available, using default 0.8.")
    default$warning <- "no mediation estimates"
    return(c(default, list(warnings = warnings)))
  }

  if (!data$is_mediation) {
    warnings <- c(warnings,
      "mo_confounding: non-mediation data, using default 0.8.")
    default$warning <- "non-mediation data"
    return(c(default, list(warnings = warnings)))
  }

  if (!data$has_mediator_instrument) {
    warnings <- c(warnings,
      "mo_confounding: no mediator instrument (Gm), using default 0.8.")
    default$warning <- "no mediator instrument"
    return(c(default, list(warnings = warnings)))
  }

  if (is.na(F_Gm) || F_Gm < 10) {
    warnings <- c(warnings,
      "mo_confounding: weak mediator instrument (F_Gm < 10), inference unreliable, using default 0.8.")
    default$warning <- "weak mediator instrument"
    return(c(default, list(warnings = warnings)))
  }

  iv2sls_nie <- estimate$NIE[estimate$method == "IV2SLS"]
  iv2sls2_nie <- estimate$NIE[estimate$method == "IV2SLS2"]

  if (length(iv2sls_nie) == 0 || length(iv2sls2_nie) == 0 ||
      all(is.na(iv2sls_nie)) || all(is.na(iv2sls2_nie))) {
    warnings <- c(warnings,
      "mo_confounding: IV2SLS or IV2SLS2 NIE missing, using default 0.8.")
    default$warning <- "NIE estimates missing"
    return(c(default, list(warnings = warnings)))
  }

  gap <- mean(abs(iv2sls_nie - iv2sls2_nie), na.rm = TRUE)
  sd_M <- if (data$is_mediation) sd(data$M[1, ], na.rm = TRUE) else 1
  if (sd_M == 0 || is.na(sd_M)) sd_M <- 1
  delta_mo_est <- min(gap / sd_M, 1.5)

  list(
    estimate = as.numeric(delta_mo_est),
    method = "IV2SLS-IV2SLS2 NIE gap",
    assumption = "valid G + Gm (both F>=10)",
    available = TRUE,
    warning = NULL,
    warnings = warnings
  )
}


#' Infer omega (NC coverage) from W--Y residual R^2 (internal)
#'
#' Regresses each NC feature on the outcome residualized on X + C and
#' takes sqrt(R^2), averaged across features. This is a composite of
#' coverage and confounder strength, not pure coverage.
#'
#' @param data An iconic_data object.
#' @param path 1 for omega_1 (W1), 2 for omega_2 (W2).
#' @param n_cores Number of parallel workers for the NC-feature loop.
#' @keywords internal
#' @noRd
.infer_omega <- function(data, path = 1, warnings, n_cores = 1) {
  default <- list(
    estimate = 0.7, method = "default (no inference)",
    assumption = "none", available = FALSE, warning = NULL)

  W_panel <- if (path == 1) data$W1 else data$W2
  path_label <- if (path == 1) "omega_1" else "omega_2"

  if (is.null(W_panel)) {
    warnings <- c(warnings,
      sprintf("%s: no negative controls available, using default 0.7.", path_label))
    default$warning <- "no negative controls"
    return(c(default, list(warnings = warnings)))
  }

  n <- data$n
  nf <- data$n_features
  X <- data$X
  cv <- data$covariates
  cnames <- if (!is.null(cv) && ncol(cv) > 0) names(cv) else character(0)
  cs <- .covar_str(cnames)

  # Residualize each outcome feature on X + C
  Y_resid <- matrix(NA_real_, nf, n)
  for (f in seq_len(nf)) {
    y <- data$Y[f, ]
    d <- .bind_covars(data.frame(y = y, X = X), cv)
    fit <- tryCatch(lm(as.formula(paste0("y ~ X", cs)), data = d),
                    error = function(e) NULL)
    Y_resid[f, ] <- if (!is.null(fit)) residuals(fit) else y
  }

  # For each NC feature, regress W on Y_resid and compute R^2
  # W_panel is q x n (features x samples); transpose to n x q
  W_mat <- t(W_panel)
  q <- ncol(W_mat)

  r2s <- unlist(.parallel_lapply(seq_len(q), function(j) {
    w <- W_mat[, j]
    # Use the outcome feature that best correlates with this NC
    best_r2 <- 0
    for (f in seq_len(nf)) {
      yr <- Y_resid[f, ]
      d <- .bind_covars(data.frame(w = w, yr = yr), cv)
      fit <- tryCatch(lm(as.formula(paste0("w ~ yr", cs)), data = d),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        r2 <- summary(fit)$r.squared
        if (!is.na(r2) && r2 > best_r2) best_r2 <- r2
      }
    }
    best_r2
  }, n_cores = n_cores, progress = paste0(path_label, " NC coverage")))

  omega_est <- sqrt(mean(r2s, na.rm = TRUE))
  if (is.na(omega_est) || omega_est > 1) omega_est <- 0.7

  list(
    estimate = as.numeric(omega_est),
    method = "sqrt(R^2) of W on Y residualized on X+C",
    assumption = "valid W (negative controls)",
    available = TRUE,
    warning = "composite: coverage x confounder strength, not pure coverage",
    warnings = warnings
  )
}


#' Infer k (number of latent confounders) via parallel analysis (internal)
#'
#' Performs Horn's parallel analysis on the correlation matrix of
#' outcomes residualized on X + C. Returns a point estimate and a
#' simple bootstrap confidence interval.
#'
#' @param data An iconic_data object.
#' @param n_cores Number of parallel workers for the permutation loop.
#' @keywords internal
#' @noRd
.infer_k <- function(data, warnings, n_cores = 1) {
  default <- list(
    estimate = 1L, ci = c(1L, 2L), method = "default (no inference)",
    assumption = "none", available = FALSE, warning = NULL)

  if (data$n_features < 5) {
    warnings <- c(warnings,
      sprintf("k: too few outcome features (%d < 5) for factor analysis, using default k=1.",
              data$n_features))
    default$warning <- "too few features for factor analysis"
    return(c(default, list(warnings = warnings)))
  }

  n <- data$n
  nf <- data$n_features
  X <- data$X
  cv <- data$covariates
  cnames <- if (!is.null(cv) && ncol(cv) > 0) names(cv) else character(0)
  cs <- .covar_str(cnames)

  # Residualize each outcome feature on X + C
  Y_resid <- matrix(NA_real_, n, nf)
  for (f in seq_len(nf)) {
    y <- data$Y[f, ]
    d <- .bind_covars(data.frame(y = y, X = X), cv)
    fit <- tryCatch(lm(as.formula(paste0("y ~ X", cs)), data = d),
                    error = function(e) NULL)
    Y_resid[, f] <- if (!is.null(fit)) residuals(fit) else y
  }

  # Correlation matrix of residualized outcomes
  cor_mat <- cor(Y_resid, use = "pairwise.complete.obs")
  if (any(is.na(cor_mat))) {
    warnings <- c(warnings,
      "k: correlation matrix has NAs (insufficient data), using default k=1.")
    default$warning <- "correlation matrix has NAs"
    return(c(default, list(warnings = warnings)))
  }

  # Eigenvalues of the correlation matrix
  eigen_vals <- eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values

  # Parallel analysis: compare to eigenvalues of random correlation matrices
  n_perms <- 100
  perm_eigen_list <- .parallel_lapply(seq_len(n_perms), function(b) {
    perm_mat <- matrix(rnorm(n * nf), n, nf)
    perm_cor <- cor(perm_mat)
    eigen(perm_cor, symmetric = TRUE, only.values = TRUE)$values
  }, n_cores = n_cores, progress = "k permutation analysis")
  perm_eigen <- do.call(rbind, perm_eigen_list)
  perm_mean <- colMeans(perm_eigen)

  # Number of eigenvalues exceeding the permutation mean
  k_est <- sum(eigen_vals > perm_mean)
  k_est <- max(1L, k_est) # at least 1 confounder

  # Simple CI: +/-1
  k_lo <- max(1L, k_est - 1L)
  k_hi <- k_est + 1L

  list(
    estimate = as.integer(k_est),
    ci = c(as.integer(k_lo), as.integer(k_hi)),
    method = "parallel analysis (Horn, 1965)",
    assumption = "factor model, >=5 outcome features",
    available = TRUE,
    warning = NULL,
    warnings = warnings
  )
}


#' Print method for iconic_confounding objects
#'
#' @param x An \code{iconic_confounding} object.
#' @param ... Unused.
#' @return Invisibly returns `x` (the `iconic_confounding` object); called for its side effect of printing a human-readable summary.
#' @export
#' @examples
#' set.seed(1)
#' data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
#'   G = rnorm(100), W = matrix(rnorm(100 * 10), 10, 100))
#' cf <- infer_confounding(data, max_infer_tasks = 5)
#' print(cf)
print.iconic_confounding <- function(x, ...) {
  cat("<iconic_confounding>\n")

  # Helper to print one parameter
  print_param <- function(label, p, is_k = FALSE) {
    if (p$available) {
      est_str <- if (is_k) {
        sprintf("%d [CI: %d, %d]", p$estimate, p$ci[1], p$ci[2])
      } else {
        sprintf("%.3f", p$estimate)
      }
      cat(sprintf(" %-16s %s (%s)\n", label, est_str, p$method))
      if (!is.null(p$warning))
        cat(sprintf(" warning: %s\n", p$warning))
    } else {
      cat(sprintf(" %-16s default (%.1f) -- %s\n", label, p$estimate,
                  ifelse(is.null(p$warning), "unavailable", p$warning)))
    }
  }

  print_param("conf_strength", x$conf_strength)
  print_param("mo_confounding", x$mo_confounding)
  print_param("omega_1", x$omega_1)
  print_param("omega_2", x$omega_2)
  print_param("k", x$k, is_k = TRUE)

  if (length(x$unavailable) > 0)
    cat("\n Unavailable:", paste(x$unavailable, collapse = ", "), "\n")

  if (length(x$warnings) > 0) {
    cat("\n Warnings:\n")
    for (w in x$warnings) cat(" ", w, "\n")
  }

  invisible(x)
}


#' Subset an iconic_data object to a random panel for confounding inference
#'
#' Draws a random subset of mediators and (when more than one) outcome
#' features, capping each at \code{max_infer_tasks}, and returns a thinned
#' \code{iconic_data} object suitable for \code{iconic_estimate()}. Used to
#' make the \code{estimate = NULL} auto-run in \code{infer_confounding()}
#' fast on large panels: the confounding-strength gaps are averages across
#' the mediator x feature grid, so a random subset is an unbiased Monte
#' Carlo estimate of the full-panel value.
#'
#' @param data An \code{iconic_data} object.
#' @param max_infer_tasks Cap on mediators and on features.
#' @return A thinned \code{iconic_data} object with an
#'   \code{inference_subset} attribute recording the subset sizes.
#' @keywords internal
.subset_data_for_inference <- function(data, max_infer_tasks = 50) {
  nm_total <- if (!is.null(data$n_mediators)) data$n_mediators else 0L
  nf_total <- if (!is.null(data$n_features)) data$n_features else 1L

  n_m_use <- min(nm_total, max_infer_tasks)
  n_f_use <- min(nf_total, max_infer_tasks)
  subsetted <- (nm_total > max_infer_tasks) || (nf_total > max_infer_tasks)

  # Nothing to subset: return the data unchanged.
  if (!subsetted) {
    attr(data, "inference_subset") <- list(
      subsetted = FALSE,
      n_mediators_used = nm_total, n_mediators_total = nm_total,
      n_features_used = nf_total, n_features_total = nf_total)
    return(data)
  }

  m_idx <- if (nm_total > 0) sort(sample(seq_len(nm_total), n_m_use)) else integer(0)
  f_idx <- if (nf_total > 1) sort(sample(seq_len(nf_total), n_f_use)) else seq_len(nf_total)

  out <- data
  # Subset the mediator panel (mediators x samples) and its names.
  if (!is.null(out$M) && length(m_idx) > 0) {
    out$M <- out$M[m_idx, , drop = FALSE]
    if (!is.null(out$mediator_names) && length(out$mediator_names) == nm_total)
      out$mediator_names <- out$mediator_names[m_idx]
    # Gm may be per-mediator (matrix with a row per mediator); subset to match.
    if (!is.null(out$Gm) && is.matrix(out$Gm) && nrow(out$Gm) == nm_total)
      out$Gm <- out$Gm[m_idx, , drop = FALSE]
    out$n_mediators <- length(m_idx)
  }
  # Subset the outcome feature panel (features x samples) when multi-feature.
  if (!is.null(out$Y) && is.matrix(out$Y) && nf_total > 1 && nrow(out$Y) == nf_total) {
    out$Y <- out$Y[f_idx, , drop = FALSE]
    if (!is.null(out$feature_names) && length(out$feature_names) == nf_total)
      out$feature_names <- out$feature_names[f_idx]
    out$n_features <- length(f_idx)
  }

  attr(out, "inference_subset") <- list(
    subsetted = TRUE,
    n_mediators_used = n_m_use, n_mediators_total = nm_total,
    n_features_used = n_f_use, n_features_total = nf_total)
  out
}
