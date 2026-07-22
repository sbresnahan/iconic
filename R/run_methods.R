# ============================================================
# Apply the estimators across the features of a dataset produced by
# generate_toy_data() or run_single_iteration().
#
# run_methods()            – internal driver used by the simulation sweeps.
# analyze_methods_robust() – exported, reference-named entry point; adds a
#                            `significant` flag and optional feature subsetting.
# analyze_methods_parallel() – same, parallelised across features.
#
# All three share .analyze_feature(), which passes the FULL W matrix
# to the matrix-bridge fit_pgc() (v0.3.1).  A scalar-bridge variant
# (fit_pgc_scalar) is exported for standalone use but is NOT included
# in the default pipeline.
# ============================================================

.methods_all <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC")

#' Estimate total effect for a single feature with specified methods (internal)
#'
#' Format-agnostic per-feature estimation: takes explicit vectors and
#' covariates (no dependency on the `dat` list format).  Called by both
#' the simulation driver (.analyze_feature) and the real-data driver
#' (iconic_estimate).
#'
#' @param Z       Numeric exposure vector (length n).
#' @param y       Numeric outcome vector (length n).
#' @param g       Numeric instrument vector (length n), or NULL.
#' @param w       Numeric NC vector (length n) or matrix (n x q), or NULL.
#'                Used for DIRECT, IV2SLS (full panel as covariates), and
#'                COCA (via W_avg scalar).
#' @param W_mat   Numeric NC matrix (n x q), or NULL.  Used for PGC
#'                (matrix bridge).  If NULL, PGC is skipped.
#' @param W_avg   Numeric vector (length n): row means of the NC panel
#'                for COCA.  If NULL but W_mat is present, computed inline.
#' @param covars  Optional data frame of covariates (n rows).
#' @param methods Character vector of methods to run.  Default: all five.
#'                Methods whose required inputs are missing are silently
#'                skipped.
#' @param feature_idx Integer or character: feature identifier for the
#'                output `feature` column.  Default 1L.
#' @return Data frame: `feature`, `method`, `beta`, `se`, `pvalue`.
#'         Returns NULL if fewer than 20 complete cases.
#' @keywords internal
.estimate_total_feature <- function(Z, y, g = NULL, w = NULL, W_mat = NULL,
                                    W_avg = NULL, covars = NULL,
                                    methods = .methods_all,
                                    feature_idx = 1L, min_f = 10) {
  na <- list(beta = NA_real_, se = NA_real_, pvalue = NA_real_)
  row <- function(m, r) data.frame(feature = feature_idx, method = m,
                                   beta = as.numeric(r$beta),
                                   se = as.numeric(r$se),
                                   pvalue = as.numeric(r$pvalue),
                                   stringsAsFactors = FALSE)

  # Determine which methods can actually run given available inputs
  can_run <- character(0)
  if ("UNADJ" %in% methods) can_run <- c(can_run, "UNADJ")
  if ("DIRECT" %in% methods && !is.null(g) && !is.null(w)) can_run <- c(can_run, "DIRECT")
  if ("COCA" %in% methods && !is.null(w)) can_run <- c(can_run, "COCA")
  if ("IV2SLS" %in% methods && !is.null(g) && !is.null(w)) can_run <- c(can_run, "IV2SLS")
  if ("PGC" %in% methods && !is.null(g) && !is.null(W_mat)) can_run <- c(can_run, "PGC")

  if (!length(can_run)) return(NULL)

  # Complete cases across all needed variables
  needed <- cbind(y, Z)
  if (!is.null(g)) needed <- cbind(needed, g)
  if (!is.null(w)) needed <- cbind(needed, as.matrix(w))
  if (!is.null(covars) && ncol(covars) > 0) needed <- cbind(needed, covars)
  ok <- stats::complete.cases(needed)
  if (sum(ok) < 20) return(NULL)

  Z_f  <- Z[ok]
  y_f  <- y[ok]
  g_f  <- if (!is.null(g)) g[ok] else NULL
  w_f  <- if (!is.null(w)) as.matrix(w)[ok, , drop = FALSE] else NULL
  cv_f <- if (!is.null(covars)) covars[ok, , drop = FALSE] else NULL
  W_mat_f <- if (!is.null(W_mat)) W_mat[ok, , drop = FALSE] else NULL
  Wa_f <- if (!is.null(W_avg)) W_avg[ok] else if (!is.null(W_mat_f)) rowMeans(W_mat_f) else NULL

  rows <- list()

  if ("UNADJ" %in% can_run) {
    r <- tryCatch({
      fit <- lm(y_f ~ Z_f); sm <- summary(fit)$coefficients
      list(beta = coef(fit)["Z_f"], se = sm["Z_f", 2], pvalue = sm["Z_f", 4])
    }, error = function(e) na)
    rows[["UNADJ"]] <- row("UNADJ", r)
  }
  if ("DIRECT" %in% can_run) {
    r <- tryCatch(fit_direct(y_f, Z_f, g_f, w_f, cv_f), error = function(e) na)
    rows[["DIRECT"]] <- row("DIRECT", r)
  }
  if ("COCA" %in% can_run) {
    r <- tryCatch(fit_coca(y_f, Z_f, Wa_f, cv_f), error = function(e) na)
    rows[["COCA"]] <- row("COCA", r)
  }
  if ("IV2SLS" %in% can_run) {
    r <- tryCatch(fit_iv2sls(y_f, Z_f, g_f, w_f, cv_f, min_f = min_f),
                  error = function(e) na)
    rows[["IV2SLS"]] <- row("IV2SLS", r)
  }
  if ("PGC" %in% can_run) {
    r <- tryCatch(fit_pgc(y_f, Z_f, g_f, W_mat_f, cv_f), error = function(e) na)
    rows[["PGC"]] <- row("PGC", r)
  }

  do.call(rbind, rows)
}

#' Analyse a single outcome feature with all estimators (internal)
#'
#' Thin wrapper around .estimate_total_feature() that extracts vectors
#' from the `dat` list format used by the simulation pipeline.
#'
#' @param dat       Dataset list (from generate_toy_data / run_single_iteration).
#' @param f         Feature (column) index into `dat$Y`, `dat$W`, `dat$G`.
#' @param W_avg     Row means of the full negative-control panel (for COCA).
#' @param W_valid   Optional: validity-screened W matrix (n x q) for the
#'                  matrix-bridge PGC.  If NULL, uses the full W panel.
#' @return Data frame of five rows (one per method) or `NULL` if too few cases.
#' @keywords internal
.analyze_feature <- function(dat, f, W_avg, W_valid = NULL) {
  Z  <- dat$Z
  cv <- dat$synthetic_data
  y  <- dat$Y[, f]
  g  <- if (!is.null(dat$G)) dat$G[, f] else NULL

  W_mat <- if (!is.null(W_valid)) W_valid else dat$W
  # v0.8.4: pass the full W panel (n x q) to all estimators.
  w <- W_mat  # full matrix (NULL when no W)

  .estimate_total_feature(
    Z = Z, y = y, g = g, w = w, W_mat = W_mat,
    W_avg = W_avg, covars = cv,
    methods = .methods_all, feature_idx = f
  )
}


#' Apply all estimators across features (internal)
#'
#' @param dat        List returned by `generate_toy_data()` / `run_single_iteration()`.
#' @param n_features Number of outcome columns to process.
#' @param W_valid    Optional: validity-screened W matrix for matrix-bridge PGC.
#' @param n_cores    Number of parallel workers. Default 1 (sequential).
#' @return Data frame with columns: feature, method, beta, se, pvalue.
#' @keywords internal
run_methods <- function(dat, n_features = ncol(dat$Y), W_valid = NULL,
                        n_cores = 1) {
  W_avg   <- rowMeans(dat$W)
  results <- .parallel_lapply(seq_len(n_features),
                    function(f) .analyze_feature(dat, f, W_avg, W_valid),
                    n_cores = n_cores)
  do.call(rbind, Filter(Negate(is.null), results))
}


#' Run all five estimators (plus UNADJ) on one synthetic dataset
#'
#' Reference-named entry point: applies UNADJ, DIRECT, COCA, IV2SLS, and PGC
#' (matrix bridge) to each tested outcome feature and returns tidy per-feature
#' results with a significance flag.  Negative controls are summarised over
#' the full W panel for COCA, so restricting `test_features` does not change
#' its control summary.  The matrix-bridge PGC uses the full W panel
#' (or a validity-screened subset if provided).
#'
#' A scalar-bridge variant ([fit_pgc_scalar()]) is exported for standalone
#' use but is not included in the default pipeline.
#'
#' @param iteration_data Dataset list from [run_single_iteration()] (or
#'   [generate_toy_data()]).
#' @param test_features  Optional integer indices of outcome features to test.
#'   Default `NULL` (all features).
#' @param alpha          Significance threshold for the `significant` flag. Default 0.05.
#' @param debug          If `TRUE`, message per-feature progress. Default `FALSE`.
#' @param n_cores        Number of parallel workers. Default 1 (sequential).
#'   Uses \code{parallel::mclapply} on Unix and a PSOCK cluster on Windows.
#'
#' @return Data frame: `feature`, `method`, `beta`, `se`, `pvalue`, `significant`.
#' @export
analyze_methods_robust <- function(iteration_data, test_features = NULL,
                                   alpha = 0.05, debug = FALSE, n_cores = 1) {
  feats <- if (is.null(test_features)) seq_len(ncol(iteration_data$Y)) else test_features
  W_avg <- rowMeans(iteration_data$W)
  res <- .parallel_lapply(feats, function(f) {
    if (debug) message("feature ", f)
    .analyze_feature(iteration_data, f, W_avg)
  }, n_cores = n_cores, progress = if (debug) "analyze_methods_robust" else NULL)
  out <- do.call(rbind, Filter(Negate(is.null), res))
  if (is.null(out)) return(out)
  out$significant <- out$pvalue < alpha
  out
}

#' Parallel version of [analyze_methods_robust()]
#'
#' Distributes the per-feature analysis across workers via
#' [`parallel::mclapply`] (falling back to sequential on Windows / one core).
#'
#' @inheritParams analyze_methods_robust
#' @param n_cores Number of workers. Default 1.
#' @return Data frame: `feature`, `method`, `beta`, `se`, `pvalue`, `significant`.
#' @export
analyze_methods_parallel <- function(iteration_data, test_features = NULL,
                                     alpha = 0.05, debug = FALSE, n_cores = 1) {
  feats <- if (is.null(test_features)) seq_len(ncol(iteration_data$Y)) else test_features
  W_avg <- rowMeans(iteration_data$W)
  res <- .parallel_lapply(feats, function(f) .analyze_feature(iteration_data, f, W_avg),
                          n_cores = n_cores, progress = "analyze_methods_parallel")
  out <- do.call(rbind, Filter(Negate(is.null), res))
  if (is.null(out)) return(out)
  out$significant <- out$pvalue < alpha
  out
}


#' Summarise simulation results across features (internal)
#'
#' @param combined   Data frame from `run_methods()`.
#' @param true_total Scalar true total causal effect.
#' @return Data frame with one row per method: mean, median, sd, bias,
#'   abs_bias, rmse, power, n.
#' @keywords internal
summarise_results <- function(combined, true_total) {
  methods <- .methods_all

  rows <- lapply(methods, function(m) {
    x <- combined$beta[combined$method == m]
    p <- combined$pvalue[combined$method == m]
    data.frame(
      method   = m,
      mean     = mean(x, na.rm = TRUE),
      median   = median(x, na.rm = TRUE),
      sd       = sd(x, na.rm = TRUE),
      bias     = mean(x, na.rm = TRUE) - true_total,
      abs_bias = abs(mean(x, na.rm = TRUE) - true_total),
      rmse     = sqrt(mean((x - true_total)^2, na.rm = TRUE)),
      power    = mean(p < 0.05, na.rm = TRUE),
      n        = sum(!is.na(x)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
