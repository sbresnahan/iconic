# ============================================================
#. The degradation-surface construction and tipping-point
# algorithm move to a technical appendix when the draft is condensed.
# No code change.
# ============================================================
# iconic_sensitivity: Degradation surface calibrated to user data.
#
# Wraps the existing gan_mediation_sensitivity() machinery with a
# 2D rho_G1 x rho_G2 sweep, producing a "degradation surface" that
# shows how each estimator's bias changes as instrument-independence
# assumptions are violated.
#
# The surface is calibrated to the user's data: n_samples, phi
# (mediator instrument strength), and mo_confounding are inferred
# from the iconic_data object and diagnosis when available.
#
# Auto-trains a GAN (or MVN fallback) from the user's
# iconic_data object when no trained_gan is supplied, so the
# synthetic texture matches the user's cohort. Also supports
# confounding = "inferred" to calibrate confounding parameters
# from the data via infer_confounding().
# ============================================================

#' Auto-train a texture model from an iconic_data object (internal)
#'
#' Converts the user's iconic_data to the load_real_input_data() format
#' and trains a GAN (or multivariate-normal fallback) on the resulting
#' tidy frame. The trained model supplies realistic covariate and
#' outcome texture to run_single_iteration(), and carries feature-level
#' residual correlation matrices for the Y, M, and W panels so the
#' simulation can inject correlated noise.
#'
#' @param data An \code{iconic_data} object.
#' @param epochs GAN training epochs. Default 100.
#' @return An \code{iconic_gan} object.
#' @keywords internal
.auto_train_gan <- function(data, epochs = 100) {
  if (!inherits(data, "iconic_data"))
    stop("data must be an iconic_data object.")

  # Convert iconic_data → load_real_input_data format.
  # X is a length-n vector; load_real_input_data expects a features x
  # samples matrix, so wrap it as 1 x n.
  X_mat <- matrix(data$X, nrow = 1L, ncol = data$n)
  Y_mat <- data$Y # already features x samples in iconic_data

  # Mediator panel: iconic_data stores M as features x samples (or NULL)
  M_mat <- data$M

  # Negative-control panel: use W (the combined panel) for correlation
  # learning. When path-specific W1/W2 are present, W is the combined
  # panel; either way, the correlation structure is what we need.
  W_mat <- data$W

  cov_df <- if (!is.null(data$covariates) && ncol(data$covariates) > 0)
    data$covariates else NULL

  input <- load_real_input_data(
    X_matrix = X_mat,
    Y_matrix = Y_mat,
    M_matrix = M_mat,
    W_matrix = W_mat,
    covariates_df = cov_df)

  train_gan_on_real_data(input$gan_training_data,
                         feature_correlations = input$feature_correlations,
                         feature_texture = input$feature_texture,
                         epochs = epochs,
                         verbose = FALSE)
}

#' Resolve the trained GAN for a simulation function (internal)
#'
#' Priority: explicit \code{trained_gan} argument > \code{data$trained_gan}
#' (attached at iconic_data() construction) > auto-train from \code{data}.
#'
#' @param trained_gan Explicit argument (may be NULL).
#' @param data An \code{iconic_data} object.
#' @param epochs Epochs for auto-training. Default 100.
#' @return A list with \code{gan} (the resolved iconic_gan) and
#' \code{source} (a string: "user-supplied GAN", "data-attached GAN",
#' or "auto-trained from data").
#' @keywords internal
.resolve_gan <- function(trained_gan, data, epochs = 100) {
  if (!is.null(trained_gan)) {
    return(list(gan = trained_gan, source = "user-supplied GAN"))
  }
  if (!is.null(data$trained_gan)) {
    return(list(gan = data$trained_gan, source = "data-attached GAN"))
  }
  # survival outcomes have no continuous Y matrix, so the GAN
  # texture model cannot be trained. Use default (NULL) texture instead.
  if (!is.null(data$outcome_type) && data$outcome_type == "survival") {
    return(list(gan = NULL, source = "default texture (survival outcome)"))
  }
  gan <- .auto_train_gan(data, epochs = epochs)
  list(gan = gan, source = "auto-trained from data")
}

#' Sensitivity (degradation) surface and effect-decomposition bias sweep
#'
#' Sweeps a 2D grid of instrument-independence violations
#' (\code{rho_G1} x \code{rho_G2}) and reports how each estimator's
#' NDE/NIE bias degrades as the genetic instruments become correlated
#' with the unmeasured confounders.
#'
#' The grid is calibrated to the user's data: sample size, mediator
#' instrument strength (\code{phi}), confounding level, and covariate
#' / outcome texture are inferred from the \code{iconic_data} object
#' and diagnosis when available.
#'
#' @aliases effect_decomposition_bias_sweep
#'
#' @section When to use:
#' Use this as an **effect-decomposition bias sweep** when you need to
#' know how robust a mediation estimate (NDE/NIE) is to violations of
#' the genetic instruments' independence from the confounders. The
#' degradation surface shows, for each estimator, the point at which
#' bias becomes material — the "tipping point" — so you can report not
#' just a point estimate but the range of violations it tolerates.
#' This is the mediation analogue of a total-effect sensitivity sweep
#' and is most informative when comparing estimators that make
#' different independence assumptions (e.g. IV2SLS2 vs PGC2Gm).
#'
#' @section Texture model:
#' When \code{trained_gan} is \code{NULL} and no GAN is attached to
#' \code{data}, a texture model is auto-trained from the user's data
#' (GAN via \code{torch} if available, otherwise a multivariate-normal
#' fallback). This ensures the synthetic covariate and outcome texture
#' matches the user's cohort. Supply \code{trained_gan} explicitly or
#' attach one via \code{\link{iconic_data}(trained_gan = ...)} to reuse
#' a pre-trained model and avoid retraining.
#'
#' @section Confounding calibration:
#' The \code{confounding} argument controls how the held-fixed
#' confounding parameters (\code{mo_confounding}, \code{omega_1},
#' \code{omega_2}) are set:
#' \itemize{
#' \item \code{"default"}: fixed defaults (0.8, 0.7, 0.7).
#' \item \code{"inferred"}: calls \code{\link{infer_confounding}()} to
#' estimate parameters from the data. Parameters that cannot be
#' inferred fall back to defaults with warnings.
#' \item \code{"manual"}: use the explicitly supplied arguments.
#' }
#'
#' @param data An \code{iconic_data} object.
#' @param diagnosis Optional \code{iconic_diagnosis} from
#' \code{\link{iconic_diagnose}()}. Used to infer \code{phi} and
#' \code{mo_confounding} when not explicitly supplied.
#' @param trained_gan Optional \code{iconic_gan} from
#' \code{\link{train_gan_on_real_data}()}. If \code{NULL}, a texture
#' model is auto-trained from \code{data}.
#' @param confounding Confounding parameter source: \code{"default"}
#' (fixed defaults), \code{"inferred"} (data-calibrated via
#' \code{\link{infer_confounding}()}), or \code{"manual"} (use
#' explicitly supplied arguments). Default \code{"default"}.
#' @param gan_epochs Epochs for auto-trained GAN. Default 100.
#' @param rho_G1_grid Values of rho_G1 (G correlation with conf_XM).
#' Default \code{c(0, 0.1, 0.2, 0.3, 0.5)}.
#' @param rho_G2_grid Values of rho_G2 (Gm correlation with conf_MY).
#' Default \code{c(0, 0.1, 0.2, 0.3, 0.5)}.
#' @param n_iter Replicates per grid cell. Default 30.
#' @param n_samples Samples per replicate. If NULL, uses \code{data$n}.
#' @param n_features Features per replicate. Default 10.
#' @param phi Mediator instrument strength. If NULL, inferred
#' from diagnosis (F_Gm) or defaults to 0.8.
#' @param mo_confounding M-O confounding strength. Default 0.8.
#' Used when \code{confounding = "default"} or \code{"manual"}.
#' @param lambda_XM Optional per-path confounder loading vector (X->M path).
#'   NULL (default) = shared loadings.
#' @param lambda_MY Optional per-path confounder loading vector (M->Y path).
#'   NULL (default) = shared loadings.
#' @param omega_1,omega_2 NC coverage of each path's confounder composite
#'   (conf_XM / conf_MY). Default \code{c(0.3, 0.7, 1.0)}, swept jointly with
#'   the rho grid. When the two vectors are identical (the default), the sweep
#'   is taken on the diagonal (\code{omega_1 == omega_2}); supply distinct
#'   vectors to cross the full \code{omega_1 x omega_2} grid.
#' Used when \code{confounding = "default"} or \code{"manual"}.
#' @param bias_threshold Absolute bias threshold for tipping-point
#' annotation. Default 0.10.
#' @param base_seed Base RNG seed. Default 700.
#' @param n_cores Number of parallel workers for simulation replicates.
#' Default 1 (sequential). Uses \code{parallel::mclapply} on Unix
#' and a PSOCK cluster on Windows.
#' @param verbose Logical: print progress messages during the sweep.
#' Default \code{FALSE} (quiet).
#' @param outcome_type \code{NULL} (inherit from \code{data}, default) or
#' \code{"continuous"} / \code{"survival"}. When survival, the
#' sensitivity sweep uses the Cox / RMST survival mediation drivers.
#' @param effect_scale \code{"loghr"} (default) or \code{"rmst"}. Only
#' used when \code{outcome_type = "survival"}.
#' @param surv_h0 Baseline hazard for survival DGP. Default 0.1.
#' @param surv_event_frac Target fraction of observed events. Default 0.6.
#' @param surv_censor_rate Explicit censoring rate. Default NULL.
#'
#' @section Defaults:
#' \tabular{lll}{
#' \strong{Parameter} \tab \strong{Default} \tab \strong{Source} \cr
#' \code{confounding} \tab "default" \tab Use DGP defaults below \cr
#' \code{gan_epochs} \tab 100 \tab Texture-model training budget \cr
#' \code{rho_G1_grid} \tab c(0,0.1,0.2,0.3,0.5) \tab Independence-violation sweep \cr
#' \code{rho_G2_grid} \tab c(0,0.1,0.2,0.3,0.5) \tab Independence-violation sweep \cr
#' \code{n_iter} \tab 30 \tab Replicates per grid cell \cr
#' \code{n_features} \tab 10 \tab Features per replicate \cr
#' \code{mo_confounding} \tab 0.8 \tab Simulation calibration (delta_mo) \cr
#' \code{lambda_XM}, \code{lambda_MY} \tab shared \tab Per-path confounder loadings \cr
#' \code{omega_1, omega_2} \tab c(0.3,0.7,1.0) \tab NC coverage (swept on the diagonal) \cr
#' \code{bias_threshold} \tab 0.10 \tab Tipping-point threshold \cr
#' }
#'
#' @return An \code{iconic_sensitivity} S3 object: a named list with
#' \code{$surface} (data frame: rho_G1, rho_G2, method, NDE_bias,
#' NIE_bias, NDE_rmse, NIE_rmse, NDE_type1, NIE_type1, tipped),
#' \code{$grid}, \code{$tipping_points}, \code{$summary},
#' \code{$texture_source} (how the GAN was obtained), and
#' \code{$inferred_confounding} (when \code{confounding = "inferred"}).
#' @export
#'
#' @examples
#' \dontrun{
#' data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100*10), 10, 100),
#' M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
#' W = matrix(rnorm(100*10), 10, 100))
#' sens <- iconic_sensitivity(data, n_iter = 10)
#' print(sens)
#' }
iconic_sensitivity <- function(data, diagnosis = NULL,
                               trained_gan = NULL,
                               confounding = c("default", "inferred", "manual"),
                               gan_epochs = 100,
                               rho_G1_grid = c(0, 0.1, 0.2, 0.3, 0.5),
                               rho_G2_grid = c(0, 0.1, 0.2, 0.3, 0.5),
                               n_iter = 30,
                               n_samples = NULL,
                               n_features = 10,
                               phi = NULL,
                               mo_confounding = 0.8,
                               lambda_XM = NULL,
                               lambda_MY = NULL,
                               omega_1 = c(0.3, 0.7, 1.0),
                               omega_2 = c(0.3, 0.7, 1.0),
                               bias_threshold = 0.10,
                               base_seed = 700,
                               n_cores = 1,
                               verbose = FALSE,
                               outcome_type = NULL,
                               effect_scale = c("loghr", "rmst"),
                               surv_h0 = 0.1,
                               surv_event_frac = 0.6,
                               surv_censor_rate = NULL) {
  if (!inherits(data, "iconic_data"))
    stop("data must be an iconic_data object from iconic_data().")
  if (!data$is_mediation)
    stop("iconic_sensitivity requires mediation data (supply M). ",
         "Use gan_sensitivity() for total-effect-only sensitivity.")

  confounding <- match.arg(confounding)
  effect_scale <- match.arg(effect_scale)

  # resolve outcome_type — inherit from data if not explicitly set.
  if (is.null(outcome_type)) {
    outcome_type <- if (!is.null(data$outcome_type)) data$outcome_type else "continuous"
  }
  outcome_type <- match.arg(outcome_type, c("continuous", "survival"))

  # ── Resolve the texture model ──
  gan_res <- .resolve_gan(trained_gan, data, epochs = gan_epochs)
  gan <- gan_res$gan
  texture_source <- gan_res$source

  # ── Resolve confounding parameters ──
  inferred_conf <- NULL
  if (confounding == "inferred") {
    # Need estimates for the gap-based inference
    est_for_inf <- iconic_estimate(data, diagnosis = diagnosis)
    inferred_conf <- infer_confounding(data, diagnosis = diagnosis,
                                       estimate = est_for_inf)
    # Use inferred values where available; fall back to defaults otherwise
    if (inferred_conf$mo_confounding$available)
      mo_confounding <- inferred_conf$mo_confounding$estimate
    if (inferred_conf$omega_1$available)
      omega_1 <- inferred_conf$omega_1$estimate
    if (inferred_conf$omega_2$available)
      omega_2 <- inferred_conf$omega_2$estimate
  }
  # "default" and "manual" use the explicitly supplied arguments as before

  # ── Build the omega sweep ──
  # omega_1/omega_2 may be vectors; the sweep crosses rho_G1 x rho_G2 x
  # omega_1 x omega_2. When confounding = "inferred" supplied a scalar
  # omega, the sweep reduces to that single value.
  omega_1_grid <- sort(unique(omega_1))
  omega_2_grid <- sort(unique(omega_2))
  omega_swept <- length(omega_1_grid) > 1 || length(omega_2_grid) > 1

  # Calibrate to user data
  if (is.null(n_samples)) n_samples <- data$n
  if (is.null(phi)) {
    if (!is.null(diagnosis) && !all(is.na(diagnosis$instrument_strength$F_Gm))) {
      # Strong Gm -> phi = 0.8; weak -> phi = 0.3; none -> phi = 0
      F_gm <- diagnosis$instrument_strength$F_Gm
      # F_Gm can be a vector (panel); use median as scalar summary
      if (length(F_gm) > 1)
        F_gm <- diagnosis$instrument_strength$F_Gm_median
      phi <- if (is.na(F_gm)) 0 else if (F_gm >= 100) 0.8 else if (F_gm >= 10) 0.5 else 0.3
    } else {
      phi <- if (data$has_mediator_instrument) 0.8 else 0
    }
  }

  # Build the sweep: rho_G1 x rho_G2 x omega_1 x omega_2. When the two omega
  # vectors are identical (the default), sweep coverage on the diagonal
  # (omega_1 == omega_2), matching the manuscript's degradation-surface
  # design; distinct omega vectors cross the full omega_1 x omega_2 grid.
  grid <- expand.grid(rho_G1 = rho_G1_grid, rho_G2 = rho_G2_grid,
                      omega_1 = omega_1_grid, omega_2 = omega_2_grid,
                      KEEP.OUT.ATTRS = FALSE)
  if (identical(omega_1_grid, omega_2_grid) && omega_swept)
    grid <- grid[grid$omega_1 == grid$omega_2, , drop = FALSE]

  true_NDE <- 0.10 # beta_X default
  true_NIE <- 0.15 # alpha_M * beta_M default

  n_grid <- nrow(grid)
  smry <- lapply(seq_len(n_grid), function(gi) {
    r1 <- grid$rho_G1[gi]
    r2 <- grid$rho_G2[gi]
    o1 <- grid$omega_1[gi]
    o2 <- grid$omega_2[gi]
    if (n_grid > 1 && isTRUE(verbose))
      message("Sensitivity grid cell ", gi, "/", n_grid,
              " (rho_G1=", r1, ", rho_G2=", r2,
              if (omega_swept) paste0(", omega_1=", o1, ", omega_2=", o2) else "", ")")

    worker <- function(i) {
      dat <- run_single_iteration(
        trained_gan = gan,
        n_synthetic_samples = n_samples, n_features = n_features,
        mo_confounding = mo_confounding, phi = phi,
        rho_G1 = r1, rho_G2 = r2, rho_pop = 0,
        lambda_XM = lambda_XM, lambda_MY = lambda_MY, omega_1 = o1, omega_2 = o2,
        seed = base_seed + gi * 1000L + i,
        outcome_type = outcome_type, surv_h0 = surv_h0,
        surv_event_frac = surv_event_frac, surv_censor_rate = surv_censor_rate)
      if (outcome_type == "survival") {
        res <- .run_surv_methods(dat, effect_scale = effect_scale,
                                 is_mediation = TRUE)
      } else {
        res <- run_mediation_methods(dat, n_features)
      }
      res$iter <- i
      res
    }

    combined <- do.call(rbind, .parallel_lapply(seq_len(n_iter), worker,
                                                n_cores = n_cores,
                                                progress = if (isTRUE(verbose)) " Replicates" else NULL))
    s <- summarise_mediation_results(combined, true_NDE, true_NIE)
    s$rho_G1 <- r1
    s$rho_G2 <- r2
    s$omega_1 <- o1
    s$omega_2 <- o2
    s
  })

  surface <- do.call(rbind, smry)

  # Tipping-point annotation: flag cells where |bias| exceeds threshold
  surface$tipped_NDE <- abs(surface$NDE_bias) > bias_threshold
  surface$tipped_NIE <- abs(surface$NIE_bias) > bias_threshold
  surface$tipped <- surface$tipped_NDE | surface$tipped_NIE

  # Reorder columns
  front <- c("rho_G1", "rho_G2", "omega_1", "omega_2", "method")
  surface <- surface[, c(front, setdiff(names(surface), front))]

  # Find tipping points: first rho_G2 (along each rho_G1 row) where a
  # method tips, for each method
  tipping <- .find_tipping_points(surface, rho_G1_grid, rho_G2_grid,
                                  bias_threshold)

  # Summary
  summary_txt <- .build_sensitivity_summary(surface, tipping, bias_threshold,
                                            n_iter)

  obj <- list(
    surface = surface,
    grid = grid,
    tipping_points = tipping,
    bias_threshold = bias_threshold,
    n_iter = n_iter,
    n_samples = n_samples,
    phi = phi,
    mo_confounding = mo_confounding,
    omega_1 = omega_1_grid,
    omega_2 = omega_2_grid,
    omega_swept = omega_swept,
    texture_source = texture_source,
    inferred_confounding = inferred_conf,
    summary = summary_txt,
    # scenario manifest for reader orientation.
    manifest = scenario_manifest(
      list(beta_X = true_NDE, alpha_M = 0.50, beta_M = 0.30,
           n_mediators = 1, n = n_samples, n_features = n_features,
           lambda_XM = lambda_XM, lambda_MY = lambda_MY, feat_cor = 0,
           conf_str = NA, mo_confounding = mo_confounding, phi = phi,
           w_signal = omega_1_grid[1], rho_G1 = NA, rho_G2 = NA, rho_pop = NA),
      rho_G1_grid = rho_G1_grid, rho_G2_grid = rho_G2_grid,
      omega_1_grid = if (omega_swept) omega_1_grid else NULL,
      omega_2_grid = if (omega_swept) omega_2_grid else NULL,
      mo_confounding_grid = if (length(unique(mo_confounding)) > 1) mo_confounding else NULL,
      phi_grid = if (length(unique(phi)) > 1) phi else NULL
    )
  )
  class(obj) <- "iconic_sensitivity"
  if (isTRUE(verbose))
    message("iconic_sensitivity complete. Call summary() or print() on the result for the full sensitivity summary.")
  obj
}


#' Safe max of absolute values, returning NA when all values are NA/NaN (internal)
#' @keywords internal
.safe_max_abs <- function(x) {
  x <- abs(x[is.finite(x)])
  if (length(x) == 0) NA_real_ else max(x)
}

#' Find tipping points from the degradation surface (internal)
#'
#' For each method, finds the first rho_G2 value (at rho_G1 = 0) where
#' |NDE_bias| or |NIE_bias| exceeds the threshold. This is the point
#' at which the estimator's assumptions are violated enough to produce
#' materially biased estimates.
#' @keywords internal
.find_tipping_points <- function(surface, rho_G1_grid, rho_G2_grid,
                                 threshold) {
  methods <- unique(surface$method)
  # When omega is swept, restrict tipping-point detection to the first
  # omega cell (the reference coverage) so the summary stays 2D.
  if ("omega_1" %in% names(surface) && "omega_2" %in% names(surface)) {
    o1_ref <- min(surface$omega_1, na.rm = TRUE)
    o2_ref <- min(surface$omega_2, na.rm = TRUE)
    surface <- surface[surface$omega_1 == o1_ref & surface$omega_2 == o2_ref, ]
  }
  rows <- lapply(methods, function(m) {
    sub <- surface[surface$method == m, ]
    # Look along the rho_G1 = 0 edge (the "pure Gm violation" axis)
    edge <- sub[sub$rho_G1 == 0, ]
    edge <- edge[order(edge$rho_G2), ]

    tip_NDE <- NA_real_
    tip_NIE <- NA_real_
    for (i in seq_len(nrow(edge))) {
      # guard against NaN/NA bias (e.g. COCA for survival outcomes).
      if (is.na(tip_NDE) && is.finite(edge$NDE_bias[i]) &&
          abs(edge$NDE_bias[i]) > threshold)
        tip_NDE <- edge$rho_G2[i]
      if (is.na(tip_NIE) && is.finite(edge$NIE_bias[i]) &&
          abs(edge$NIE_bias[i]) > threshold)
        tip_NIE <- edge$rho_G2[i]
    }

    # Also check along rho_G2 = 0 edge (the "pure G violation" axis)
    edge2 <- sub[sub$rho_G2 == 0, ]
    edge2 <- edge2[order(edge2$rho_G1), ]

    tip_NDE_G1 <- NA_real_
    for (i in seq_len(nrow(edge2))) {
      if (is.na(tip_NDE_G1) && is.finite(edge2$NDE_bias[i]) &&
          abs(edge2$NDE_bias[i]) > threshold)
        tip_NDE_G1 <- edge2$rho_G1[i]
    }

    data.frame(
      method = m,
      tip_rho_G2_NDE = tip_NDE,
      tip_rho_G2_NIE = tip_NIE,
      tip_rho_G1_NDE = tip_NDE_G1,
      max_NDE_bias = .safe_max_abs(sub$NDE_bias),
      max_NIE_bias = .safe_max_abs(sub$NIE_bias),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


#' Build sensitivity summary (internal)
#' @keywords internal
.build_sensitivity_summary <- function(surface, tipping, threshold, n_iter) {
  lines <- character(0)
  lines <- c(lines,
    sprintf(" Grid: %d x %d (rho_G1 x rho_G2), %d reps per cell",
            length(unique(surface$rho_G1)),
            length(unique(surface$rho_G2)),
            n_iter),
    sprintf(" Bias threshold for tipping: %.2f", threshold))

  # Per-method summary
  methods <- unique(surface$method)
  for (m in methods) {
    sub <- surface[surface$method == m, ]
    tip <- tipping[tipping$method == m, ]
    max_nde <- .safe_max_abs(sub$NDE_bias)
    max_nie <- .safe_max_abs(sub$NIE_bias)
    n_tipped <- sum(sub$tipped, na.rm = TRUE)

    tip_str <- ""
    if (!is.na(tip$tip_rho_G2_NDE))
      tip_str <- paste0(tip_str, sprintf(" NDE tips at rho_G2=%.1f", tip$tip_rho_G2_NDE))
    if (!is.na(tip$tip_rho_G1_NDE))
      tip_str <- paste0(tip_str, sprintf(" NDE tips at rho_G1=%.1f", tip$tip_rho_G1_NDE))
    if (nchar(tip_str) == 0)
      tip_str <- " no tipping within grid"

    lines <- c(lines, sprintf(" %s: max|NDE bias|=%.3f, max|NIE bias|=%.3f,%s",
                              m, max_nde, max_nie, tip_str))
  }

  paste(lines, collapse = "\n")
}

# Null-coalescing operator (internal)
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Print method for iconic_sensitivity objects
#'
#' @param x An \code{iconic_sensitivity} object.
#' @param ... Unused.
#' @return Invisibly returns `x` (the `iconic_sensitivity` object); called for its side effect of printing a human-readable summary.
#' @export
print.iconic_sensitivity <- function(x, ...) {
  cat("<iconic_sensitivity>\n")
  cat(x$summary, "\n")
  if (!is.null(x$texture_source))
    cat(" Texture:", x$texture_source, "\n")
  if (!is.null(x$inferred_confounding))
    cat(" Confounding: inferred from data\n")
  invisible(x)
}

#' Summary method for iconic_sensitivity objects
#'
#' Prints the full sensitivity summary (same as \code{print()}).
#' @param object An \code{iconic_sensitivity} object.
#' @param ... Unused.
#' @return Invisibly returns \code{object}.
#' @export
summary.iconic_sensitivity <- function(object, ...) {
  print.iconic_sensitivity(object, ...)
  invisible(object)
}
