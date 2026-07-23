# ============================================================
# TODO(v1.0): supplement — move technical sensitivity detail to appendix
# (JYH #478). The degradation-surface construction and tipping-point
# algorithm move to a technical appendix when the draft is condensed.
# No code change in v0.9.2.
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
# v0.7.0: Auto-trains a GAN (or MVN fallback) from the user's
# iconic_data object when no trained_gan is supplied, so the
# synthetic texture matches the user's cohort. Also supports
# confounding = "inferred" to calibrate confounding parameters
# from the data via infer_confounding().
# ============================================================

#' Auto-train a texture model from an iconic_data object (internal)
#'
#' Converts the user's iconic_data to the load_real_input_data() format
#' and trains a GAN (or multivariate-normal fallback) on the resulting
#' tidy frame.  The trained model supplies realistic covariate and
#' outcome texture to run_single_iteration(), and carries feature-level
#' residual correlation matrices for the Y, M, and W panels so the
#' simulation can inject correlated noise.
#'
#' @param data   An \code{iconic_data} object.
#' @param epochs GAN training epochs. Default 100.
#' @return An \code{iconic_gan} object.
#' @keywords internal
.auto_train_gan <- function(data, epochs = 100) {
  if (!inherits(data, "iconic_data"))
    stop("data must be an iconic_data object.")

  # Convert iconic_data → load_real_input_data format.
  # Z is a length-n vector; load_real_input_data expects a features x
  # samples matrix, so wrap it as 1 x n.
  Z_mat <- matrix(data$Z, nrow = 1L, ncol = data$n)
  Y_mat <- data$Y  # already features x samples in iconic_data

  # Mediator panel: iconic_data stores M as features x samples (or NULL)
  M_mat <- data$M

  # Negative-control panel: use W (the combined panel) for correlation
  # learning.  When path-specific W1/W2 are present, W is the combined
  # panel; either way, the correlation structure is what we need.
  W_mat <- data$W

  cov_df <- if (!is.null(data$covariates) && ncol(data$covariates) > 0)
    data$covariates else NULL

  input <- load_real_input_data(
    Z_matrix     = Z_mat,
    Y_matrix     = Y_mat,
    M_matrix     = M_mat,
    W_matrix     = W_mat,
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
#' @param data        An \code{iconic_data} object.
#' @param epochs      Epochs for auto-training. Default 100.
#' @return A list with \code{gan} (the resolved iconic_gan) and
#'   \code{source} (a string: "user-supplied GAN", "data-attached GAN",
#'   or "auto-trained from data").
#' @keywords internal
.resolve_gan <- function(trained_gan, data, epochs = 100) {
  if (!is.null(trained_gan)) {
    return(list(gan = trained_gan, source = "user-supplied GAN"))
  }
  if (!is.null(data$trained_gan)) {
    return(list(gan = data$trained_gan, source = "data-attached GAN"))
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
#' @section Texture model (v0.7.0):
#' When \code{trained_gan} is \code{NULL} and no GAN is attached to
#' \code{data}, a texture model is auto-trained from the user's data
#' (GAN via \code{torch} if available, otherwise a multivariate-normal
#' fallback).  This ensures the synthetic covariate and outcome texture
#' matches the user's cohort.  Supply \code{trained_gan} explicitly or
#' attach one via \code{\link{iconic_data}(trained_gan = ...)} to reuse
#' a pre-trained model and avoid retraining.
#'
#' @section Confounding calibration (v0.7.0):
#' The \code{confounding} argument controls how the held-fixed
#' confounding parameters (\code{mo_confounding}, \code{omega_1},
#' \code{omega_2}) are set:
#' \itemize{
#'   \item \code{"default"}: fixed defaults (0.8, 0.7, 0.7).
#'   \item \code{"inferred"}: calls \code{\link{infer_confounding}()} to
#'     estimate parameters from the data.  Parameters that cannot be
#'     inferred fall back to defaults with warnings.
#'   \item \code{"manual"}: use the explicitly supplied arguments.
#' }
#'
#' @param data      An \code{iconic_data} object.
#' @param diagnosis Optional \code{iconic_diagnosis} from
#'   \code{\link{iconic_diagnose}()}.  Used to infer \code{phi} and
#'   \code{mo_confounding} when not explicitly supplied.
#' @param trained_gan Optional \code{iconic_gan} from
#'   \code{\link{train_gan_on_real_data}()}.  If \code{NULL}, a texture
#'   model is auto-trained from \code{data} (v0.7.0).
#' @param confounding Confounding parameter source: \code{"default"}
#'   (fixed defaults), \code{"inferred"} (data-calibrated via
#'   \code{\link{infer_confounding}()}), or \code{"manual"} (use
#'   explicitly supplied arguments).  Default \code{"default"}.
#' @param gan_epochs Epochs for auto-trained GAN. Default 100.
#' @param rho_G1_grid Values of rho_G1 (G correlation with U_XM).
#'   Default \code{c(0, 0.1, 0.2, 0.3, 0.5)}.
#' @param rho_G2_grid Values of rho_G2 (Gm correlation with U_MY).
#'   Default \code{c(0, 0.1, 0.2, 0.3, 0.5)}.
#' @param n_iter     Replicates per grid cell. Default 30.
#' @param n_samples  Samples per replicate. If NULL, uses \code{data$n}.
#' @param n_features Features per replicate. Default 10.
#' @param phi        Mediator instrument strength. If NULL, inferred
#'   from diagnosis (F_Gm) or defaults to 0.8.
#' @param mo_confounding M-O confounding strength. Default 0.8.
#'   Used when \code{confounding = "default"} or \code{"manual"}.
#' @param separate_U Use separate confounders for Z->M and M->Y paths.
#'   Default TRUE (matches v0.5.0 DGP).
#' @param omega_1,omega_2 NC coverage of U_XM / U_MY. Default 0.7.
#'   Used when \code{confounding = "default"} or \code{"manual"}.
#' @param bias_threshold Absolute bias threshold for tipping-point
#'   annotation. Default 0.10.
#' @param base_seed  Base RNG seed. Default 700.
#' @param n_cores    Number of parallel workers for simulation replicates.
#'   Default 1 (sequential).  Uses \code{parallel::mclapply} on Unix
#'   and a PSOCK cluster on Windows.
#'
#' @section Defaults:
#' \tabular{lll}{
#'   \strong{Parameter} \tab \strong{Default} \tab \strong{Source} \cr
#'   \code{confounding} \tab "default" \tab Use DGP defaults below \cr
#'   \code{gan_epochs} \tab 100 \tab Texture-model training budget \cr
#'   \code{rho_G1_grid} \tab c(0,0.1,0.2,0.3,0.5) \tab Independence-violation sweep \cr
#'   \code{rho_G2_grid} \tab c(0,0.1,0.2,0.3,0.5) \tab Independence-violation sweep \cr
#'   \code{n_iter} \tab 30 \tab Replicates per grid cell \cr
#'   \code{n_features} \tab 10 \tab Features per replicate \cr
#'   \code{mo_confounding} \tab 0.8 \tab Simulation calibration (delta_mo) \cr
#'   \code{separate_U} \tab TRUE \tab Path-specific confounders \cr
#'   \code{omega_1, omega_2} \tab 0.7 \tab NC coverage (simulation calibration) \cr
#'   \code{bias_threshold} \tab 0.10 \tab Tipping-point threshold \cr
#' }
#'
#' @return An \code{iconic_sensitivity} S3 object: a named list with
#'   \code{$surface} (data frame: rho_G1, rho_G2, method, NDE_bias,
#'   NIE_bias, NDE_rmse, NIE_rmse, NDE_type1, NIE_type1, tipped),
#'   \code{$grid}, \code{$tipping_points}, \code{$summary},
#'   \code{$texture_source} (how the GAN was obtained), and
#'   \code{$inferred_confounding} (when \code{confounding = "inferred"}).
#' @export
#'
#' @examples
#' \dontrun{
#' data <- iconic_data(Z = rnorm(100), Y = matrix(rnorm(100*10), 10, 100),
#'                     M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
#'                     W = matrix(rnorm(100*10), 10, 100))
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
                               separate_U = TRUE,
                               omega_1 = 0.7, omega_2 = 0.7,
                               bias_threshold = 0.10,
                               base_seed = 700,
                               n_cores = 1) {
  if (!inherits(data, "iconic_data"))
    stop("data must be an iconic_data object from iconic_data().")
  if (!data$is_mediation)
    stop("iconic_sensitivity requires mediation data (supply M). ",
         "Use gan_sensitivity() for total-effect-only sensitivity.")

  confounding <- match.arg(confounding)

  # ── Resolve the texture model (v0.7.0) ──
  gan_res <- .resolve_gan(trained_gan, data, epochs = gan_epochs)
  gan <- gan_res$gan
  texture_source <- gan_res$source

  # ── Resolve confounding parameters (v0.7.0) ──
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

  # Build the 2D sweep
  grid <- expand.grid(rho_G1 = rho_G1_grid, rho_G2 = rho_G2_grid,
                      KEEP.OUT.ATTRS = FALSE)

  true_NDE <- 0.10  # beta_Z default
  true_NIE <- 0.15  # alpha_M * beta_M default

  n_grid <- nrow(grid)
  smry <- lapply(seq_len(n_grid), function(gi) {
    r1 <- grid$rho_G1[gi]
    r2 <- grid$rho_G2[gi]
    if (n_grid > 1)
      message("Sensitivity grid cell ", gi, "/", n_grid,
              " (rho_G1=", r1, ", rho_G2=", r2, ")")

    worker <- function(i) {
      dat <- run_single_iteration(
        trained_gan = gan,
        n_synthetic_samples = n_samples, n_features = n_features,
        mo_confounding = mo_confounding, phi = phi,
        rho_G1 = r1, rho_G2 = r2, rho_pop = 0,
        separate_U = separate_U, omega_1 = omega_1, omega_2 = omega_2,
        seed = base_seed + gi * 1000L + i)
      res <- run_mediation_methods(dat, n_features)
      res$iter <- i
      res
    }

    combined <- do.call(rbind, .parallel_lapply(seq_len(n_iter), worker,
                                                n_cores = n_cores,
                                                progress = "  Replicates"))
    s <- summarise_mediation_results(combined, true_NDE, true_NIE)
    s$rho_G1 <- r1
    s$rho_G2 <- r2
    s
  })

  surface <- do.call(rbind, smry)

  # Tipping-point annotation: flag cells where |bias| exceeds threshold
  surface$tipped_NDE <- abs(surface$NDE_bias) > bias_threshold
  surface$tipped_NIE <- abs(surface$NIE_bias) > bias_threshold
  surface$tipped <- surface$tipped_NDE | surface$tipped_NIE

  # Reorder columns
  front <- c("rho_G1", "rho_G2", "method")
  surface <- surface[, c(front, setdiff(names(surface), front))]

  # Find tipping points: first rho_G2 (along each rho_G1 row) where a
  # method tips, for each method
  tipping <- .find_tipping_points(surface, rho_G1_grid, rho_G2_grid,
                                  bias_threshold)

  # Summary
  summary_txt <- .build_sensitivity_summary(surface, tipping, bias_threshold,
                                            n_iter)

  obj <- list(
    surface             = surface,
    grid                = grid,
    tipping_points      = tipping,
    bias_threshold      = bias_threshold,
    n_iter              = n_iter,
    n_samples           = n_samples,
    phi                 = phi,
    mo_confounding      = mo_confounding,
    omega_1             = omega_1,
    omega_2             = omega_2,
    texture_source      = texture_source,
    inferred_confounding = inferred_conf,
    summary             = summary_txt,
    # v0.9.2 (JYH #543, #582): scenario manifest for reader orientation.
    manifest            = scenario_manifest(
      list(beta_Z = true_NDE, alpha_M = 0.50, beta_M = 0.30,
           n_mediators = 1, n = n_samples, n_features = n_features,
           separate_U = separate_U, feat_cor = 0,
           conf_str = NA, mo_confounding = mo_confounding, phi = phi,
           w_signal = omega_1, rho_G1 = NA, rho_G2 = NA, rho_pop = NA),
      rho_G1_grid = rho_G1_grid, rho_G2_grid = rho_G2_grid,
      mo_confounding_grid = if (length(unique(mo_confounding)) > 1) mo_confounding else NULL,
      phi_grid = if (length(unique(phi)) > 1) phi else NULL
    )
  )
  class(obj) <- c("iconic_sensitivity", "list")
  obj
}


#' Find tipping points from the degradation surface (internal)
#'
#' For each method, finds the first rho_G2 value (at rho_G1 = 0) where
#' |NDE_bias| or |NIE_bias| exceeds the threshold.  This is the point
#' at which the estimator's assumptions are violated enough to produce
#' materially biased estimates.
#' @keywords internal
.find_tipping_points <- function(surface, rho_G1_grid, rho_G2_grid,
                                 threshold) {
  methods <- unique(surface$method)
  rows <- lapply(methods, function(m) {
    sub <- surface[surface$method == m, ]
    # Look along the rho_G1 = 0 edge (the "pure Gm violation" axis)
    edge <- sub[sub$rho_G1 == 0, ]
    edge <- edge[order(edge$rho_G2), ]

    tip_NDE <- NA_real_
    tip_NIE <- NA_real_
    for (i in seq_len(nrow(edge))) {
      if (is.na(tip_NDE) && abs(edge$NDE_bias[i]) > threshold)
        tip_NDE <- edge$rho_G2[i]
      if (is.na(tip_NIE) && abs(edge$NIE_bias[i]) > threshold)
        tip_NIE <- edge$rho_G2[i]
    }

    # Also check along rho_G2 = 0 edge (the "pure G violation" axis)
    edge2 <- sub[sub$rho_G2 == 0, ]
    edge2 <- edge2[order(edge2$rho_G1), ]

    tip_NDE_G1 <- NA_real_
    for (i in seq_len(nrow(edge2))) {
      if (is.na(tip_NDE_G1) && abs(edge2$NDE_bias[i]) > threshold)
        tip_NDE_G1 <- edge2$rho_G1[i]
    }

    data.frame(
      method = m,
      tip_rho_G2_NDE = tip_NDE,
      tip_rho_G2_NIE = tip_NIE,
      tip_rho_G1_NDE = tip_NDE_G1,
      max_NDE_bias = max(abs(sub$NDE_bias), na.rm = TRUE),
      max_NIE_bias = max(abs(sub$NIE_bias), na.rm = TRUE),
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
    sprintf("  Grid: %d x %d (rho_G1 x rho_G2), %d reps per cell",
            length(unique(surface$rho_G1)),
            length(unique(surface$rho_G2)),
            n_iter),
    sprintf("  Bias threshold for tipping: %.2f", threshold))

  # Per-method summary
  methods <- unique(surface$method)
  for (m in methods) {
    sub <- surface[surface$method == m, ]
    tip <- tipping[tipping$method == m, ]
    max_nde <- max(abs(sub$NDE_bias), na.rm = TRUE)
    max_nie <- max(abs(sub$NIE_bias), na.rm = TRUE)
    n_tipped <- sum(sub$tipped, na.rm = TRUE)

    tip_str <- ""
    if (!is.na(tip$tip_rho_G2_NDE))
      tip_str <- paste0(tip_str, sprintf(" NDE tips at rho_G2=%.1f", tip$tip_rho_G2_NDE))
    if (!is.na(tip$tip_rho_G1_NDE))
      tip_str <- paste0(tip_str, sprintf(" NDE tips at rho_G1=%.1f", tip$tip_rho_G1_NDE))
    if (nchar(tip_str) == 0)
      tip_str <- " no tipping within grid"

    lines <- c(lines, sprintf("  %s: max|NDE bias|=%.3f, max|NIE bias|=%.3f,%s",
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
#' @export
print.iconic_sensitivity <- function(x, ...) {
  cat("<iconic_sensitivity>\n")
  cat(x$summary, "\n")
  if (!is.null(x$texture_source))
    cat("  Texture:", x$texture_source, "\n")
  if (!is.null(x$inferred_confounding))
    cat("  Confounding: inferred from data\n")
  invisible(x)
}
