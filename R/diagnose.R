# ============================================================
# iconic_diagnose: Diagnostic eligibility layer for the
# model selection workflow.
#
# Runs instrument-strength checks, negative-control validity
# screens, and path-completeness assessments on an iconic_data
# object, then returns an eligibility report for all 8 estimators.
#
# Adapts the existing nc_validity_screen(), nc_independence_check(),
# nc_independence_check_gm(), and nc_completeness_check() functions
# from nc_diagnostics.R to operate on iconic_data objects.
# ============================================================

#' Diagnose data and determine estimator eligibility
#'
#' Runs a battery of diagnostic checks on the user's data and returns
#' an eligibility report for all eight ICONIC estimators. The checks
#' include:
#'
#' \itemize{
#' \item \strong{Instrument strength}: first-stage partial F-statistics
#' for G (exposure instrument) and Gm (mediator instrument), with
#' the Stock-Yogo weak-instrument threshold (F >= 10).
#' \item \strong{NC validity (A1)}: screens each negative-control
#' feature for association with the exposure Z
#' (\code{\link{nc_validity_screen}()}).
#' \item \strong{NC independence (A2)}: screens each NC for association
#' with the instrument G (\code{\link{nc_independence_check}()}).
#' \item \strong{NC independence (A2')}: screens each NC for association
#' with the mediator instrument Gm
#' (\code{\link{nc_independence_check_gm}()}).
#' \item \strong{Path completeness}: checks whether the valid NC panel
#' has enough features to span the confounder subspace
#' (\code{\link{nc_completeness_check}()}).
#' }
#'
#' @param data An \code{iconic_data} object from \code{\link{iconic_data}()}.
#' @param fdr_level FDR level for NC screens. Default 0.10.
#' @param min_f Minimum partial F for instrument strength. Default 10.
#' Used as the scalar threshold when \code{g_threshold}/\code{gm_threshold}
#' are NULL (legacy behavior). When thresholds are supplied, \code{min_f}
#' is unused for the panel decision (the threshold's \code{R} governs).
#' @param k Number of latent confounders assumed for the
#' completeness check. Default 1 (a single confounder
#' is the typical case in mediation settings).
#' @param g_threshold Optional list \code{list(E = 0.5, R = 10)} controlling
#' G-dependent method eligibility via the panel distribution. A method is
#' eligible if at least fraction \code{E} of instruments have F_G >= \code{R}.
#' Default NULL: legacy scalar behavior (median F_G vs \code{min_f}).
#' @param gm_threshold Optional list \code{list(E = 0.5, R = 10)} controlling
#' Gm-dependent method eligibility (IV2SLS2, PGC2Gm) via the panel
#' distribution. A method is eligible if at least fraction \code{E} of
#' mediators have F_Gm >= \code{R}. Default NULL: legacy scalar behavior
#' (median F_Gm vs \code{min_f}).
#' @param n_cores Number of parallel workers for NC validity screens and
#' panel instrument-strength computation. Default 1 (sequential). Uses
#' \code{parallel::mclapply} on Unix and a PSOCK cluster on Windows.
#' @param allow_no_proxy Logical: when \code{TRUE}
#' (default), proceed with a message if no instruments and no NCs are
#' supplied (only UNADJ is eligible). When \code{FALSE}, error instead.
#'
#' @section Defaults:
#' \tabular{lll}{
#' \strong{Parameter} \tab \strong{Default} \tab \strong{Source} \cr
#' \code{fdr_level} \tab 0.10 \tab BH-FDR level for NC screens \cr
#' \code{min_f} \tab 10 \tab Stock-Yogo weak-instrument threshold \cr
#' \code{k} \tab 1 \tab Single-confounder assumption (typical mediation) \cr
#' \code{g_threshold} \tab NULL \tab Legacy scalar (median F_G vs min_f) \cr
#' \code{gm_threshold} \tab NULL \tab Legacy scalar (median F_Gm vs min_f) \cr
#' \code{n_cores} \tab 1 \tab Sequential \cr
#' \code{allow_no_proxy} \tab TRUE \tab Proceed with UNADJ-only when no IV/NC \cr
#' }
#'
#' @return An \code{iconic_diagnosis} S3 object: a named list with
#' \code{$instrument_strength}, \code{$nc_validity},
#' \code{$nc_independence}, \code{$nc_independence_gm},
#' \code{$completeness}, \code{$eligibility}, and \code{$summary}.
#' \code{$instrument_strength$F_G} and \code{$instrument_strength$F_Gm}
#' are numeric vectors containing the full panel distributions (one entry
#' per instrument / mediator).
#' @export
#'
#' @examples
#' set.seed(1)
#' data <- iconic_data(Z = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
#' G = rnorm(100), W = matrix(rnorm(100 * 10), 10, 100))
#' diag <- iconic_diagnose(data)
#' print(diag)
iconic_diagnose <- function(data, fdr_level = 0.10, min_f = 10, k = 1,
                            g_threshold = NULL, gm_threshold = NULL,
                            n_cores = 1, allow_no_proxy = TRUE) {
  if (!inherits(data, "iconic_data"))
    stop("data must be an iconic_data object from iconic_data().")

  # explicit control over proceed-without-proxy behavior.
  # When allow_no_proxy = FALSE, error if neither instruments nor NCs are
  # present; when TRUE (default), proceed with a message noting only UNADJ
  # is eligible.
  has_iv <- !is.null(data$G) || !is.null(data$Gm)
  has_nc <- !is.null(data$W)
  if (!has_iv && !has_nc) {
    if (!allow_no_proxy) {
      stop("No genetic instruments (G/Gm) and no negative controls (W) ",
           "supplied. Set allow_no_proxy = TRUE to proceed with UNADJ-only ",
           "estimation, or supply instruments/NCs for the core estimators.")
    }
    message("No instruments or negative controls supplied: only UNADJ ",
            "estimation is eligible. Set allow_no_proxy = FALSE to make ",
            "this an error.")
  }

  # ── Instrument strength ──
  inst <- .check_instrument_strength(data, min_f = min_f,
                                     g_threshold = g_threshold,
                                     gm_threshold = gm_threshold,
                                     n_cores = n_cores)

  # ── NC validity screens ──
  nc_val <- NULL
  nc_ind <- NULL
  nc_ind_gm <- NULL
  completeness <- NULL

  if (data$has_nc) {
    nc_dat <- .to_nc_dat(data) # bridge to nc_diagnostics format
    nc_val <- nc_validity_screen(nc_dat, fdr_level = fdr_level, n_cores = n_cores)
    if (data$has_instrument)
      nc_ind <- nc_independence_check(nc_dat, fdr_level = fdr_level, n_cores = n_cores)
    if (data$has_mediator_instrument)
      nc_ind_gm <- nc_independence_check_gm(nc_dat, fdr_level = fdr_level, n_cores = n_cores)

    # Completeness
    n_valid <- .count_valid_ncs(nc_val, nc_ind)
    # COCA-specific valid count (A2 exempt per).
    # COCA's identifying assumption does not involve the instrument, so it
    # should not be gated on A2. Because the FDR-based A1 screen flags valid
    # NCs that share U with Z (the "always significant via U" problem), the
    # COCA count uses a magnitude-based A1 with a high threshold that flags
    # only gross Z-dependence (downstream-of-Z violations), not the intended
    # confounder-sharing signal.
    nc_val_coca <- nc_validity_screen(nc_dat, fdr_level = fdr_level,
                                      n_cores = n_cores,
                                      criterion = "magnitude",
                                      magnitude_threshold = 0.75)
    n_valid_coca <- .count_valid_ncs(nc_val_coca, nc_ind, for_estimator = "COCA")
    W_for_comp <- data$W
    if (is.null(W_for_comp) && data$has_path_nc)
      W_for_comp <- data$W1
    completeness <- .assess_completeness(n_valid, k, ncol(W_for_comp))
    completeness$n_valid_coca <- n_valid_coca
  }

  # ── Eligibility ──
  eligibility <- .assess_eligibility(data, inst, completeness, min_f)

  # ── Summary ──
  summary_txt <- .build_diagnosis_summary(data, inst, nc_val, nc_ind,
                                          nc_ind_gm, completeness,
                                          eligibility)

  obj <- list(
    instrument_strength = inst,
    nc_validity = nc_val,
    nc_independence = nc_ind,
    nc_independence_gm = nc_ind_gm,
    completeness = completeness,
    eligibility = eligibility,
    summary = summary_txt,
    min_f = min_f,
    k = k,
    fdr_level = fdr_level
  )
  class(obj) <- c("iconic_diagnosis", "list")
  obj
}


#' Check instrument strength via first-stage partial F (internal)
#'
#' Computes partial F-statistics for G (in Z ~ G + W + covars) and
#' Gm (in M ~ Z_hat + Gm + W + covars), matching the first-stage
#' regressions used by fit_iv2sls() and fit_iv2sls_mediation2().
#'
#' @param data iconic_data object.
#' @param min_f Minimum acceptable partial F. Default 10.
#' @return List: F_G, F_Gm, weak_G (logical), weak_Gm (logical).
#' @keywords internal
.check_instrument_strength <- function(data, min_f = 10, g_threshold = NULL,
                                       gm_threshold = NULL, n_cores = 1) {
  n <- data$n
  cv <- data$covariates
  cnames <- if (!is.null(cv) && ncol(cv) > 0) names(cv) else character(0)
  cs <- .covar_str(cnames)
  Z <- data$Z
  G <- data$G

  # Helpers for panel-distribution summaries
  panel_summary <- function(F_vec, R) {
    F_vec <- F_vec[!is.na(F_vec)]
    if (!length(F_vec)) {
      return(list(min = NA_real_, median = NA_real_, mean = NA_real_,
                  p_pass = NA_real_, n_pass = 0L, n_total = 0L))
    }
    n_pass <- sum(F_vec >= R)
    list(
      min = min(F_vec),
      median = median(F_vec),
      mean = mean(F_vec),
      p_pass = n_pass / length(F_vec),
      n_pass = n_pass,
      n_total = length(F_vec)
    )
  }
  # Resolve threshold decision: NULL -> legacy scalar (median vs min_f);
  # list(E, R) -> panel fraction passing R vs E.
  threshold_decision <- function(F_vec, threshold, min_f) {
    s <- panel_summary(F_vec, if (!is.null(threshold)) threshold$R else min_f)
    if (is.null(threshold)) {
      weak <- is.na(s$median) || s$median < min_f
    } else {
      weak <- is.na(s$p_pass) || s$p_pass < threshold$E
    }
    list(weak = weak, summary = s)
  }

  # ── G (exposure instrument) ──
  F_G_vec <- NA_real_
  weak_G <- TRUE
  fs <- NULL

  if (data$has_instrument) {
    # First stage: Z ~ G + W + covars (use W_avg as scalar W)
    # Fall back to W1 if W is absent but path-specific NCs are supplied.
    W_for_fs <- data$W
    if (is.null(W_for_fs) && data$has_path_nc)
      W_for_fs <- data$W1
    w_avg <- if (!is.null(W_for_fs)) colMeans(W_for_fs) else NULL

    # G may be a scalar vector (common) or a matrix (multiple instruments).
    # For scalar G, F_G_vec is length-1. For matrix G, compute F for each
    # instrument column (panel distribution).
    G_mat <- if (is.matrix(G)) G else matrix(G, ncol = 1)
    F_G_vec <- .parallel_lapply(seq_len(ncol(G_mat)), function(j) {
      d_fs <- .bind_covars(data.frame(Z = Z, g = G_mat[, j], w = w_avg), cv)
      fit <- tryCatch(lm(as.formula(paste0("Z ~ g + w", cs)), data = d_fs),
                      error = function(e) NULL)
      if (is.null(fit)) NA_real_ else .partial_F(fit, "g")
    }, n_cores = n_cores)
    F_G_vec <- as.numeric(F_G_vec)
    if (is.null(dimnames(G_mat)[[2]]) || dimnames(G_mat)[[2]][1] == "")
      names(F_G_vec) <- paste0("G", seq_along(F_G_vec))
    else
      names(F_G_vec) <- dimnames(G_mat)[[2]]

    # Keep fs (first instrument's fit) for Z_hat used in the Gm stage
    d_fs0 <- .bind_covars(data.frame(Z = Z, g = G_mat[, 1], w = w_avg), cv)
    fs <- tryCatch(lm(as.formula(paste0("Z ~ g + w", cs)), data = d_fs0),
                   error = function(e) NULL)

    g_dec <- threshold_decision(F_G_vec, g_threshold, min_f)
    weak_G <- g_dec$weak
  } else {
    g_dec <- list(summary = panel_summary(NA_real_, min_f))
  }

  # ── Gm (mediator instrument) ──
  F_Gm_vec <- NA_real_
  weak_Gm <- TRUE

  if (data$has_mediator_instrument && data$is_mediation) {
    # Need Z_hat from the G first stage, then M ~ Z_hat + Gm + W + covars
    if (!is.null(fs)) {
      Z_hat <- fitted(fs)
    } else {
      # Fallback: regress Z on G alone
      d_fs0 <- .bind_covars(data.frame(Z = Z, g = G), cv)
      fs0 <- tryCatch(lm(as.formula(paste0("Z ~ g", cs)), data = d_fs0),
                      error = function(e) NULL)
      Z_hat <- if (!is.null(fs0)) fitted(fs0) else Z
    }

    W_for_ms <- data$W
    if (is.null(W_for_ms) && data$has_path_nc)
      W_for_ms <- data$W1
    w_avg <- if (!is.null(W_for_ms)) colMeans(W_for_ms) else NULL

    # Gm: n_mediators x n matrix (rows = mediators). Compute F_Gm for each.
    Gm_mat <- if (is.matrix(data$Gm)) data$Gm else matrix(data$Gm, nrow = 1)
    M_mat <- if (is.matrix(data$M)) data$M else matrix(data$M, nrow = 1)
    nm <- nrow(Gm_mat)
    F_Gm_vec <- .parallel_lapply(seq_len(nm), function(m) {
      d_ms <- .bind_covars(data.frame(M = M_mat[m, ], Z_hat = Z_hat,
                                      gm = Gm_mat[m, ], w = w_avg), cv)
      fit <- tryCatch(lm(as.formula(paste0("M ~ Z_hat + gm + w", cs)),
                         data = d_ms), error = function(e) NULL)
      if (is.null(fit)) NA_real_ else .partial_F(fit, "gm")
    }, n_cores = n_cores)
    F_Gm_vec <- as.numeric(F_Gm_vec)
    mnames <- data$mediator_names
    if (!is.null(mnames) && length(mnames) == nm)
      names(F_Gm_vec) <- mnames
    else
      names(F_Gm_vec) <- paste0("M", seq_along(F_Gm_vec))

    gm_dec <- threshold_decision(F_Gm_vec, gm_threshold, min_f)
    weak_Gm <- gm_dec$weak
  } else {
    gm_dec <- list(summary = panel_summary(NA_real_, min_f))
  }

  gs <- g_dec$summary
  gms <- gm_dec$summary

  list(
    F_G = F_G_vec,
    F_G_panel = F_G_vec, # alias for clarity
    F_G_min = gs$min,
    F_G_median = gs$median,
    F_G_mean = gs$mean,
    F_G_p_pass = gs$p_pass,
    F_G_n_pass = gs$n_pass,
    F_G_n_total = gs$n_total,
    F_Gm = F_Gm_vec,
    F_Gm_panel = F_Gm_vec, # alias for clarity
    F_Gm_min = gms$min,
    F_Gm_median = gms$median,
    F_Gm_mean = gms$mean,
    F_Gm_p_pass = gms$p_pass,
    F_Gm_n_pass = gms$n_pass,
    F_Gm_n_total = gms$n_total,
    weak_G = weak_G,
    weak_Gm = weak_Gm,
    min_f = min_f,
    g_threshold = g_threshold,
    gm_threshold = gm_threshold
  )
}


#' Bridge iconic_data to nc_diagnostics dat format (internal)
#'
#' The existing nc_validity_screen / nc_independence_check functions
#' expect a `dat` list with $W (n x q), $Z, $G (n x p), $Gm, and
#' $synthetic_data. This helper converts an iconic_data object to
#' that format.
#' @keywords internal
.to_nc_dat <- function(data) {
  # Use W if available; fall back to W1 when only path-specific NCs supplied
  W_mat <- data$W
  if (is.null(W_mat) && data$has_path_nc)
    W_mat <- data$W1
  W <- if (!is.null(W_mat)) t(W_mat) else NULL # q x n -> n x q
  list(
    W = W,
    Z = data$Z,
    G = if (!is.null(data$G)) matrix(data$G, ncol = 1) else NULL,
    Gm = if (!is.null(data$Gm)) {
      if (is.matrix(data$Gm)) {
        if (ncol(data$Gm) == data$n) data$Gm[1, ] else as.numeric(data$Gm[1, ])
      } else as.numeric(data$Gm)
    } else NULL,
    synthetic_data = data$covariates
  )
}


#' Count valid NCs from screen results (internal)
#'
#' Uses only the A2 (instrument-independence) screen as a hard gate,
#' since A1 (W ~ Z | C) cannot distinguish "W downstream of Z" (true
#' violation) from "W shares a cause with Z" (intended NC behavior).
#' A1 is reported but not used for eligibility.
#'
#' `for_estimator` controls which screen is used as the gate,
#' implementing the COCA A2 exemption. COCA's identifying
#' assumption does not involve the instrument, so COCA uses A1 only;
#' IV2SLS/PGC/IV2SLS2/PGC2/PGC2Gm use A2 (the legacy behavior).
#' @keywords internal
.count_valid_ncs <- function(nc_val, nc_ind, for_estimator = NULL) {
  # COCA does not require A2 (instrument-independence): its identifying
  # assumption is the negative-control outcome structure, not the instrument.
  # When for_estimator == "COCA", gate on A1 only.
  if (!is.null(for_estimator) && for_estimator == "COCA") {
    if (!is.null(nc_val))
      return(sum(!nc_val$significant, na.rm = TRUE))
    # No A1 screen — fall back to A2 if available, else 0.
    if (!is.null(nc_ind))
      return(sum(!nc_ind$significant, na.rm = TRUE))
    return(0)
  }
  # Default (IV2SLS, PGC, IV2SLS2, PGC2, PGC2Gm): A2 is the meaningful gate.
  if (!is.null(nc_ind)) {
    return(sum(!nc_ind$significant, na.rm = TRUE))
  }
  # No A2 screen (no instrument) — fall back to A1
  if (!is.null(nc_val))
    return(sum(!nc_val$significant, na.rm = TRUE))
  0
}


#' Assess completeness (internal)
#' @keywords internal
.assess_completeness <- function(n_valid, k, dim_W) {
  status <- if (n_valid > k) {
    "satisfied"
  } else if (n_valid == k) {
    "borderline"
  } else {
    "under-identified"
  }
  list(
    n_valid_controls = n_valid,
    k = k,
    dim_W = dim_W,
    completeness = status
  )
}


#' Assess estimator eligibility from diagnostics (internal)
#'
#' Applies the eligibility rules:
#' UNADJ: always eligible
#' DIRECT: requires G + W
#' COCA: requires W, A1 passed (at least some valid NCs); A2 NOT required
#' IV2SLS: requires G, F_G >= min_f, W
#' PGC: requires G, F_G >= min_f, completeness satisfied
#' IV2SLS2: requires G + Gm + W, F_G >= min_f, F_Gm >= min_f
#' PGC2: requires G + W1 + W2, F_G >= min_f, path completeness
#' PGC2Gm: requires G + Gm + W1 + W2, F_G >= min_f, path completeness
#'
#' COCA no longer requires the A2 (instrument-independence) screen,
#' only A1 + completeness. The eligibility table gains an
#' `a2_required` column documenting which estimators need A2.
#' @keywords internal
.assess_eligibility <- function(data, inst, completeness, min_f) {
  has_G <- data$has_instrument
  has_Gm <- data$has_mediator_instrument
  has_W <- data$has_nc
  has_W1W2 <- data$has_path_nc
  is_med <- data$is_mediation

  F_G_ok <- !inst$weak_G
  F_Gm_ok <- !inst$weak_Gm

  # Completeness
  comp_ok <- if (!is.null(completeness))
    completeness$completeness != "under-identified" else FALSE
  has_valid_nc <- if (!is.null(completeness))
    completeness$n_valid_controls > 0 else has_W
  # COCA uses the A1-only valid count (A2 exempt).
  has_valid_nc_coca <- if (!is.null(completeness))
    completeness$n_valid_coca > 0 else has_W

  # Build eligibility table
  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
               "IV2SLS2", "PGC2", "PGC2Gm")

  eligible <- c(
    TRUE, # UNADJ
    has_G && has_W, # DIRECT
    has_W && has_valid_nc_coca, # COCA (A2 exempt)
    has_G && has_W && F_G_ok, # IV2SLS
    has_G && has_W && F_G_ok && comp_ok, # PGC
    is_med && has_G && has_Gm && has_W && F_G_ok && F_Gm_ok, # IV2SLS2
    is_med && has_G && has_W1W2 && F_G_ok && comp_ok, # PGC2
    is_med && has_G && has_Gm && has_W1W2 && F_G_ok && comp_ok # PGC2Gm
  )

  # document which estimators require A2 (instrument-independence).
  a2_required <- c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)

  # Build reason strings
  # F_G / F_Gm are now vectors (panel distributions). For concise reason
  # strings, report the median (legacy scalar) or the panel decision when
  # a threshold was supplied.
  fmt_F_G <- function() {
    if (is.null(inst$g_threshold)) {
      sprintf("F_G=%.1f", inst$F_G_median)
    } else {
      sprintf("%.0f%% of %d transcripts F_G>=%g (E threshold: %.0f%%)",
              100 * inst$F_G_p_pass, inst$F_G_n_total,
              inst$g_threshold$R, 100 * inst$g_threshold$E)
    }
  }
  fmt_F_Gm <- function() {
    if (is.null(inst$gm_threshold)) {
      sprintf("F_Gm=%.1f", inst$F_Gm_median)
    } else {
      sprintf("%.0f%% of %d transcripts F_Gm>=%g (E threshold: %.0f%%)",
              100 * inst$F_Gm_p_pass, inst$F_Gm_n_total,
              inst$gm_threshold$R, 100 * inst$gm_threshold$E)
    }
  }

  reasons <- character(8)
  reasons[1] <- "always eligible"
  reasons[2] <- if (eligible[2]) "G + W present" else
    paste0("requires G + W (have: G=", has_G, ", W=", has_W, ")")
  reasons[3] <- if (eligible[3]) "W present, valid NCs available (A2 exempt)" else
    paste0("requires W with valid NCs (have W=", has_W, ", valid(A1)=", has_valid_nc_coca, ")")
  reasons[4] <- if (eligible[4]) "G + W present, F_G ok" else
    paste0("requires G + W + F_G>=", min_f, " (", fmt_F_G(), ")")
  reasons[5] <- if (eligible[5]) "G + W present, F_G ok, completeness ok" else
    paste0("requires completeness (status: ",
           ifelse(is.null(completeness), "unknown", completeness$completeness), ")")
  reasons[6] <- if (eligible[6]) "G + Gm + W present, both F ok" else
    paste0("requires G + Gm + W + both F ok (", fmt_F_G(), ", ", fmt_F_Gm(), ")")
  reasons[7] <- if (eligible[7]) "G + W1/W2 present, F_G ok, completeness ok" else
    paste0("requires G + W1/W2 + completeness (have W1/W2=", has_W1W2, ")")
  reasons[8] <- if (eligible[8]) "G + Gm + W1/W2 present, F_G ok, completeness ok" else
    paste0("requires G + Gm + W1/W2 + completeness (have Gm=", has_Gm, ", W1/W2=", has_W1W2, ")")

  data.frame(
    estimator = methods,
    eligible = eligible,
    a2_required = a2_required,
    reason = reasons,
    stringsAsFactors = FALSE
  )
}


#' Build human-readable diagnosis summary (internal)
#' @keywords internal
.build_diagnosis_summary <- function(data, inst, nc_val, nc_ind,
                                     nc_ind_gm, completeness,
                                     eligibility) {
  lines <- character(0)

  # Instrument strength
  if (data$has_instrument) {
    if (is.null(inst$g_threshold)) {
      lines <- c(lines, sprintf(" G (exposure instrument): partial F = %.1f %s",
                                inst$F_G_median,
                                ifelse(inst$weak_G, "(WEAK)", "(ok)")))
    } else {
      lines <- c(lines, sprintf(
        " G (exposure instrument): median F = %.1f, %.0f%% of %d transcripts >= %g %s (E threshold: %.0f%%)",
        inst$F_G_median, 100 * inst$F_G_p_pass, inst$F_G_n_total,
        inst$g_threshold$R,
        ifelse(inst$weak_G, "(BELOW E)", "(ok)"),
        100 * inst$g_threshold$E))
    }
  }
  if (data$has_mediator_instrument) {
    if (is.null(inst$gm_threshold)) {
      lines <- c(lines, sprintf(" Gm (mediator instrument): partial F = %.1f %s",
                                inst$F_Gm_median,
                                ifelse(inst$weak_Gm, "(WEAK)", "(ok)")))
    } else {
      lines <- c(lines, sprintf(
        " Gm (mediator instrument): median F = %.1f, %.0f%% of %d transcripts >= %g %s (E threshold: %.0f%%)",
        inst$F_Gm_median, 100 * inst$F_Gm_p_pass, inst$F_Gm_n_total,
        inst$gm_threshold$R,
        ifelse(inst$weak_Gm, "(BELOW E)", "(ok)"),
        100 * inst$gm_threshold$E))
    }
  }

  # NC validity
  if (!is.null(nc_val)) {
    n_drop <- sum(nc_val$significant, na.rm = TRUE)
    n_total <- nrow(nc_val)
    lines <- c(lines, sprintf(" NC validity (A1): %d/%d controls valid (%d flagged)",
                              n_total - n_drop, n_total, n_drop))
  }
  if (!is.null(nc_ind)) {
    n_drop <- sum(nc_ind$significant, na.rm = TRUE)
    n_total <- nrow(nc_ind)
    lines <- c(lines, sprintf(" NC independence (A2): %d/%d controls valid (%d flagged)",
                              n_total - n_drop, n_total, n_drop))
  }
  if (!is.null(nc_ind_gm)) {
    n_drop <- sum(nc_ind_gm$significant, na.rm = TRUE)
    n_total <- nrow(nc_ind_gm)
    lines <- c(lines, sprintf(" NC independence (A2'): %d/%d controls valid (%d flagged)",
                              n_total - n_drop, n_total, n_drop))
  }

  # Completeness
  if (!is.null(completeness)) {
    lines <- c(lines, sprintf(" Completeness: %d valid NCs vs k=%d -> %s",
                              completeness$n_valid_controls,
                              completeness$k,
                              completeness$completeness))
  }

  # Eligibility
  n_elig <- sum(eligibility$eligible)
  lines <- c(lines, "", sprintf(" Eligible estimators: %d/%d", n_elig, nrow(eligibility)))
  elig_methods <- eligibility$estimator[eligibility$eligible]
  if (length(elig_methods))
    lines <- c(lines, paste(" ", paste(elig_methods, collapse = ", ")))

  paste(lines, collapse = "\n")
}


#' Print method for iconic_diagnosis objects
#'
#' @param x An \code{iconic_diagnosis} object.
#' @param ... Unused.
#' @return Invisibly returns `x` (the `iconic_diagnosis` object); called for its side effect of printing a human-readable summary.
#' @export
print.iconic_diagnosis <- function(x, ...) {
  cat("<iconic_diagnosis>\n")
  cat("Diagnostic summary:\n")
  cat(x$summary, "\n")
  invisible(x)
}
