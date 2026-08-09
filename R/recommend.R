# ============================================================
# iconic_recommend: Model recommendation layer for the
# model selection workflow.
#
# Ranks eligible estimators directly on per-estimand robustness
# (NDE and NIE separately) from a sensitivity analysis. Returns a
# transparent recommendation with rationale. No tier system: every
# estimator is ranked on its measured robustness, not on an a priori
# identification-strength class.
# ============================================================

# What each estimator requires (for rationale)
# COCA requirement updated to reflect A2 exemption.
.estimator_requirements <- c(
  UNADJ = "no assumptions (naive OLS)",
  DIRECT = "G + W as covariates (no causal identification)",
  COCA = "valid NCs (A1), completeness (A2 not required)",
  IV2SLS = "valid G (F>=10), exclusion restriction",
  PGC = "valid G (F>=10), valid NCs (A2), completeness",
  IV2SLS2 = "valid G + Gm (F>=10), exclusion for both; optional path-specific NCs (W1/W2)",
  PGC2 = "valid G (F>=10), path-specific NCs (W1/W2), completeness",
  PGC2Gm = "valid G + Gm (F>=10), path-specific NCs (W1/W2), completeness"
)


#' Recommend the best causal estimator for the user's data
#'
#' Ranks all eligible estimators directly on per-estimand robustness
#' (NDE and NIE separately) from a sensitivity analysis. Returns the
#' top-ranked estimator with a transparent rationale. There is no tier
#' system: every estimator is ranked on its measured robustness to
#' assumption violations, not on an a priori identification-strength
#' class.
#'
#' When \code{sensitivity} is supplied (from
#' \code{\link{iconic_sensitivity}()}), eligible estimators are ranked by
#' robustness: the estimator whose estimates degrade least across the
#' assumption-violation surface ranks higher. NDE and NIE robustness are
#' computed separately, so a mediation estimator is ranked on the
#' estimand of interest rather than a pooled maximum.
#'
#' @param data An \code{iconic_data} object.
#' @param diagnosis Optional \code{iconic_diagnosis} from
#' \code{\link{iconic_diagnose}()}. If \code{NULL}, auto-eligibility
#' is computed from the data.
#' @param estimate Optional estimate data frame from
#' \code{\link{iconic_estimate}()}. Used to report point estimates
#' alongside the recommendation.
#' @param sensitivity Optional \code{iconic_sensitivity} from
#' \code{\link{iconic_sensitivity}()}. Used for robustness ranking. When
#' \code{NULL} and \code{auto_sensitivity = TRUE} (the default), the
#' sensitivity suite is run automatically so the recommendation is
#' robustness-based out of the box.
#' @param criterion Character: \code{"combined"} (default),
#' \code{"minimax_bias"}, or \code{"ci_coverage"}. Controls how robustness
#' is ranked. \code{"combined"} blends bias and CI coverage;
#' \code{"minimax_bias"} ranks by worst-case bias across the violation grid;
#' \code{"ci_coverage"} ranks by CI coverage of the true effect. When
#' \code{sensitivity} is supplied, also populates \code{$per_scenario}
#' (best estimator at the origin vs at violation cells).
#' @param auto_sensitivity Logical: when \code{TRUE} (default) and
#' \code{sensitivity = NULL}, run \code{\link{iconic_sensitivity}()}
#' internally to obtain the robustness surface. Requires the torch backend;
#' when torch is unavailable the function falls back to eligibility-only
#' ranking with a message. Set \code{FALSE} to skip the auto-run.
#' @param rho_G1_grid,rho_G2_grid Instrument-exogeneity violation grid used
#' for the auto-run sensitivity suite. Default \code{c(0, 0.1, 0.2, 0.3, 0.5)}.
#' @param omega_1,omega_2 Negative-control coverage grid for the auto-run.
#' Default \code{c(0.3, 0.7, 1.0)}, swept on the diagonal
#' (\code{omega_1 == omega_2}).
#' @param n_iter_sens Replicates per grid cell for the auto-run. Default 30.
#' @param gan_epochs GAN training epochs for the auto-run. Default 100.
#' @param n_cores Cores for the auto-run sensitivity sweep. Default 1.
#' @param verbose Logical: print progress messages. Default \code{FALSE}
#' (quiet). Also silences the auto-run sensitivity suite.
#'
#' @return An \code{iconic_recommendation} S3 object: a named list with
#' \code{$ranking} (data frame: estimator, eligible, rank, per-estimand
#' robustness scores, rationale), \code{$recommended} (top estimator for
#' NDE), \code{$recommended_NIE} (top estimator for NIE, when it
#' differs), \code{$per_scenario} (when \code{sensitivity} is supplied),
#' and \code{$summary}.
#' @export
#'
#' @examples
#' \donttest{
#' data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100*10), 10, 100),
#' G = rnorm(100), W = matrix(rnorm(100*10), 10, 100))
#' diag <- iconic_diagnose(data)
#' est <- iconic_estimate(data, diagnosis = diag)
#' rec <- iconic_recommend(data, diagnosis = diag, estimate = est)
#' print(rec)
#' }
iconic_recommend <- function(data, diagnosis = NULL, estimate = NULL,
                             sensitivity = NULL,
                             criterion = c("combined","minimax_bias","ci_coverage"),
                             auto_sensitivity = TRUE,
                             rho_G1_grid = c(0, 0.1, 0.2, 0.3, 0.5),
                             rho_G2_grid = c(0, 0.1, 0.2, 0.3, 0.5),
                             omega_1 = c(0.3, 0.7, 1.0),
                             omega_2 = c(0.3, 0.7, 1.0),
                             n_iter_sens = 30,
                             gan_epochs = 100,
                             n_cores = 1,
                             verbose = FALSE) {
  criterion <- match.arg(criterion)
  if (!inherits(data, "iconic_data"))
    stop("data must be an iconic_data object from iconic_data().")

  # Auto-run the sensitivity suite when no surface is supplied, so the
  # recommendation is robustness-based by default (max |bias| per estimand
  # over the joint rho x omega grid). Requires the torch backend to train the
  # texture model; without it, fall back to eligibility-only ranking.
  if (is.null(sensitivity) && isTRUE(auto_sensitivity)) {
    if (isTRUE(check_torch_setup())) {
      if (isTRUE(verbose))
        message("iconic_recommend: no sensitivity surface supplied; running ",
                "iconic_sensitivity() (rho x omega grid) for robustness ranking.")
      sensitivity <- tryCatch(
        iconic_sensitivity(data, diagnosis = diagnosis,
                           rho_G1_grid = rho_G1_grid, rho_G2_grid = rho_G2_grid,
                           omega_1 = omega_1, omega_2 = omega_2,
                           n_iter = n_iter_sens, gan_epochs = gan_epochs,
                           n_cores = n_cores, verbose = verbose),
        error = function(e) {
          warning("auto-sensitivity failed (", e$message,
                  "); falling back to eligibility-only ranking.")
          NULL
        })
    } else {
      if (isTRUE(verbose))
        message("iconic_recommend: torch backend unavailable; skipping ",
                "auto-sensitivity and using eligibility-only ranking.")
    }
  }

  # Get eligibility
  if (!is.null(diagnosis)) {
    elig <- diagnosis$eligibility
  } else {
    diag_auto <- iconic_diagnose(data)
    elig <- diag_auto$eligibility
  }

  # Build ranking table
  all_methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
                   "IV2SLS2", "PGC2", "PGC2Gm")
  eligible <- elig$eligible[match(all_methods, elig$estimator)]

  ranking <- data.frame(
    estimator = all_methods,
    eligible = eligible,
    stringsAsFactors = FALSE
  )

  # Add point estimates if available
  if (!is.null(estimate)) {
    if (data$is_mediation && "NDE" %in% names(estimate)) {
      # Mediation: use mean NDE across features/mediators per method
      est_summary <- aggregate(cbind(NDE, NIE) ~ method,
                               data = estimate, FUN = mean, na.rm = TRUE)
      ranking$mean_NDE <- est_summary$NDE[match(ranking$estimator, est_summary$method)]
      ranking$mean_NIE <- est_summary$NIE[match(ranking$estimator, est_summary$method)]
    } else if ("beta" %in% names(estimate)) {
      # Total effect: mean beta per method
      est_summary <- aggregate(beta ~ method,
                               data = estimate, FUN = mean, na.rm = TRUE)
      ranking$mean_beta <- est_summary$beta[match(ranking$estimator, est_summary$method)]
    }
  }

  # Add sensitivity-based robustness if available (per estimand)
  if (!is.null(sensitivity)) {
    robustness <- .extract_robustness(sensitivity, criterion = criterion)
    ranking$max_bias_NDE <- robustness$max_bias_NDE[match(ranking$estimator, robustness$method)]
    ranking$max_bias_NIE <- robustness$max_bias_NIE[match(ranking$estimator, robustness$method)]
    ranking$coverage_dist_NDE <- robustness$coverage_dist_NDE[match(ranking$estimator, robustness$method)]
    ranking$coverage_dist_NIE <- robustness$coverage_dist_NIE[match(ranking$estimator, robustness$method)]
    ranking$robustness_NDE <- robustness$score_NDE[match(ranking$estimator, robustness$method)]
    ranking$robustness_NIE <- robustness$score_NIE[match(ranking$estimator, robustness$method)]
  }

  # Rank: eligible first, then by per-estimand robustness (NDE primary,
  # NIE secondary). robustness scores may be absent when no sensitivity
  # analysis is supplied; then rank on eligibility alone (stable order).
  if ("robustness_NDE" %in% names(ranking)) {
    ranking$robustness_NDE[is.na(ranking$robustness_NDE)] <- -Inf
    ranking$robustness_NIE[is.na(ranking$robustness_NIE)] <- -Inf
    ranking <- ranking[order(!ranking$eligible,
                             -ranking$robustness_NDE,
                             -ranking$robustness_NIE), ]
  } else {
    # No sensitivity surface: there is no valid data-driven robustness proxy
    # (the truth is unknown, and deviation from the confounded UNADJ estimate
    # is not evidence of bias). The only defensible tie-break is structural:
    # UNADJ and DIRECT assume no unmeasured confounding -- the very problem
    # this package exists to solve -- so they must never be surface as the
    # top recommendation merely because they appear first in the method list.
    # Demote them below the instrument/NC-based estimators; among the latter,
    # keep the stable (eligibility) order. This is a last-resort ordering used
    # only when the caller supplied no sensitivity analysis.
    is_naive <- ranking$estimator %in% c("UNADJ", "DIRECT")
    ranking <- ranking[order(!ranking$eligible, is_naive), ]
  }
  ranking$rank <- seq_len(nrow(ranking))

  # Build rationale for each
  ranking$rationale <- vapply(seq_len(nrow(ranking)), function(i) {
    m <- ranking$estimator[i]
    e <- ranking$eligible[i]
    req <- .estimator_requirements[m]
    if (!e) {
      paste0("ineligible -- ", elig$reason[elig$estimator == m])
    } else {
      paste0("requires: ", req)
    }
  }, character(1))

  # Recommended estimator (top eligible by NDE robustness)
  elig_ranking <- ranking[ranking$eligible, ]
  recommended <- if (nrow(elig_ranking) > 0) elig_ranking$estimator[1] else NA
  # Per-estimand recommendation: the top eligible estimator for NIE may
  # differ from NDE when their robustness profiles differ.
  recommended_NIE <- if (nrow(elig_ranking) > 0 && "robustness_NIE" %in% names(ranking)) {
    elig_nie <- elig_ranking[order(-elig_ranking$robustness_NIE), ]
    elig_nie$estimator[1]
  } else recommended

  # per-scenario top estimator.
  per_scenario <- if (!is.null(sensitivity))
    .extract_per_scenario(sensitivity, criterion = criterion) else NULL

  # Summary
  summary_txt <- .build_recommendation_summary(ranking, recommended,
                                               recommended_NIE, data)

  obj <- list(
    ranking = ranking,
    recommended = recommended,
    recommended_NIE = recommended_NIE,
    per_scenario = per_scenario,
    criterion = criterion,
    summary = summary_txt
  )
  class(obj) <- c("iconic_recommendation", "list")
  if (isTRUE(verbose))
    message("iconic_recommend complete. Call summary() or print() on the result for the full recommendation.")
  obj
}


#' Extract robustness scores from sensitivity analysis (internal)
#'
#' Computes a robustness score for each estimator from the sensitivity
#' surface. adds a `criterion` argument:
#' - "minimax_bias" (legacy): lower maximum absolute bias across the
#' violation grid means higher robustness.
#' - "ci_coverage": distance to nominal 95% CI coverage (|coverage - 0.95|),
#' lower is better. Requires NDE_coverage/NIE_coverage columns.
#' - "combined" (default): normalized max|bias| + normalized |coverage-0.95|,
#' lower is better. Implements the "second-best bias with nominal coverage"
#' tradeoff.
#' @keywords internal
.extract_robustness <- function(sensitivity, criterion = c("combined","minimax_bias","ci_coverage")) {
  criterion <- match.arg(criterion)
  empty <- data.frame(method = character(0), max_bias_NDE = numeric(0),
                      max_bias_NIE = numeric(0), coverage_dist_NDE = numeric(0),
                      coverage_dist_NIE = numeric(0), score_NDE = numeric(0),
                      score_NIE = numeric(0), stringsAsFactors = FALSE)
  if (is.null(sensitivity) || is.null(sensitivity$surface))
    return(empty)

  surface <- sensitivity$surface

  # Determine which bias columns exist
  is_mediation <- "NDE_bias" %in% names(surface)
  if (!is_mediation && !"bias" %in% names(surface)) return(empty)

  has_cov <- "NDE_coverage" %in% names(surface) && "NIE_coverage" %in% names(surface)
  methods <- unique(surface$method)

  # Per-estimand max absolute bias and coverage distance.
  # For total-effect surfaces (single "bias" column), NDE = the total
  # effect and NIE is NA.
  per_estimand <- function(bias_col, cov_col) {
    mb <- vapply(methods, function(m) {
      sub <- surface[surface$method == m, ]
      max(abs(sub[[bias_col]]), na.rm = TRUE)
    }, numeric(1))
    cd <- if (has_cov && !is.null(cov_col)) {
      vapply(methods, function(m) {
        sub <- surface[surface$method == m, ]
        mean(abs(sub[[cov_col]] - 0.95), na.rm = TRUE)
      }, numeric(1))
    } else rep(NA_real_, length(methods))
    list(max_bias = mb, coverage_dist = cd)
  }

  if (is_mediation) {
    nde <- per_estimand("NDE_bias", if (has_cov) "NDE_coverage" else NULL)
    nie <- per_estimand("NIE_bias", if (has_cov) "NIE_coverage" else NULL)
  } else {
    nde <- per_estimand("bias", NULL)
    nie <- list(max_bias = rep(NA_real_, length(methods)),
                coverage_dist = rep(NA_real_, length(methods)))
  }

  # Score by criterion, computed per estimand (higher = better).
  # Normalize over finite values only: an estimator with an all-NA bias
  # column (e.g. COCA on survival outcomes, where it is structurally
  # incompatible) yields max_bias = -Inf, which would otherwise poison the
  # min-max scaling for every estimator (min = -Inf -> NaN throughout).
  # Non-finite scores are set to -Inf (worst) so such estimators rank last.
  score_one <- function(mb, cd) {
    switch(criterion,
      minimax_bias = {
        s <- 1 / (1 + mb)
        s[!is.finite(s)] <- -Inf
        s
      },
      ci_coverage = {
        if (!has_cov || all(is.na(cd))) rep(NA_real_, length(methods))
        else {
          s <- 1 / (1 + cd)
          s[!is.finite(s)] <- -Inf
          s
        }
      },
      combined = {
        norm <- function(v) {
          vf <- v[is.finite(v)]
          if (length(vf) < 2 || diff(range(vf)) == 0)
            return(ifelse(is.finite(v), 0, NA_real_))
          out <- (v - min(vf)) / (diff(range(vf)) + 1e-8)
          out[!is.finite(out)] <- NA_real_
          out
        }
        n_bias <- norm(mb)
        s <- if (has_cov && !all(is.na(cd))) {
          n_cov <- norm(cd)
          # if exactly one of n_bias/n_cov is NA, treat the NA as worst (1)
          nb <- ifelse(is.na(n_bias), 1, n_bias)
          nc <- ifelse(is.na(n_cov), 1, n_cov)
          1 / (1 + nb + nc)
        } else {
          1 / (1 + ifelse(is.na(n_bias), 1, n_bias))
        }
        # estimators with non-finite bias get the worst score
        s[!is.finite(mb)] <- -Inf
        s
      })
  }

  data.frame(method = methods,
             max_bias_NDE = nde$max_bias, max_bias_NIE = nie$max_bias,
             coverage_dist_NDE = nde$coverage_dist, coverage_dist_NIE = nie$coverage_dist,
             score_NDE = score_one(nde$max_bias, nde$coverage_dist),
             score_NIE = score_one(nie$max_bias, nie$coverage_dist),
             stringsAsFactors = FALSE)
}


#' Extract per-scenario top estimator from sensitivity surface (internal)
#'
#' For each cell of the sensitivity surface, identifies the
#' top-ranked estimator by the chosen criterion. Returns a data frame
#' with one row per cell (
#' insufficient; some scenarios favor IV, others PGC).
#' @keywords internal
.extract_per_scenario <- function(sensitivity, criterion = c("combined","minimax_bias","ci_coverage")) {
  criterion <- match.arg(criterion)
  if (is.null(sensitivity) || is.null(sensitivity$surface)) return(NULL)
  surface <- sensitivity$surface

  # Identify the cell-defining columns (rho_G1/rho_G2/omega for the
  # degradation surface; conf_strength/coverage/k for the total-effect
  # surface).
  cell_cols <- intersect(c("rho_G1","rho_G2","omega_1","omega_2",
                           "conf_strength","coverage","k"),
                         names(surface))
  if (length(cell_cols) == 0) return(NULL)

  bias_cols <- if ("NDE_bias" %in% names(surface)) c("NDE_bias","NIE_bias")
               else if ("bias" %in% names(surface)) "bias" else return(NULL)
  has_cov <- "NDE_coverage" %in% names(surface) && "NIE_coverage" %in% names(surface)
  cov_cols <- if (has_cov) c("NDE_coverage","NIE_coverage") else NULL

  cells <- unique(surface[, cell_cols, drop = FALSE])
  rows <- vector("list", nrow(cells))
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, , drop = FALSE]
    sub <- surface[do.call(paste, surface[, cell_cols, drop=FALSE]) ==
                   do.call(paste, cell), , drop = FALSE]
    if (nrow(sub) == 0) next
    methods <- unique(sub$method)
    methods <- methods[!is.na(methods)]
    if (length(methods) == 0) next
    # Per-method metric for this cell.
    mb <- vapply(methods, function(m) {
      mm <- sub[sub$method == m, ]
      v <- abs(unlist(mm[, bias_cols, drop = FALSE]))
      v <- v[is.finite(v)]
      if (length(v) == 0) NA_real_ else max(v)
    }, numeric(1))
    if (has_cov) {
      cd <- vapply(methods, function(m) {
        mm <- sub[sub$method == m, ]
        v <- abs(unlist(mm[, cov_cols, drop = FALSE]) - 0.95)
        v <- v[is.finite(v)]
        if (length(v) == 0) NA_real_ else mean(v)
      }, numeric(1))
    } else cd <- rep(NA_real_, length(methods))
    # Skip degenerate cells where every method returned a non-finite metric
    # (e.g. all estimators failed in this scenario); which.min() on an all-NA
    # vector returns integer(0) and breaks the cbind below.
    if (all(!is.finite(mb))) next
    # Pick top by criterion (lower distance = better). Normalise within the
    # cell using only finite values so a subset of failed estimators (NA)
    # does not poison the scaling.
    mb_f <- mb[is.finite(mb)]; cd_f <- cd[is.finite(cd)]
    norm <- function(v, vf) if (length(vf) > 1 && diff(range(vf)) > 0)
      (v - min(vf)) / (diff(range(vf)) + 1e-8) else rep(0, length(v))
    dist <- switch(criterion,
      minimax_bias = mb,
      ci_coverage = if (has_cov) cd else mb,
      combined = {
        nmb <- norm(mb, mb_f)
        if (has_cov) nmb + norm(cd, cd_f) else nmb
      })
    # Treat non-finite distances as worst so they are never picked as top.
    dist[!is.finite(dist)] <- Inf
    wi <- which.min(dist)
    if (length(wi) == 0 || !is.finite(dist[wi])) next
    top <- methods[wi]
    rows[[i]] <- cbind(cell, data.frame(top_estimator = top,
                                        top_bias = mb[wi],
                                        top_coverage_dist = cd[wi],
                                        stringsAsFactors = FALSE))
  }
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}


#' Build recommendation summary (internal)
#' @keywords internal
.build_recommendation_summary <- function(ranking, recommended,
                                          recommended_NIE, data) {
  lines <- character(0)

  if (is.na(recommended)) {
    lines <- c(lines, " No eligible estimators beyond UNADJ.",
               " Consider supplying instruments (G) or negative controls (W).")
  } else {
    lines <- c(lines,
      sprintf(" Recommended (NDE): %s", recommended))
    if (!is.null(recommended_NIE) && !is.na(recommended_NIE) &&
        recommended_NIE != recommended)
      lines <- c(lines,
        sprintf(" Recommended (NIE): %s", recommended_NIE))

    # Report point estimate if available
    if ("mean_beta" %in% names(ranking)) {
      val <- ranking$mean_beta[ranking$estimator == recommended]
      if (!is.na(val))
        lines <- c(lines, sprintf(" Estimated total effect (mean): %.4f", val))
    }
    if ("mean_NDE" %in% names(ranking)) {
      nde <- ranking$mean_NDE[ranking$estimator == recommended]
      nie <- ranking$mean_NIE[ranking$estimator == recommended]
      if (!is.na(nde))
        lines <- c(lines, sprintf(" Estimated NDE (mean): %.4f", nde))
      if (!is.na(nie))
        lines <- c(lines, sprintf(" Estimated NIE (mean): %.4f", nie))
    }

    # Robustness info (per estimand)
    if ("max_bias_NDE" %in% names(ranking)) {
      mb_nde <- ranking$max_bias_NDE[ranking$estimator == recommended]
      mb_nie <- ranking$max_bias_NIE[ranking$estimator == recommended]
      if (!is.na(mb_nde))
        lines <- c(lines, sprintf(" Max |NDE bias| across sensitivity surface: %.4f", mb_nde))
      if (!is.na(mb_nie))
        lines <- c(lines, sprintf(" Max |NIE bias| across sensitivity surface: %.4f", mb_nie))
    }

    # List eligible alternatives
    elig <- ranking[ranking$eligible, ]
    if (nrow(elig) > 1) {
      alts <- elig$estimator[-1]
      lines <- c(lines, "",
                 sprintf(" Eligible alternatives: %s",
                         paste(alts, collapse = ", ")))
    }
  }

  paste(lines, collapse = "\n")
}


#' Print method for iconic_recommendation objects
#'
#' @param x An \code{iconic_recommendation} object.
#' @param ... Unused.
#' @return Invisibly returns `x` (the `iconic_recommendation` object); called for its side effect of printing a human-readable summary.
#' @export
print.iconic_recommendation <- function(x, ...) {
  cat("<iconic_recommendation>\n")
  cat(x$summary, "\n")
  cat("\n Full ranking:\n")
  e <- x$ranking[x$ranking$eligible, ]
  if (nrow(e) > 0) {
    has_rob <- "robustness_NDE" %in% names(e)
    for (i in seq_len(nrow(e))) {
      rob <- if (has_rob && is.finite(e$robustness_NDE[i]))
        sprintf(" NDE score %.3f, NIE score %.3f;",
                e$robustness_NDE[i], e$robustness_NIE[i]) else ""
      cat(sprintf(" %d. %s --%s %s\n",
                  e$rank[i], e$estimator[i], rob, e$rationale[i]))
    }
  }
  inelig <- x$ranking[!x$ranking$eligible, ]
  if (nrow(inelig) > 0) {
    cat("\n Ineligible:\n")
    for (i in seq_len(nrow(inelig))) {
      cat(sprintf(" %s -- %s\n", inelig$estimator[i], inelig$rationale[i]))
    }
  }
  invisible(x)
}

#' Summary method for iconic_recommendation objects
#'
#' Prints the full recommendation summary (same as \code{print()}).
#' @param object An \code{iconic_recommendation} object.
#' @param ... Unused.
#' @return Invisibly returns \code{object}.
#' @export
summary.iconic_recommendation <- function(object, ...) {
  print.iconic_recommendation(object, ...)
  invisible(object)
}
