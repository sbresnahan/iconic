# ============================================================
# iconic_recommend: Model recommendation layer for the v0.6.0
# model selection workflow.
#
# Ranks eligible estimators by identification strength (tier) and
# optionally by robustness from a sensitivity analysis.  Returns a
# transparent recommendation with rationale.
#
# Tier system:
#   A (identified):     IV2SLS2, PGC2Gm  -- two instruments + path NCs
#   B (negative-control): PGC, PGC2, COCA -- NC-based confounding correction
#   C (instrument-based): IV2SLS          -- single instrument, no NC
#   D (bias reference):   UNADJ, DIRECT   -- naive, no causal identification
# ============================================================

# Tier assignments (internal)
.estimator_tiers <- c(
  UNADJ   = "D",
  DIRECT  = "D",
  COCA    = "B",
  IV2SLS  = "C",
  PGC     = "B",
  IV2SLS2 = "A",
  PGC2    = "B",
  PGC2Gm  = "A"
)

# Human-readable tier descriptions
# v0.9.2: Tier A label changed to "identified under stated assumptions" to
# acknowledge that no tier is guaranteed identified without its assumptions
# holding (JYH #414).
.tier_labels <- c(
  A = "identified under stated assumptions (two instruments + path NCs)",
  B = "negative-control based",
  C = "instrument-based (single IV)",
  D = "bias reference (naive)"
)

# What each estimator requires (for rationale)
# v0.9.2: COCA requirement updated to reflect A2 exemption (JYH #509).
.estimator_requirements <- c(
  UNADJ   = "no assumptions (naive OLS)",
  DIRECT  = "G + W as covariates (no causal identification)",
  COCA    = "valid NCs (A1), completeness (A2 not required)",
  IV2SLS  = "valid G (F>=10), exclusion restriction",
  PGC     = "valid G (F>=10), valid NCs (A2), completeness",
  IV2SLS2 = "valid G + Gm (F>=10), exclusion for both, valid NCs",
  PGC2    = "valid G (F>=10), path-specific NCs (W1/W2), completeness",
  PGC2Gm  = "valid G + Gm (F>=10), path-specific NCs (W1/W2), completeness"
)


#' Recommend the best causal estimator for the user's data
#'
#' Ranks all eligible estimators by identification strength (tier) and,
#' when a sensitivity analysis is supplied, by robustness to assumption
#' violations.  Returns the top-ranked estimator with a transparent
#' rationale.
#'
#' The tier system reflects the strength of identification:
#' \itemize{
#'   \item \strong{Tier A (identified)}: IV2SLS2, PGC2Gm — use two
#'     instruments and path-specific negative controls, providing point
#'     identification under mediator-outcome confounding.
#'   \item \strong{Tier B (negative-control)}: PGC, PGC2, COCA — use
#'     negative controls to proxy unmeasured confounding.
#'   \item \strong{Tier C (instrument-based)}: IV2SLS — single
#'     instrument without NC correction.
#'   \item \strong{Tier D (bias reference)}: UNADJ, DIRECT — naive
#'     adjustment with no causal identification.
#' }
#'
#' When \code{sensitivity} is supplied (from
#' \code{\link{iconic_sensitivity}()}), estimators within the same tier
#' are further ranked by robustness: the estimator whose estimates
#' degrade least across the assumption-violation surface ranks higher.
#'
#' @param data       An \code{iconic_data} object.
#' @param diagnosis  Optional \code{iconic_diagnosis} from
#'   \code{\link{iconic_diagnose}()}.  If \code{NULL}, auto-eligibility
#'   is computed from the data.
#' @param estimate   Optional estimate data frame from
#'   \code{\link{iconic_estimate}()}.  Used to report point estimates
#'   alongside the recommendation.
#' @param sensitivity Optional \code{iconic_sensitivity} from
#'   \code{\link{iconic_sensitivity}()}.  Used for robustness ranking
#'   within tiers.
#' @param criterion Character (v0.9.2, JYH #544): \code{"combined"} (default),
#'   \code{"minimax_bias"}, or \code{"ci_coverage"}. Controls how robustness
#'   is ranked within tiers. \code{"combined"} blends bias and CI coverage;
#'   \code{"minimax_bias"} ranks by worst-case bias across the violation grid;
#'   \code{"ci_coverage"} ranks by CI coverage of the true effect. When
#'   \code{sensitivity} is supplied, also populates \code{$per_scenario}
#'   (best estimator at the origin vs at violation cells).
#'
#' @return An \code{iconic_recommendation} S3 object: a named list with
#'   \code{$ranking} (data frame: estimator, tier, eligible, rank,
#'   rationale), \code{$recommended} (name of top estimator),
#'   \code{$recommended_tier}, \code{$per_scenario} (v0.9.2, when
#'   \code{sensitivity} is supplied), and \code{$summary}.
#' @export
#'
#' @examples
#' \dontrun{
#' data <- iconic_data(Z = rnorm(100), Y = matrix(rnorm(100*10), 10, 100),
#'                     G = rnorm(100), W = matrix(rnorm(100*10), 10, 100))
#' diag <- iconic_diagnose(data)
#' est  <- iconic_estimate(data, diagnosis = diag)
#' rec  <- iconic_recommend(data, diagnosis = diag, estimate = est)
#' print(rec)
#' }
iconic_recommend <- function(data, diagnosis = NULL, estimate = NULL,
                             sensitivity = NULL,
                             criterion = c("combined","minimax_bias","ci_coverage")) {
  criterion <- match.arg(criterion)
  if (!inherits(data, "iconic_data"))
    stop("data must be an iconic_data object from iconic_data().")

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
  tiers <- .estimator_tiers[all_methods]
  eligible <- elig$eligible[match(all_methods, elig$estimator)]

  ranking <- data.frame(
    estimator = all_methods,
    tier      = tiers,
    eligible  = eligible,
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

  # Add sensitivity-based robustness if available
  if (!is.null(sensitivity)) {
    robustness <- .extract_robustness(sensitivity, criterion = criterion)
    ranking$max_bias <- robustness$max_bias[match(ranking$estimator, robustness$method)]
    ranking$coverage_dist <- robustness$coverage_dist[match(ranking$estimator, robustness$method)]
    ranking$robustness_score <- robustness$score[match(ranking$estimator, robustness$method)]
  }

  # Rank: eligible first, then by tier (A > B > C > D), then by robustness
  # (if available).  robustness_score may be absent when no sensitivity
  # analysis is supplied.
  if ("robustness_score" %in% names(ranking)) {
    ranking$robustness_score[is.na(ranking$robustness_score)] <- -Inf
    ranking <- ranking[order(!ranking$eligible,
                             factor(ranking$tier, levels = c("A", "B", "C", "D")),
                             -ranking$robustness_score), ]
  } else {
    ranking <- ranking[order(!ranking$eligible,
                             factor(ranking$tier, levels = c("A", "B", "C", "D"))), ]
  }
  ranking$rank <- seq_len(nrow(ranking))

  # Build rationale for each
  ranking$rationale <- vapply(seq_len(nrow(ranking)), function(i) {
    m <- ranking$estimator[i]
    t <- ranking$tier[i]
    e <- ranking$eligible[i]
    req <- .estimator_requirements[m]
    if (!e) {
      paste0("ineligible — ", elig$reason[elig$estimator == m])
    } else {
      paste0("tier ", t, " (", .tier_labels[t], "); requires: ", req)
    }
  }, character(1))

  # Recommended estimator
  elig_ranking <- ranking[ranking$eligible, ]
  recommended <- if (nrow(elig_ranking) > 0) elig_ranking$estimator[1] else NA
  recommended_tier <- if (!is.na(recommended)) as.character(.estimator_tiers[recommended]) else NA

  # v0.9.2: per-scenario top estimator (JYH #544).
  per_scenario <- if (!is.null(sensitivity))
    .extract_per_scenario(sensitivity, criterion = criterion) else NULL

  # Summary
  summary_txt <- .build_recommendation_summary(ranking, recommended,
                                               recommended_tier, data)

  obj <- list(
    ranking          = ranking,
    recommended      = recommended,
    recommended_tier = recommended_tier,
    per_scenario     = per_scenario,
    criterion        = criterion,
    summary          = summary_txt
  )
  class(obj) <- c("iconic_recommendation", "list")
  obj
}


#' Extract robustness scores from sensitivity analysis (internal)
#'
#' Computes a robustness score for each estimator from the sensitivity
#' surface. v0.9.2 adds a `criterion` argument:
#'   - "minimax_bias" (legacy): lower maximum absolute bias across the
#'     violation grid means higher robustness.
#'   - "ci_coverage": distance to nominal 95% CI coverage (|coverage - 0.95|),
#'     lower is better. Requires NDE_coverage/NIE_coverage columns.
#'   - "combined" (default): normalized max|bias| + normalized |coverage-0.95|,
#'     lower is better. Implements the "second-best bias with nominal coverage"
#'     tradeoff (JYH #480).
#' @keywords internal
.extract_robustness <- function(sensitivity, criterion = c("combined","minimax_bias","ci_coverage")) {
  criterion <- match.arg(criterion)
  if (is.null(sensitivity) || is.null(sensitivity$surface))
    return(data.frame(method = character(0), max_bias = numeric(0),
                      coverage_dist = numeric(0), score = numeric(0)))

  surface <- sensitivity$surface

  # Determine which bias columns exist
  if ("NDE_bias" %in% names(surface)) {
    bias_cols <- c("NDE_bias", "NIE_bias")
  } else if ("bias" %in% names(surface)) {
    bias_cols <- "bias"
  } else {
    return(data.frame(method = character(0), max_bias = numeric(0),
                      coverage_dist = numeric(0), score = numeric(0)))
  }

  # Coverage columns (v0.9.2). May be absent if the surface was generated
  # by a pre-v0.9.2 sensitivity run.
  has_cov <- "NDE_coverage" %in% names(surface) && "NIE_coverage" %in% names(surface)
  cov_cols <- if (has_cov) c("NDE_coverage", "NIE_coverage") else NULL

  methods <- unique(surface$method)

  # Max absolute bias per method across the surface
  max_bias <- vapply(methods, function(m) {
    sub <- surface[surface$method == m, ]
    max(abs(unlist(sub[, bias_cols, drop = FALSE])), na.rm = TRUE)
  }, numeric(1))

  # Coverage distance: mean |coverage - 0.95| across cells (0 = nominal).
  coverage_dist <- if (has_cov) {
    vapply(methods, function(m) {
      sub <- surface[surface$method == m, ]
      mean(abs(unlist(sub[, cov_cols, drop = FALSE]) - 0.95), na.rm = TRUE)
    }, numeric(1))
  } else rep(NA_real_, length(methods))

  # Score by criterion.
  score <- switch(criterion,
    minimax_bias = 1 / (1 + max_bias),
    ci_coverage  = {
      if (!has_cov) return(data.frame(method = methods, max_bias = max_bias,
                                       coverage_dist = coverage_dist,
                                       score = rep(NA_real_, length(methods))))
      1 / (1 + coverage_dist)
    },
    combined = {
      # Normalize each component to [0,1] across estimators, then sum (lower better).
      # Convert to a "higher = better" score via 1/(1+combined_distance).
      n_bias <- (max_bias - min(max_bias, na.rm = TRUE)) /
                (max(max_bias, na.rm = TRUE) - min(max_bias, na.rm = TRUE) + 1e-8)
      if (has_cov) {
        n_cov <- (coverage_dist - min(coverage_dist, na.rm = TRUE)) /
                 (max(coverage_dist, na.rm = TRUE) - min(coverage_dist, na.rm = TRUE) + 1e-8)
        combined_dist <- n_bias + n_cov
      } else {
        combined_dist <- n_bias
      }
      1 / (1 + combined_dist)
    }
  )

  data.frame(method = methods, max_bias = max_bias,
             coverage_dist = coverage_dist, score = score,
             stringsAsFactors = FALSE)
}


#' Extract per-scenario top estimator from sensitivity surface (internal)
#'
#' v0.9.2: For each cell of the sensitivity surface, identifies the
#' top-ranked estimator by the chosen criterion. Returns a data frame
#' with one row per cell (JYH #544: a single global recommendation is
#' insufficient; some scenarios favor IV, others PGC).
#' @keywords internal
.extract_per_scenario <- function(sensitivity, criterion = c("combined","minimax_bias","ci_coverage")) {
  criterion <- match.arg(criterion)
  if (is.null(sensitivity) || is.null(sensitivity$surface)) return(NULL)
  surface <- sensitivity$surface

  # Identify the cell-defining columns (rho_G1/rho_G2 for the degradation
  # surface; conf_strength/coverage/k for the total-effect surface).
  cell_cols <- intersect(c("rho_G1","rho_G2","conf_strength","coverage","k"),
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
    # Per-method metric for this cell.
    mb <- vapply(methods, function(m) {
      mm <- sub[sub$method == m, ]
      max(abs(unlist(mm[, bias_cols, drop = FALSE])), na.rm = TRUE)
    }, numeric(1))
    if (has_cov) {
      cd <- vapply(methods, function(m) {
        mm <- sub[sub$method == m, ]
        mean(abs(unlist(mm[, cov_cols, drop = FALSE]) - 0.95), na.rm = TRUE)
      }, numeric(1))
    } else cd <- rep(NA_real_, length(methods))
    # Pick top by criterion (lower distance = better).
    dist <- switch(criterion,
      minimax_bias = mb,
      ci_coverage  = if (has_cov) cd else mb,
      combined     = {
        nmb <- (mb - min(mb)) / (max(mb) - min(mb) + 1e-8)
        if (has_cov) {
          ncd <- (cd - min(cd)) / (max(cd) - min(cd) + 1e-8)
          nmb + ncd
        } else nmb
      })
    top <- methods[which.min(dist)]
    rows[[i]] <- cbind(cell, data.frame(top_estimator = top,
                                        top_bias = mb[which.min(dist)],
                                        top_coverage_dist = cd[which.min(dist)],
                                        stringsAsFactors = FALSE))
  }
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}


#' Build recommendation summary (internal)
#' @keywords internal
.build_recommendation_summary <- function(ranking, recommended,
                                          recommended_tier, data) {
  lines <- character(0)

  if (is.na(recommended)) {
    lines <- c(lines, "  No eligible estimators beyond UNADJ.",
               "  Consider supplying instruments (G) or negative controls (W).")
  } else {
    lines <- c(lines,
      sprintf("  Recommended: %s (tier %s — %s)",
              recommended, recommended_tier, .tier_labels[recommended_tier]))

    # Report point estimate if available
    if ("mean_beta" %in% names(ranking)) {
      val <- ranking$mean_beta[ranking$estimator == recommended]
      if (!is.na(val))
        lines <- c(lines, sprintf("  Estimated total effect (mean): %.4f", val))
    }
    if ("mean_NDE" %in% names(ranking)) {
      nde <- ranking$mean_NDE[ranking$estimator == recommended]
      nie <- ranking$mean_NIE[ranking$estimator == recommended]
      if (!is.na(nde))
        lines <- c(lines, sprintf("  Estimated NDE (mean): %.4f", nde))
      if (!is.na(nie))
        lines <- c(lines, sprintf("  Estimated NIE (mean): %.4f", nie))
    }

    # Robustness info
    if ("max_bias" %in% names(ranking)) {
      mb <- ranking$max_bias[ranking$estimator == recommended]
      if (!is.na(mb))
        lines <- c(lines, sprintf("  Max |bias| across sensitivity surface: %.4f", mb))
    }

    # List eligible alternatives
    elig <- ranking[ranking$eligible, ]
    if (nrow(elig) > 1) {
      alts <- elig$estimator[-1]
      lines <- c(lines, "",
                 sprintf("  Eligible alternatives: %s",
                         paste(alts, collapse = ", ")))
    }
  }

  paste(lines, collapse = "\n")
}


#' Print method for iconic_recommendation objects
#'
#' @param x An \code{iconic_recommendation} object.
#' @param ... Unused.
#' @export
print.iconic_recommendation <- function(x, ...) {
  cat("<iconic_recommendation>\n")
  cat(x$summary, "\n")
  cat("\n  Full ranking:\n")
  e <- x$ranking[x$ranking$eligible, ]
  if (nrow(e) > 0) {
    for (i in seq_len(nrow(e))) {
      cat(sprintf("    %d. %s [tier %s] — %s\n",
                  e$rank[i], e$estimator[i], e$tier[i], e$rationale[i]))
    }
  }
  inelig <- x$ranking[!x$ranking$eligible, ]
  if (nrow(inelig) > 0) {
    cat("\n  Ineligible:\n")
    for (i in seq_len(nrow(inelig))) {
      cat(sprintf("    %s — %s\n", inelig$estimator[i], inelig$rationale[i]))
    }
  }
  invisible(x)
}
