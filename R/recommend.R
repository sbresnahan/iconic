# ============================================================
# iconic_recommend: Model recommendation layer for the
# model selection workflow.
#
# Ranks eligible estimators on a data-driven composite robustness score
# from a sensitivity analysis: the worst-estimand (min of NDE and NIE)
# robustness, discounted by a confidence multiplier derived from the graded
# diagnostic verdict for the assumptions each estimator depends on. Returns
# a transparent recommendation with rationale. Every estimator is ranked on
# its measured robustness, not on an a priori identification-strength class.
# ============================================================

# What each estimator requires (for rationale). Threshold-aware: the
# instrument-strength label interpolates the actual min_f used in the
# diagnosis rather than a hardcoded value, so the rationale matches the
# gate that was actually applied.
# COCA requires valid negative controls (A1) but not completeness (A2).
.estimator_requirements <- function(min_f = 10) {
  f <- sprintf("F>=%s", min_f)
  c(
    UNADJ = "no assumptions (naive OLS)",
    DIRECT = "G + W as covariates (no causal identification)",
    COCA = "valid NCs (A1), completeness (A2 not required)",
    IV2SLS = sprintf("valid G (%s), exclusion restriction", f),
    PGC = sprintf("valid G (%s), valid NCs (A2), completeness", f),
    IV2SLS2 = sprintf("valid G + Gm (%s), exclusion for both; optional path-specific NCs (W1/W2)", f),
    PGC2 = sprintf("valid G (%s), path-specific NCs (W1/W2), completeness", f),
    PGC2Gm = sprintf("valid G + Gm (%s), path-specific NCs (W1/W2), completeness", f)
  )
}


#' Recommend the best causal estimator for the user's data
#'
#' Ranks all eligible estimators directly on per-estimand robustness
#' (NDE and NIE separately) from a sensitivity analysis. Returns the
#' top-ranked estimator with a transparent rationale. Every estimator is
#' ranked on its measured robustness to assumption violations, not on an
#' a priori identification-strength class.
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
#' @param completeness_penalty Named numeric vector with values in
#' \eqn{[0, 1]} mapping the
#' graded completeness verdict to a confidence multiplier applied to
#' bridge-dependent estimators (DIRECT, COCA, PGC, PGC2, PGC2Gm). Default
#' \code{c(satisfied = 1.0, borderline = 0.7, "weak-capture" = 0.5,
#' "under-identified" = 0)}. Instrument-only estimators (IV2SLS, IV2SLS2)
#' are not discounted. Override to sensitivity-check the discount.
#' @param min_f Instrument-strength gate (first-stage F threshold) used when
#' \code{diagnosis} is \code{NULL} and \code{iconic_diagnose()} is run
#' internally. Default 10. Also used for the requirement labels when the
#' supplied \code{diagnosis} predates the stored \code{min_f} field.
#' @param g_threshold,gm_threshold Optional instrument-selection thresholds
#' passed to \code{\link{iconic_diagnose}()} when \code{diagnosis} is
#' \code{NULL}. Default \code{NULL} (use the \code{iconic_diagnose}
#' defaults).
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
#' robustness scores, \code{composite}, \code{confidence_mult},
#' \code{final_score}, rationale), \code{$recommended} (top eligible,
#' non-naive estimator by the data-driven composite score),
#' \code{$recommended_NDE} and \code{$recommended_NIE} (top estimator per
#' estimand), \code{$per_scenario} (when \code{sensitivity} is supplied),
#' \code{$completeness_penalty}, and \code{$summary}.
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
                             completeness_penalty = c("satisfied" = 1.0,
                                                      "borderline" = 0.7,
                                                      "weak-capture" = 0.5,
                                                      "under-identified" = 0),
                             min_f = 10,
                             g_threshold = NULL,
                             gm_threshold = NULL,
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

  # Get eligibility (and keep the diagnosis object so the composite
  # confidence multiplier can read the graded completeness verdict).
  # When the caller does not supply a diagnosis, auto-diagnose using the
  # SAME instrument-strength thresholds the caller specified -- otherwise a
  # non-default min_f would be silently discarded here.
  if (!is.null(diagnosis)) {
    diag_obj <- diagnosis
  } else {
    diag_obj <- iconic_diagnose(data, min_f = min_f,
                                g_threshold = g_threshold,
                                gm_threshold = gm_threshold)
  }
  elig <- diag_obj$eligibility
  # Use the threshold actually applied in the diagnosis for the rationale
  # labels (falls back to the min_f argument if the diagnosis predates the
  # stored field).
  min_f_used <- if (!is.null(diag_obj$min_f)) diag_obj$min_f else min_f

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

  # Rank: eligible first, then by a data-driven composite robustness score.
  # The composite is the WORST-estimand robustness (min of NDE and NIE), so an
  # estimator is only as trustworthy as its weakest estimand, multiplied by a
  # confidence factor derived from the graded diagnostic verdict for the
  # assumptions that estimator depends on (e.g. path completeness for the
  # bridge estimators). robustness scores may be absent when no sensitivity
  # analysis is supplied; then rank on eligibility alone (stable order).
  if ("robustness_NDE" %in% names(ranking)) {
    ranking$robustness_NDE[is.na(ranking$robustness_NDE)] <- -Inf
    ranking$robustness_NIE[is.na(ranking$robustness_NIE)] <- -Inf
    comp <- .composite_robustness(ranking, diag_obj, completeness_penalty)
    ranking$composite <- comp$composite
    ranking$confidence_mult <- comp$confidence_mult
    ranking$final_score <- comp$final_score
    # Structural-naive estimators (UNADJ, DIRECT) assume no unmeasured
    # confounding -- the very problem this package exists to solve -- so they
    # must never out-rank an eligible instrument/NC-based estimator merely
    # because they dodge the completeness discount. Demote them below the
    # IV/NC-based estimators; among the latter, rank by final_score.
    is_naive <- ranking$estimator %in% c("UNADJ", "DIRECT")
    ranking <- ranking[order(!ranking$eligible,
                             is_naive,
                             -ranking$final_score), ]
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

  # Build rationale for each (threshold-aware requirement labels)
  reqs <- .estimator_requirements(min_f_used)
  ranking$rationale <- vapply(seq_len(nrow(ranking)), function(i) {
    m <- ranking$estimator[i]
    e <- ranking$eligible[i]
    req <- reqs[m]
    if (!e) {
      paste0("ineligible -- ", elig$reason[elig$estimator == m])
    } else {
      paste0("requires: ", req)
    }
  }, character(1))

  # Recommended estimator: the top eligible, non-naive estimator by the
  # data-driven composite robustness score (worst-estimand x confidence
  # multiplier). The ranking is already sorted eligible-first, naive-last,
  # final_score descending, so the first eligible non-naive row is the
  # recommendation. When no sensitivity surface is available there is no
  # composite; fall back to the top eligible non-naive estimator.
  elig_ranking <- ranking[ranking$eligible, ]
  elig_nonnaive <- elig_ranking[!elig_ranking$estimator %in% c("UNADJ", "DIRECT"), ]
  recommended <- if (nrow(elig_nonnaive) > 0) elig_nonnaive$estimator[1] else
    (if (nrow(elig_ranking) > 0) elig_ranking$estimator[1] else NA)
  # Per-estimand recommendations: the top eligible estimator for each
  # estimand may differ when their robustness profiles differ.
  recommended_NDE <- if (nrow(elig_ranking) > 0 && "robustness_NDE" %in% names(ranking)) {
    elig_nde <- elig_ranking[order(-elig_ranking$robustness_NDE), ]
    elig_nde$estimator[1]
  } else recommended
  recommended_NIE <- if (nrow(elig_ranking) > 0 && "robustness_NIE" %in% names(ranking)) {
    elig_nie <- elig_ranking[order(-elig_ranking$robustness_NIE), ]
    elig_nie$estimator[1]
  } else recommended

  # per-scenario top estimator.
  per_scenario <- if (!is.null(sensitivity))
    .extract_per_scenario(sensitivity, criterion = criterion) else NULL

  # Summary
  summary_txt <- .build_recommendation_summary(ranking, recommended,
                                               recommended_NDE, recommended_NIE,
                                               data)

  obj <- list(
    ranking = ranking,
    recommended = recommended,
    recommended_NDE = recommended_NDE,
    recommended_NIE = recommended_NIE,
    per_scenario = per_scenario,
    criterion = criterion,
    completeness_penalty = completeness_penalty,
    summary = summary_txt
  )
  class(obj) <- c("iconic_recommendation", "list")
  if (isTRUE(verbose))
    message("iconic_recommend complete. Call summary() or print() on the result for the full recommendation.")
  obj
}


#' Composite data-driven robustness score (internal)
#'
#' Combines per-estimand robustness into a single score and discounts each
#' estimator by the graded strength of the diagnostic assumptions it depends
#' on. Three steps:
#' \enumerate{
#'   \item \strong{Worst-estimand composite}: \code{pmin(robustness_NDE,
#'   robustness_NIE)}. An estimator is only as trustworthy as its weakest
#'   estimand, so a strong NDE cannot compensate for a weak NIE (the usual
#'   object of a mediation analysis).
#'   \item \strong{Confidence multiplier}: each estimator's composite is
#'   multiplied by a factor in \eqn{(0, 1]} read from the graded completeness
#'   verdict (\code{diagnosis$completeness$completeness}). Only estimators
#'   whose identification routes through the negative-control bridge
#'   (DIRECT, COCA, PGC, PGC2, PGC2Gm) are discounted; instrument-only
#'   estimators (IV2SLS, IV2SLS2) and UNADJ keep their full score. When no
#'   completeness object exists (no NC panel), bridge-dependent estimators
#'   are already ineligible, so the multiplier is never misapplied.
#'   \item \strong{Final score}: \code{composite * confidence_mult}.
#' }
#'
#' @param ranking Data frame with \code{estimator}, \code{robustness_NDE},
#'   \code{robustness_NIE} columns.
#' @param diagnosis An \code{iconic_diagnosis} object (may carry a NULL
#'   \code{$completeness}).
#' @param completeness_penalty Named numeric vector mapping completeness
#'   verdicts to multipliers in \eqn{[0, 1]}.
#' @return Data frame with \code{composite}, \code{confidence_mult}, and
#'   \code{final_score} aligned to \code{ranking$estimator}.
#' @keywords internal
.composite_robustness <- function(ranking, diagnosis, completeness_penalty) {
  n <- nrow(ranking)
  # Worst-estimand composite. Treat non-finite (all-NA bias, e.g. COCA on
  # survival) as worst so such estimators rank last.
  nde <- ranking$robustness_NDE
  nie <- ranking$robustness_NIE
  composite <- pmin(nde, nie)
  composite[!is.finite(composite)] <- -Inf

  # Confidence multiplier from the graded completeness verdict.
  bridge_dep <- ranking$estimator %in% c("DIRECT", "COCA", "PGC", "PGC2", "PGC2Gm")
  verdict <- if (!is.null(diagnosis) && !is.null(diagnosis$completeness) &&
                 !is.null(diagnosis$completeness$completeness))
    diagnosis$completeness$completeness else NA_character_
  # Default penalties; honour user overrides via completeness_penalty.
  pen <- c("satisfied" = 1.0, "borderline" = 0.7,
           "weak-capture" = 0.5, "under-identified" = 0)
  if (!is.null(completeness_penalty)) {
    nm <- names(completeness_penalty)
    pen[nm] <- completeness_penalty
  }
  mult <- rep(1.0, n)
  if (!is.na(verdict) && verdict %in% names(pen)) {
    mult[bridge_dep] <- pen[[verdict]]
  } else if (!is.na(verdict)) {
    # Unrecognised verdict: be conservative for bridge-dependent estimators.
    mult[bridge_dep] <- 0.5
  }
  # If there is no completeness verdict at all, bridge-dependent estimators
  # are already ineligible; leave their multiplier at 1 (moot).

  final_score <- composite * mult
  final_score[!is.finite(final_score)] <- -Inf

  data.frame(composite = composite,
             confidence_mult = mult,
             final_score = final_score,
             stringsAsFactors = FALSE)
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
                                          recommended_NDE, recommended_NIE,
                                          data) {
  lines <- character(0)

  if (is.na(recommended)) {
    lines <- c(lines, " No eligible estimators beyond UNADJ.",
               " Consider supplying instruments (G) or negative controls (W).")
  } else {
    lines <- c(lines,
      sprintf(" Recommended (composite): %s", recommended))
    # Per-estimand detail, shown when a sensitivity surface was used.
    if ("robustness_NDE" %in% names(ranking)) {
      if (!is.null(recommended_NDE) && !is.na(recommended_NDE))
        lines <- c(lines,
          sprintf("   (by NDE robustness: %s)", recommended_NDE))
      if (!is.null(recommended_NIE) && !is.na(recommended_NIE))
        lines <- c(lines,
          sprintf("   (by NIE robustness: %s)", recommended_NIE))
    }

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
    has_fin <- "final_score" %in% names(e)
    for (i in seq_len(nrow(e))) {
      rob <- if (has_rob && is.finite(e$robustness_NDE[i]))
        sprintf(" NDE %.3f, NIE %.3f;",
                e$robustness_NDE[i], e$robustness_NIE[i]) else ""
      fin <- if (has_fin && is.finite(e$final_score[i]))
        sprintf(" composite %.3f (x%.2f) = %.3f;",
                e$composite[i], e$confidence_mult[i], e$final_score[i]) else ""
      cat(sprintf(" %d. %s --%s%s %s\n",
                  e$rank[i], e$estimator[i], rob, fin, e$rationale[i]))
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
