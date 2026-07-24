# ── Survival DGP helpers (v0.9.4) ───────────────────────────────────────
# Converts a continuous linear predictor (the Y matrix that generate_toy_data
# / run_single_iteration already build) into a time-to-event outcome via an
# exponential proportional-hazards model with independent exponential
# censoring.  The log-hazard ratio for Z (and M) equals the coefficient in
# the linear predictor, so true_total / true_NDE / true_NIE are valid on the
# Cox log-HR scale.
#
# T  = -log(U) / (h0 * exp(eta))      U ~ Uniform(0,1)
# C  ~ Exp(rate = censor_rate)
# observed time   = min(T, C)
# observed event  = 1[T <= C]
#
# The baseline hazard h0 and censor_rate are chosen so that the median
# follow-up and event fraction are realistic (~50-60% events).

#' Convert a linear predictor to survival data (internal, v0.9.4)
#'
#' Given an n-vector (or n x p matrix) of linear predictors `eta` and a
#' target event fraction, returns a list with \code{surv_time} and
#' \code{surv_event} vectors.  When \code{eta} is a matrix, the first
#' column is used (survival outcomes are scalar — n_features = 1).
#'
#' @param eta numeric vector or matrix of linear predictors.
#' @param h0 baseline hazard (default 0.1).
#' @param event_frac target event fraction; \code{censor_rate} is tuned
#'   so that approximately this fraction of subjects are observed events.
#' @param censor_rate optional explicit censoring rate; if NULL it is
#'   solved from \code{event_frac}.
#' @return list(surv_time, surv_event, true_h0, true_censor_rate)
#' @keywords internal
.linpred_to_surv <- function(eta, h0 = 0.1, event_frac = 0.6,
                             censor_rate = NULL) {
  if (is.matrix(eta)) eta <- eta[, 1]
  n <- length(eta)

  # Center eta so the baseline hazard corresponds to the mean subject.
  eta_c <- eta - mean(eta)

  # True event times under exponential PH
  U <- stats::runif(n)
  T_true <- -log(U) / (h0 * exp(eta_c))

  # Censoring: solve for censor_rate that yields ~event_frac events.
  # P(event) = h0*exp(eta) / (h0*exp(eta) + censor_rate)  (in expectation).
  # Use the mean hazard as the reference for tuning.
  if (is.null(censor_rate)) {
    mean_haz <- h0 * mean(exp(eta_c))
    # P(event) = mean_haz / (mean_haz + cr) = event_frac
    # => cr = mean_haz * (1 - event_frac) / event_frac
    censor_rate <- mean_haz * (1 - event_frac) / event_frac
  }
  C <- stats::rexp(n, rate = censor_rate)

  surv_time  <- pmin(T_true, C)
  surv_event <- as.integer(T_true <= C)

  list(surv_time = surv_time, surv_event = surv_event,
       true_h0 = h0, true_censor_rate = censor_rate)
}
