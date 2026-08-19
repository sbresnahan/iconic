# ── Binary DGP helpers ───────────────────────────────────────
# Converts a continuous linear predictor (the Y matrix that
# generate_toy_data / run_single_iteration already build) into a
# binary (0/1) outcome via a logistic (Bernoulli) model. The
# coefficient on X (and M) in the linear predictor equals the
# conditional log-odds ratio, so true_total / true_NDE / true_NIE
# are valid on the logistic log-OR scale — the same conditional-scale
# convention used for the Cox log-HR truth in the survival DGP.
#
# P(Y = 1) = plogis(e0 + eta),  Y ~ Bernoulli(P)
#
# The intercept e0 is solved by root-finding so that the marginal
# event probability mean(P) matches the target prevalence `prev`.

#' Convert a linear predictor to a binary outcome
#'
#' Given an n-vector (or n x p matrix) of linear predictors \code{eta}
#' and a target prevalence, returns a list with the 0/1 outcome vector
#' \code{y_bin}. When \code{eta} is a matrix, the first column is used
#' (binary outcomes are scalar — \code{n_features = 1}).
#'
#' @param eta numeric vector or matrix of linear predictors.
#' @param prev target prevalence (marginal event probability); the
#' intercept is solved so that approximately this fraction of subjects
#' are cases. Default 0.5.
#' @return list(y_bin, p, intercept, prev)
#' @keywords internal
.linpred_to_binary <- function(eta, prev = 0.5) {
  if (is.matrix(eta)) eta <- eta[, 1]
  n <- length(eta)
  if (!is.numeric(prev) || length(prev) != 1L || !is.finite(prev) ||
      prev <= 0 || prev >= 1)
    stop("bin_prev must be a single number strictly between 0 and 1.",
         call. = FALSE)

  # Center eta so the intercept corresponds to the mean subject.
  eta_c <- eta - mean(eta)

  # Solve the intercept e0 so that mean(plogis(e0 + eta_c)) = prev.
  # plogis is monotone in e0, so the root is unique and bracketed.
  f <- function(e0) mean(stats::plogis(e0 + eta_c)) - prev
  e0 <- stats::uniroot(f, interval = c(-50, 50), tol = 1e-10)$root

  p <- stats::plogis(e0 + eta_c)
  y_bin <- as.numeric(stats::rbinom(n, 1L, p))

  list(y_bin = y_bin, p = p, intercept = e0, prev = mean(y_bin))
}
