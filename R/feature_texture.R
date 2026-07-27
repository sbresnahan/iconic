# ============================================================
#. The copula texture model documentation moves to a
# supplementary methods section when the draft is condensed. No code
# change.
# ============================================================
# Feature-level texture model for the mediator (M) panel.
#
# Learns the full joint distribution of the mediator panel using a
# Gaussian copula with flexible marginal distributions. Unlike the
# sample-level GAN (which learns per-sample summaries), this model
# captures the per-feature marginal shapes (skewness, tails, bounded
# support) and their cross-feature dependence structure.
#
# Marginal fitting:
# - "empirical": empirical CDF (step function from observed data)
# - "parametric": best fit among {normal, log-normal, gamma, beta} by AIC
# - "auto" (default): fit parametric families, run KS test; use the
# best parametric fit if it passes KS at p > 0.05, otherwise fall
# back to empirical CDF. The user can override with marginal_method.
#
# Dependence:
# Gaussian copula — the p x p correlation matrix of the normal scores
# (rank-transformed to uniform, then qnorm). Near-PD corrected.
#
# Sampling:
# Draw n samples from MVN(0, copula_cor), then transform each feature
# through its inverse marginal CDF. The result is centered and scaled
# so the texture is zero-mean noise; the structural signal provides the
# mean in run_single_iteration().
# ============================================================


#' Fit a parametric distribution to a data vector (internal)
#'
#' Attempts to fit normal, log-normal, gamma, and beta distributions
#' (where applicable) by maximum likelihood, and returns the best fit
#' by AIC along with the KS test p-value.
#'
#' @param x Numeric vector (the feature's observed values).
#' @return A list with \code{family}, \code{params}, \code{aic}, \code{ks_p},
#' or NULL if no parametric fit succeeds.
#' @keywords internal
.fit_parametric_marginal <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 5) return(NULL)

  fits <- list()

  # Normal
  fits$normal <- tryCatch({
    mu <- mean(x); sigma <- stats::sd(x)
    if (sigma <= 0) sigma <- 1e-8
    loglik <- sum(stats::dnorm(x, mu, sigma, log = TRUE))
    aic <- -2 * loglik + 2 * 2
    ks_p <- tryCatch(stats::ks.test(x, "pnorm", mu, sigma)$p.value, error = function(e) 0)
    list(family = "normal", params = list(mu = mu, sigma = sigma),
         aic = aic, ks_p = ks_p)
  }, error = function(e) NULL)

  # Log-normal (only if all positive)
  if (all(x > 0)) {
    fits$lognormal <- tryCatch({
      lx <- log(x)
      mu <- mean(lx); sigma <- stats::sd(lx)
      if (sigma <= 0) sigma <- 1e-8
      loglik <- sum(stats::dlnorm(x, mu, sigma, log = TRUE))
      aic <- -2 * loglik + 2 * 2
      ks_p <- tryCatch(stats::ks.test(x, "plnorm", mu, sigma)$p.value, error = function(e) 0)
      list(family = "lognormal", params = list(mu = mu, sigma = sigma),
           aic = aic, ks_p = ks_p)
    }, error = function(e) NULL)
  }

  # Gamma (only if all positive)
  if (all(x > 0)) {
    fits$gamma <- tryCatch({
      fit <- MASS::fitdistr(x, "gamma")
      shape <- fit$estimate["shape"]; rate <- fit$estimate["rate"]
      loglik <- as.numeric(fit$loglik)
      aic <- -2 * loglik + 2 * 2
      ks_p <- tryCatch(
        stats::ks.test(x, "pgamma", shape = shape, rate = rate)$p.value,
        error = function(e) 0)
      list(family = "gamma",
           params = list(shape = as.numeric(shape), rate = as.numeric(rate)),
           aic = aic, ks_p = ks_p)
    }, error = function(e) NULL)
  }

  # Beta (only if all in (0, 1))
  if (all(x > 0 & x < 1)) {
    fits$beta <- tryCatch({
      fit <- MASS::fitdistr(x, dbeta, start = list(shape1 = 2, shape2 = 2),
                            lower = c(0.01, 0.01))
      s1 <- fit$estimate["shape1"]; s2 <- fit$estimate["shape2"]
      loglik <- as.numeric(fit$loglik)
      aic <- -2 * loglik + 2 * 2
      ks_p <- tryCatch(
        stats::ks.test(x, "pbeta", shape1 = s1, shape2 = s2)$p.value,
        error = function(e) 0)
      list(family = "beta",
           params = list(shape1 = as.numeric(s1), shape2 = as.numeric(s2)),
           aic = aic, ks_p = ks_p)
    }, error = function(e) NULL)
  }

  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (!length(fits)) return(NULL)

  # Select best by AIC
  aics <- vapply(fits, function(f) f$aic, numeric(1))
  best <- fits[[which.min(aics)]]
  best
}


#' Inverse CDF for a fitted marginal (internal)
#'
#' Returns the quantile function for a parametric or empirical marginal.
#'
#' @param marginal A marginal specification list from .fit_parametric_marginal
#' or an empirical CDF object.
#' @param p Vector of probabilities in (0, 1).
#' @return Numeric vector of quantiles.
#' @keywords internal
.marginal_quantile <- function(marginal, p) {
  # Clamp to avoid Inf at boundaries
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)

  if (marginal$type == "parametric") {
    fam <- marginal$family
    prm <- marginal$params
    switch(fam,
      normal = stats::qnorm(p, prm$mu, prm$sigma),
      lognormal = stats::qlnorm(p, prm$mu, prm$sigma),
      gamma = stats::qgamma(p, shape = prm$shape, rate = prm$rate),
      beta = stats::qbeta(p, shape1 = prm$shape1, shape2 = prm$shape2),
      # Fallback: normal
      stats::qnorm(p)
    )
  } else if (marginal$type == "empirical") {
    # Empirical quantile via inverse of the ecdf
    # marginal$sorted_values are the sorted observed data
    n <- length(marginal$sorted_values)
    idx <- p * (n + 1)
    lo <- floor(idx); hi <- ceiling(idx)
    lo <- pmax(pmin(lo, n), 1)
    hi <- pmax(pmin(hi, n), 1)
    frac <- idx - lo
    marginal$sorted_values[lo] * (1 - frac) + marginal$sorted_values[hi] * frac
  } else {
    stats::qnorm(p)
  }
}


#' CDF for a fitted marginal (internal)
#'
#' Returns the cumulative distribution function for a parametric or
#' empirical marginal.
#'
#' @param marginal A marginal specification list.
#' @param q Numeric vector of quantiles.
#' @return Numeric vector of probabilities.
#' @keywords internal
.marginal_cdf <- function(marginal, q) {
  if (marginal$type == "parametric") {
    fam <- marginal$family
    prm <- marginal$params
    switch(fam,
      normal = stats::pnorm(q, prm$mu, prm$sigma),
      lognormal = stats::plnorm(q, prm$mu, prm$sigma),
      gamma = stats::pgamma(q, shape = prm$shape, rate = prm$rate),
      beta = stats::pbeta(q, shape1 = prm$shape1, shape2 = prm$shape2),
      stats::pnorm(q)
    )
  } else if (marginal$type == "empirical") {
    # Empirical CDF: proportion of observed values <= q
    sv <- marginal$sorted_values
    n <- length(sv)
    vapply(q, function(qi) sum(sv <= qi) / n, numeric(1))
  } else {
    stats::pnorm(q)
  }
}


#' Near-PD correction for a correlation matrix (internal)
#'
#' Clips eigenvalues to a small positive floor and renormalizes to unit
#' diagonal. Same approach as .residual_correlation() in load_data.R.
#'
#' @param cor_mat A symmetric matrix (intended to be a correlation matrix).
#' @return A valid correlation matrix (symmetric, positive-definite, unit diagonal).
#' @keywords internal
.make_pd_cor <- function(cor_mat) {
  if (is.null(cor_mat)) return(NULL)
  cor_mat[is.na(cor_mat)] <- 0
  diag(cor_mat) <- 1
  eig <- eigen(cor_mat, symmetric = TRUE)
  eig_vals <- pmax(eig$values, 1e-6)
  cor_mat <- eig$vectors %*% diag(eig_vals) %*% t(eig$vectors)
  dd <- sqrt(diag(cor_mat))
  cor_mat <- cor_mat / outer(dd, dd)
  diag(cor_mat) <- 1
  cor_mat
}


#' Train a feature-level texture model for the mediator panel
#'
#' Learns the joint distribution of a features x samples mediator matrix
#' using a Gaussian copula with flexible marginal distributions. For each
#' feature, the marginal is fitted as either an empirical CDF (default) or
#' the best parametric family (normal, log-normal, gamma, beta) selected by
#' AIC with a Kolmogorov-Smirnov goodness-of-fit check. Cross-feature
#' dependence is captured by the Gaussian copula correlation matrix.
#'
#' The resulting model can be sampled via [sample_feature_texture()] to
#' produce realistic synthetic feature vectors that preserve the marginal
#' shapes and correlation structure of the user's mediator panel.
#'
#' @param M_matrix Mediator panel, `features x samples` (e.g. placental
#' transcript expression, one row per transcript, one column per sample).
#' @param marginal_method Marginal fitting method: \code{"auto"} (default,
#' uses parametric if KS test passes at p > 0.05, otherwise empirical),
#' \code{"empirical"} (always use empirical CDF), or \code{"parametric"}
#' (always use best parametric fit by AIC).
#'
#' @return An \code{iconic_feature_texture} S3 object: a named list with
#' \code{marginals} (list of per-feature marginal specs), \code{copula_cor}
#' (p x p correlation matrix), \code{n_features}, \code{n_samples},
#' \code{marginal_method}, and \code{marginal_types} (summary of how many
#' features used each method).
#' @export
#'
#' @examples
#' \donttest{
#' M <- matrix(rnorm(30 * 200), 30, 200) # 30 transcripts, 200 samples
#' ft <- train_feature_texture(M)
#' draws <- sample_feature_texture(ft, 100) # 30 x 100 synthetic draws
#' }
train_feature_texture <- function(M_matrix, marginal_method = "auto") {
  M_matrix <- as.matrix(M_matrix)
  p <- nrow(M_matrix)
  n <- ncol(M_matrix)

  if (p < 1 || n < 5)
    stop("M_matrix must have at least 1 feature and 5 samples.")

  # ── Fit marginals for each feature ──
  marginals <- vector("list", p)
  marginal_types <- character(p)

  for (f in seq_len(p)) {
    x <- as.numeric(M_matrix[f, ])
    x <- x[is.finite(x)]

    if (marginal_method == "empirical") {
      marginals[[f]] <- list(type = "empirical",
                             sorted_values = sort(x))
      marginal_types[f] <- "empirical"
      next
    }

    if (marginal_method == "parametric") {
      fit <- .fit_parametric_marginal(x)
      if (!is.null(fit)) {
        marginals[[f]] <- list(type = "parametric",
                               family = fit$family,
                               params = fit$params,
                               ks_p = fit$ks_p)
        marginal_types[f] <- fit$family
      } else {
        marginals[[f]] <- list(type = "empirical",
                               sorted_values = sort(x))
        marginal_types[f] <- "empirical"
      }
      next
    }

    # "auto": try parametric, check KS, fall back to empirical
    fit <- .fit_parametric_marginal(x)
    if (!is.null(fit) && fit$ks_p > 0.05) {
      marginals[[f]] <- list(type = "parametric",
                             family = fit$family,
                             params = fit$params,
                             ks_p = fit$ks_p)
      marginal_types[f] <- fit$family
    } else {
      marginals[[f]] <- list(type = "empirical",
                             sorted_values = sort(x))
      marginal_types[f] <- "empirical"
    }
  }

  # ── Gaussian copula: correlation of normal scores ──
  # Transform each feature to uniform via its marginal CDF, then to normal.
  normal_scores <- matrix(NA_real_, p, n)
  for (f in seq_len(p)) {
    u <- .marginal_cdf(marginals[[f]], as.numeric(M_matrix[f, ]))
    # Clamp to avoid Inf at exact 0 or 1
    u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
    normal_scores[f, ] <- stats::qnorm(u)
  }

  copula_cor <- tryCatch(
    cor(t(normal_scores), use = "pairwise.complete.obs"),
    error = function(e) NULL
  )
  copula_cor <- .make_pd_cor(copula_cor)
  if (is.null(copula_cor)) {
    # Degenerate: use identity (independent features)
    copula_cor <- diag(p)
  }

  structure(
    list(
      marginals = marginals,
      copula_cor = copula_cor,
      n_features = p,
      n_samples = n,
      marginal_method = marginal_method,
      marginal_types = marginal_types
    ),
    class = "iconic_feature_texture"
  )
}


#' Draw synthetic feature vectors from a trained feature texture model
#'
#' Samples from the Gaussian copula: draws n_samples from
#' MVN(0, copula_cor), then transforms each feature through its inverse
#' marginal CDF. The result is centered and scaled to zero mean and unit
#' variance per feature, so the texture acts as noise that the structural
#' signal in [run_single_iteration()] provides the mean for.
#'
#' @param feature_texture An \code{iconic_feature_texture} object from
#' [train_feature_texture()].
#' @param n_samples Number of synthetic samples to draw.
#' @param n_features Target number of features. If NULL, uses the number
#' of features in the training data. If larger, additional features are
#' drawn by sampling existing columns with replacement and adding
#' independent noise. If smaller, the first n_features are used.
#'
#' @return A \code{n_features x n_samples} matrix of synthetic texture
#' values, centered and scaled per feature.
#' @export
#'
#' @examples
#' \donttest{
#' M <- matrix(rnorm(30 * 200), 30, 200) # 30 transcripts, 200 samples
#' ft <- train_feature_texture(M)
#' draws <- sample_feature_texture(ft, 500, n_features = 20)
#' }
sample_feature_texture <- function(feature_texture, n_samples, n_features = NULL) {
  stopifnot(inherits(feature_texture, "iconic_feature_texture"))

  p_train <- feature_texture$n_features
  p <- if (is.null(n_features)) p_train else n_features

  # Draw from the Gaussian copula
  # If p <= p_train: use the first p features of the copula correlation
  # If p > p_train: draw p_train, then expand
  p_use <- min(p, p_train)
  cor_sub <- feature_texture$copula_cor[seq_len(p_use), seq_len(p_use), drop = FALSE]

  z <- MASS::mvrnorm(n_samples, mu = rep(0, p_use), Sigma = cor_sub)
  if (n_samples == 1) z <- matrix(z, nrow = 1)

  # Transform through inverse marginal CDFs
  draws <- matrix(NA_real_, p_use, n_samples)
  for (f in seq_len(p_use)) {
    u <- stats::pnorm(z[, f])
    draws[f, ] <- .marginal_quantile(feature_texture$marginals[[f]], u)
  }

  # Expand if p > p_train: sample additional features from existing ones
  if (p > p_train) {
    extra <- matrix(NA_real_, p - p_train, n_samples)
    for (f in seq_len(p - p_train)) {
      src <- sample.int(p_train, 1)
      extra[f, ] <- draws[src, ] + rnorm(n_samples, 0, 0.1)
    }
    draws <- rbind(draws, extra)
  }

  # Center and scale per feature (texture = zero-mean noise)
  for (f in seq_len(p)) {
    mu <- mean(draws[f, ], na.rm = TRUE)
    sd_f <- stats::sd(draws[f, ], na.rm = TRUE)
    if (is.na(sd_f) || sd_f == 0) sd_f <- 1
    draws[f, ] <- (draws[f, ] - mu) / sd_f
  }

  draws
}


#' Print method for iconic_feature_texture objects
#'
#' @param x An `iconic_feature_texture` object.
#' @param ... Unused.
#' @return Invisibly returns `x` (the `iconic_feature_texture` object); called for its side effect of printing a human-readable summary.
#' @export
print.iconic_feature_texture <- function(x, ...) {
  cat("<iconic_feature_texture>\n")
  cat(" features :", x$n_features, "\n")
  cat(" samples :", x$n_samples, "\n")
  cat(" method :", x$marginal_method, "\n")

  # Summarise marginal types
  type_tab <- table(x$marginal_types)
  type_str <- paste(sprintf("%s (%d)", names(type_tab), as.integer(type_tab)),
                    collapse = ", ")
  cat(" marginals :", type_str, "\n")

  # Copula correlation summary
  cor_vals <- x$copula_cor[lower.tri(x$copula_cor)]
  cat(sprintf(" copula cor: mean=%.3f, range=[%.3f, %.3f]\n",
              mean(cor_vals, na.rm = TRUE),
              min(cor_vals, na.rm = TRUE),
              max(cor_vals, na.rm = TRUE)))

  invisible(x)
}
