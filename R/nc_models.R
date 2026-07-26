# ============================================================
# Pluggable negative-control (NC) mechanisms.
#
# The package does not hard-code how negative controls arise. A
# negative-control model is any function with the signature
#
# function(U, covariates, params) -> W
#
# where
# U : n x k matrix of latent confounders (columns are confounders)
# covariates : data frame of n observed covariates (may have 0 columns)
# params : named list of model-specific settings; always contains
# `n_features` (how many control columns to return)
#
# and the return value is an n x n_features numeric matrix of negative
# controls. Ship-with built-ins: `nc_proxy` (direct confounder proxy) and
# `nc_cpg` (CpG-predicted expression, the SCENIC case). Users may register
# their own or pass a function directly to run_single_iteration().
#
# The scientific point of the interface: negative controls are only valid
# when they span the SAME confounder subspace that drives the outcome. Each
# built-in exposes a `coverage`/`captured` knob controlling exactly that, so
# the sensitivity analysis can sweep NC validity from perfect to broken.
#
# nc_proxy() gains a `mode` parameter. mode = "shared" (default,
# backward-compatible) gives every column the same rowMeans(U[, captured])
# signal. mode = "distinct" assigns each column to a DIFFERENT confounder,
# so the W matrix has genuine dimensional structure — this is what makes
# the proximal completeness condition (dim(W_valid) >= k) binding for the
# matrix-bridge PGC estimator.
# ============================================================


#' Direct-proxy negative-control model
#'
#' Each control is a coverage-weighted mixture of the captured confounders plus
#' noise. With one confounder and `coverage = omega` this reduces to
#' the classic `w_signal` proxy in [generate_toy_data()]. Setting `captured` to a
#' strict subset of the confounders models controls that miss part of the
#' confounder subspace (invalidating the negative-control assumption).
#'
#' When `noise_cor` is supplied (a `p x p` correlation matrix), the noise
#' component is drawn from a multivariate normal with that correlation
#' structure, so the negative controls retain realistic cross-feature
#' correlations conditional on the confounder. When `noise_cor` is `NULL`
#' (default), the noise is independent across features.
#'
#' @param U `n x k` confounder matrix.
#' @param covariates Covariate data frame (unused; kept for the NC contract).
#' @param params List with `n_features`, and optionally `coverage`
#' (scalar in `[0, 1]`, default 0.7), `captured` (integer confounder indices
#' the controls see, default all), `noise_sd` (default 0.3), `MMCon`
#' (loading multiplier, default 1), `mode` (default `"shared"`), and
#' `noise_cor` (a `p x p` correlation matrix for correlated noise, or
#' `NULL` for independent noise).
#'
#' @section Modes:
#' \describe{
#' \item{`"shared"`}{All columns carry the same `rowMeans(U[, captured])`
#' signal. This is the original behaviour (backward-compatible) and is
#' numerically stable, but the W matrix has only one effective dimension
#' regardless of `n_features`, so the proximal completeness condition is
#' never binding for the matrix-bridge PGC.}
#' \item{`"distinct"`}{Column `f` captures confounder
#' `captured[((f - 1) \%\% length(captured)) + 1]`. Different columns
#' therefore carry signals from different confounders, giving the W
#' matrix genuine dimensional structure. This is the mode to use when
#' benchmarking the completeness cliff: the matrix-bridge PGC is
#' identified only when `n_features >= k`.}
#' }
#'
#' @return `n x n_features` matrix of negative controls.
#' @export
nc_proxy <- function(U, covariates, params) {
  n <- nrow(U); k <- ncol(U)
  p <- params$n_features
  cov <- if (is.null(params$coverage)) 0.7 else params$coverage
  cap <- if (is.null(params$captured)) seq_len(k) else params$captured
  sd0 <- if (is.null(params$noise_sd)) 0.3 else params$noise_sd
  mm <- if (is.null(params$MMCon)) 1 else params$MMCon
  mode <- if (is.null(params$mode)) "shared" else params$mode
  noise_cor <- params$noise_cor # NULL or p x p correlation matrix

  cap <- cap[cap >= 1 & cap <= k]

  # ── Generate the noise component ──
  # When noise_cor is a valid p x p correlation matrix, draw correlated
  # noise via MASS::mvrnorm. Otherwise use independent rnorm (backward
  # compatible). The noise has two parts: (1 - cov) * noise (the
  # "uncaptured" component) and rnorm(n, 0, sd0) (idiosyncratic). We
  # apply the correlation to the combined noise term.
  noise_mat <- .generate_nc_noise(n, p, noise_cor, sd0)

  W <- matrix(NA_real_, n, p)

  if (mode == "distinct" && length(cap) > 0) {
    # Each column captures a DIFFERENT confounder (cycling through `cap`).
    # This gives W genuine dimensional structure so that the matrix-bridge
    # PGC requires ncol(W) >= k for identification.
    for (f in seq_len(p)) {
      conf_idx <- cap[((f - 1) %% length(cap)) + 1]
      signal <- U[, conf_idx]
      W[, f] <- mm * cov * signal +
                (1 - cov) * noise_mat[, f] +
                noise_mat[, f + p]
    }
  } else {
    # Original "shared" mode: all columns carry the same averaged signal.
    signal <- if (length(cap)) rowMeans(U[, cap, drop = FALSE]) else rep(0, n)
    for (f in seq_len(p)) {
      W[, f] <- mm * cov * signal +
                (1 - cov) * noise_mat[, f] +
                noise_mat[, f + p]
    }
  }

  scale(W)
}


#' Generate correlated noise for negative-control models (internal)
#'
#' When `noise_cor` is a valid `p x p` correlation matrix, draws `n` samples
#' from `MVN(0, noise_cor)` and scales by `sd0`. Returns an `n x (2*p)`
#' matrix: columns `1:p` are the "uncaptured" noise component and columns
#' `(p+1):(2*p)` are the idiosyncratic noise component (both with the same
#' correlation structure). When `noise_cor` is NULL, uses independent
#' `rnorm` (backward compatible).
#'
#' @param n Number of samples.
#' @param p Number of features.
#' @param noise_cor A `p x p` correlation matrix, or NULL.
#' @param sd0 Idiosyncratic noise SD.
#' @return An `n x (2*p)` numeric matrix.
#' @keywords internal
.generate_nc_noise <- function(n, p, noise_cor, sd0) {
  use_cor <- !is.null(noise_cor) &&
             is.matrix(noise_cor) &&
             nrow(noise_cor) == p && ncol(noise_cor) == p

  if (use_cor) {
    # Draw 2*p correlated noise columns: first p for the (1-cov) component,
    # next p for the idiosyncratic component.
    # Use sd0 for idiosyncratic, unit variance for the (1-cov) component.
    Sigma1 <- noise_cor # unit variance
    Sigma2 <- noise_cor * (sd0^2) # scaled by sd0^2
    noise1 <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma1)
    noise2 <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma2)
    if (n == 1) { noise1 <- matrix(noise1, nrow = 1); noise2 <- matrix(noise2, nrow = 1) }
    cbind(noise1, noise2)
  } else {
    # Independent noise (backward compatible)
    noise1 <- matrix(rnorm(n * p), n, p)
    noise2 <- matrix(rnorm(n * p, 0, sd0), n, p)
    cbind(noise1, noise2)
  }
}


#' CpG-predicted-expression negative-control model (SCENIC case)
#'
#' Simulates spatially correlated CpG methylation whose signal is partly driven
#' by the captured confounders, then forms each negative control as a linear
#' prediction from the methylation sites ("CpG-predicted expression"). The
#' controls therefore carry confounder information only to the extent the
#' methylation does, mediated through a realistic spatial methylation layer.
#'
#' When `noise_cor` is supplied (a `p x p` correlation matrix), the
#' idiosyncratic noise added to each control is drawn from a multivariate
#' normal with that correlation structure, so the controls retain realistic
#' cross-feature correlations conditional on the confounder.
#'
#' @param U `n x k` confounder matrix.
#' @param covariates Covariate data frame (unused; kept for the NC contract).
#' @param params List with `n_features`, and optionally `coverage`
#' (confounder->methylation strength, default 0.7), `captured` (confounder
#' indices, default all), `n_cpg` (methylation sites, default 60), `rho`
#' (AR(1) spatial correlation across sites, default 0.6), `MMCpG`
#' (methylation-confounding multiplier, default 1), `MMCon` (default 1),
#' and `noise_cor` (a `p x p` correlation matrix for correlated
#' idiosyncratic noise, or `NULL` for independent noise).
#'
#' @return `n x n_features` matrix of CpG-predicted negative controls.
#' @export
nc_cpg <- function(U, covariates, params) {
  n <- nrow(U); k <- ncol(U)
  p <- params$n_features
  cov <- if (is.null(params$coverage)) 0.7 else params$coverage
  cap <- if (is.null(params$captured)) seq_len(k) else params$captured
  ncpg<- if (is.null(params$n_cpg)) 60 else params$n_cpg
  rho <- if (is.null(params$rho)) 0.6 else params$rho
  mmc <- if (is.null(params$MMCpG)) 1 else params$MMCpG
  mm <- if (is.null(params$MMCon)) 1 else params$MMCon
  noise_cor <- params$noise_cor # NULL or p x p correlation matrix

  cap <- cap[cap >= 1 & cap <= k]
  conf <- if (length(cap)) rowMeans(U[, cap, drop = FALSE]) else rep(0, n)

  # Spatially correlated methylation: AR(1) across neighbouring CpG sites,
  # with a confounder-driven component shared across all sites.
  eps <- matrix(rnorm(n * ncpg), n, ncpg)
  meth <- matrix(NA_real_, n, ncpg)
  meth[, 1] <- eps[, 1]
  for (s in 2:ncpg)
    meth[, s] <- rho * meth[, s - 1] + sqrt(1 - rho^2) * eps[, s]
  meth <- meth + mmc * cov * matrix(rep(conf, ncpg), n, ncpg)

  # Each negative control is a sparse linear prediction from the CpG sites,
  # plus idiosyncratic noise (optionally correlated across features).
  noise_mat <- .generate_nc_noise(n, p, noise_cor, 0.3)

  W <- matrix(NA_real_, n, p)
  for (f in seq_len(p)) {
    wts <- rnorm(ncpg) * rbinom(ncpg, 1, 0.3) # sparse predictor weights
    if (all(wts == 0)) wts[sample.int(ncpg, 1)] <- 1
    W[, f] <- mm * as.numeric(meth %*% wts) + noise_mat[, f + p]
  }
  scale(W)
}


# Registry of built-in negative-control models by name.
.iconic_nc_registry <- list(proxy = nc_proxy, cpg = nc_cpg)

#' Resolve a negative-control model argument to a function (internal)
#'
#' Accepts a function (used as-is), or a registered name (`"proxy"`, `"cpg"`).
#' @keywords internal
.resolve_nc_model <- function(nc_model) {
  if (is.function(nc_model)) return(nc_model)
  if (is.character(nc_model) && length(nc_model) == 1 &&
      nc_model %in% names(.iconic_nc_registry))
    return(.iconic_nc_registry[[nc_model]])
  stop("`nc_model` must be a function or one of: ",
       paste(names(.iconic_nc_registry), collapse = ", "))
}

#' Validate that a negative-control model honours the interface (internal)
#'
#' Calls the model on tiny inputs and checks the returned shape.
#' @keywords internal
.validate_nc_model <- function(nc_model, k = 2, n = 10, p = 3) {
  f <- .resolve_nc_model(nc_model)
  U <- matrix(rnorm(n * k), n, k)
  W <- f(U, data.frame(row.names = seq_len(n)), list(n_features = p))
  if (!is.matrix(W) || nrow(W) != n || ncol(W) != p)
    stop("negative-control model must return an n x n_features numeric matrix.")
  invisible(TRUE)
}

#' List the built-in negative-control models
#' @return Character vector of registered NC model names.
#' @export
list_nc_models <- function() names(.iconic_nc_registry)


#' Simulate a single genetic instrument
#'
#' Draws one standardised genetic instrument for the exposure, as an additive
#' score over `n_snps` biallelic SNPs (dosages `Binomial(2, maf)`). Mirrors the
#' single-instrument (PRS) design used throughout the SCENIC framework.
#'
#' @param n Number of samples.
#' @param n_snps SNPs contributing to the score. Default 20.
#' @param maf Minor-allele frequency. Default 0.3.
#' @param seed Optional RNG seed.
#'
#' @return List with `G` (standardised instrument, length `n`), `dosages`
#' (`n x n_snps`), and `maf`.
#' @export
simulate_single_genetic_instrument <- function(n, n_snps = 20, maf = 0.3,
                                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dosages <- matrix(rbinom(n * n_snps, 2, maf), n, n_snps)
  score <- rowSums(dosages)
  list(G = as.numeric(scale(score)), dosages = dosages, maf = maf)
}
