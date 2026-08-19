# ============================================================
# Internal data-generating process for the SCENIC toy simulation.
#
# The mo_confounding parameter (default 0) controls whether the
# unmeasured confounder U1 also affects the mediator M, creating
# mediator-outcome (M-O) confounding. When mo_confounding = 0 the
# original DGP is preserved exactly (backward compatible).
#
# the phi parameter (default 0) adds a mediator-specific
# genetic instrument Gm (e.g. a cis-eQTL for a mediator transcript).
# When phi > 0, Gm affects M but is independent
# of U and has no direct path to Y, making the 2-stage MR mediation
# estimator (fit_iv2sls_mediation2) point-identified under M-O
# confounding. When phi = 0 (default), no Gm is generated and the
# original DGP is preserved exactly.
#
# imperfect instrument independence. The generates
# G and Gm as pure noise, completely independent of U -- the best
# possible case for an instrument. adds parameters that let the
# instruments be correlated with the confounders (rho_G1, rho_G2) via a
# shared population-structure factor (rho_pop), draws separate
# confounder loadings for the X->M and M->Y paths (lambda_XM/lambda_MY), and generates
# path-specific negative controls W1 and W2 with independent coverage
# (omega_1, omega_2). This sets up the "tipping point" simulation: as
# instrument exogeneity is violated, IV2SLS2 degrades while the
# negative-control-augmented PGC-2 estimator (fit_pgc_mediation2) may
# become preferable. All new parameters default to values that
# reproduce the exactly.
# ============================================================

#' Generate one synthetic dataset (internal)
#'
#' When \code{mo_confounding > 0}, the unmeasured confounder U1 also
#' affects the mediator M, creating mediator-outcome confounding — the
#' key extension for the mediation simulation. When \code{mo_confounding = 0}
#' (default), the original DGP is preserved exactly.
#'
#' When \code{phi > 0}, a mediator-specific genetic instrument
#' \code{Gm} is generated: \eqn{Gm \sim \mathcal{N}(0,1)}, independent of
#' \code{G}, \code{U}, and all other variables, and the mediator equation
#' becomes \eqn{M = \alpha_M X + \delta_{mo} 0.5 U_1 + \phi Gm + \varepsilon_M}.
#' \code{Gm} is a valid instrument for \code{M}: it moves \code{M} but is
#' independent of the confounders and has no direct path to \code{Y}. This
#' enables the 2-stage MR mediation estimator (\code{\link{fit_iv2sls_mediation2}})
#' to point-identify NDE and NIE even under M-O confounding. When
#' \code{phi = 0} (default), no \code{Gm} is generated and the original DGP
#' is preserved exactly.
#'
#' @section imperfect instrument independence:
#' The generates instruments as pure noise, independent of the
#' confounders — the best possible case. adds parameters that
#' violate this independence, modelling realistic genomic structure:
#' \describe{
#' \item{\code{rho_G1}}{Correlation of the exposure instrument G1 with the
#' X->M confounder conf_XM (instrument exogeneity violation). Default 0.}
#' \item{\code{rho_G2}}{Correlation of the mediator instrument G2 with the
#' M->Y confounder conf_MY (instrument exogeneity violation). Default 0.}
#' \item{\code{rho_pop}}{Shared population-structure factor P that induces
#' correlation between G1 and G2 (linkage / stratification). Default 0.}
#' \item{\code{lambda_XM}, \code{lambda_MY}}{Per-path confounder loading vectors; distinct values draw conf_XM and conf_MY as
#' independent confounders for the X->M and M->Y paths respectively.
#' If \code{FALSE} (default), conf_XM = conf_MY = U1 (backward compatible).}
#' \item{\code{omega_1}}{Coverage of conf_XM by the path-specific negative
#' controls W1. \code{NULL} (default) uses \code{w_signal}.}
#' \item{\code{omega_2}}{Coverage of conf_MY by the path-specific negative
#' controls W2. \code{NULL} (default) uses \code{w_signal}.}
#' }
#' When all are at their defaults, the DGP is identical
#' to the default. When any is non-default, the output additionally includes
#' \code{G1}, \code{G2}, \code{W1}, \code{W2}, \code{conf_XM}, \code{conf_MY},
#' and (when \code{rho_pop > 0}) \code{P}.
#'
#' @param n Sample size. Default 500.
#' @param n_features Number of outcome and negative-control features. Default 20.
#' @param n_mediators Number of independent mediators. When > 1,
#' each mediator M_m has its own genetic instrument Gm_m
#' and contributes additively to Y. M and Gm are returned
#' as \code{n_mediators x n} matrices. Default 1 (single
#' mediator, backward compatible).
#' @param beta_X Direct effect of X on Y (true NDE). Default 0.10.
#' @param alpha_M Effect of X on mediator M. Default 0.50.
#' @param beta_M Effect of M on Y (per-mediator NIE = alpha_M * beta_M;
#' total NIE = n_mediators * alpha_M * beta_M). Default 0.30.
#' @param conf_str Confounding strength delta. Default 0.80.
#' @param w_signal Proxy quality omega (0 = noise, 1 = perfect U proxy). Default 0.70.
#' @param mo_confounding Strength of U1 -> M (mediator-outcome confounding).
#' 0 = no M-O confounding (original DGP). Default 0.
#' @param pleio Strength of a direct G -> Y path (horizontal pleiotropy),
#' violating the exclusion restriction. 0 = no pleiotropy
#' (valid instrument, original DGP). Default 0.
#' @param phi Strength of the mediator instrument Gm -> M.
#' 0 = no mediator instrument (original DGP, no Gm generated).
#' Default 0. When > 0, a valid instrument for M is generated,
#' enabling point identification of NDE/NIE via
#' \code{\link{fit_iv2sls_mediation2}}.
#' @param gamma_G Strength of the exposure instrument G -> X. Default 0.6.
#' @param rho_G1 Correlation of G1 with conf_XM. Default 0.
#' @param rho_G2 Correlation of G2 with conf_MY. Default 0.
#' @param rho_pop Shared population structure inducing G1-G2 correlation
#'   . Default 0.
#' @param lambda_XM Optional per-path loading vector for the X->M confounder
#'   composite. NULL (default) = shared loadings.
#' @param lambda_MY Optional per-path loading vector for the M->Y confounder
#'   composite. NULL (default) = shared loadings.
#' @param omega_1 Coverage of conf_XM by W1. NULL = use w_signal.
#' @param omega_2 Coverage of conf_MY by W2. NULL = use w_signal.
#' @param feat_cor Within-module feature correlation. When > 0,
#' the idiosyncratic noise in the Y and W panels is drawn
#' from a multivariate normal with a block-diagonal
#' correlation matrix modelling co-expression modules.
#' 0 = independent noise (backward compatible). Default 0.
#' Sensitivity sweeps in the ICONIC preprint show estimator
#' performance is largely insensitive to \code{feat_cor}, so it
#' is retained for completeness rather than as a primary
#' stress axis.
#' @param u_strength Numeric scalar: scales the single
#' confounder's effect on X, M, and Y. Default NULL → 1
#' (backward compatible). See \code{\link{run_single_iteration}}
#' for the length-k vector version (multiple confounders).
#' @param w_coverage_profile A list with optional \code{w1} and \code{w2}
#' numeric vectors: per-control coverage
#' of conf_XM / conf_MY. Default NULL → scalar omega applied
#' uniformly (backward compatible).
#' @param outcome_type \code{"continuous"} (default), \code{"survival"},
#'   or \code{"binary"}.
#'   When \code{"survival"}, the continuous
#' linear predictor Y is converted to \code{surv_time}
#' and \code{surv_event} via an exponential PH model.
#'   When \code{"binary"}, the linear predictor is converted to a 0/1
#' outcome \code{y_bin} via a logistic (Bernoulli) model; the true
#' effects are then on the conditional log-odds-ratio scale.
#' @param surv_h0 Baseline hazard for the survival DGP.
#' Default 0.1.
#' @param surv_event_frac Target fraction of observed events.
#' Default 0.6. The censoring rate is solved
#' internally from this.
#' @param surv_censor_rate Explicit censoring rate. Default NULL
#' → solved from \code{surv_event_frac}.
#' @param bin_prev Target prevalence (marginal event probability) for
#' the binary DGP. Default 0.5. The logistic intercept is solved
#' internally from this.
#' @param seed Optional integer RNG seed for reproducibility.
#' @param separate_U Defunct. Passing a value errors with a message pointing
#'   to the replacement per-path loading vectors
#'   \code{lambda_XM} / \code{lambda_MY}. Retained in the signature only to
#'   catch and redirect old calls.
#' @param beta_Z Defunct. Renamed to \code{beta_X}; passing a value
#'   errors with a message pointing to \code{beta_X}. Retained in the
#'   signature only to catch and redirect old calls.
#'
#' @return A named list with elements X, G (n x n_features matrix), Y, W, U1, M,
#' synthetic_data, true_total, true_NDE, true_NIE. When \code{phi > 0}, also
#' includes \code{Gm} (numeric vector, length n, or \code{n_mediators x n}
#' matrix when \code{n_mediators > 1}). When any parameter
#' is non-default, also includes \code{G1}, \code{G2}, \code{W1}, \code{W2},
#' \code{conf_XM}, \code{conf_MY}, and (when \code{rho_pop > 0}) \code{P}.
#' When \code{n_mediators > 1}, \code{M} is an \code{n_mediators x n} matrix.
#' @export
#'
#' @examples
#' # Basic total-effect data
#' dat <- generate_toy_data(n = 200, n_features = 10, seed = 42)
#'
#' # Mediation with mediator instrument and path-specific NCs
#' dat <- generate_toy_data(n = 200, n_features = 10, phi = 0.8,
#' mo_confounding = 0.8, lambda_XM = c(1, 0), lambda_MY = c(0, 1),
#' omega_1 = 0.7, omega_2 = 0.7, seed = 42)
#' idata <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
#' G = dat$G[, 1], Gm = dat$Gm,
#' W1 = dat$W1, W2 = dat$W2)
generate_toy_data <- function(n = 500,
                              n_features = 20,
                              n_mediators = 1,
                              beta_X = 0.10,
                              alpha_M = 0.50,
                              beta_M = 0.30,
                              conf_str = 0.80,
                              w_signal = 0.70,
                              mo_confounding = 0,
                              pleio = 0,
                              phi = 0,
                              gamma_G = 0.6,
                              rho_G1 = 0,
                              rho_G2 = 0,
                              rho_pop = 0,
                              lambda_XM = NULL,
                              lambda_MY = NULL,
                              omega_1 = NULL,
                              omega_2 = NULL,
                              feat_cor = 0,
                              separate_U = NULL,
                              # per-confounder strength and
                              # per-control coverage for heterogeneous U/W.
                              u_strength = NULL,
                              w_coverage_profile = NULL,
                              # survival / time-to-event outcome.
                              # When "survival", the continuous linear predictor
                              # Y is converted to (surv_time, surv_event) via
                              # an exponential PH model with censoring. The
                              # true_total / true_NDE / true_NIE are then on
                              # the Cox log-HR scale.
                              # When "binary", Y is converted to a 0/1 outcome
                              # via a logistic (Bernoulli) model; the true
                              # effects are on the conditional log-OR scale.
                              outcome_type = c("continuous", "survival", "binary"),
                              surv_h0 = 0.1,
                              surv_event_frac = 0.6,
                              surv_censor_rate = NULL,
                              bin_prev = 0.5,
                              seed = NULL,
                              beta_Z = NULL) {
  if (!is.null(seed)) withr::local_seed(seed)
  outcome_type <- match.arg(outcome_type)

  # Deprecated-argument trap: the exposure effect was renamed beta_Z -> beta_X.
  if (!is.null(beta_Z))
    stop("argument `beta_Z` was renamed to `beta_X`; please use `beta_X = ...`.",
         call. = FALSE)

  # Deprecated-argument trap: the separate_U toggle was removed and
  # replaced by per-path confounder loadings lambda_XM / lambda_MY.
  if (!is.null(separate_U))
    stop("argument `separate_U` was removed. Confounding is now ",
         "specified by per-path loading vectors `lambda_XM` and `lambda_MY`. ",
         "Shared loadings (the default) reproduce the old separate_U = FALSE; ",
         "lambda_XM = c(1, 0), lambda_MY = c(0, 1) reproduce the old ",
         "separate_U = TRUE.", call. = FALSE)

  # ── Determine whether the is active ──
  v05_active <- rho_G1 != 0 || rho_G2 != 0 || rho_pop != 0 ||
                !is.null(lambda_XM) || !is.null(lambda_MY) ||
                !is.null(omega_1) || !is.null(omega_2)

  # u_strength scales the single confounder's effect.
  # generate_toy_data has k=1, so u_strength is a scalar (default 1).
  # (run_single_iteration supports k>1 with a length-k u_strength vector.)
  u_str <- if (is.null(u_strength)) 1 else as.numeric(rep_len(u_strength, 1))

  # w_coverage_profile — per-control coverage vectors.
  wcp <- w_coverage_profile
  if (!is.null(wcp)) {
    if (!is.list(wcp)) stop("w_coverage_profile must be a list with 'w1'/'w2'")
    if (!is.null(wcp$w1)) wcp$w1 <- rep_len(as.numeric(wcp$w1), n_features)
    if (!is.null(wcp$w2)) wcp$w2 <- rep_len(as.numeric(wcp$w2), n_features)
  }

  # ── Confounders ──
  # generate_toy_data is the k=1 toy generator: two independent latent
  # confounders U1 and U2 are drawn, and per-path scalar loadings select the
  # confounder composite driving each backdoor path. Default (both lambda
  # NULL) shares U1 across both paths (the old separate_U = FALSE). Setting
  # lambda_XM = c(1, 0), lambda_MY = c(0, 1) over the (U1, U2) space recovers
  # the old separate_U = TRUE (independent per-path confounders).
  U1 <- rnorm(n)
  U2 <- rnorm(n)
  Umat <- cbind(U1, U2)

  .resolve_lambda_toy <- function(lam) {
    if (is.null(lam)) return(c(1, 0))  # shared: load on U1 only
    lam <- rep_len(as.numeric(lam), 2)
    nrm <- sqrt(sum(lam^2))
    if (nrm == 0) stop("lambda vectors must have nonzero norm.", call. = FALSE)
    lam / nrm
  }
  lambda_XM <- .resolve_lambda_toy(lambda_XM)
  lambda_MY <- .resolve_lambda_toy(lambda_MY)

  # Per-path confounder composites (length-n).
  conf_XM <- as.numeric(Umat %*% lambda_XM)
  conf_MY <- as.numeric(Umat %*% lambda_MY)

  # ── Instruments ──
  # G ~ N(0,1), Gm ~ N(0,1), both pure noise.
  # G1 = rho_G1*conf_XM + rho_pop*P + sqrt(1-rho_G1^2-rho_pop^2)*N(0,1)
  # G2 = rho_G2*conf_MY + rho_pop*P + sqrt(1-rho_G2^2-rho_pop^2)*N(0,1)
  # The sqrt term keeps Var(G) = 1 when rho terms are within valid range.
  if (v05_active) {
    P <- rnorm(n) # shared population structure

    # Variance budget for G1: rho_G1^2 + rho_pop^2 must be <= 1
    resid_var_G1 <- max(1 - rho_G1^2 - rho_pop^2, 0)
    G1 <- rho_G1 * conf_XM + rho_pop * P + sqrt(resid_var_G1) * rnorm(n)
    G1 <- as.numeric(scale(G1))

    # G2 (mediator instrument): only generated when phi > 0
    # when n_mediators > 1, generate one instrument per mediator
    # as an n_mediators x n matrix (each an independent draw with
    # the same rho_G2/rho_pop structure).
    if (phi > 0) {
      resid_var_G2 <- max(1 - rho_G2^2 - rho_pop^2, 0)
      if (n_mediators > 1) {
        G2 <- matrix(NA_real_, n_mediators, n)
        for (mm in seq_len(n_mediators)) {
          G2[mm, ] <- rho_G2 * conf_MY + rho_pop * P +
                      sqrt(resid_var_G2) * rnorm(n)
          G2[mm, ] <- as.numeric(scale(G2[mm, ]))
        }
      } else {
        G2 <- rho_G2 * conf_MY + rho_pop * P + sqrt(resid_var_G2) * rnorm(n)
        G2 <- as.numeric(scale(G2))
      }
    } else {
      G2 <- NULL
    }

    G <- G1 # G is the exposure instrument (alias for backward compat)
    Gm <- G2 # Gm is the mediator instrument (alias for backward compat)
  } else {
    # path: pure-noise instruments
    G <- rnorm(n)
    if (phi > 0) {
      if (n_mediators > 1) {
        Gm <- matrix(rnorm(n_mediators * n), n_mediators, n)
      } else {
        Gm <- rnorm(n)
      }
    } else {
      Gm <- NULL
    }
    P <- NULL
    G1 <- G
    G2 <- Gm
  }

  # ── Exposure ──
  # X = scale(0.6*G + conf_str*conf_XM + eps)
  # conf_XM = U1 in backward-compatible mode.
  # u_str scales the confounder's effect.
  X_raw <- gamma_G * G + u_str * conf_str * conf_XM + rnorm(n, 0, 0.5)
  X <- as.numeric(scale(X_raw))

  # ── Mediator ──
  # M = alpha_M*X + mo_confounding*0.5*conf_MY + phi*G2 + eps
  # conf_MY = U1 in backward-compatible mode.
  # u_str scales the M<-conf_MY confounding.
  # when n_mediators > 1, generate M as an n_mediators x n matrix,
  # each mediator with its own Gm instrument and independent noise.
  if (n_mediators > 1) {
    M <- matrix(NA_real_, n_mediators, n)
    for (mm in seq_len(n_mediators)) {
      M[mm, ] <- alpha_M * X
      if (mo_confounding > 0) M[mm, ] <- M[mm, ] + u_str * mo_confounding * 0.5 * conf_MY
      if (phi > 0) {
        gm_mm <- if (is.matrix(Gm)) Gm[mm, ] else Gm
        M[mm, ] <- M[mm, ] + phi * gm_mm
      }
      M[mm, ] <- M[mm, ] + rnorm(n, 0, 0.05)
    }
  } else {
    M <- alpha_M * X
    if (mo_confounding > 0) M <- M + u_str * mo_confounding * 0.5 * conf_MY
    if (phi > 0) M <- M + phi * Gm
    M <- M + rnorm(n, 0, 0.05)
  }

  # ── Negative controls ──
  # single W panel, W_f = w_signal*U1 + (1-w_signal)*U2 + eps.
  # path-specific W1 (captures conf_XM) and W2 (captures conf_MY),
  # with independent coverage omega_1, omega_2.
  # when feat_cor > 0, the idiosyncratic noise is drawn from a
  # multivariate normal with a block-diagonal correlation matrix
  # (co-expression modules), so W features are correlated conditional
  # on U.
  om1 <- if (is.null(omega_1)) w_signal else omega_1
  om2 <- if (is.null(omega_2)) w_signal else omega_2
  # per-control coverage vectors (length n_features).
  # NULL w_coverage_profile → scalar applied uniformly (legacy behavior).
  om1_vec <- if (!is.null(wcp$w1)) wcp$w1 else rep(om1, n_features)
  om2_vec <- if (!is.null(wcp$w2)) wcp$w2 else rep(om2, n_features)
  ws_vec <- if (!is.null(wcp$w1)) wcp$w1 else rep(w_signal, n_features)

  # Build block-diagonal correlation matrix for noise
  cor_mat <- if (feat_cor > 0) .block_diag_cor(n_features, feat_cor) else NULL

  if (v05_active) {
    W1 <- matrix(NA_real_, n, n_features)
    W2 <- matrix(NA_real_, n, n_features)
    if (!is.null(cor_mat)) {
      # Correlated noise: draw n x p from MVN(0, Sigma * 0.3^2)
      noise1 <- MASS::mvrnorm(n, rep(0, n_features), cor_mat * 0.3^2)
      noise2 <- MASS::mvrnorm(n, rep(0, n_features), cor_mat * 0.3^2)
      for (f in seq_len(n_features)) {
        # The non-covered portion of each NC must be fresh noise, NOT the
        # other path's confounder composite. With path-specific loadings
        # (conf_XM = U1, conf_MY = U2), using U2 here would contaminate W1
        # with the M->Y confounder and bias the PGC purge. (In the legacy
        # separate_U DGP, U2 was independent of the fresh U_MY, so this was
        # harmless; with the loading-vector model it is not.)
        W1[, f] <- om1_vec[f] * conf_XM + (1 - om1_vec[f]) * rnorm(n) + noise1[, f]
        W2[, f] <- om2_vec[f] * conf_MY + (1 - om2_vec[f]) * rnorm(n) + noise2[, f]
      }
    } else {
      for (f in seq_len(n_features)) {
        W1[, f] <- om1_vec[f] * conf_XM + (1 - om1_vec[f]) * rnorm(n) + rnorm(n, 0, 0.3)
        W2[, f] <- om2_vec[f] * conf_MY + (1 - om2_vec[f]) * rnorm(n) + rnorm(n, 0, 0.3)
      }
    }
    W1 <- scale(W1)
    W2 <- scale(W2)
    # Combined W for backward-compatible estimators. When the two paths share
    # the same confounder composite (identical loadings), W1 = W2 = W exactly;
    # otherwise average the two panels so W carries both path composites.
    W <- if (identical(lambda_XM, lambda_MY)) W1 else scale((W1 + W2) / 2)
  } else {
    W <- matrix(NA_real_, n, n_features)
    if (!is.null(cor_mat)) {
      noise_w <- MASS::mvrnorm(n, rep(0, n_features), cor_mat * 0.3^2)
      for (f in seq_len(n_features))
        W[, f] <- ws_vec[f] * U1 + (1 - ws_vec[f]) * U2 + noise_w[, f]
    } else {
      for (f in seq_len(n_features))
        W[, f] <- ws_vec[f] * U1 + (1 - ws_vec[f]) * U2 + rnorm(n, 0, 0.3)
    }
    W <- scale(W)
    W1 <- W
    W2 <- W
  }

  # ── Outcome ──
  # Y_f = beta_M*M + beta_X*X + gamma_f*conf_MY + pleio*G1 + eps
  # conf_MY = U1 in backward-compatible mode.
  # when feat_cor > 0, the idiosyncratic noise is drawn from a
  # multivariate normal with a block-diagonal correlation matrix,
  # so Y features are correlated conditional on U.
  Y <- matrix(NA_real_, n, n_features)
  # when n_mediators > 1, M is a matrix (n_mediators x n);
  # sum across mediators so each contributes to Y.
  M_eff <- if (is.matrix(M)) colSums(M) else M
  for (f in seq_len(n_features)) {
    gamma_f <- runif(1, 0.4, 0.8) * conf_str
    Y[, f] <- beta_M * M_eff + beta_X * X + u_str * gamma_f * conf_MY + pleio * G
  }
  if (!is.null(cor_mat)) {
    Y_noise <- MASS::mvrnorm(n, rep(0, n_features), cor_mat * 0.2^2)
    Y <- Y + Y_noise
  } else {
    for (f in seq_len(n_features))
      Y[, f] <- Y[, f] + rnorm(n, 0, 0.2)
  }

  # survival outcome — convert the linear predictor to time-to-event.
  # The first column of Y carries the causal signal (all columns share the
  # same beta_X/beta_M/M_eff; they differ only in the confounder loading and
  # noise, which are nuisance for the survival DGP). The resulting
  # surv_time / surv_event are scalar vectors; true_total / true_NDE /
  # true_NIE are on the Cox log-HR scale.
  surv <- NULL
  if (outcome_type == "survival") {
    surv <- .linpred_to_surv(Y[, 1], h0 = surv_h0,
                             event_frac = surv_event_frac,
                             censor_rate = surv_censor_rate)
  }

  # binary outcome — convert the linear predictor to a 0/1 outcome via a
  # logistic (Bernoulli) model. As for survival, the first column of Y
  # carries the causal signal; true_total / true_NDE / true_NIE are on
  # the conditional log-OR scale.
  bin <- NULL
  if (outcome_type == "binary") {
    bin <- .linpred_to_binary(Y[, 1], prev = bin_prev)
  }

  out <- list(
    X = X,
    G = matrix(rep(G, n_features), n, n_features),
    Y = Y,
    W = W,
    U1 = U1,
    M = M,
    synthetic_data = data.frame(
      fetal_sex = rbinom(n, 1, 0.5),
      gestational_age = rnorm(n)
    ),
    true_total = beta_X + n_mediators * alpha_M * beta_M,
    true_NDE = beta_X,
    true_NIE = n_mediators * alpha_M * beta_M
  )
  if (outcome_type == "survival") {
    out$surv_time <- surv$surv_time
    out$surv_event <- surv$surv_event
    out$outcome_type <- "survival"
  } else if (outcome_type == "binary") {
    out$y_bin <- bin$y_bin
    out$outcome_type <- "binary"
  } else {
    out$outcome_type <- "continuous"
  }
  if (phi > 0) out$Gm <- Gm
  out$n_mediators <- n_mediators

  # record the DGP parameters used, so scenario_manifest() can
  # recover the modifiable-parameter scalars and feat_cor from a dat
  # object without requiring the caller to pass them again.
  out$dgp_params <- list(
    n = n, n_features = n_features, n_mediators = n_mediators,
    beta_X = beta_X, alpha_M = alpha_M, beta_M = beta_M,
    conf_str = conf_str, w_signal = w_signal,
    mo_confounding = mo_confounding, pleio = pleio, phi = phi,
    gamma_G = gamma_G, rho_G1 = rho_G1, rho_G2 = rho_G2,
    rho_pop = rho_pop,
    lambda_XM = lambda_XM, lambda_MY = lambda_MY,
    omega_1 = omega_1, omega_2 = omega_2, feat_cor = feat_cor,
    u_strength = u_str,
    w_coverage_profile = wcp,
    outcome_type = outcome_type,
    surv_h0 = surv_h0, surv_event_frac = surv_event_frac,
    surv_censor_rate = surv_censor_rate,
    bin_prev = bin_prev
  )

  # (only when the is active)
  if (v05_active) {
    out$G1 <- G1
    out$G2 <- G2
    out$W1 <- W1
    out$W2 <- W2
    out$conf_XM <- conf_XM
    out$conf_MY <- conf_MY
    if (rho_pop > 0) out$P <- P
  }

  out
}


# ── Scenario manifest ─────────────────────────
# Surfaces the ground-truth estimands and the modifiable/fixed parameter
# split so the manuscript can render an orientation table preceding the
# simulation results. Accepts either a `dat` list returned by
# `generate_toy_data()` (estimands + fixed params read from the object)
# or a bare parameter list (estimands recomputed from beta_X/alpha_M/beta_M).

#' Scenario manifest: truth and parameter ranges for a simulation
#'
#' Returns a reader-orienting summary of a simulation scenario: the
#' ground-truth estimands (NDE, NIE, total effect), the modifiable
#' parameters with their swept ranges, and the fixed parameters with
#' their values. Intended to be rendered as a table preceding simulation
#' results in the manuscript (, #582: "state the truth and the
#' parameter ranges up front to orient the reader").
#'
#' `dat_or_params` may be either:
#' \itemize{
#' \item a list returned by [generate_toy_data()], in which case the
#' estimands and fixed parameters are read from the object
#' (`true_NDE`, `true_NIE`, `true_total`, `n`, `n_features`,
#' `n_mediators`, `lambda_XM`, `lambda_MY`, `feat_cor`); or
#' \item a bare named list of DGP parameters
#' (`beta_X`, `alpha_M`, `beta_M`, `n_mediators`, `n`,
#' `n_features`, `lambda_XM`, `lambda_MY`, `feat_cor`, ...), in which case
#' the estimands are recomputed as
#' `NDE = beta_X`, `NIE = n_mediators * alpha_M * beta_M`,
#' `total = NDE + NIE`.
#' }
#'
#' The `*_grid` arguments record the swept ranges for the modifiable
#' parameters. When a grid argument is omitted or NULL, the
#' corresponding modifiable parameter is reported with its scalar value
#' only (no range).
#'
#' @param dat_or_params A `dat` list from [generate_toy_data()] or a
#' named list of DGP parameters.
#' @param conf_grid Swept values of `conf_str` (confounding strength).
#' Default NULL.
#' @param coverage_grid Swept values of NC coverage (`w_signal` /
#' `omega`). Default NULL.
#' @param mo_confounding_grid Swept values of mediator-outcome
#' confounding strength. Default NULL.
#' @param phi_grid Swept values of the mediator-instrument strength
#' `phi`. Default NULL.
#' @param rho_G1_grid Swept values of the G1-conf_XM correlation. Default NULL.
#' @param rho_G2_grid Swept values of the G2-conf_MY correlation. Default NULL.
#' @param rho_pop_grid Swept values of the population-stratification
#' correlation. Default NULL.
#' @param omega_1_grid Swept values of `omega_1` (coverage of conf_XM by W1).
#' Default NULL.
#' @param omega_2_grid Swept values of `omega_2` (coverage of conf_MY by W2).
#' Default NULL.
#' @return A named list with elements:
#' \describe{
#' \item{`estimands`}{Named numeric vector: `NDE`, `NIE`, `total`.}
#' \item{`modifiable_parameters`}{Data frame with columns
#' `parameter`, `value`, `swept_range`. `value` is the scalar
#' used (or NA when only a grid was supplied);
#' `swept_range` is a comma-separated string of grid values,
#' or NA when the parameter was held fixed.}
#' \item{`fixed_parameters`}{Data frame with columns `parameter`,
#' `value`.}
#' }
#' @export
#' @examples
#' dat <- generate_toy_data(n = 200, phi = 0.8, lambda_XM = c(1, 0), lambda_MY = c(0, 1),
#' omega_1 = 0.7, omega_2 = 0.7, seed = 42)
#' scenario_manifest(dat, conf_grid = c(0.3, 0.8),
#' coverage_grid = c(0.3, 0.7, 1))
scenario_manifest <- function(dat_or_params,
                              conf_grid = NULL,
                              coverage_grid = NULL,
                              mo_confounding_grid = NULL,
                              phi_grid = NULL,
                              rho_G1_grid = NULL,
                              rho_G2_grid = NULL,
                              rho_pop_grid = NULL,
                              omega_1_grid = NULL,
                              omega_2_grid = NULL) {
  p <- dat_or_params

  # Detect whether `p` is a generate_toy_data() output (has true_NDE) or a
  # bare parameter list.
  is_dat <- !is.null(p$true_NDE) && !is.null(p$true_NIE)

  if (is_dat) {
    nde <- p$true_NDE
    nie <- p$true_NIE
    tot <- if (!is.null(p$true_total)) p$true_total else nde + nie
    n <- attr(p$Y, "n"); if (is.null(n)) n <- nrow(p$Y)
    nf <- ncol(p$Y)
    nm <- if (!is.null(p$n_mediators)) p$n_mediators else 1L
    # Confounder structure: shared vs path-specific loadings.
    dp0 <- p$dgp_params
    if (!is.null(dp0) && !is.null(dp0$lambda_XM)) {
      conf_struct <- if (identical(dp0$lambda_XM, dp0$lambda_MY))
        "shared loadings" else "path-specific loadings"
    } else {
      conf_struct <- if (!is.null(p$conf_XM)) "path-specific loadings" else "shared loadings"
    }
    # recover scalar modifiable params + feat_cor from dgp_params
    # when generate_toy_data() recorded them.
    dp <- p$dgp_params
    if (!is.null(dp)) {
      fcor <- dp$feat_cor
      conf_val <- dp$conf_str
      mo_val <- dp$mo_confounding
      phi_val <- dp$phi
      cov_val <- if (!is.null(dp$omega_1)) dp$omega_1 else dp$w_signal
      rho1 <- dp$rho_G1; rho2 <- dp$rho_G2; rhop <- dp$rho_pop
    } else {
      fcor <- attr(p$Y, "feat_cor"); if (is.null(fcor)) fcor <- 0
      conf_val <- NA_real_; mo_val <- NA_real_; phi_val <- NA_real_
      cov_val <- NA_real_; rho1 <- NA_real_; rho2 <- NA_real_; rhop <- NA_real_
    }
  } else {
    beta_X <- if (!is.null(p$beta_X)) p$beta_X else 0.10
    alpha_M <- if (!is.null(p$alpha_M)) p$alpha_M else 0.50
    beta_M <- if (!is.null(p$beta_M)) p$beta_M else 0.30
    nm <- if (!is.null(p$n_mediators)) p$n_mediators else 1L
    nde <- beta_X
    nie <- nm * alpha_M * beta_M
    tot <- nde + nie
    n <- if (!is.null(p$n)) p$n else 500
    nf <- if (!is.null(p$n_features)) p$n_features else 20
    conf_struct <- if (!is.null(p$lambda_XM) && !is.null(p$lambda_MY)) {
      if (identical(p$lambda_XM, p$lambda_MY)) "shared loadings" else "path-specific loadings"
    } else "shared loadings"
    fcor <- if (!is.null(p$feat_cor)) p$feat_cor else 0
    conf_val <- if (!is.null(p$conf_str)) p$conf_str else NA_real_
    mo_val <- if (!is.null(p$mo_confounding)) p$mo_confounding else NA_real_
    phi_val <- if (!is.null(p$phi)) p$phi else NA_real_
    cov_val <- if (!is.null(p$w_signal)) p$w_signal else
                if (!is.null(p$omega_1)) p$omega_1 else NA_real_
    rho1 <- if (!is.null(p$rho_G1)) p$rho_G1 else NA_real_
    rho2 <- if (!is.null(p$rho_G2)) p$rho_G2 else NA_real_
    rhop <- if (!is.null(p$rho_pop)) p$rho_pop else NA_real_
  }

  fmt_range <- function(g) {
    if (is.null(g) || length(g) == 0) NA_character_
    else paste(format(g, trim = TRUE, scientific = FALSE), collapse = ", ")
  }

  mod <- data.frame(
    parameter = c("conf_str", "w_signal / omega",
                  "mo_confounding", "phi",
                  "rho_G1", "rho_G2", "rho_pop",
                  "omega_1", "omega_2"),
    value = c(conf_val, cov_val, mo_val, phi_val, rho1, rho2, rhop,
              cov_val, cov_val),
    swept_range = c(fmt_range(conf_grid), fmt_range(coverage_grid),
                    fmt_range(mo_confounding_grid), fmt_range(phi_grid),
                    fmt_range(rho_G1_grid), fmt_range(rho_G2_grid),
                    fmt_range(rho_pop_grid),
                    fmt_range(omega_1_grid), fmt_range(omega_2_grid)),
    stringsAsFactors = FALSE
  )

  fixed <- data.frame(
    parameter = c("n", "n_features", "n_mediators",
                  "confounder_structure", "feat_cor"),
    value = c(n, nf, nm, conf_struct, fcor),
    stringsAsFactors = FALSE
  )

  list(
    estimands = c(NDE = nde, NIE = nie, total = tot),
    modifiable_parameters = mod,
    fixed_parameters = fixed
  )
}


#' Build a block-diagonal correlation matrix (internal)
#'
#' Creates a \code{p x p} block-diagonal correlation matrix modelling
#' co-expression modules: features within the same block have pairwise
#' correlation \code{rho}, features in different blocks are uncorrelated.
#' The number of blocks is \code{ceiling(sqrt(p))}, with block sizes as
#' equal as possible.
#'
#' @param p Number of features.
#' @param rho Within-block pairwise correlation (off-diagonal).
#' @return A \code{p x p} numeric matrix with unit diagonal.
#' @keywords internal
.block_diag_cor <- function(p, rho) {
  if (p < 2 || rho == 0) return(NULL)
  n_modules <- ceiling(sqrt(p))
  # Distribute features as evenly as possible across modules
  base_size <- p %/% n_modules
  remainder <- p %% n_modules
  sizes <- c(rep(base_size + 1, remainder), rep(base_size, n_modules - remainder))
  sizes <- sizes[sizes > 0]

  # Build the full matrix block by block
  Sigma <- matrix(0, p, p)
  idx <- 1
  for (sz in sizes) {
    if (sz > 0) {
      block <- matrix(rho, sz, sz)
      diag(block) <- 1
      Sigma[idx:(idx + sz - 1), idx:(idx + sz - 1)] <- block
      idx <- idx + sz
    }
  }
  Sigma
}
