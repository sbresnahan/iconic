# ============================================================
# Internal data-generating process for the SCENIC toy simulation.
#
# The mo_confounding parameter (default 0) controls whether the
# unmeasured confounder U1 also affects the mediator M, creating
# mediator-outcome (M-O) confounding.  When mo_confounding = 0 the
# original DGP is preserved exactly (backward compatible).
#
# v0.4.0: the phi parameter (default 0) adds a mediator-specific
# genetic instrument Gm (e.g. an eQTL from fetal genotype -> placental
# isoform expression).  When phi > 0, Gm affects M but is independent
# of U and has no direct path to Y, making the 2-stage MR mediation
# estimator (fit_iv2sls_mediation2) point-identified under M-O
# confounding.  When phi = 0 (default), no Gm is generated and the
# original DGP is preserved exactly.
#
# v0.5.0: imperfect instrument independence.  The v0.4.0 DGP generates
# G and Gm as pure noise, completely independent of U -- the best
# possible case for an instrument.  v0.5.0 adds parameters that let the
# instruments be correlated with the confounders (rho_G1, rho_G2) via a
# shared population-structure factor (rho_pop), draws separate
# confounders for the Z->M and M->Y paths (separate_U), and generates
# path-specific negative controls W1 and W2 with independent coverage
# (omega_1, omega_2).  This sets up the "tipping point" simulation: as
# instrument exogeneity is violated, IV2SLS2 degrades while the
# negative-control-augmented PGC-2 estimator (fit_pgc_mediation2) may
# become preferable.  All new parameters default to values that
# reproduce the v0.4.0 DGP exactly.
# ============================================================

#' Generate one synthetic dataset (internal)
#'
#' When \code{mo_confounding > 0}, the unmeasured confounder U1 also
#' affects the mediator M, creating mediator-outcome confounding — the
#' key extension for the mediation simulation.  When \code{mo_confounding = 0}
#' (default), the original DGP is preserved exactly.
#'
#' When \code{phi > 0} (v0.4.0), a mediator-specific genetic instrument
#' \code{Gm} is generated: \eqn{Gm \sim \mathcal{N}(0,1)}, independent of
#' \code{G}, \code{U}, and all other variables, and the mediator equation
#' becomes \eqn{M = \alpha_M Z + \delta_{mo} 0.5 U_1 + \phi Gm + \varepsilon_M}.
#' \code{Gm} is a valid instrument for \code{M}: it moves \code{M} but is
#' independent of the confounders and has no direct path to \code{Y}.  This
#' enables the 2-stage MR mediation estimator (\code{\link{fit_iv2sls_mediation2}})
#' to point-identify NDE and NIE even under M-O confounding.  When
#' \code{phi = 0} (default), no \code{Gm} is generated and the original DGP
#' is preserved exactly.
#'
#' @section v0.5.0 — imperfect instrument independence:
#' The v0.4.0 DGP generates instruments as pure noise, independent of the
#' confounders — the best possible case.  v0.5.0 adds parameters that
#' violate this independence, modelling realistic genomic structure:
#' \describe{
#'   \item{\code{rho_G1}}{Correlation of the exposure instrument G1 with the
#'     Z->M confounder U_XM (instrument exogeneity violation).  Default 0.}
#'   \item{\code{rho_G2}}{Correlation of the mediator instrument G2 with the
#'     M->Y confounder U_MY (instrument exogeneity violation).  Default 0.}
#'   \item{\code{rho_pop}}{Shared population-structure factor P that induces
#'     correlation between G1 and G2 (linkage / stratification).  Default 0.}
#'   \item{\code{separate_U}}{If \code{TRUE}, draw U_XM and U_MY as
#'     independent confounders for the Z->M and M->Y paths respectively.
#'     If \code{FALSE} (default), U_XM = U_MY = U1 (backward compatible).}
#'   \item{\code{omega_1}}{Coverage of U_XM by the path-specific negative
#'     controls W1.  \code{NULL} (default) uses \code{w_signal}.}
#'   \item{\code{omega_2}}{Coverage of U_MY by the path-specific negative
#'     controls W2.  \code{NULL} (default) uses \code{w_signal}.}
#' }
#' When all v0.5.0 parameters are at their defaults, the DGP is identical
#' to v0.4.0.  When any is non-default, the output additionally includes
#' \code{G1}, \code{G2}, \code{W1}, \code{W2}, \code{U_XM}, \code{U_MY},
#' and (when \code{rho_pop > 0}) \code{P}.
#'
#' @param n              Sample size. Default 500.
#' @param n_features     Number of outcome and negative-control features. Default 20.
#' @param n_mediators    Number of independent mediators (v0.8.4). When > 1,
#'                       each mediator M_m has its own genetic instrument Gm_m
#'                       and contributes additively to Y. M and Gm are returned
#'                       as \code{n_mediators x n} matrices. Default 1 (single
#'                       mediator, backward compatible).
#' @param beta_Z         Direct effect of Z on Y (true NDE). Default 0.10.
#' @param alpha_M        Effect of Z on mediator M. Default 0.50.
#' @param beta_M         Effect of M on Y (per-mediator NIE = alpha_M * beta_M;
#'                       total NIE = n_mediators * alpha_M * beta_M). Default 0.30.
#' @param conf_str       Confounding strength delta. Default 0.80.
#' @param w_signal       Proxy quality omega (0 = noise, 1 = perfect U proxy). Default 0.70.
#' @param mo_confounding Strength of U1 -> M (mediator-outcome confounding).
#'                       0 = no M-O confounding (original DGP). Default 0.
#' @param pleio          Strength of a direct G -> Y path (horizontal pleiotropy),
#'                       violating the exclusion restriction. 0 = no pleiotropy
#'                       (valid instrument, original DGP). Default 0.
#' @param phi            Strength of the mediator instrument Gm -> M (v0.4.0).
#'                       0 = no mediator instrument (original DGP, no Gm generated).
#'                       Default 0. When > 0, a valid instrument for M is generated,
#'                       enabling point identification of NDE/NIE via
#'                       \code{\link{fit_iv2sls_mediation2}}.
#' @param gamma_G        Strength of the exposure instrument G -> Z. Default 0.6.
#' @param rho_G1         Correlation of G1 with U_XM (v0.5.0). Default 0.
#' @param rho_G2         Correlation of G2 with U_MY (v0.5.0). Default 0.
#' @param rho_pop        Shared population structure inducing G1-G2 correlation
#'                       (v0.5.0). Default 0.
#' @param separate_U     If TRUE, draw separate confounders U_XM and U_MY for
#'                       the Z->M and M->Y paths (v0.5.0). Default FALSE.
#' @param omega_1        Coverage of U_XM by W1 (v0.5.0). NULL = use w_signal.
#' @param omega_2        Coverage of U_MY by W2 (v0.5.0). NULL = use w_signal.
#' @param feat_cor       Within-module feature correlation (v0.8.1).  When > 0,
#'                       the idiosyncratic noise in the Y and W panels is drawn
#'                       from a multivariate normal with a block-diagonal
#'                       correlation matrix modelling co-expression modules.
#'                       0 = independent noise (backward compatible). Default 0.
#' @param u_strength     Numeric scalar (v0.9.2, JYH #864): scales the single
#'                       confounder's effect on Z, M, and Y. Default NULL → 1
#'                       (backward compatible). See \code{\link{run_single_iteration}}
#'                       for the length-k vector version (multiple confounders).
#' @param w_coverage_profile A list with optional \code{w1} and \code{w2}
#'                       numeric vectors (v0.9.2, JYH #864): per-control coverage
#'                       of U_XM / U_MY. Default NULL → scalar omega applied
#'                       uniformly (backward compatible).
#' @param seed           Optional integer RNG seed for reproducibility.
#'
#' @return A named list with elements Z, G (n x n_features matrix), Y, W, U1, M,
#'   synthetic_data, true_total, true_NDE, true_NIE.  When \code{phi > 0}, also
#'   includes \code{Gm} (numeric vector, length n, or \code{n_mediators x n}
#'   matrix when \code{n_mediators > 1}).  When any v0.5.0 parameter
#'   is non-default, also includes \code{G1}, \code{G2}, \code{W1}, \code{W2},
#'   \code{U_XM}, \code{U_MY}, and (when \code{rho_pop > 0}) \code{P}.
#'   When \code{n_mediators > 1}, \code{M} is an \code{n_mediators x n} matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic total-effect data
#' dat <- generate_toy_data(n = 200, n_features = 10, seed = 42)
#'
#' # Mediation with mediator instrument and path-specific NCs
#' dat <- generate_toy_data(n = 200, n_features = 10, phi = 0.8,
#'                          mo_confounding = 0.8, separate_U = TRUE,
#'                          omega_1 = 0.7, omega_2 = 0.7, seed = 42)
#' idata <- iconic_data(Z = dat$Z, Y = dat$Y, M = dat$M,
#'                      G = dat$G[, 1], Gm = dat$Gm,
#'                      W1 = dat$W1, W2 = dat$W2)
#' }
generate_toy_data <- function(n              = 500,
                              n_features     = 20,
                              n_mediators    = 1,
                              beta_Z         = 0.10,
                              alpha_M        = 0.50,
                              beta_M         = 0.30,
                              conf_str       = 0.80,
                              w_signal       = 0.70,
                              mo_confounding = 0,
                              pleio          = 0,
                              phi            = 0,
                              gamma_G        = 0.6,
                              rho_G1         = 0,
                              rho_G2         = 0,
                              rho_pop        = 0,
                              separate_U     = FALSE,
                              omega_1        = NULL,
                              omega_2        = NULL,
                              feat_cor       = 0,
                              # v0.9.2 (JYH #864): per-confounder strength and
                              # per-control coverage for heterogeneous U/W.
                              u_strength     = NULL,
                              w_coverage_profile = NULL,
                              seed           = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # ── Determine whether the v0.5.0 DGP is active ──
  v05_active <- rho_G1 != 0 || rho_G2 != 0 || rho_pop != 0 ||
                separate_U || !is.null(omega_1) || !is.null(omega_2)

  # v0.9.2 (JYH #864): u_strength scales the single confounder's effect.
  # generate_toy_data has k=1, so u_strength is a scalar (default 1).
  # (run_single_iteration supports k>1 with a length-k u_strength vector.)
  u_str <- if (is.null(u_strength)) 1 else as.numeric(rep_len(u_strength, 1))

  # v0.9.2 (JYH #864): w_coverage_profile — per-control coverage vectors.
  wcp <- w_coverage_profile
  if (!is.null(wcp)) {
    if (!is.list(wcp)) stop("w_coverage_profile must be a list with 'w1'/'w2'")
    if (!is.null(wcp$w1)) wcp$w1 <- rep_len(as.numeric(wcp$w1), n_features)
    if (!is.null(wcp$w2)) wcp$w2 <- rep_len(as.numeric(wcp$w2), n_features)
  }

  # ── Confounders ──
  # v0.4.0: single U1 confounds both paths.
  # v0.5.0 (separate_U = TRUE): U_XM confounds Z->M, U_MY confounds M->Y.
  U1 <- rnorm(n)
  U2 <- rnorm(n)

  if (v05_active) {
    if (separate_U) {
      U_XM <- rnorm(n)
      U_MY <- rnorm(n)
    } else {
      U_XM <- U1
      U_MY <- U1
    }
  } else {
    U_XM <- U1
    U_MY <- U1
  }

  # ── Instruments ──
  # v0.4.0: G ~ N(0,1), Gm ~ N(0,1), both pure noise.
  # v0.5.0: G1 = rho_G1*U_XM + rho_pop*P + sqrt(1-rho_G1^2-rho_pop^2)*N(0,1)
  #         G2 = rho_G2*U_MY + rho_pop*P + sqrt(1-rho_G2^2-rho_pop^2)*N(0,1)
  # The sqrt term keeps Var(G) = 1 when rho terms are within valid range.
  if (v05_active) {
    P <- rnorm(n)   # shared population structure

    # Variance budget for G1: rho_G1^2 + rho_pop^2 must be <= 1
    resid_var_G1 <- max(1 - rho_G1^2 - rho_pop^2, 0)
    G1 <- rho_G1 * U_XM + rho_pop * P + sqrt(resid_var_G1) * rnorm(n)
    G1 <- as.numeric(scale(G1))

    # G2 (mediator instrument): only generated when phi > 0
    # v0.8.4: when n_mediators > 1, generate one instrument per mediator
    #         as an n_mediators x n matrix (each an independent draw with
    #         the same rho_G2/rho_pop structure).
    if (phi > 0) {
      resid_var_G2 <- max(1 - rho_G2^2 - rho_pop^2, 0)
      if (n_mediators > 1) {
        G2 <- matrix(NA_real_, n_mediators, n)
        for (mm in seq_len(n_mediators)) {
          G2[mm, ] <- rho_G2 * U_MY + rho_pop * P +
                      sqrt(resid_var_G2) * rnorm(n)
          G2[mm, ] <- as.numeric(scale(G2[mm, ]))
        }
      } else {
        G2 <- rho_G2 * U_MY + rho_pop * P + sqrt(resid_var_G2) * rnorm(n)
        G2 <- as.numeric(scale(G2))
      }
    } else {
      G2 <- NULL
    }

    G  <- G1          # G is the exposure instrument (alias for backward compat)
    Gm <- G2           # Gm is the mediator instrument (alias for backward compat)
  } else {
    # v0.4.0 path: pure-noise instruments
    G  <- rnorm(n)
    if (phi > 0) {
      if (n_mediators > 1) {
        Gm <- matrix(rnorm(n_mediators * n), n_mediators, n)
      } else {
        Gm <- rnorm(n)
      }
    } else {
      Gm <- NULL
    }
    P  <- NULL
    G1 <- G
    G2 <- Gm
  }

  # ── Exposure ──
  # Z = scale(0.6*G + conf_str*U_XM + eps)
  # U_XM = U1 in backward-compatible mode.
  # v0.9.2 (JYH #864): u_str scales the confounder's effect.
  Z_raw <- gamma_G * G + u_str * conf_str * U_XM + rnorm(n, 0, 0.5)
  Z <- as.numeric(scale(Z_raw))

  # ── Mediator ──
  # M = alpha_M*Z + mo_confounding*0.5*U_MY + phi*G2 + eps
  # U_MY = U1 in backward-compatible mode.
  # v0.9.2 (JYH #864): u_str scales the M<-U_MY confounding.
  # v0.8.4: when n_mediators > 1, generate M as an n_mediators x n matrix,
  #         each mediator with its own Gm instrument and independent noise.
  if (n_mediators > 1) {
    M <- matrix(NA_real_, n_mediators, n)
    for (mm in seq_len(n_mediators)) {
      M[mm, ] <- alpha_M * Z
      if (mo_confounding > 0) M[mm, ] <- M[mm, ] + u_str * mo_confounding * 0.5 * U_MY
      if (phi > 0) {
        gm_mm <- if (is.matrix(Gm)) Gm[mm, ] else Gm
        M[mm, ] <- M[mm, ] + phi * gm_mm
      }
      M[mm, ] <- M[mm, ] + rnorm(n, 0, 0.05)
    }
  } else {
    M <- alpha_M * Z
    if (mo_confounding > 0) M <- M + u_str * mo_confounding * 0.5 * U_MY
    if (phi > 0)            M <- M + phi * Gm
    M <- M + rnorm(n, 0, 0.05)
  }

  # ── Negative controls ──
  # v0.4.0: single W panel, W_f = w_signal*U1 + (1-w_signal)*U2 + eps.
  # v0.5.0: path-specific W1 (captures U_XM) and W2 (captures U_MY),
  #         with independent coverage omega_1, omega_2.
  # v0.8.1: when feat_cor > 0, the idiosyncratic noise is drawn from a
  #         multivariate normal with a block-diagonal correlation matrix
  #         (co-expression modules), so W features are correlated conditional
  #         on U.
  om1 <- if (is.null(omega_1)) w_signal else omega_1
  om2 <- if (is.null(omega_2)) w_signal else omega_2
  # v0.9.2 (JYH #864): per-control coverage vectors (length n_features).
  # NULL w_coverage_profile → scalar applied uniformly (legacy behavior).
  om1_vec <- if (!is.null(wcp$w1)) wcp$w1 else rep(om1, n_features)
  om2_vec <- if (!is.null(wcp$w2)) wcp$w2 else rep(om2, n_features)
  ws_vec  <- if (!is.null(wcp$w1)) wcp$w1 else rep(w_signal, n_features)

  # Build block-diagonal correlation matrix for noise (v0.8.1)
  cor_mat <- if (feat_cor > 0) .block_diag_cor(n_features, feat_cor) else NULL

  if (v05_active) {
    W1 <- matrix(NA_real_, n, n_features)
    W2 <- matrix(NA_real_, n, n_features)
    if (!is.null(cor_mat)) {
      # Correlated noise: draw n x p from MVN(0, Sigma * 0.3^2)
      noise1 <- MASS::mvrnorm(n, rep(0, n_features), cor_mat * 0.3^2)
      noise2 <- MASS::mvrnorm(n, rep(0, n_features), cor_mat * 0.3^2)
      for (f in seq_len(n_features)) {
        W1[, f] <- om1_vec[f] * U_XM + (1 - om1_vec[f]) * U2 + noise1[, f]
        W2[, f] <- om2_vec[f] * U_MY + (1 - om2_vec[f]) * rnorm(n) + noise2[, f]
      }
    } else {
      for (f in seq_len(n_features)) {
        W1[, f] <- om1_vec[f] * U_XM + (1 - om1_vec[f]) * U2 + rnorm(n, 0, 0.3)
        W2[, f] <- om2_vec[f] * U_MY + (1 - om2_vec[f]) * rnorm(n) + rnorm(n, 0, 0.3)
      }
    }
    W1 <- scale(W1)
    W2 <- scale(W2)
    # Combined W for backward-compatible estimators: average of W1 and W2
    # (when separate_U = FALSE, W1 = W2 = W, so this is exact).
    W <- if (separate_U) scale((W1 + W2) / 2) else W1
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
  # Y_f = beta_M*M + beta_Z*Z + gamma_f*U_MY + pleio*G1 + eps
  # U_MY = U1 in backward-compatible mode.
  # v0.8.1: when feat_cor > 0, the idiosyncratic noise is drawn from a
  #         multivariate normal with a block-diagonal correlation matrix,
  #         so Y features are correlated conditional on U.
  Y <- matrix(NA_real_, n, n_features)
  # v0.8.4: when n_mediators > 1, M is a matrix (n_mediators x n);
  # sum across mediators so each contributes to Y.
  M_eff <- if (is.matrix(M)) colSums(M) else M
  for (f in seq_len(n_features)) {
    gamma_f <- runif(1, 0.4, 0.8) * conf_str
    Y[, f]  <- beta_M * M_eff + beta_Z * Z + u_str * gamma_f * U_MY + pleio * G
  }
  if (!is.null(cor_mat)) {
    Y_noise <- MASS::mvrnorm(n, rep(0, n_features), cor_mat * 0.2^2)
    Y <- Y + Y_noise
  } else {
    for (f in seq_len(n_features))
      Y[, f] <- Y[, f] + rnorm(n, 0, 0.2)
  }

  out <- list(
    Z = Z,
    G = matrix(rep(G, n_features), n, n_features),
    Y = Y,
    W = W,
    U1 = U1,
    M = M,
    synthetic_data = data.frame(
      fetal_sex = rbinom(n, 1, 0.5),
      gestational_age = rnorm(n)
    ),
    true_total = beta_Z + n_mediators * alpha_M * beta_M,
    true_NDE   = beta_Z,
    true_NIE   = n_mediators * alpha_M * beta_M
  )
  if (phi > 0) out$Gm <- Gm
  out$n_mediators <- n_mediators

  # v0.9.2: record the DGP parameters used, so scenario_manifest() can
  # recover the modifiable-parameter scalars and feat_cor from a dat
  # object without requiring the caller to pass them again (JYH #543).
  out$dgp_params <- list(
    n = n, n_features = n_features, n_mediators = n_mediators,
    beta_Z = beta_Z, alpha_M = alpha_M, beta_M = beta_M,
    conf_str = conf_str, w_signal = w_signal,
    mo_confounding = mo_confounding, pleio = pleio, phi = phi,
    gamma_G = gamma_G, rho_G1 = rho_G1, rho_G2 = rho_G2,
    rho_pop = rho_pop, separate_U = separate_U,
    omega_1 = omega_1, omega_2 = omega_2, feat_cor = feat_cor,
    # v0.9.2 (JYH #864)
    u_strength = u_str,
    w_coverage_profile = wcp
  )

  # v0.5.0 extras (only when the v0.5.0 DGP is active)
  if (v05_active) {
    out$G1   <- G1
    out$G2   <- G2
    out$W1   <- W1
    out$W2   <- W2
    out$U_XM <- U_XM
    out$U_MY <- U_MY
    if (rho_pop > 0) out$P <- P
  }

  out
}


# ── Scenario manifest (v0.9.2, JYH #543 / #582) ─────────────────────────
# Surfaces the ground-truth estimands and the modifiable/fixed parameter
# split so the manuscript can render an orientation table preceding the
# simulation results. Accepts either a `dat` list returned by
# `generate_toy_data()` (estimands + fixed params read from the object)
# or a bare parameter list (estimands recomputed from beta_Z/alpha_M/beta_M).

#' Scenario manifest: truth and parameter ranges for a simulation
#'
#' Returns a reader-orienting summary of a simulation scenario: the
#' ground-truth estimands (NDE, NIE, total effect), the modifiable
#' parameters with their swept ranges, and the fixed parameters with
#' their values. Intended to be rendered as a table preceding simulation
#' results in the manuscript (JYH #543, #582: "state the truth and the
#' parameter ranges up front to orient the reader").
#'
#' `dat_or_params` may be either:
#' \itemize{
#'   \item a list returned by [generate_toy_data()], in which case the
#'         estimands and fixed parameters are read from the object
#'         (`true_NDE`, `true_NIE`, `true_total`, `n`, `n_features`,
#'         `n_mediators`, `separate_U`, `feat_cor`); or
#'   \item a bare named list of DGP parameters
#'         (`beta_Z`, `alpha_M`, `beta_M`, `n_mediators`, `n`,
#'         `n_features`, `separate_U`, `feat_cor`, ...), in which case
#'         the estimands are recomputed as
#'         `NDE = beta_Z`, `NIE = n_mediators * alpha_M * beta_M`,
#'         `total = NDE + NIE`.
#' }
#'
#' The `*_grid` arguments record the swept ranges for the modifiable
#' parameters. When a grid argument is omitted or NULL, the
#' corresponding modifiable parameter is reported with its scalar value
#' only (no range).
#'
#' @param dat_or_params A `dat` list from [generate_toy_data()] or a
#'   named list of DGP parameters.
#' @param conf_grid    Swept values of `conf_str` (confounding strength).
#'   Default NULL.
#' @param coverage_grid Swept values of NC coverage (`w_signal` /
#'   `omega`). Default NULL.
#' @param mo_confounding_grid Swept values of mediator-outcome
#'   confounding strength. Default NULL.
#' @param phi_grid     Swept values of the mediator-instrument strength
#'   `phi`. Default NULL.
#' @param rho_G1_grid  Swept values of the G1-U_XM correlation. Default NULL.
#' @param rho_G2_grid  Swept values of the G2-U_MY correlation. Default NULL.
#' @param rho_pop_grid Swept values of the population-stratification
#'   correlation. Default NULL.
#' @return A named list with elements:
#'   \describe{
#'     \item{`estimands`}{Named numeric vector: `NDE`, `NIE`, `total`.}
#'     \item{`modifiable_parameters`}{Data frame with columns
#'           `parameter`, `value`, `swept_range`. `value` is the scalar
#'           used (or NA when only a grid was supplied);
#'           `swept_range` is a comma-separated string of grid values,
#'           or NA when the parameter was held fixed.}
#'     \item{`fixed_parameters`}{Data frame with columns `parameter`,
#'           `value`.}
#'   }
#' @export
#' @examples
#' dat <- generate_toy_data(n = 200, phi = 0.8, separate_U = TRUE,
#'                          omega_1 = 0.7, omega_2 = 0.7, seed = 42)
#' scenario_manifest(dat, conf_grid = c(0.3, 0.8),
#'                   coverage_grid = c(0.3, 0.7, 1))
scenario_manifest <- function(dat_or_params,
                              conf_grid = NULL,
                              coverage_grid = NULL,
                              mo_confounding_grid = NULL,
                              phi_grid = NULL,
                              rho_G1_grid = NULL,
                              rho_G2_grid = NULL,
                              rho_pop_grid = NULL) {
  p <- dat_or_params

  # Detect whether `p` is a generate_toy_data() output (has true_NDE) or a
  # bare parameter list.
  is_dat <- !is.null(p$true_NDE) && !is.null(p$true_NIE)

  if (is_dat) {
    nde  <- p$true_NDE
    nie  <- p$true_NIE
    tot  <- if (!is.null(p$true_total)) p$true_total else nde + nie
    n    <- attr(p$Y, "n"); if (is.null(n)) n <- nrow(p$Y)
    nf   <- ncol(p$Y)
    nm   <- if (!is.null(p$n_mediators)) p$n_mediators else 1L
    sepU <- !is.null(p$U_XM)
    # v0.9.2: recover scalar modifiable params + feat_cor from dgp_params
    # when generate_toy_data() recorded them.
    dp <- p$dgp_params
    if (!is.null(dp)) {
      fcor <- dp$feat_cor
      conf_val <- dp$conf_str
      mo_val   <- dp$mo_confounding
      phi_val  <- dp$phi
      cov_val  <- if (!is.null(dp$omega_1)) dp$omega_1 else dp$w_signal
      rho1 <- dp$rho_G1; rho2 <- dp$rho_G2; rhop <- dp$rho_pop
    } else {
      fcor <- attr(p$Y, "feat_cor"); if (is.null(fcor)) fcor <- 0
      conf_val <- NA_real_; mo_val <- NA_real_; phi_val <- NA_real_
      cov_val  <- NA_real_; rho1 <- NA_real_; rho2 <- NA_real_; rhop <- NA_real_
    }
  } else {
    beta_Z  <- if (!is.null(p$beta_Z))  p$beta_Z  else 0.10
    alpha_M <- if (!is.null(p$alpha_M)) p$alpha_M else 0.50
    beta_M  <- if (!is.null(p$beta_M))  p$beta_M  else 0.30
    nm      <- if (!is.null(p$n_mediators)) p$n_mediators else 1L
    nde <- beta_Z
    nie <- nm * alpha_M * beta_M
    tot <- nde + nie
    n    <- if (!is.null(p$n))  p$n  else 500
    nf   <- if (!is.null(p$n_features)) p$n_features else 20
    sepU <- if (!is.null(p$separate_U)) p$separate_U else FALSE
    fcor <- if (!is.null(p$feat_cor))   p$feat_cor   else 0
    conf_val <- if (!is.null(p$conf_str))       p$conf_str       else NA_real_
    mo_val   <- if (!is.null(p$mo_confounding)) p$mo_confounding else NA_real_
    phi_val  <- if (!is.null(p$phi))            p$phi            else NA_real_
    cov_val  <- if (!is.null(p$w_signal))       p$w_signal       else
                if (!is.null(p$omega_1))        p$omega_1        else NA_real_
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
                  "rho_G1", "rho_G2", "rho_pop"),
    value = c(conf_val, cov_val, mo_val, phi_val, rho1, rho2, rhop),
    swept_range = c(fmt_range(conf_grid), fmt_range(coverage_grid),
                    fmt_range(mo_confounding_grid), fmt_range(phi_grid),
                    fmt_range(rho_G1_grid), fmt_range(rho_G2_grid),
                    fmt_range(rho_pop_grid)),
    stringsAsFactors = FALSE
  )

  fixed <- data.frame(
    parameter = c("n", "n_features", "n_mediators",
                  "separate_U", "feat_cor"),
    value = c(n, nf, nm, sepU, fcor),
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
#' @param p   Number of features.
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
