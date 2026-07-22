# ============================================================
# Generalised structural generator: one synthetic dataset.
#
# This is the multi-confounder, pluggable-negative-control analogue of
# generate_toy_data(). It keeps the causal skeleton (so the total effect is
# known exactly) while (a) drawing realistic covariate / outcome texture from
# a trained generator, (b) supporting k latent confounders, and (c) sourcing
# the negative controls from a swappable nc_model. The return shape matches
# generate_toy_data() so run_methods()/analyze_methods_robust() consume it
# unchanged.
#
# The mo_confounding parameter (default 0) optionally adds U1 -> M
# (mediator-outcome confounding) for the mediation extension.  When 0,
# the original DGP is preserved exactly (backward compatible).
#
# v0.4.0: the phi parameter (default 0) adds a mediator-specific genetic
# instrument Gm (e.g. an eQTL from fetal genotype -> placental isoform
# expression).  When phi > 0, Gm affects M but is independent of U and has
# no direct path to Y, making the 2-stage MR mediation estimator
# (fit_iv2sls_mediation2) point-identified under M-O confounding.  When
# phi = 0 (default), no Gm is generated and the original DGP is preserved.
#
# v0.5.0: imperfect instrument independence.  Adds rho_G1, rho_G2, rho_pop,
# separate_U, omega_1, omega_2 — see generate_toy_data() for the full
# description.  When separate_U = TRUE, the confounder matrix U is split:
# the first confounder drives U_XM (Z->M path), a second independent
# confounder drives U_MY (M->Y path).  Path-specific negative controls W1
# and W2 are generated from the pluggable nc_model with independent coverage.
# All new parameters default to v0.4.0-reproducing values.
#
# Identification note: the generator draws texture INDEPENDENTLY of the latent
# confounders U and never feeds a shared texture factor into both Z and Y, so
# the only Z<->Y backdoor is through U. That keeps true_total a closed-form
# function of the causal knobs.
# ============================================================

#' Generate one synthetic dataset under the generalised SCM
#'
#' @param trained_gan  Optional `iconic_gan` from [train_gan_on_real_data()],
#'   supplying realistic covariate/outcome texture. If `NULL`, default synthetic
#'   covariates are used.
#' @param n_synthetic_samples Sample size. Default 500.
#' @param n_features   Number of outcome (and control) features. Default 20.
#' @param n_confounders Number of latent confounders `k`. Default 1
#'   (backward-compatible single-confounder model).
#' @param n_mediators  Number of independent mediators (v0.8.4). When > 1,
#'   each mediator has its own genetic instrument Gm and contributes
#'   additively to Y. Default 1 (single mediator, backward compatible).
#' @param beta_Z       Direct effect of Z on Y (true NDE). Default 0.10.
#' @param alpha_M      Effect of Z on the mediator M. Default 0.50.
#' @param beta_M       Effect of M on Y. Default 0.30.
#' @param effect_size  Optional shortcut: if non-`NULL`, sets a pure direct total
#'   effect (`beta_Z = effect_size`, `alpha_M = beta_M = 0`). Use `0` for null
#'   (Type I error) simulations. Default `NULL` (use the mediation parameters).
#' @param conf_strength Overall confounding strength (analogue of `conf_str`);
#'   scales the confounder loadings into Z and Y. Default 0.80.
#' @param coverage     How well the negative controls span the confounder
#'   subspace, in `[0, 1]` (passed to `nc_model`). Default 1.
#' @param captured     Integer indices of the confounders the negative controls
#'   see (passed to `nc_model`). Default all `k`.
#' @param nc_model     Negative-control model: a function `(U, covariates,
#'   params) -> W`, or a registered name (`"proxy"`, `"cpg"`). Default `"proxy".
#' @param nc_params    Extra named parameters forwarded to `nc_model`.
#' @param mo_confounding Strength of U1 -> M (mediator-outcome confounding).
#'   0 = no M-O confounding (original DGP). Default 0. When > 0, the first
#'   confounder U\[,1\] also affects M, creating the M-O confounding that the
#'   mediation estimators are benchmarked against.
#' @param pleio      Strength of a direct G -> Y path (horizontal pleiotropy),
#'   violating the exclusion restriction. 0 = no pleiotropy (valid instrument,
#'   original DGP). Default 0. When > 0, G affects Y directly in addition to
#'   through Z, allowing benchmarking of IV/2SLS under instrument invalidity.
#' @param phi      Strength of the mediator instrument Gm -> M (v0.4.0).
#'   0 = no mediator instrument (original DGP, no Gm generated). Default 0.
#'   When > 0, a valid instrument for M is generated (independent of U and G,
#'   no direct path to Y), enabling point identification of NDE/NIE via
#'   [fit_iv2sls_mediation2()] even under M-O confounding.
#' @param gamma_G  Strength of the exposure instrument G -> Z. Default 0.6.
#' @param rho_G1   Correlation of G1 with U_XM (v0.5.0, instrument exogeneity
#'   violation). Default 0.
#' @param rho_G2   Correlation of G2 with U_MY (v0.5.0, instrument exogeneity
#'   violation). Default 0.
#' @param rho_pop  Shared population structure inducing G1-G2 correlation
#'   (v0.5.0). Default 0.
#' @param separate_U If TRUE, draw separate confounders for the Z->M and M->Y
#'   paths (v0.5.0). Default FALSE.
#' @param omega_1  Coverage of U_XM by W1 (v0.5.0). NULL = use `coverage`.
#' @param omega_2  Coverage of U_MY by W2 (v0.5.0). NULL = use `coverage`.
#' @param feat_cor Within-module feature correlation (v0.8.1). When > 0 and
#'   the GAN does not provide feature_correlations, a block-diagonal
#'   correlation matrix is used for the NC and outcome noise. GAN-learned
#'   correlations take precedence. Default 0.
#' @param u_strength Numeric vector (v0.9.2, JYH #864): per-confounder strength
#'   profile (length k). Default NULL → rep(1, k) (equal strength, backward
#'   compatible). Recycled to length k and normalized so the total confounding
#'   budget is unchanged.
#' @param w_coverage_profile A list with optional \code{w1} and \code{w2}
#'   numeric vectors (v0.9.2, JYH #864): per-control coverage of U_XM / U_MY
#'   (length n_features). Default NULL → scalar omega applied uniformly.
#' @param MMExp,MMOut,MMCon,MMCpG Per-pathway confounding multipliers (exposure,
#'   outcome, controls, methylation). Default 1.
#' @param seed         Optional RNG seed.
#'
#' @return A named list matching [generate_toy_data()] — `Z`, `G` (n x
#'   n_features), `Y`, `W`, `U1`, `M`, `synthetic_data`, `true_total`,
#'   `true_NDE`, `true_NIE` — plus `U` (full n x k confounder matrix),
#'   `genetic_instrument`, `successful_features`, `failed_features`, and
#'   `params`. When `phi > 0`, also includes `Gm` (numeric vector, length n,
#'   or `n_mediators x n` matrix when `n_mediators > 1`).
#'   When any v0.5.0 parameter is non-default, also includes `G1`, `G2`,
#'   `W1`, `W2`, `U_XM`, `U_MY`, and (when `rho_pop > 0`) `P`.
#'   When `n_mediators > 1`, `M` is an `n_mediators x n` matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' gan <- train_gan_on_real_data(load_real_input_data(example = TRUE)$gan_training_data)
#' dat <- run_single_iteration(gan, n_features = 10, n_confounders = 2,
#'                             nc_model = "cpg", coverage = 0.5)
#' analyze_methods_robust(dat)
#' }
run_single_iteration <- function(trained_gan          = NULL,
                                 n_synthetic_samples   = 500,
                                 n_features            = 20,
                                 n_confounders         = 1,
                                 n_mediators           = 1,
                                 beta_Z                = 0.10,
                                 alpha_M               = 0.50,
                                 beta_M                = 0.30,
                                 effect_size           = NULL,
                                 conf_strength         = 0.80,
                                 coverage              = 1,
                                 captured              = NULL,
                                 nc_model              = "proxy",
                                 nc_params             = list(),
                                 mo_confounding        = 0,
                                 pleio                 = 0,
                                 phi                   = 0,
                                 gamma_G               = 0.6,
                                 rho_G1                = 0,
                                 rho_G2                = 0,
                                 rho_pop               = 0,
                                 separate_U            = FALSE,
                                 omega_1               = NULL,
                                 omega_2               = NULL,
                                 feat_cor              = 0,
                                 # v0.9.2 (JYH #864): per-confounder strength and
                                 # per-control coverage profiles for heterogeneous
                                 # U/W characterization in the A3 simulation.
                                 u_strength            = NULL,
                                 w_coverage_profile    = NULL,
                                 MMExp = 1, MMOut = 1, MMCon = 1, MMCpG = 1,
                                 seed                  = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(effect_size)) { beta_Z <- effect_size; alpha_M <- 0; beta_M <- 0 }
  true_total <- beta_Z + n_mediators * alpha_M * beta_M
  true_NDE   <- beta_Z
  true_NIE   <- n_mediators * alpha_M * beta_M

  n  <- n_synthetic_samples
  p  <- n_features
  k  <- max(1L, as.integer(n_confounders))
  nc <- .resolve_nc_model(nc_model)

  # v0.9.2 (JYH #864): per-confounder strength profile. Default rep(1, k)
  # preserves the legacy "equal-strength confounders" behavior. When
  # supplied, u_strength[j] scales confounder j's contribution to both Z
  # and Y, so some confounders can be stronger than others. Recycled to
  # length k and re-normalized so the overall confounding budget is
  # unchanged (sum(u_strength) == k); this keeps conf_strength comparable
  # across specifications.
  # NOTE: u_strength is resolved AFTER the k-bump below (v05_active &&
  # separate_U && k < 2 → k <- 2L), so its length matches the actual k.

  # v0.9.2 (JYH #864): per-control coverage profile. A list with optional
  # `w1` and `w2` numeric vectors (length = n_features) giving per-control
  # coverage of U_XM / U_MY. NULL → scalar omega applied uniformly (legacy).
  wcp <- w_coverage_profile
  if (!is.null(wcp)) {
    if (!is.list(wcp)) stop("w_coverage_profile must be a list with 'w1'/'w2'")
    if (!is.null(wcp$w1)) wcp$w1 <- rep_len(as.numeric(wcp$w1), p)
    if (!is.null(wcp$w2)) wcp$w2 <- rep_len(as.numeric(wcp$w2), p)
  }

  # ── v0.5.0: is the imperfect-independence DGP active? ──
  v05_active <- rho_G1 != 0 || rho_G2 != 0 || rho_pop != 0 ||
                separate_U || !is.null(omega_1) || !is.null(omega_2)

  # Latent confounders (columns are independent standardised confounders).
  # v0.5.0 (separate_U = TRUE): ensure at least 2 confounders so U_XM and
  # U_MY can be drawn as independent columns.
  if (v05_active && separate_U && k < 2) k <- 2L

  # v0.9.2 (JYH #864): resolve u_strength AFTER the k-bump so its length
  # matches the actual number of confounders.
  if (is.null(u_strength)) {
    u_strength <- rep(1, k)
  } else {
    u_strength <- rep_len(as.numeric(u_strength), k)
    u_strength <- u_strength * k / sum(u_strength)
  }

  U <- matrix(rnorm(n * k), n, k)

  # v0.5.0: path-specific confounders.
  # U_XM confounds Z -> M (first confounder).
  # U_MY confounds M -> Y (second confounder when separate_U, else first).
  U_XM <- U[, 1]
  U_MY <- if (separate_U) U[, 2] else U[, 1]

  # Realistic covariate / outcome texture, drawn independently of U.
  cov_df <- .iteration_covariates(trained_gan, n)
  exo_Y  <- .iteration_outcome_texture(trained_gan, n)   # exogenous, Y-only

  # Feature-level copula texture for the mediator panel (v0.9.0).
  # When the GAN carries a feature_texture (iconic_feature_texture), draw
  # realistic p-dimensional mediator feature vectors that preserve the
  # marginal distributions and cross-feature correlation structure of the
  # user's mediator panel.  When absent, fall back to scalar texture.
  exo_M_features <- .iteration_mediator_texture_features(trained_gan, n, p)

  # Feature-level residual correlation matrices from the trained texture
  # model.  When present, these inject realistic cross-feature correlations
  # into the negative-control noise and outcome noise, so that conditional
  # on U the simulated features are not mutually independent.
  # v0.8.1: when the GAN does not provide feature_correlations, fall back
  #         to a synthetic block-diagonal correlation matrix from feat_cor.
  gan_feat_cor <- trained_gan$feature_correlations %||% list()
  noise_cor_W <- gan_feat_cor$W   # p x p correlation matrix for NC noise
  noise_cor_Y <- gan_feat_cor$Y   # p x p correlation matrix for outcome noise

  # Synthetic block-diagonal fallback (v0.8.1)
  if (is.null(noise_cor_W) && feat_cor > 0)
    noise_cor_W <- .block_diag_cor(p, feat_cor)
  if (is.null(noise_cor_Y) && feat_cor > 0)
    noise_cor_Y <- .block_diag_cor(p, feat_cor)

  # ── Instruments ──
  # v0.4.0: G from simulate_single_genetic_instrument(), Gm ~ N(0,1).
  # v0.5.0: G1 = rho_G1*U_XM + rho_pop*P + sqrt(1-rho_G1^2-rho_pop^2)*N(0,1)
  #         G2 = rho_G2*U_MY + rho_pop*P + sqrt(1-rho_G2^2-rho_pop^2)*N(0,1)
  if (v05_active) {
    P <- rnorm(n)   # shared population structure
    resid_var_G1 <- max(1 - rho_G1^2 - rho_pop^2, 0)
    G1 <- as.numeric(scale(rho_G1 * U_XM + rho_pop * P +
                           sqrt(resid_var_G1) * rnorm(n)))
    if (phi > 0) {
      resid_var_G2 <- max(1 - rho_G2^2 - rho_pop^2, 0)
      # v0.8.4: per-mediator Gm when n_mediators > 1
      if (n_mediators > 1) {
        G2 <- matrix(NA_real_, n_mediators, n)
        for (mm in seq_len(n_mediators)) {
          G2[mm, ] <- as.numeric(scale(rho_G2 * U_MY + rho_pop * P +
                                       sqrt(resid_var_G2) * rnorm(n)))
        }
      } else {
        G2 <- as.numeric(scale(rho_G2 * U_MY + rho_pop * P +
                               sqrt(resid_var_G2) * rnorm(n)))
      }
    } else {
      G2 <- NULL
    }
    G  <- G1
    Gm <- G2
    gi <- list(G = G, dosages = NULL, maf = NULL)
  } else {
    gi <- simulate_single_genetic_instrument(n)
    G  <- gi$G
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

  # Structural exposure (scaled to unit variance BEFORE building M and Y,
  # mirroring generate_toy_data()).  U_XM = U[,1] in backward-compatible mode.
  # v0.9.2 (JYH #864): u_strength scales each confounder's contribution.
  aZ <- conf_strength * MMExp * u_strength / sqrt(k)       # confounder -> exposure
  Z  <- as.numeric(scale(gamma_G * G + as.numeric(U %*% aZ) + rnorm(n, 0, 0.5)))

  # Mediator: Z -> M, optionally U_MY -> M (mediator-outcome confounding),
  # optionally Gm -> M (mediator instrument, v0.4.0), plus exogenous mediator
  # texture from the copula model (v0.9.0: realistic feature-level marginals
  # and cross-feature correlations) or scalar fallback.
  # U_MY = U[,1] in backward-compatible mode.
  # v0.8.4: per-mediator M and Gm when n_mediators > 1
  # v0.9.0: when exo_M_features (copula) is available, each mediator gets a
  #         realistic feature-level texture draw; otherwise scalar fallback.
  if (n_mediators > 1) {
    M <- matrix(NA_real_, n_mediators, n)
    for (mm in seq_len(n_mediators)) {
      M[mm, ] <- alpha_M * Z
      if (mo_confounding > 0) M[mm, ] <- M[mm, ] + mo_confounding * 0.5 * U_MY
      if (phi > 0) {
        gm_mm <- if (is.matrix(Gm)) Gm[mm, ] else Gm
        M[mm, ] <- M[mm, ] + phi * gm_mm
      }
      if (!is.null(exo_M_features)) {
        tex_mm <- if (is.matrix(exo_M_features)) exo_M_features[mm, ] else exo_M_features
        M[mm, ] <- M[mm, ] + tex_mm + rnorm(n, 0, 0.02)
      } else {
        M[mm, ] <- M[mm, ] + 0.3 * rnorm(n) + rnorm(n, 0, 0.05)
      }
    }
  } else {
    M <- alpha_M * Z
    if (mo_confounding > 0) M <- M + mo_confounding * 0.5 * U_MY
    if (phi > 0)            M <- M + phi * Gm
    if (!is.null(exo_M_features)) {
      tex <- if (is.matrix(exo_M_features)) exo_M_features[1, ] else exo_M_features
      M <- M + tex + rnorm(n, 0, 0.02)
    } else {
      M <- M + 0.3 * rnorm(n) + rnorm(n, 0, 0.05)
    }
  }

  # Outcome panel: causal path + per-feature confounder loadings + shared
  # exogenous texture (adds realistic cross-feature correlation without a Z
  # backdoor) + idiosyncratic noise.
  # v0.5.0: outcome confounding uses U_MY (M->Y path confounder).
  # When noise_cor_Y is available, the idiosyncratic noise is drawn from a
  # multivariate normal with that correlation structure, so the outcome
  # features retain realistic cross-feature correlations conditional on U.
  Y <- matrix(NA_real_, n, p)
  gY_norm <- numeric(p)
  # v0.8.4: when n_mediators > 1, M is a matrix (n_mediators x n);
  # sum across mediators so each contributes to Y.
  M_eff <- if (is.matrix(M)) colSums(M) else M
  for (f in seq_len(p)) {
    # v0.9.2 (JYH #864): u_strength scales each confounder's outcome loading.
    gY      <- conf_strength * MMOut * u_strength * runif(k, 0.4, 0.8) / sqrt(k)
    gY_norm[f] <- sqrt(sum(gY^2))
    Y[, f]  <- beta_M * M_eff + beta_Z * Z + as.numeric(U %*% gY) +
               pleio * G + 0.4 * exo_Y
  }
  # Add idiosyncratic noise (optionally correlated across features)
  Y_noise <- .generate_outcome_noise(n, p, noise_cor_Y)
  Y <- Y + Y_noise

  # ── Negative controls ──
  # v0.4.0: single W panel from the pluggable nc_model.
  # v0.5.0: path-specific W1 (captures U_XM) and W2 (captures U_MY),
  #         with independent coverage omega_1, omega_2.
  # v0.9.2 (JYH #864): w_coverage_profile supplies per-control coverage
  #   vectors (length p) for W1 and/or W2, allowing heterogeneous proxy
  #   strength (e.g. one strong W_2 capturing most of U_MY, others weak).
  om1 <- if (is.null(omega_1)) coverage else omega_1
  om2 <- if (is.null(omega_2)) coverage else omega_2
  # Per-control coverage vectors (recycled to p); NULL → scalar (legacy).
  om1_vec <- if (!is.null(wcp$w1)) wcp$w1 else rep(om1, p)
  om2_vec <- if (!is.null(wcp$w2)) wcp$w2 else rep(om2, p)

  # Helper: build an n x p NC panel where column f has coverage om_vec[f].
  # Calls the nc_model once per distinct coverage value and assigns columns,
  # preserving the nc_model's texture within each coverage stratum.
  .nc_panel_per_coverage <- function(U_col, om_vec) {
    U_mat <- matrix(U_col, n, 1)
    out <- matrix(NA_real_, n, p)
    for (ov in unique(om_vec)) {
      cols <- which(om_vec == ov)
      nc_p <- utils::modifyList(
        list(n_features = length(cols), coverage = ov, captured = 1,
             MMCon = MMCon, MMCpG = MMCpG,
             noise_cor = .match_cor(noise_cor_W, length(cols))), nc_params)
      Wc <- nc(U_mat, cov_df, nc_p)
      if (!is.matrix(Wc)) Wc <- as.matrix(Wc)
      if (ncol(Wc) != length(cols))
        Wc <- Wc[, rep(seq_len(ncol(Wc)), length.out = length(cols)), drop = FALSE]
      out[, cols] <- Wc
    }
    out
  }

  if (v05_active) {
    W1 <- .nc_panel_per_coverage(U_XM, om1_vec)
    W2 <- .nc_panel_per_coverage(U_MY, om2_vec)
    # Combined W for backward-compatible estimators
    W <- if (separate_U) scale((W1 + W2) / 2) else W1
  } else {
    nc_p <- utils::modifyList(
      list(n_features = p, coverage = coverage, captured = captured,
           MMCon = MMCon, MMCpG = MMCpG,
           noise_cor = .match_cor(noise_cor_W, p)), nc_params)
    W <- nc(U, cov_df, nc_p)
    if (!is.matrix(W)) W <- as.matrix(W)
    if (ncol(W) != p) W <- W[, rep(seq_len(ncol(W)), length.out = p), drop = FALSE]
    W1 <- W
    W2 <- W
  }

  # Features with strong vs weak confounding (reference-compatible bookkeeping).
  thr <- stats::median(gY_norm)
  successful <- which(gY_norm >= thr)
  failed     <- which(gY_norm <  thr)

  out <- list(
    Z                   = Z,
    G                   = matrix(rep(G, p), n, p),
    Y                   = Y,
    W                   = W,
    U1                  = U[, 1],
    M                   = M,
    synthetic_data      = cov_df,
    true_total          = true_total,
    true_NDE            = true_NDE,
    true_NIE            = true_NIE,
    U                   = U,
    genetic_instrument  = gi,
    successful_features = successful,
    failed_features     = failed,
    params              = list(n = n, n_features = p, n_confounders = k,
                               n_mediators = n_mediators,
                               beta_Z = beta_Z, alpha_M = alpha_M, beta_M = beta_M,
                               conf_strength = conf_strength, coverage = coverage,
                               mo_confounding = mo_confounding,
                               pleio = pleio,
                               phi = phi,
                               gamma_G = gamma_G,
                               rho_G1 = rho_G1, rho_G2 = rho_G2,
                               rho_pop = rho_pop, separate_U = separate_U,
                               omega_1 = omega_1, omega_2 = omega_2,
                               nc_model = if (is.character(nc_model)) nc_model else "custom",
                               # v0.9.2 (JYH #864): record the effective
                               # per-confounder strength and per-control
                               # coverage vectors used.
                               u_strength = u_strength,
                               w_coverage_profile = if (!is.null(wcp)) wcp else list(w1 = om1_vec, w2 = om2_vec))
  )
  if (phi > 0) out$Gm <- Gm
  out$n_mediators <- n_mediators

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


#' Covariates for one iteration: from texture if available, else default (internal)
#' @keywords internal
.iteration_covariates <- function(trained_gan, n) {
  default <- data.frame(fetal_sex = rbinom(n, 1, 0.5),
                        gestational_age = rnorm(n))
  if (is.null(trained_gan)) return(default)

  tex  <- tryCatch(sample_texture(trained_gan, n), error = function(e) NULL)
  if (is.null(tex)) return(default)
  drop <- c("exposure_level", "outcome_level", "mediator_level")
  cv   <- tex[, setdiff(names(tex), drop), drop = FALSE]
  if (!ncol(cv)) return(default)
  cv
}

#' Exogenous, outcome-only texture component for one iteration (internal)
#'
#' Returns a centred n-vector taken from the generator's outcome channel. It is
#' added to every outcome feature to inject a realistic shared component, but is
#' never fed into Z, so it cannot open a confounding backdoor.
#' @keywords internal
.iteration_outcome_texture <- function(trained_gan, n) {
  if (is.null(trained_gan)) return(rnorm(n))
  tex <- tryCatch(sample_texture(trained_gan, n), error = function(e) NULL)
  if (is.null(tex) || !"outcome_level" %in% names(tex)) return(rnorm(n))
  as.numeric(scale(tex$outcome_level))
}

#' Exogenous, mediator-only texture component for one iteration (internal)
#'
#' Returns a p x n matrix of realistic mediator feature vectors drawn from
#' the feature-level copula texture model (v0.9.0).  Each row is a feature
#' (transcript), each column a sample.  The draws preserve the marginal
#' distributions and cross-feature correlation structure of the user's
#' mediator panel.  The texture is drawn independently of U and never fed
#' into Z, so it cannot open a confounding backdoor.
#'
#' When no feature_texture is available (trained_gan is NULL or lacks
#' $feature_texture), returns NULL and the caller falls back to scalar
#' noise.
#'
#' @param trained_gan An iconic_gan object (may be NULL).
#' @param n           Number of synthetic samples.
#' @param p           Number of features (target n_features for this iteration).
#' @return A p x n matrix, or NULL.
#' @keywords internal
.iteration_mediator_texture_features <- function(trained_gan, n, p) {
  if (is.null(trained_gan)) return(NULL)
  ft <- trained_gan$feature_texture %||% NULL
  if (is.null(ft)) return(NULL)
  tryCatch(sample_feature_texture(ft, n_samples = n, n_features = p),
           error = function(e) NULL)
}

#' Match a correlation matrix to the target feature count (internal)
#'
#' If `cor_mat` is a valid `p_target x p_target` matrix, return it as-is.
#' If it's a different size, return NULL (the simulation falls back to
#' independent noise).  This guards against dimension mismatches when the
#' simulation uses a different `n_features` than the training data.
#'
#' @param cor_mat   A correlation matrix or NULL.
#' @param p_target  Target number of features.
#' @return A `p_target x p_target` correlation matrix, or NULL.
#' @keywords internal
.match_cor <- function(cor_mat, p_target) {
  if (is.null(cor_mat)) return(NULL)
  if (!is.matrix(cor_mat)) return(NULL)
  if (nrow(cor_mat) != p_target || ncol(cor_mat) != p_target) return(NULL)
  cor_mat
}

#' Generate outcome noise, optionally correlated across features (internal)
#'
#' When `noise_cor` is a valid `p x p` correlation matrix, draws `n` samples
#' from `MVN(0, 0.2^2 * noise_cor)`.  Otherwise uses independent
#' `rnorm(n, 0, 0.2)` (backward compatible).
#'
#' @param n         Number of samples.
#' @param p         Number of features.
#' @param noise_cor A `p x p` correlation matrix, or NULL.
#' @return An `n x p` numeric matrix of outcome noise.
#' @keywords internal
.generate_outcome_noise <- function(n, p, noise_cor) {
  if (!is.null(noise_cor) && is.matrix(noise_cor) &&
      nrow(noise_cor) == p && ncol(noise_cor) == p) {
    Sigma <- noise_cor * (0.2^2)
    noise <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    if (n == 1) noise <- matrix(noise, nrow = 1)
    noise
  } else {
    matrix(rnorm(n * p, 0, 0.2), n, p)
  }
}
