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
# (mediator-outcome confounding) for the mediation extension. When 0,
# the original DGP is preserved exactly (backward compatible).
#
# the phi parameter (default 0) adds a mediator-specific genetic
# instrument Gm (e.g. a cis-eQTL for a mediator transcript). When phi > 0,
# Gm affects M but is independent of U and has
# no direct path to Y, making the 2-stage MR mediation estimator
# (fit_iv2sls_mediation2) point-identified under M-O confounding. When
# phi = 0 (default), no Gm is generated and the original DGP is preserved.
#
# imperfect instrument independence. Adds rho_G1, rho_G2, rho_pop,
# lambda_XM, lambda_MY, omega_1, omega_2 — see generate_toy_data() for the full
# description. Confounding is structured by per-path loading vectors over a
# k-dimensional confounder space: the X->M path confounds on
# conf_XM = U %*% lambda_XM and the M->Y path on conf_MY = U %*% lambda_MY.
# Default shared loadings reproduce the single-confounder behavior; distinct
# loadings give many-confounders-per-path structure. Path-specific negative
# controls W1 and W2 are generated from the pluggable nc_model with independent
# coverage of each path's composite. All new parameters default to
# default-reproducing values.
#
# Identification note: the generator draws texture INDEPENDENTLY of the latent
# confounders U and never feeds a shared texture factor into both X and Y, so
# the only X<->Y backdoor is through U. That keeps true_total a closed-form
# function of the causal knobs.
# ============================================================

#' Generate one synthetic dataset under the generalised SCM
#'
#' @param trained_gan Optional `iconic_gan` from [train_gan_on_real_data()],
#' supplying realistic covariate/outcome texture. If `NULL`, default synthetic
#' covariates are used.
#' @param n_synthetic_samples Sample size. Default 500.
#' @param n_features Number of outcome (and control) features. Default 20.
#' @param n_confounders Number of latent confounders `k`. Default 1
#' (backward-compatible single-confounder model).
#' @param n_mediators Number of independent mediators. When > 1,
#' each mediator has its own genetic instrument Gm and contributes
#' additively to Y. Default 1 (single mediator, backward compatible).
#' @param beta_X Direct effect of X on Y (true NDE). Default 0.10.
#' @param alpha_M Effect of X on the mediator M. Default 0.50.
#' @param beta_M Effect of M on Y. Default 0.30.
#' @param effect_size Optional shortcut: if non-`NULL`, sets a pure direct total
#' effect (`beta_X = effect_size`, `alpha_M = beta_M = 0`). Use `0` for null
#' (Type I error) simulations. Default `NULL` (use the mediation parameters).
#' @param conf_strength Overall confounding strength (analogue of `conf_str`);
#' scales the confounder loadings into X and Y. Default 0.80.
#' @param coverage How well the negative controls span the confounder
#' subspace, in `[0, 1]` (passed to `nc_model`). Default 1.
#' @param captured Integer indices of the confounders the negative controls
#' see (passed to `nc_model`). Default all `k`.
#' @param nc_model Negative-control model: a function `(U, covariates,
#' params) -> W`, or a registered name (`"proxy"`, `"cpg"`). Default `"proxy".
#' @param nc_params Extra named parameters forwarded to `nc_model`.
#' @param mo_confounding Strength of U1 -> M (mediator-outcome confounding).
#' 0 = no M-O confounding (original DGP). Default 0. When > 0, the first
#' confounder U\[,1\] also affects M, creating the M-O confounding that the
#' mediation estimators are benchmarked against.
#' @param pleio Strength of a direct G -> Y path (horizontal pleiotropy),
#' violating the exclusion restriction. 0 = no pleiotropy (valid instrument,
#' original DGP). Default 0. When > 0, G affects Y directly in addition to
#' through X, allowing benchmarking of IV/2SLS under instrument invalidity.
#' @param phi Strength of the mediator instrument Gm -> M.
#' 0 = no mediator instrument (original DGP, no Gm generated). Default 0.
#' When > 0, a valid instrument for M is generated (independent of U and G,
#' no direct path to Y), enabling point identification of NDE/NIE via
#' [fit_iv2sls_mediation2()] even under M-O confounding.
#' @param gamma_G Strength of the exposure instrument G -> X. Default 0.6.
#' @param rho_G1 Correlation of G1 with conf_XM (instrument exogeneity
#'   violation). Default 0.
#' @param rho_G2 Correlation of G2 with conf_MY (instrument exogeneity
#'   violation). Default 0.
#' @param rho_pop Shared population structure inducing G1-G2 correlation
#'   Default 0.
#' @param lambda_XM Optional length-k loading vector giving each confounder's
#'   weight on the X->M backdoor path. NULL (default) = shared loadings.
#' @param lambda_MY Optional length-k loading vector giving each confounder's
#'   weight on the M->Y backdoor path. NULL (default) = shared loadings.
#' @param omega_1 Coverage of conf_XM by W1. NULL = use `coverage`.
#' @param omega_2 Coverage of conf_MY by W2. NULL = use `coverage`.
#' @param feat_cor Within-module feature correlation. When > 0 and
#' the GAN does not provide feature_correlations, a block-diagonal
#' correlation matrix is used for the NC and outcome noise. GAN-learned
#' correlations take precedence. Default 0.
#' @param u_strength Numeric vector: per-confounder strength
#' profile (length k). Default NULL → rep(1, k) (equal strength, backward
#' compatible). Recycled to length k and normalized so the total confounding
#' budget is unchanged.
#' @param w_coverage_profile A list with optional \code{w1} and \code{w2}
#' numeric vectors: per-control coverage of conf_XM / conf_MY
#' (length n_features). Default NULL → scalar omega applied uniformly.
#' @param MMExp,MMOut,MMCon,MMCpG Per-pathway confounding multipliers (exposure,
#' outcome, controls, methylation). Default 1.
#' @param outcome_type \code{"continuous"} (default), \code{"survival"},
#'   or \code{"binary"}.
#'   When \code{"survival"}, the linear predictor is converted to
#' \code{surv_time} and \code{surv_event} via an exponential PH model.
#'   When \code{"binary"}, the linear predictor is converted to a 0/1
#' outcome \code{y_bin} via a logistic (Bernoulli) model; true effects
#' are on the conditional log-odds-ratio scale.
#' @param surv_h0 Baseline hazard for the survival DGP. Default 0.1.
#' @param surv_event_frac Target fraction of observed events. Default 0.6.
#' @param surv_censor_rate Explicit censoring rate. Default NULL.
#' @param bin_prev Target prevalence (marginal event probability) for the
#' binary DGP. Default 0.5.
#' @param seed Optional RNG seed.
#' @param separate_U Defunct. Passing a value errors with a
#'   message pointing to the replacement per-path loading vectors
#'   \code{lambda_XM} / \code{lambda_MY}. Retained in the signature only to
#'   catch and redirect old calls.
#'
#' @return A named list matching [generate_toy_data()] — `X`, `G` (n x
#' n_features), `Y`, `W`, `U1`, `M`, `synthetic_data`, `true_total`,
#' `true_NDE`, `true_NIE` — plus `U` (full n x k confounder matrix),
#' `genetic_instrument`, `successful_features`, `failed_features`, and
#' `params`. When `phi > 0`, also includes `Gm` (numeric vector, length n,
#' or `n_mediators x n` matrix when `n_mediators > 1`).
#' When any parameter is non-default, also includes `G1`, `G2`,
#' `W1`, `W2`, `conf_XM`, `conf_MY`, and (when `rho_pop > 0`) `P`.
#' When `n_mediators > 1`, `M` is an `n_mediators x n` matrix.
#' @export
#'
#' @examples
#' dat <- run_single_iteration(NULL, n_synthetic_samples = 100,
#'     n_features = 5, n_confounders = 1, seed = 1)
#' analyze_methods_robust(dat)
run_single_iteration <- function(trained_gan = NULL,
                                 n_synthetic_samples = 500,
                                 n_features = 20,
                                 n_confounders = 1,
                                 n_mediators = 1,
                                 beta_X = 0.10,
                                 alpha_M = 0.50,
                                 beta_M = 0.30,
                                 effect_size = NULL,
                                 conf_strength = 0.80,
                                 coverage = 1,
                                 captured = NULL,
                                 nc_model = "proxy",
                                 nc_params = list(),
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
                                 # per-control coverage profiles for heterogeneous
                                 # U/W characterization in the A3 simulation.
                                 u_strength = NULL,
                                 w_coverage_profile = NULL,
                                 MMExp = 1, MMOut = 1, MMCon = 1, MMCpG = 1,
                                 # survival / time-to-event outcome.
                                 outcome_type = c("continuous", "survival", "binary"),
                                 surv_h0 = 0.1,
                                 surv_event_frac = 0.6,
                                 surv_censor_rate = NULL,
                                 bin_prev = 0.5,
                                 seed = NULL) {
    if (!is.null(seed)) withr::local_seed(seed)
    outcome_type <- match.arg(outcome_type)
    .rsi_check_defunct(separate_U)

    if (!is.null(effect_size)) {
        beta_X <- effect_size
        alpha_M <- 0
        beta_M <- 0
    }
    true_NDE <- beta_X
    true_NIE <- n_mediators * alpha_M * beta_M
    true_total <- true_NDE + true_NIE

    n <- n_synthetic_samples
    p <- n_features
    k <- max(1L, as.integer(n_confounders))
    nc <- .resolve_nc_model(nc_model)
    wcp <- .rsi_coverage_profile(w_coverage_profile, p)

    # Is the imperfect-independence DGP active?
    v05_active <- rho_G1 != 0 || rho_G2 != 0 || rho_pop != 0 ||
        !is.null(lambda_XM) || !is.null(lambda_MY) ||
        !is.null(omega_1) || !is.null(omega_2)

    # Latent confounders and per-path composites (conf_XM, conf_MY).
    u_strength <- .rsi_u_strength(u_strength, k)
    U <- matrix(rnorm(n * k), n, k)
    lambda_XM <- .rsi_resolve_lambda(lambda_XM, k, u_strength)
    lambda_MY <- .rsi_resolve_lambda(lambda_MY, k, u_strength)
    conf_XM <- as.numeric(U %*% lambda_XM)
    conf_MY <- as.numeric(U %*% lambda_MY)

    # Realistic covariate / outcome / mediator texture, independent of U.
    cov_df <- .iteration_covariates(trained_gan, n)
    exo_Y <- .iteration_outcome_texture(trained_gan, n)
    exo_M_features <- .iteration_mediator_texture_features(trained_gan, n, p)
    noise_cor <- .rsi_noise_cor(trained_gan, p, feat_cor)

    inst <- .rsi_instruments(v05_active, n, n_mediators, phi, rho_G1, rho_G2,
        rho_pop, conf_XM, conf_MY)

    # Structural exposure, scaled to unit variance before building M and Y
    # (mirrors generate_toy_data()).
    aX <- conf_strength * MMExp * u_strength / sqrt(k)
    X <- as.numeric(scale(gamma_G * inst$G + as.numeric(U %*% aX) +
        rnorm(n, 0, 0.5)))

    M <- .rsi_mediator(X, n, n_mediators, alpha_M, mo_confounding, conf_MY,
        phi, inst$Gm, exo_M_features)
    yy <- .rsi_outcome(X, M, U, inst$G, exo_Y, n, p, k, beta_X, beta_M,
        conf_strength, MMOut, u_strength, pleio, noise_cor$Y)

    ncp <- .rsi_negative_controls(
        v05_active, nc, U, conf_XM, conf_MY, cov_df, n, p, coverage, captured,
        omega_1, omega_2, wcp, lambda_XM, lambda_MY, MMCon, MMCpG,
        noise_cor$W, nc_params
    )

    # Features with strong vs weak confounding (reference-compatible
    # bookkeeping).
    thr <- stats::median(yy$gY_norm)

    out <- list(
        X = X, G = matrix(rep(inst$G, p), n, p), Y = yy$Y, W = ncp$W,
        U1 = U[, 1], M = M, synthetic_data = cov_df,
        true_total = true_total, true_NDE = true_NDE, true_NIE = true_NIE,
        U = U, genetic_instrument = inst$gi,
        successful_features = which(yy$gY_norm >= thr),
        failed_features = which(yy$gY_norm < thr),
        params = list(
            n = n, n_features = p, n_confounders = k, n_mediators = n_mediators,
            beta_X = beta_X, alpha_M = alpha_M, beta_M = beta_M,
            conf_strength = conf_strength, coverage = coverage,
            mo_confounding = mo_confounding, pleio = pleio, phi = phi,
            gamma_G = gamma_G, rho_G1 = rho_G1, rho_G2 = rho_G2, rho_pop = rho_pop,
            lambda_XM = lambda_XM, lambda_MY = lambda_MY,
            omega_1 = omega_1, omega_2 = omega_2,
            nc_model = if (is.character(nc_model)) nc_model else "custom",
            u_strength = u_strength,
            w_coverage_profile = if (!is.null(wcp)) wcp else
                list(w1 = ncp$om1_vec, w2 = ncp$om2_vec),
            outcome_type = outcome_type,
            surv_h0 = surv_h0, surv_event_frac = surv_event_frac,
            surv_censor_rate = surv_censor_rate, bin_prev = bin_prev
        )
    )
    out <- .rsi_add_outcome(out, outcome_type, yy$Y[, 1], surv_h0,
        surv_event_frac, surv_censor_rate, bin_prev)
    if (phi > 0) out$Gm <- inst$Gm
    out$n_mediators <- n_mediators
    if (v05_active) {
        out$G1 <- inst$G1
        out$G2 <- inst$G2
        out$W1 <- ncp$W1
        out$W2 <- ncp$W2
        out$conf_XM <- conf_XM
        out$conf_MY <- conf_MY
        if (rho_pop > 0) out$P <- inst$P
    }
    out
}

# ------------------------------------------------------------
# Internal helpers for run_single_iteration(). Each one is a single stage
# of the generator; they are called in a fixed order so the random-number
# stream (and hence every seeded dataset) matches the original
# single-function implementation.
# ------------------------------------------------------------

#' Error on the defunct separate_U argument (internal)
#' @keywords internal
#' @noRd
.rsi_check_defunct <- function(separate_U) {
    # The separate_U toggle was replaced by per-path confounder loading
    # vectors lambda_XM / lambda_MY.
    if (!is.null(separate_U))
        stop("argument `separate_U` was removed. Confounding is now ",
            "specified by per-path loading vectors `lambda_XM` and `lambda_MY` ",
            "over a k-dimensional confounder space. Shared loadings (the default) ",
            "reproduce the old separate_U = FALSE; lambda_XM = e1, lambda_MY = e2 ",
            "at k = 2 reproduce the old separate_U = TRUE.", call. = FALSE)
    invisible(NULL)
}

#' Validate and recycle a per-control coverage profile (internal)
#'
#' A list with optional `w1` / `w2` numeric vectors (recycled to length p)
#' giving per-control coverage of conf_XM / conf_MY. NULL means the scalar
#' omega is applied uniformly (legacy).
#' @keywords internal
#' @noRd
.rsi_coverage_profile <- function(wcp, p) {
    if (is.null(wcp)) return(NULL)
    if (!is.list(wcp)) stop("w_coverage_profile must be a list with 'w1'/'w2'")
    if (!is.null(wcp$w1)) wcp$w1 <- rep_len(as.numeric(wcp$w1), p)
    if (!is.null(wcp$w2)) wcp$w2 <- rep_len(as.numeric(wcp$w2), p)
    wcp
}

#' Resolve the per-confounder strength profile to length k (internal)
#'
#' Default rep(1, k) gives equal-strength confounders. A supplied profile is
#' recycled to length k and renormalised so sum(u_strength) == k, keeping the
#' overall confounding budget (and conf_strength) comparable.
#' @keywords internal
#' @noRd
.rsi_u_strength <- function(u_strength, k) {
    if (is.null(u_strength)) return(rep(1, k))
    u_strength <- rep_len(as.numeric(u_strength), k)
    u_strength * k / sum(u_strength)
}

#' Resolve a per-path confounder loading vector (internal)
#'
#' NULL gives shared loadings (normalised u_strength), so the same composite
#' drives both backdoor paths. lambda_XM = e1, lambda_MY = e2 at k = 2
#' recovers fully separate path confounders.
#' @keywords internal
#' @noRd
.rsi_resolve_lambda <- function(lam, k, u_strength) {
    if (is.null(lam)) return(u_strength / sqrt(sum(u_strength^2)))
    lam <- rep_len(as.numeric(lam), k)
    nrm <- sqrt(sum(lam^2))
    if (nrm == 0) stop("lambda vectors must have nonzero norm.", call. = FALSE)
    lam / nrm
}

#' Feature-level noise correlation matrices for W and Y (internal)
#'
#' Uses the texture model's learned residual correlations when present, so
#' simulated features stay correlated conditional on U; otherwise falls back
#' to a synthetic block-diagonal matrix when feat_cor > 0.
#' @keywords internal
#' @noRd
.rsi_noise_cor <- function(trained_gan, p, feat_cor) {
    gan_feat_cor <- trained_gan$feature_correlations %||% list()
    out <- list(W = gan_feat_cor$W, Y = gan_feat_cor$Y)
    if (is.null(out$W) && feat_cor > 0) out$W <- .block_diag_cor(p, feat_cor)
    if (is.null(out$Y) && feat_cor > 0) out$Y <- .block_diag_cor(p, feat_cor)
    out
}

#' Exposure (G) and mediator (Gm) instruments for one iteration (internal)
#'
#' Default DGP: G from simulate_single_genetic_instrument(), Gm ~ N(0, 1).
#' Imperfect-independence DGP:
#'   G1 = rho_G1*conf_XM + rho_pop*P + sqrt(1-rho_G1^2-rho_pop^2)*N(0,1)
#'   G2 = rho_G2*conf_MY + rho_pop*P + sqrt(1-rho_G2^2-rho_pop^2)*N(0,1)
#' with one Gm row per mediator when n_mediators > 1.
#' @keywords internal
#' @noRd
.rsi_instruments <- function(v05_active, n, n_mediators, phi, rho_G1, rho_G2,
                             rho_pop, conf_XM, conf_MY) {
    if (!v05_active) {
        gi <- simulate_single_genetic_instrument(n)
        Gm <- NULL
        if (phi > 0) {
            Gm <- if (n_mediators > 1) {
                matrix(rnorm(n_mediators * n), n_mediators, n)
            } else {
                rnorm(n)
            }
        }
        return(list(G = gi$G, Gm = Gm, G1 = gi$G, G2 = Gm, P = NULL, gi = gi))
    }
    P <- rnorm(n) # shared population structure
    resid_var_G1 <- max(1 - rho_G1^2 - rho_pop^2, 0)
    G1 <- as.numeric(scale(rho_G1 * conf_XM + rho_pop * P +
        sqrt(resid_var_G1) * rnorm(n)))
    G2 <- NULL
    if (phi > 0) {
        resid_var_G2 <- max(1 - rho_G2^2 - rho_pop^2, 0)
        draw_g2 <- function() {
            as.numeric(scale(rho_G2 * conf_MY + rho_pop * P +
                sqrt(resid_var_G2) * rnorm(n)))
        }
        if (n_mediators > 1) {
            G2 <- matrix(NA_real_, n_mediators, n)
            for (mm in seq_len(n_mediators)) G2[mm, ] <- draw_g2()
        } else {
            G2 <- draw_g2()
        }
    }
    list(G = G1, Gm = G2, G1 = G1, G2 = G2, P = P,
        gi = list(G = G1, dosages = NULL, maf = NULL))
}

#' Simulate the mediator(s) for one iteration (internal)
#'
#' M = alpha_M * X, plus optional conf_MY -> M (mediator-outcome
#' confounding), optional Gm -> M (mediator instrument), and exogenous
#' texture: a feature-level copula draw when available, otherwise a scalar
#' fallback. Returns a vector, or an n_mediators x n matrix.
#' @keywords internal
#' @noRd
.rsi_mediator <- function(X, n, n_mediators, alpha_M, mo_confounding, conf_MY,
                          phi, Gm, exo_M_features) {
    one_mediator <- function(gm, tex) {
        m <- alpha_M * X
        if (mo_confounding > 0) m <- m + mo_confounding * 0.5 * conf_MY
        if (phi > 0) m <- m + phi * gm
        if (!is.null(exo_M_features)) {
            m + tex + rnorm(n, 0, 0.02)
        } else {
            m + 0.3 * rnorm(n) + rnorm(n, 0, 0.05)
        }
    }
    if (n_mediators == 1) {
        tex <- if (is.matrix(exo_M_features)) exo_M_features[1, ] else
            exo_M_features
        return(one_mediator(Gm, tex))
    }
    M <- matrix(NA_real_, n_mediators, n)
    for (mm in seq_len(n_mediators)) {
        gm <- if (is.matrix(Gm)) Gm[mm, ] else Gm
        tex <- if (is.matrix(exo_M_features)) exo_M_features[mm, ] else
            exo_M_features
        M[mm, ] <- one_mediator(gm, tex)
    }
    M
}

#' Simulate the outcome panel for one iteration (internal)
#'
#' Each feature gets the causal path, its own confounder loadings (scaled by
#' u_strength), optional pleiotropy, a shared exogenous texture term (no X
#' backdoor) and idiosyncratic noise, optionally correlated across features.
#' Multiple mediators contribute additively. Returns Y and the per-feature
#' confounding norm used for feature bookkeeping.
#' @keywords internal
#' @noRd
.rsi_outcome <- function(X, M, U, G, exo_Y, n, p, k, beta_X, beta_M,
                         conf_strength, MMOut, u_strength, pleio,
                         noise_cor_Y) {
    Y <- matrix(NA_real_, n, p)
    gY_norm <- numeric(p)
    M_eff <- if (is.matrix(M)) colSums(M) else M
    for (f in seq_len(p)) {
        gY <- conf_strength * MMOut * u_strength * runif(k, 0.4, 0.8) / sqrt(k)
        gY_norm[f] <- sqrt(sum(gY^2))
        Y[, f] <- beta_M * M_eff + beta_X * X + as.numeric(U %*% gY) +
            pleio * G + 0.4 * exo_Y
    }
    list(Y = Y + .generate_outcome_noise(n, p, noise_cor_Y), gY_norm = gY_norm)
}

#' Build an n x p NC panel with per-column coverage (internal)
#'
#' Calls the nc_model once per distinct coverage value and assigns those
#' columns, preserving the model's texture within each coverage stratum.
#' @keywords internal
#' @noRd
.rsi_nc_panel_per_coverage <- function(nc, U_col, om_vec, cov_df, n, p,
                                       MMCon, MMCpG, noise_cor_W, nc_params) {
    U_mat <- matrix(U_col, n, 1)
    out <- matrix(NA_real_, n, p)
    for (ov in unique(om_vec)) {
        cols <- which(om_vec == ov)
        nc_p <- utils::modifyList(
            list(n_features = length(cols), coverage = ov, captured = 1,
                MMCon = MMCon, MMCpG = MMCpG,
                noise_cor = .match_cor(noise_cor_W, length(cols))), nc_params)
        Wc <- .rsi_fit_nc_columns(nc(U_mat, cov_df, nc_p), length(cols))
        out[, cols] <- Wc
    }
    out
}

#' Coerce an nc_model result to a matrix with exactly `p` columns (internal)
#' @keywords internal
#' @noRd
.rsi_fit_nc_columns <- function(W, p) {
    if (!is.matrix(W)) W <- as.matrix(W)
    if (ncol(W) != p) W <- W[, rep(seq_len(ncol(W)), length.out = p),
        drop = FALSE]
    W
}

#' Negative-control panels for one iteration (internal)
#'
#' Default DGP: a single W panel from the pluggable nc_model. Imperfect-
#' independence DGP: path-specific W1 (captures conf_XM) and W2 (captures
#' conf_MY) with coverage omega_1 / omega_2 or per-control profiles; the
#' combined W is W1 when both paths share loadings, else the scaled average.
#' @keywords internal
#' @noRd
.rsi_negative_controls <- function(v05_active, nc, U, conf_XM, conf_MY,
                                   cov_df, n, p, coverage, captured,
                                   omega_1, omega_2, wcp, lambda_XM,
                                   lambda_MY, MMCon, MMCpG, noise_cor_W,
                                   nc_params) {
    om1 <- if (is.null(omega_1)) coverage else omega_1
    om2 <- if (is.null(omega_2)) coverage else omega_2
    om1_vec <- if (!is.null(wcp$w1)) wcp$w1 else rep(om1, p)
    om2_vec <- if (!is.null(wcp$w2)) wcp$w2 else rep(om2, p)

    if (v05_active) {
        panel <- function(U_col, om_vec) {
            .rsi_nc_panel_per_coverage(nc, U_col, om_vec, cov_df, n, p, MMCon,
                MMCpG, noise_cor_W, nc_params)
        }
        W1 <- panel(conf_XM, om1_vec)
        W2 <- panel(conf_MY, om2_vec)
        W <- if (identical(lambda_XM, lambda_MY)) W1 else scale((W1 + W2) / 2)
    } else {
        nc_p <- utils::modifyList(
            list(n_features = p, coverage = coverage, captured = captured,
                MMCon = MMCon, MMCpG = MMCpG,
                noise_cor = .match_cor(noise_cor_W, p)), nc_params)
        W <- .rsi_fit_nc_columns(nc(U, cov_df, nc_p), p)
        W1 <- W
        W2 <- W
    }
    list(W = W, W1 = W1, W2 = W2, om1_vec = om1_vec, om2_vec = om2_vec)
}

#' Attach the survival or binary outcome to a simulated dataset (internal)
#'
#' Survival: the first outcome column's linear predictor is converted to
#' surv_time / surv_event (true effects on the Cox log-HR scale). Binary: it
#' is converted to a 0/1 outcome via a logistic model (log-OR scale).
#' @keywords internal
#' @noRd
.rsi_add_outcome <- function(out, outcome_type, linpred, surv_h0,
                             surv_event_frac, surv_censor_rate, bin_prev) {
    if (outcome_type == "survival") {
        surv <- .linpred_to_surv(linpred, h0 = surv_h0,
            event_frac = surv_event_frac,
            censor_rate = surv_censor_rate)
        out$surv_time <- surv$surv_time
        out$surv_event <- surv$surv_event
    } else if (outcome_type == "binary") {
        out$y_bin <- .linpred_to_binary(linpred, prev = bin_prev)$y_bin
    }
    out$outcome_type <- outcome_type
    out
}


#' Covariates for one iteration: from texture if available, else default (internal)
#' @keywords internal
#' @noRd
.iteration_covariates <- function(trained_gan, n) {
    default <- data.frame(fetal_sex = rbinom(n, 1, 0.5),
        gestational_age = rnorm(n))
    if (is.null(trained_gan)) return(default)

    tex <- tryCatch(sample_texture(trained_gan, n), error = function(e) NULL)
    if (is.null(tex)) return(default)
    drop <- c("exposure_level", "outcome_level", "mediator_level")
    cv <- tex[, setdiff(names(tex), drop), drop = FALSE]
    if (!ncol(cv)) return(default)
    cv
}

#' Exogenous, outcome-only texture component for one iteration (internal)
#'
#' Returns a centred n-vector taken from the generator's outcome channel. It is
#' added to every outcome feature to inject a realistic shared component, but is
#' never fed into X, so it cannot open a confounding backdoor.
#' @keywords internal
#' @noRd
.iteration_outcome_texture <- function(trained_gan, n) {
    if (is.null(trained_gan)) return(rnorm(n))
    tex <- tryCatch(sample_texture(trained_gan, n), error = function(e) NULL)
    if (is.null(tex) || !"outcome_level" %in% names(tex)) return(rnorm(n))
    as.numeric(scale(tex$outcome_level))
}

#' Exogenous, mediator-only texture component for one iteration (internal)
#'
#' Returns a p x n matrix of realistic mediator feature vectors drawn from
#' the feature-level copula texture model. Each row is a feature
#' (transcript), each column a sample. The draws preserve the marginal
#' distributions and cross-feature correlation structure of the user's
#' mediator panel. The texture is drawn independently of U and never fed
#' into X, so it cannot open a confounding backdoor.
#'
#' When no feature_texture is available (trained_gan is NULL or lacks
#' $feature_texture), returns NULL and the caller falls back to scalar
#' noise.
#'
#' @param trained_gan An iconic_gan object (may be NULL).
#' @param n Number of synthetic samples.
#' @param p Number of features (target n_features for this iteration).
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
#' independent noise). This guards against dimension mismatches when the
#' simulation uses a different `n_features` than the training data.
#'
#' @param cor_mat A correlation matrix or NULL.
#' @param p_target Target number of features.
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
#' from `MVN(0, 0.2^2 * noise_cor)`. Otherwise uses independent
#' `rnorm(n, 0, 0.2)` (backward compatible).
#'
#' @param n Number of samples.
#' @param p Number of features.
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
