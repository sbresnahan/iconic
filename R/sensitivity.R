# ============================================================
# Sensitivity analysis on generator-realistic data.
#
# gan_sensitivity()   – benchmark all estimators across a grid of confounding
#                       scenarios (confounding strength x negative-control
#                       coverage x number of confounders).
# recommend_estimator() – summarise which method is preferred, per scenario and
#                       robustly across scenarios.
# nc_validity_check() – probe specifically whether the negative controls hold as
#                       coverage drops / confounders multiply, including the
#                       proximal-identification check (#controls vs #confounders).
# gan_mediation_sensitivity() – mediation analogue: NDE/NIE bias and Type I
#                       error across confounding scenarios with M-O confounding.
#
# v0.3.1: PGC now uses the matrix bridge. The scalar-bridge variant
# (fit_pgc_scalar) is a standalone exported function, not included in
# the default sensitivity pipeline.
#
# v0.4.0: gan_mediation_sensitivity() gains a phi parameter for the
# mediator instrument (Gm).  When phi > 0, the 2-stage MR estimator
# (IV2SLS2) is included in the results.
#
# v0.9.4: survival / time-to-event outcomes.  When outcome_type = "survival",
# the simulation generates (surv_time, surv_event) via run_single_iteration()
# and estimates via iconic_estimate() with the survival drivers, instead of
# run_methods() / run_mediation_methods().  The summary functions are
# format-compatible because iconic_estimate() returns the same column names.
# ============================================================

# v0.9.4: Run survival estimators on a dat list from run_single_iteration().
# Converts dat to an iconic_data object and calls iconic_estimate().
# Returns a data frame with the same columns as run_methods() (total-effect)
# or run_mediation_methods() (mediation), so summarise_results() /
# summarise_mediation_results() consume it unchanged.
.run_surv_methods <- function(dat, effect_scale = "loghr", tau = NULL,
                              is_mediation = FALSE) {
  G_vec <- if (!is.null(dat$G)) {
    if (is.matrix(dat$G)) dat$G[, 1] else dat$G
  } else NULL
  Gm_vec <- if (!is.null(dat$Gm)) {
    if (is.matrix(dat$Gm)) dat$Gm[1, ] else dat$Gm
  } else NULL

  if (is_mediation) {
    sdat <- iconic_data(
      Z = dat$Z, outcome_type = "survival",
      surv_time = dat$surv_time, surv_event = dat$surv_event,
      M = dat$M, G = G_vec, Gm = Gm_vec,
      W = t(dat$W))
  } else {
    sdat <- iconic_data(
      Z = dat$Z, outcome_type = "survival",
      surv_time = dat$surv_time, surv_event = dat$surv_event,
      G = G_vec, W = t(dat$W))
  }
  iconic_estimate(sdat, effect_scale = effect_scale, tau = tau)
}

#' Benchmark estimators across confounding scenarios on synthetic data
#'
#' For each scenario in the grid, generates `n_iter` synthetic datasets with
#' [run_single_iteration()] (using the trained generator for texture), runs every
#' estimator, and summarises bias / RMSE / power. This is the multi-confounder,
#' negative-control-aware generalisation of the package's parameter sweeps.
#'
#' @param trained_gan   An `iconic_gan` (or `NULL` to use default texture).
#' @param conf_grid     Confounding-strength values to sweep. Default `c(0.2, 0.5, 0.8)`.
#' @param coverage_grid Negative-control coverage values in `[0,1]`. Default `c(0.3, 0.7, 1)`.
#' @param k_grid        Numbers of latent confounders to sweep. Default `1`.
#' @param nc_model      Negative-control model (function or name). Default `"proxy"`.
#' @param n_iter        Replicates per scenario. Default 50.
#' @param n_samples     Samples per replicate. Default 500.
#' @param n_features    Features per replicate. Default 20.
#' @param beta_Z,alpha_M,beta_M Causal paths (ground truth). Defaults 0.10 / 0.50 / 0.30.
#' @param effect_size   Optional pure-direct total effect override (see
#'   [run_single_iteration()]).
#' @param base_seed     Base RNG seed. Default 700.
#' @param n_cores       Parallel workers across replicates. Default 1.
#' @param outcome_type  \code{"continuous"} (default) or \code{"survival"}
#'   (v0.9.4).  When survival, the DGP generates time-to-event outcomes and
#'   estimation uses the Cox / RMST survival drivers via
#'   [iconic_estimate()].
#' @param effect_scale  \code{"loghr"} (default) or \code{"rmst"}.  Only
#'   used when \code{outcome_type = "survival"}.
#' @param surv_h0       Baseline hazard for the survival DGP (v0.9.4). Default 0.1.
#' @param surv_event_frac Target fraction of observed events (v0.9.4). Default 0.6.
#' @param surv_censor_rate Explicit censoring rate (v0.9.4). Default NULL.
#'
#' @return A list with `summary` (one row per scenario x method, with
#'   `conf_strength`, `coverage`, `k`, `true_total` and the columns from
#'   `summarise_results`) and `grid` (the scenario grid).
#' @export
gan_sensitivity <- function(trained_gan  = NULL,
                            conf_grid     = c(0.2, 0.5, 0.8),
                            coverage_grid = c(0.3, 0.7, 1),
                            k_grid        = 1,
                            nc_model      = "proxy",
                            n_iter        = 50,
                            n_samples     = 500,
                            n_features    = 20,
                            beta_Z = 0.10, alpha_M = 0.50, beta_M = 0.30,
                            effect_size   = NULL,
                            base_seed     = 700,
                            n_cores       = 1,
                            outcome_type  = c("continuous", "survival"),
                            effect_scale  = c("loghr", "rmst"),
                            surv_h0       = 0.1,
                            surv_event_frac = 0.6,
                            surv_censor_rate = NULL) {
  outcome_type <- match.arg(outcome_type)
  effect_scale <- match.arg(effect_scale)

  grid <- expand.grid(conf_strength = conf_grid, coverage = coverage_grid,
                      k = k_grid, KEEP.OUT.ATTRS = FALSE)

  true_total <- if (is.null(effect_size)) beta_Z + alpha_M * beta_M else effect_size

  smry <- lapply(seq_len(nrow(grid)), function(gi) {
    cs <- grid$conf_strength[gi]; cov <- grid$coverage[gi]; kk <- grid$k[gi]

    worker <- function(i) {
      dat <- run_single_iteration(
        trained_gan, n_synthetic_samples = n_samples, n_features = n_features,
        n_confounders = kk, beta_Z = beta_Z, alpha_M = alpha_M, beta_M = beta_M,
        effect_size = effect_size, conf_strength = cs, coverage = cov,
        nc_model = nc_model, seed = base_seed + gi * 1000L + i,
        outcome_type = outcome_type, surv_h0 = surv_h0,
        surv_event_frac = surv_event_frac, surv_censor_rate = surv_censor_rate)
      if (outcome_type == "survival") {
        res <- .run_surv_methods(dat, effect_scale = effect_scale,
                                 is_mediation = FALSE)
      } else {
        res <- run_methods(dat, n_features)
      }
      res$iter <- i
      res
    }

    combined <- do.call(rbind, .parallel_lapply(seq_len(n_iter), worker, n_cores))
    s <- summarise_results(combined, true_total)
    s$conf_strength <- cs; s$coverage <- cov; s$k <- kk; s$true_total <- true_total
    s
  })

  summary <- do.call(rbind, smry)
  front   <- c("conf_strength", "coverage", "k", "true_total", "method")
  summary <- summary[, c(front, setdiff(names(summary), front))]
  list(summary = summary, grid = grid)
}


#' Recommend the preferred estimator from a sensitivity sweep
#'
#' Ranks methods by RMSE within each scenario, and identifies an overall pick
#' that is robust across scenarios (smallest worst-case RMSE). UNADJ is excluded
#' from recommendations (it is a bias reference).
#'
#' @param sens   Object returned by [gan_sensitivity()].
#' @param exclude Methods to exclude from the recommendation. Default
#'   `c("UNADJ")`.
#'
#' @return A list with `per_scenario` (best method + RMSE per scenario),
#'   `worst_case` (max RMSE per method across scenarios), and `overall` (the
#'   method minimising worst-case RMSE).
#' @export
recommend_estimator <- function(sens, exclude = c("UNADJ")) {
  s <- sens$summary
  s <- s[!s$method %in% exclude, , drop = FALSE]

  scen <- interaction(s$conf_strength, s$coverage, s$k, drop = TRUE)
  per_scenario <- do.call(rbind, lapply(split(s, scen), function(d) {
    d <- d[is.finite(d$rmse), , drop = FALSE]
    if (!nrow(d)) return(NULL)
    best <- d[which.min(d$rmse), ]
    data.frame(conf_strength = best$conf_strength, coverage = best$coverage,
               k = best$k, best_method = best$method, rmse = best$rmse,
               stringsAsFactors = FALSE)
  }))

  worst <- vapply(split(s$rmse, s$method), function(v) max(v, na.rm = TRUE), numeric(1))
  worst <- worst[is.finite(worst)]
  overall <- if (length(worst)) names(worst)[which.min(worst)] else NA_character_

  list(
    per_scenario = per_scenario,
    worst_case   = data.frame(method = names(worst), worst_rmse = as.numeric(worst),
                              row.names = NULL, stringsAsFactors = FALSE),
    overall      = overall
  )
}


#' Check negative-control validity across confounding scenarios
#'
#' Focuses on whether the negative-control estimators (COCA, PGC) hold as
#' the controls' coverage of the confounder subspace drops and as the number
#' of latent confounders grows. Reports their bias/RMSE alongside IV2SLS,
#' and flags scenarios where identification is not credible — in particular
#' when the number of latent confounders exceeds the effective number of
#' valid controls (the proximal-inference completeness condition).
#'
#' The matrix-bridge PGC is the estimator for which the completeness condition
#' is binding: when `k > n_valid_controls`, PGC is under-identified while
#' IV2SLS remains unbiased (it does not depend on NC completeness).
#'
#' @param trained_gan   An `iconic_gan` (or `NULL`).
#' @param coverage_grid Coverage values to sweep. Default `c(0.2, 0.5, 0.8, 1)`.
#' @param k_grid        Numbers of confounders. Default `c(1, 2, 3)`.
#' @param conf_strength Fixed confounding strength. Default 0.8.
#' @param n_valid_controls Number of *distinct valid* controls the design
#'   provides (for the identifiability flag). Default 1.
#' @param n_iter,n_samples,n_features,nc_model,base_seed,n_cores As in
#'   [gan_sensitivity()].
#'
#' @return A list with `summary` (COCA/PGC/IV2SLS bias & RMSE per scenario,
#'   with an `identified` flag) and `verdict` (short per-scenario diagnosis).
#' @export
nc_validity_check <- function(trained_gan   = NULL,
                              coverage_grid  = c(0.2, 0.5, 0.8, 1),
                              k_grid         = c(1, 2, 3),
                              conf_strength  = 0.8,
                              n_valid_controls = 1,
                              n_iter         = 50,
                              n_samples      = 500,
                              n_features     = 20,
                              nc_model       = "proxy",
                              base_seed      = 800,
                              n_cores        = 1) {

  sens <- gan_sensitivity(trained_gan, conf_grid = conf_strength,
                          coverage_grid = coverage_grid, k_grid = k_grid,
                          nc_model = nc_model, n_iter = n_iter,
                          n_samples = n_samples, n_features = n_features,
                          base_seed = base_seed, n_cores = n_cores)

  s <- sens$summary
  keep <- s$method %in% c("COCA", "PGC", "IV2SLS")
  s <- s[keep, c("coverage", "k", "method", "bias", "abs_bias", "rmse", "power"), drop = FALSE]

  # Proximal identifiability: need at least as many valid controls as confounders.
  s$identified <- s$k <= n_valid_controls

  # Best-case (highest-coverage) negative-control bias for each k, used as the
  # reference for judging whether controls *degrade* as coverage drops.
  # Only COCA and matrix-bridge PGC are NC-dependent; IV2SLS is
  # instrument-based and not affected by completeness.
  nc_bias <- function(d) mean(d$abs_bias[d$method %in% c("COCA", "PGC")], na.rm = TRUE)
  base_by_k <- vapply(split(s, s$k), function(dk) {
    top <- dk[dk$coverage == max(dk$coverage), , drop = FALSE]
    nc_bias(top)
  }, numeric(1))

  verdict <- do.call(rbind, lapply(split(s, interaction(s$coverage, s$k, drop = TRUE)),
    function(d) {
      kk       <- d$k[1]
      here     <- nc_bias(d)
      baseline <- base_by_k[[as.character(kk)]]
      degrading <- is.finite(here) && is.finite(baseline) && here > 1.5 * baseline + 1e-8
      data.frame(
        coverage = d$coverage[1], k = kk, identified = d$identified[1],
        nc_abs_bias = round(here, 3),
        diagnosis = if (!d$identified[1])
                      "under-identified: #confounders > #valid controls"
                    else if (degrading)
                      "negative controls degrading (low coverage)"
                    else "negative controls holding",
        stringsAsFactors = FALSE)
    }))

  if (any(!s$identified))
    warning("Some scenarios have more confounders than valid negative controls; ",
            "the matrix-bridge PGC and COCA are not point-identified there ",
            "(proximal completeness fails). IV2SLS remains valid as it does ",
            "not depend on NC completeness.")

  list(summary = s, verdict = verdict)
}


# ============================================================
# Mediation sensitivity analysis: extends gan_sensitivity() to the
# NDE/NIE setting with mediator-outcome confounding.
# ============================================================

#' Benchmark mediation estimators across confounding scenarios
#'
#' For each scenario in the grid, generates `n_iter` synthetic datasets
#' with mediator-outcome confounding (via `mo_confounding`), runs every
#' mediation estimator, and summarises NDE/NIE bias / RMSE / Type I error.
#' This is the mediation analogue of [gan_sensitivity()].
#'
#' When `phi > 0` (v0.4.0), a mediator-specific genetic instrument (Gm) is
#' generated and the 2-stage MR estimator (IV2SLS2) is included in the
#' results, enabling point identification of NDE/NIE under M-O confounding.
#'
#' @param trained_gan   An `iconic_gan` (or `NULL` to use default texture).
#' @param conf_grid     Confounding-strength values to sweep. Default `c(0.2, 0.5, 0.8)`.
#' @param coverage_grid Negative-control coverage values in `[0,1]`. Default `c(0.3, 0.7, 1)`.
#' @param k_grid        Numbers of latent confounders to sweep. Default `1`.
#' @param mo_confounding Strength of U1 -> M (mediator-outcome confounding). Default 0.80.
#' @param phi           Strength of the mediator instrument Gm -> M (v0.4.0).
#'                      0 = no mediator instrument (five estimators). > 0 =
#'                      generates Gm and includes the 2-stage MR estimator
#'                      (IV2SLS2). Default 0.
#' @param rho_G1        Correlation of G1 with U_XM (v0.5.0). Default 0.
#' @param rho_G2        Correlation of G2 with U_MY (v0.5.0). Default 0.
#' @param rho_pop       Shared population structure (v0.5.0). Default 0.
#' @param separate_U    Draw separate confounders for Z->M and M->Y (v0.5.0). Default FALSE.
#' @param omega_1       Coverage of U_XM by W1 (v0.5.0). NULL = use `coverage`.
#' @param omega_2       Coverage of U_MY by W2 (v0.5.0). NULL = use `coverage`.
#' @param nc_model      Negative-control model (function or name). Default `"proxy"`.
#' @param n_iter        Replicates per scenario. Default 50.
#' @param n_samples     Samples per replicate. Default 500.
#' @param n_features    Features per replicate. Default 20.
#' @param beta_Z,alpha_M,beta_M Causal paths (ground truth). Defaults 0.10 / 0.50 / 0.30.
#' @param base_seed     Base RNG seed. Default 750.
#' @param n_cores       Parallel workers across replicates. Default 1.
#' @param outcome_type  \code{"continuous"} (default) or \code{"survival"}
#'   (v0.9.4).  When survival, the DGP generates time-to-event outcomes and
#'   estimation uses the Cox / RMST survival mediation drivers via
#'   [iconic_estimate()].
#' @param effect_scale  \code{"loghr"} (default) or \code{"rmst"}.  Only
#'   used when \code{outcome_type = "survival"}.
#' @param surv_h0       Baseline hazard for survival DGP (v0.9.4).  See
#'   [run_single_iteration()].
#' @param surv_event_frac Target event fraction for survival DGP (v0.9.4).
#' @param surv_censor_rate Censoring rate for survival DGP (v0.9.4).
#'
#' @return A list with `summary` (one row per scenario x method, with
#'   `conf_strength`, `coverage`, `k`, `mo_confounding`, `phi`, `true_NDE`,
#'   `true_NIE` and NDE/NIE bias/RMSE/Type I columns) and `grid`.
#' @export
#'
#' @examples
#' \dontrun{
#' sens <- gan_mediation_sensitivity(NULL, conf_grid = c(0.3, 0.8),
#'         mo_confounding = 0.8, n_iter = 20)
#' head(sens$summary)
#' }
gan_mediation_sensitivity <- function(trained_gan    = NULL,
                                      conf_grid      = c(0.2, 0.5, 0.8),
                                      coverage_grid  = c(0.3, 0.7, 1),
                                      k_grid         = 1,
                                      mo_confounding = 0.80,
                                      phi            = 0,
                                      rho_G1         = 0,
                                      rho_G2         = 0,
                                      rho_pop        = 0,
                                      separate_U     = FALSE,
                                      omega_1        = NULL,
                                      omega_2        = NULL,
                                      nc_model       = "proxy",
                                      n_iter         = 50,
                                      n_samples      = 500,
                                      n_features     = 20,
                                      beta_Z = 0.10, alpha_M = 0.50, beta_M = 0.30,
                                      base_seed      = 750,
                                      n_cores        = 1,
                                      outcome_type   = c("continuous", "survival"),
                                      effect_scale   = c("loghr", "rmst"),
                                      surv_h0        = 0.1,
                                      surv_event_frac  = 0.6,
                                      surv_censor_rate = NULL) {
  outcome_type <- match.arg(outcome_type)
  effect_scale <- match.arg(effect_scale)

  grid <- expand.grid(conf_strength = conf_grid, coverage = coverage_grid,
                      k = k_grid, KEEP.OUT.ATTRS = FALSE)

  true_NDE <- beta_Z
  true_NIE <- alpha_M * beta_M

  smry <- lapply(seq_len(nrow(grid)), function(gi) {
    cs <- grid$conf_strength[gi]; cov <- grid$coverage[gi]; kk <- grid$k[gi]

    worker <- function(i) {
      dat <- run_single_iteration(
        trained_gan, n_synthetic_samples = n_samples, n_features = n_features,
        n_confounders = kk, beta_Z = beta_Z, alpha_M = alpha_M, beta_M = beta_M,
        conf_strength = cs, coverage = cov, nc_model = nc_model,
        mo_confounding = mo_confounding, phi = phi,
        rho_G1 = rho_G1, rho_G2 = rho_G2, rho_pop = rho_pop,
        separate_U = separate_U, omega_1 = omega_1, omega_2 = omega_2,
        seed = base_seed + gi * 1000L + i,
        outcome_type = outcome_type, surv_h0 = surv_h0,
        surv_event_frac = surv_event_frac, surv_censor_rate = surv_censor_rate)
      if (outcome_type == "survival") {
        res <- .run_surv_methods(dat, effect_scale = effect_scale,
                                 is_mediation = TRUE)
      } else {
        res <- run_mediation_methods(dat, n_features)
      }
      res$iter <- i
      res
    }

    combined <- do.call(rbind, .parallel_lapply(seq_len(n_iter), worker, n_cores))
    s <- summarise_mediation_results(combined, true_NDE, true_NIE)
    s$conf_strength  <- cs
    s$coverage       <- cov
    s$k              <- kk
    s$mo_confounding <- mo_confounding
    s$phi            <- phi
    s$rho_G1         <- rho_G1
    s$rho_G2         <- rho_G2
    s$rho_pop        <- rho_pop
    s$separate_U     <- separate_U
    s$true_NDE       <- true_NDE
    s$true_NIE       <- true_NIE
    s
  })

  summary <- do.call(rbind, smry)
  front   <- c("conf_strength", "coverage", "k", "mo_confounding", "phi",
               "rho_G1", "rho_G2", "rho_pop", "separate_U",
               "true_NDE", "true_NIE", "method")
  summary <- summary[, c(front, setdiff(names(summary), front))]

  # v0.9.2 (JYH #543, #582): attach a scenario manifest so the manuscript
  # can render the truth + parameter ranges as an orientation table.
  manifest <- scenario_manifest(
    list(beta_Z = beta_Z, alpha_M = alpha_M, beta_M = beta_M,
         n_mediators = 1, n = n_samples, n_features = n_features,
         separate_U = separate_U, feat_cor = 0,
         conf_str = NA, mo_confounding = mo_confounding, phi = phi,
         w_signal = NA, rho_G1 = rho_G1, rho_G2 = rho_G2,
         rho_pop = rho_pop),
    conf_grid = conf_grid, coverage_grid = coverage_grid,
    mo_confounding_grid = if (length(unique(mo_confounding)) > 1) mo_confounding else NULL,
    phi_grid = if (length(unique(phi)) > 1) phi else NULL,
    rho_G1_grid = if (length(unique(rho_G1)) > 1) rho_G1 else NULL,
    rho_G2_grid = if (length(unique(rho_G2)) > 1) rho_G2 else NULL,
    rho_pop_grid = if (length(unique(rho_pop)) > 1) rho_pop else NULL
  )

  list(summary = summary, grid = grid, manifest = manifest)
}
