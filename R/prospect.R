# ============================================================
# iconic_prospect: Prospective analysis for users without
# instruments or negative controls.
#
# When the user has X, M, Y but no instruments (G, Gm) or negative
# controls (W), ICONIC cannot apply any of its causal estimators
# beyond UNADJ. iconic_prospect() answers the question:
#
# "If I were to collect an instrument / NC, how much would my
# estimate change, and how strong would the instrument need to be?"
#
# Phase 1 -- Sensitivity surface: sweeps instrument strength (gamma_G)
# and confounding level, showing how estimates converge to the true
# effect as the instrument strengthens.
#
# Phase 2 -- Prospective simulation: for a target instrument strength,
# simulates the full DGP with synthetic instruments and NCs to show
# what the user could expect if they collected such data.
#
# Auto-trains a GAN (or MVN fallback) from the user's
# iconic_data object when no trained_gan is supplied, so the
# synthetic texture matches the user's cohort. Also supports
# confounding = "inferred" to calibrate confounding parameters
# from the data via infer_confounding(). Note: in the prospective
# setting (no instruments/NCs), most confounding parameters cannot
# be inferred and will fall back to defaults with warnings.
# ============================================================

#' Prospective bias-reduction analysis for data without instruments or negative controls
#'
#' For users who have exposure (X), mediator (M), and outcome (Y) but
#' no genetic instruments (G, Gm) or negative controls (W), this
#' function simulates what estimates they could expect if they were to
#' collect such data.
#'
#' Phase 1 sweeps instrument strength (\code{gamma_G}) to show how
#' estimates converge as the instrument strengthens. Phase 2 runs a
#' full prospective simulation at a target strength, generating
#' synthetic instruments and NCs calibrated to the user's sample size
#' and confounding level.
#'
#' @aliases bias_reduction_prospective
#'
#' @section When to use:
#' Use this as a **bias-reduction prospective** when you have an
#' observational exposure-mediator-outcome triplet but lack the genetic
#' instruments and negative controls that the core estimators require.
#' It quantifies the **relative bias improvement** you could expect by
#' collecting such data: the sweep shows how much of the naive
#' confounding bias is removed as instrument strength and NC coverage
#' increase, letting you decide whether the marginal gain justifies the
#' cost of genotyping / profiling the additional assays. It is a
#' planning tool, not an estimator -- it does not produce a causal
#' estimate from your current data, but tells you what a future
#' instrumented study would yield.
#'
#' @section Texture model:
#' When \code{trained_gan} is \code{NULL} and no GAN is attached to
#' \code{data}, a texture model is auto-trained from the user's data.
#'
#' @section Confounding calibration:
#' The \code{confounding} argument controls how the held-fixed
#' confounding parameters are set. In the prospective setting (no
#' instruments or NCs), most parameters cannot be inferred and will
#' fall back to defaults with warnings -- this is an honest limitation,
#' not a silent failure.
#'
#' @param data An \code{iconic_data} object (must have X and Y;
#' M is required for mediation prospect).
#' @param trained_gan Optional \code{iconic_gan} from
#' \code{\link{train_gan_on_real_data}()}. If \code{NULL}, a texture
#' model is auto-trained from \code{data}.
#' @param confounding Confounding parameter source: \code{"default"},
#' \code{"inferred"}, or \code{"manual"}. Default \code{"default"}.
#' @param gan_epochs Epochs for auto-trained GAN. Default 100.
#' @param gamma_G_grid Instrument strength values to sweep (Phase 1).
#' Default \code{c(0.2, 0.4, 0.6, 0.8, 1.0)}.
#' @param target_gamma_G Target instrument strength for Phase 2.
#' Default 0.6 (matching the DGP default).
#' @param n_iter Replicates per grid cell. Default 30.
#' @param n_features Features per replicate. Default 10.
#' @param mo_confounding Assumed M-O confounding strength. Default 0.8.
#' @param phi Assumed mediator instrument strength. Default 0.8.
#' @param lambda_XM Optional per-path confounder loading vector (X->M path).
#' @param lambda_MY Optional per-path confounder loading vector (M->Y path).
#' Default TRUE.
#' @param omega_1,omega_2 Assumed NC coverage. Default 0.7.
#' @param rho_G1_grid,rho_G2_grid Instrument-exogeneity violation values
#' (correlation of each instrument with its path's confounder composite) to
#' sweep in Phase 3 at the target instrument strength. Default
#' \code{c(0, 0.1, 0.2, 0.3, 0.5)}, matching \code{\link{iconic_sensitivity}}.
#' @param omega_grid_rho Negative-control coverage values swept jointly with
#' the rho grid in Phase 3, on the diagonal (\code{omega_1 == omega_2}).
#' Default \code{c(0.3, 0.7, 1.0)}, matching \code{\link{iconic_sensitivity}}.
#' This is independent of the Phase 1 \code{omega_1}/\code{omega_2} sweep.
#' @param run_rho_sweep Logical: run the Phase 3 robustness sweep
#' (default \code{TRUE}). The sweep crosses instrument-exogeneity violations
#' (rho) with negative-control coverage (omega) and feeds the resulting
#' degradation surface into \code{\link{iconic_recommend}()} so the
#' recommended estimator is chosen by robustness to both imperfect
#' instruments and weakening controls, rather than by eligibility alone.
#' Set \code{FALSE} to skip (faster, but the recommendation then falls back
#' to a single-point / eligibility ranking).
#' @param bias_threshold Tipping-point threshold. Default 0.10.
#' @param base_seed Base RNG seed. Default 500.
#' @param verbose Logical: print progress messages during the sweep.
#' Default \code{FALSE} (quiet).
#' @param allow_no_proxy Logical: when \code{TRUE}
#' (default), proceed with the prospective sweep even if the data
#' already has instruments/NCs (with a message). When \code{FALSE},
#' error if the data already has IV+NC (use iconic_estimate instead).
#' @param outcome_type \code{"continuous"} (default) or \code{"survival"}
#'   Threads through to the simulation DGP.
#' @param effect_scale \code{"loghr"} (default) or \code{"rmst"}. Only
#' used when \code{outcome_type = "survival"}.
#' @param surv_h0,surv_event_frac,surv_censor_rate Survival DGP parameters
#'   See \code{\link{generate_toy_data}}.
#'
#' @section Defaults:
#' \tabular{lll}{
#' \strong{Parameter} \tab \strong{Default} \tab \strong{Source} \cr
#' \code{confounding} \tab "default" \tab Use DGP defaults below \cr
#' \code{gan_epochs} \tab 100 \tab Texture-model training budget \cr
#' \code{gamma_G_grid} \tab c(0.2,0.4,0.6,0.8,1.0) \tab Instrument-strength sweep \cr
#' \code{target_gamma_G} \tab 0.6 \tab DGP default (gamma_G) \cr
#' \code{n_iter} \tab 30 \tab Replicates per grid cell \cr
#' \code{n_features} \tab 10 \tab Features per replicate \cr
#' \code{mo_confounding} \tab 0.8 \tab Simulation calibration (delta_mo) \cr
#' \code{phi} \tab 0.8 \tab Strong mediator instrument assumption \cr
#' \code{lambda_XM}, \code{lambda_MY} \tab shared \tab Per-path confounder loadings \cr
#' \code{omega_1, omega_2} \tab 0.7 \tab NC coverage (simulation calibration) \cr
#' \code{rho_G1_grid}, \code{rho_G2_grid} \tab c(0,0.1,0.2,0.3,0.5) \tab Phase 3 exogeneity sweep \cr
#' \code{omega_grid_rho} \tab c(0.3,0.7,1.0) \tab Phase 3 NC-coverage sweep (diagonal) \cr
#' \code{run_rho_sweep} \tab TRUE \tab Run Phase 3 robustness sweep \cr
#' \code{bias_threshold} \tab 0.10 \tab Tipping-point threshold \cr
#' \code{allow_no_proxy} \tab TRUE \tab Proceed in prospective setting \cr
#' }
#'
#' @return An \code{iconic_prospect} S3 object: a named list with
#' \code{$strength_surface} (Phase 1: gamma_G x method estimates),
#' \code{$prospective} (Phase 2: full simulation at target strength),
#' \code{$rho_surface} (Phase 3: rho_G1 x rho_G2 exogeneity-robustness
#' surface at the target strength; \code{NULL} when
#' \code{run_rho_sweep = FALSE}),
#' \code{$summary}, \code{$recommendation},
#' \code{$texture_source}, and \code{$inferred_confounding} (when
#' \code{confounding = "inferred"}).
#' @export
#'
#' @examples
#' if (check_torch_setup()) {
#'   data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
#'     M = rnorm(100))
#'   result <- iconic_prospect(data, n_iter = 2, gan_epochs = 5,
#'     gamma_G_grid = c(0.4, 0.8), run_rho_sweep = FALSE)
#'   print(result)
#' }
iconic_prospect <- function(data,
                            trained_gan = NULL,
                            confounding = c("default", "inferred", "manual"),
                            gan_epochs = 100,
                            gamma_G_grid = c(0.2, 0.4, 0.6, 0.8, 1.0),
                            target_gamma_G = 0.6,
                            n_iter = 30,
                            n_features = 10,
                            mo_confounding = 0.8,
                            phi = 0.8,
                            lambda_XM = NULL,
                            lambda_MY = NULL,
                            omega_1 = 0.7, omega_2 = 0.7,
                            rho_G1_grid = c(0, 0.1, 0.2, 0.3, 0.5),
                            rho_G2_grid = c(0, 0.1, 0.2, 0.3, 0.5),
                            omega_grid_rho = c(0.3, 0.7, 1.0),
                            run_rho_sweep = TRUE,
                            bias_threshold = 0.10,
                            base_seed = 500,
                            verbose = FALSE,
                            allow_no_proxy = TRUE,
                            outcome_type = c("continuous", "survival"),
                            effect_scale = c("loghr", "rmst"),
                            surv_h0 = 0.1,
                            surv_event_frac = 0.6,
                            surv_censor_rate = NULL) {
  if (!inherits(data, "iconic_data"))
    stop("data must be an iconic_data object from iconic_data().")
  if (!data$is_mediation)
    stop("iconic_prospect requires mediation data (supply M). ",
         "This function answers: 'if I collected instruments/NCs, ",
         "how would my NDE/NIE estimates change?'")

  outcome_type <- match.arg(outcome_type)
  effect_scale <- match.arg(effect_scale)
  # inherit outcome_type from data if not explicitly set
  if (outcome_type == "continuous" && data$outcome_type == "survival")
    outcome_type <- "survival"

  # Capture which calibration arguments the user left at their defaults
  # BEFORE any reassignment: inferred values may fill in defaults, but
  # user-supplied vectors (e.g. an omega sweep) always take precedence.
  omega_1_default <- missing(omega_1)
  omega_2_default <- missing(omega_2)
  mo_confounding_default <- missing(mo_confounding)

  # explicit control over proceed-without-proxy behavior.
  # iconic_prospect is designed for the no-IV/no-NC case, so the default
  # (TRUE) is to proceed. allow_no_proxy = FALSE makes the user acknowledge
  # the prospective setting explicitly.
  has_iv <- !is.null(data$G) || !is.null(data$Gm)
  has_nc <- !is.null(data$W)
  if (has_iv && has_nc && !allow_no_proxy) {
    stop("Data already has instruments and negative controls. ",
         "iconic_prospect is designed for the prospective (no-IV/no-NC) ",
         "setting. Use iconic_estimate() / iconic_sensitivity() instead, ",
         "or set allow_no_proxy = TRUE to run the prospective sweep anyway.")
  }
  if (!has_iv && !has_nc && isTRUE(verbose)) {
    message("No instruments or negative controls supplied: running the ",
            "prospective sweep (simulating what estimates you could expect ",
            "if you collected such data). Set allow_no_proxy = FALSE to ",
            "silence this note.")
  }

  confounding <- match.arg(confounding)

  # -- Resolve the texture model --
  gan_res <- .resolve_gan(trained_gan, data, epochs = gan_epochs)
  gan <- gan_res$gan
  texture_source <- gan_res$source

  # -- Resolve confounding parameters --
  inferred_conf <- NULL
  if (confounding == "inferred") {
    inferred_conf <- infer_confounding(data, diagnosis = NULL,
                                       estimate = NULL)
    # Use inferred values where available, but ONLY to fill in arguments the
    # user left at their defaults; explicit user-supplied values (e.g. an
    # omega sweep) always take precedence over the inferred scalars.
    if (mo_confounding_default && inferred_conf$mo_confounding$available)
      mo_confounding <- inferred_conf$mo_confounding$estimate
    if (omega_1_default && inferred_conf$omega_1$available)
      omega_1 <- inferred_conf$omega_1$estimate
    if (omega_2_default && inferred_conf$omega_2$available)
      omega_2 <- inferred_conf$omega_2$estimate
  }

  n <- data$n

  # === Phase 1: Instrument-strength surface ===
  # Sweep gamma_G to show how estimates change as the instrument strengthens.
  # At each gamma_G, generate full DGP with synthetic G, Gm, W and run all
  # estimators. The key comparison: UNADJ (no instrument) vs IV2SLS2/PGC2Gm
  # (with instrument).

  true_NDE <- 0.10
  true_NIE <- 0.15

  # ── Build the omega sweep ──
  # omega_1/omega_2 may be vectors; Phase 1 sweeps gamma_G x omega_1 x
  # omega_2. When confounding = "inferred" supplied a scalar omega, the
  # sweep reduces to that single value.
  omega_1_grid <- sort(unique(omega_1))
  omega_2_grid <- sort(unique(omega_2))
  omega_swept <- length(omega_1_grid) > 1 || length(omega_2_grid) > 1

  # Phase 1 uses the reference (first) omega cell for the prospective
  # Phase 2 simulation; the full omega grid is swept on the strength
  # surface only.
  omega_1_ref <- omega_1_grid[1]
  omega_2_ref <- omega_2_grid[1]

  grid <- expand.grid(gamma_G = gamma_G_grid,
                      omega_1 = omega_1_grid, omega_2 = omega_2_grid,
                      KEEP.OUT.ATTRS = FALSE)
  n_grid <- nrow(grid)
  strength_rows <- lapply(seq_len(n_grid), function(gi) {
    gg <- grid$gamma_G[gi]
    o1 <- grid$omega_1[gi]
    o2 <- grid$omega_2[gi]
    if (n_grid > 1 && isTRUE(verbose)) {
      omega_txt <- if (omega_swept) paste0(", omega_1=", o1, ", omega_2=", o2) else ""
      message("Phase 1: gamma_G=", gg, omega_txt,
              " (", gi, "/", n_grid, ")")
    }
    worker <- function(i) {
      dat <- run_single_iteration(
        trained_gan = gan,
        n_synthetic_samples = n, n_features = n_features,
        mo_confounding = mo_confounding, phi = phi, gamma_G = gg,
        lambda_XM = lambda_XM, lambda_MY = lambda_MY, omega_1 = o1, omega_2 = o2,
        seed = base_seed + as.integer(gg * 10000) + gi * 1000L + i,
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
    combined <- do.call(rbind, lapply(seq_len(n_iter), worker))
    s <- summarise_mediation_results(combined, true_NDE, true_NIE)
    s$gamma_G <- gg
    s$omega_1 <- o1
    s$omega_2 <- o2
    s
  })

  strength_surface <- do.call(rbind, strength_rows)
  # Reorder columns
  front <- c("gamma_G", "omega_1", "omega_2", "method")
  strength_surface <- strength_surface[, c(front, setdiff(names(strength_surface), front))]

  # === Phase 2: Prospective simulation at target strength ===
  # Full simulation at the target instrument strength, showing what the
  # user could expect if they collected an instrument of that strength.

  prospect_worker <- function(i) {
    dat <- run_single_iteration(
      trained_gan = gan,
      n_synthetic_samples = n, n_features = n_features,
      mo_confounding = mo_confounding, phi = phi, gamma_G = target_gamma_G,
      lambda_XM = lambda_XM, lambda_MY = lambda_MY, omega_1 = omega_1_ref, omega_2 = omega_2_ref,
      seed = base_seed + 99999 + i,
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
  prospect_combined <- do.call(rbind, lapply(seq_len(n_iter),
                                             prospect_worker))
  prospective <- summarise_mediation_results(prospect_combined, true_NDE, true_NIE)

  # First-stage F at target strength (from one representative dataset)
  rep_dat <- run_single_iteration(
    trained_gan = gan,
    n_synthetic_samples = n, n_features = n_features,
    mo_confounding = mo_confounding, phi = phi, gamma_G = target_gamma_G,
    lambda_XM = lambda_XM, lambda_MY = lambda_MY, omega_1 = omega_1_ref, omega_2 = omega_2_ref,
    seed = base_seed + 77777,
    outcome_type = outcome_type, surv_h0 = surv_h0,
    surv_event_frac = surv_event_frac, surv_censor_rate = surv_censor_rate)
  if (outcome_type == "survival") {
    rep_idata <- iconic_data(
      X = rep_dat$X, outcome_type = "survival",
      surv_time = rep_dat$surv_time, surv_event = rep_dat$surv_event,
      M = rep_dat$M, W = t(rep_dat$W), W1 = t(rep_dat$W1), W2 = t(rep_dat$W2),
      G = rep_dat$G[, 1], Gm = rep_dat$Gm,
      covariates = rep_dat$synthetic_data)
  } else {
    rep_idata <- iconic_data(
      X = rep_dat$X, Y = t(rep_dat$Y), M = rep_dat$M,
      W = t(rep_dat$W), W1 = t(rep_dat$W1), W2 = t(rep_dat$W2),
      G = rep_dat$G[, 1], Gm = rep_dat$Gm,
      covariates = rep_dat$synthetic_data)
  }
  rep_diag <- iconic_diagnose(rep_idata)

  # === Phase 3: Exogeneity-robustness sweep at target strength ===
  # Phases 1-2 hold the instruments perfectly exogenous (rho_G1 = rho_G2 = 0).
  # But the planning question "which estimator should I collect instruments
  # for?" depends on robustness to *imperfect* instruments, since any real
  # instrument will be somewhat pleiotropic / confounded. Phase 3 sweeps the
  # instrument-exogeneity violations (rho_G1 x rho_G2) jointly with
  # negative-control coverage (omega) at the target strength, and feeds the
  # resulting degradation surface into iconic_recommend() so the
  # recommendation is robustness-based (per-estimand max|bias| + coverage
  # distance) rather than a single-point or eligibility-only ranking. This is
  # what prevents a high-bias estimator such as UNADJ from being "recommended"
  # merely because it is eligible. The omega sweep matches iconic_sensitivity:
  # when omega_1/omega_2 are identical vectors, coverage is swept on the
  # diagonal (omega_1 == omega_2); distinct vectors cross the full grid.
  rho_surface <- NULL
  if (isTRUE(run_rho_sweep)) {
    # Phase 3 sweeps NC coverage on the diagonal (omega_1 == omega_2) over
    # omega_grid_rho, matching iconic_sensitivity's default and the
    # manuscript's degradation-surface design. This is independent of the
    # Phase 1 omega_1/omega_2 sweep (which maps bias across coverage on the
    # strength surface); Phase 3 varies coverage jointly with exogeneity.
    omega_rho_grid <- sort(unique(omega_grid_rho))
    rho_grid <- expand.grid(rho_G1 = rho_G1_grid, rho_G2 = rho_G2_grid,
                            omega_1 = omega_rho_grid, omega_2 = omega_rho_grid,
                            KEEP.OUT.ATTRS = FALSE)
    rho_grid <- rho_grid[rho_grid$omega_1 == rho_grid$omega_2, , drop = FALSE]
    n_rho <- nrow(rho_grid)
    rho_rows <- lapply(seq_len(n_rho), function(gi) {
      r1 <- rho_grid$rho_G1[gi]
      r2 <- rho_grid$rho_G2[gi]
      o1 <- rho_grid$omega_1[gi]
      o2 <- rho_grid$omega_2[gi]
      if (n_rho > 1 && isTRUE(verbose))
        message("Phase 3: rho_G1=", r1, ", rho_G2=", r2,
                ", omega_1=", o1, ", omega_2=", o2, " (", gi, "/", n_rho, ")")
      worker <- function(i) {
        dat <- run_single_iteration(
          trained_gan = gan,
          n_synthetic_samples = n, n_features = n_features,
          mo_confounding = mo_confounding, phi = phi, gamma_G = target_gamma_G,
          rho_G1 = r1, rho_G2 = r2, rho_pop = 0,
          lambda_XM = lambda_XM, lambda_MY = lambda_MY,
          omega_1 = o1, omega_2 = o2,
          seed = base_seed + 55555 + gi * 1000L + i,
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
      combined <- do.call(rbind, lapply(seq_len(n_iter), worker))
      s <- summarise_mediation_results(combined, true_NDE, true_NIE)
      s$rho_G1 <- r1
      s$rho_G2 <- r2
      s$omega_1 <- o1
      s$omega_2 <- o2
      s
    })
    rho_surface <- do.call(rbind, rho_rows)
    # Tipping-point annotation, matching iconic_sensitivity's schema.
    rho_surface$tipped_NDE <- abs(rho_surface$NDE_bias) > bias_threshold
    rho_surface$tipped_NIE <- abs(rho_surface$NIE_bias) > bias_threshold
    rho_surface$tipped <- rho_surface$tipped_NDE | rho_surface$tipped_NIE
    front <- c("rho_G1", "rho_G2", "omega_1", "omega_2", "method")
    rho_surface <- rho_surface[, c(front, setdiff(names(rho_surface), front))]
  }

  # === Summary ===
  summary_txt <- .build_prospect_summary(strength_surface, prospective,
                                         target_gamma_G, rep_diag, n,
                                         mo_confounding, phi,
                                         rho_surface = rho_surface,
                                         bias_threshold = bias_threshold)

  # Recommendation: which estimator would be best at the target strength?
  # When the Phase 3 rho sweep ran, pass its surface as `sensitivity` so the
  # ranking is robustness-based (per-estimand max|bias| + coverage distance)
  # rather than eligibility-only. .extract_robustness() reads $surface.
  rec <- iconic_recommend(rep_idata, diagnosis = rep_diag,
                          sensitivity = if (!is.null(rho_surface))
                            list(surface = rho_surface) else NULL,
                          verbose = verbose)

  # Conditional recommendation: the best estimator under each collection
  # scenario (which instruments / negative controls the user might collect).
  # An estimator is only recommendable for a scenario if that scenario
  # supplies the data it requires -- e.g. COCA uses only negative controls,
  # so it never appears under an instrument-collection scenario.
  rec_by_scenario <- .recommend_by_scenario(rec)

  obj <- list(
    strength_surface = strength_surface,
    prospective = prospective,
    rho_surface = rho_surface,
    recommendation_by_scenario = rec_by_scenario,
    target_gamma_G = target_gamma_G,
    instrument_F = rep_diag$instrument_strength,
    recommendation = rec,
    n_iter = n_iter,
    n_samples = n,
    mo_confounding = mo_confounding,
    phi = phi,
    omega_1 = omega_1_grid,
    omega_2 = omega_2_grid,
    omega_swept = omega_swept,
    rho_G1_grid = if (isTRUE(run_rho_sweep)) rho_G1_grid else NULL,
    rho_G2_grid = if (isTRUE(run_rho_sweep)) rho_G2_grid else NULL,
    texture_source = texture_source,
    inferred_confounding = inferred_conf,
    summary = summary_txt
  )
  class(obj) <- c("iconic_prospect", "list")
  if (isTRUE(verbose))
    message("iconic_prospect complete. Call summary() or print() on the result for the full prospective summary.")
  obj
}


#' Best estimator under each collection scenario (internal)
#'
#' Maps each instrument/NC collection scenario to the estimators it makes
#' available, then picks the most robust eligible one from the
#' \code{iconic_recommend} ranking (robustness = per-estimand max|bias| over
#' the Phase 3 rho x omega surface). An estimator appears only under scenarios
#' that supply its required data:
#' \itemize{
#'   \item UNADJ: none (always available)
#'   \item COCA: W (negative controls only; no instrument)
#'   \item DIRECT, IV2SLS, PGC: G + W
#'   \item IV2SLS2: G + Gm (optional path-specific W1/W2 augmentation)
#'   \item PGC2: G + W1 + W2
#'   \item PGC2Gm: G + Gm + W1 + W2
#' }
#' @param rec An \code{iconic_recommendation} object.
#' @return A data frame with columns \code{scenario}, \code{estimator}, and
#'   \code{robustness_NDE}.
#' @keywords internal
.recommend_by_scenario <- function(rec) {
  ranking <- rec$ranking
  if (is.null(ranking) || !"estimator" %in% names(ranking)) return(NULL)

  # Robustness score used to order within a scenario (higher = better).
  # Fall back to eligibility order when no sensitivity surface was used.
  has_rob <- "robustness_NDE" %in% names(ranking)
  score <- if (has_rob) ranking$robustness_NDE else rep(NA_real_, nrow(ranking))

  # Scenario -> estimators it makes available. The IV estimators (IV2SLS,
  # IV2SLS2) are identified by the instrument(s) alone (classic 2SLS); NC
  # augmentation is optional, not a requirement. IV2SLS uses a single pooled
  # W; IV2SLS2 uses optional path-specific panels (W1 in stage 1, W2 in
  # stages 2-3) and runs plain 2-stage MR when they are absent or identical
  # (so under a pooled-only "G1 + Gm + W" scenario IV2SLS2 is unaugmented).
  # The proximal bridge estimators (PGC, PGC2, PGC2Gm) and COCA/DIRECT
  # require W.
  scenarios <- list(
    "G1 only"                    = c("UNADJ", "IV2SLS"),
    "Gm only"                    = c("UNADJ"),
    "G1 + Gm"                    = c("UNADJ", "IV2SLS", "IV2SLS2"),
    "W only"                     = c("UNADJ", "COCA"),
    "W1 + W2"                    = c("UNADJ", "COCA"),
    "G1 + W"                     = c("UNADJ", "COCA", "DIRECT", "IV2SLS", "PGC"),
    "G1 + Gm + W"                = c("UNADJ", "COCA", "DIRECT", "IV2SLS", "PGC", "IV2SLS2"),
    "G1 + W1 + W2"               = c("UNADJ", "COCA", "DIRECT", "IV2SLS", "PGC", "PGC2"),
    "G1 + Gm + W1 + W2 (full)"   = c("UNADJ", "COCA", "DIRECT", "IV2SLS", "PGC", "IV2SLS2", "PGC2", "PGC2Gm")
  )

  rows <- lapply(names(scenarios), function(sc) {
    avail <- scenarios[[sc]]
    sub <- ranking[ranking$estimator %in% avail, , drop = FALSE]
    sub <- sub[!is.na(sub$eligible) & sub$eligible, , drop = FALSE]
    if (nrow(sub) == 0) {
      return(data.frame(scenario = sc, estimator = NA_character_,
                        robustness_NDE = NA_real_, stringsAsFactors = FALSE))
    }
    if (has_rob) {
      ord <- order(-sub$robustness_NDE)
      sub <- sub[ord, , drop = FALSE]
    }
    data.frame(scenario = sc,
               estimator = sub$estimator[1],
               robustness_NDE = if (has_rob) sub$robustness_NDE[1] else NA_real_,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Build prospect summary (internal)
#' @keywords internal
#' @noRd
.build_prospect_summary <- function(strength_surface, prospective,
                                    target_gamma_G, rep_diag, n,
                                    mo_confounding, phi,
                                    rho_surface = NULL,
                                    bias_threshold = 0.10) {
  lines <- character(0)
  lines <- c(lines,
    sprintf(" Sample size: %d (calibrated to user data)", n),
    sprintf(" Assumed M-O confounding: %.1f, mediator instrument phi: %.1f",
            mo_confounding, phi),
    sprintf(" Target instrument strength (gamma_G): %.1f", target_gamma_G))

  # Instrument strength at target
  F_G_val <- rep_diag$instrument_strength$F_G
  F_Gm_val <- rep_diag$instrument_strength$F_Gm
  if (length(F_G_val) > 1) F_G_val <- rep_diag$instrument_strength$F_G_median
  if (length(F_Gm_val) > 1) F_Gm_val <- rep_diag$instrument_strength$F_Gm_median
  if (!is.na(F_G_val))
    lines <- c(lines, sprintf(" First-stage F (G) at target: %.1f", F_G_val))
  if (!is.na(F_Gm_val))
    lines <- c(lines, sprintf(" First-stage F (Gm) at target: %.1f", F_Gm_val))

  # Phase 1: how UNADJ vs best estimator change with instrument strength
  lines <- c(lines, "", " Phase 1 -- Instrument-strength surface (NDE bias):")
  gammas <- sort(unique(strength_surface$gamma_G))
  for (gg in gammas) {
    sub <- strength_surface[strength_surface$gamma_G == gg, ]
    unadj_bias <- sub$NDE_bias[sub$method == "UNADJ"]
    iv2sls2_bias <- if ("IV2SLS2" %in% sub$method) sub$NDE_bias[sub$method == "IV2SLS2"] else NA
    pgc2gm_bias <- if ("PGC2Gm" %in% sub$method) sub$NDE_bias[sub$method == "PGC2Gm"] else NA
    lines <- c(lines, sprintf(" gamma_G=%.1f: UNADJ bias=%+.3f, IV2SLS2 bias=%+.3f, PGC2Gm bias=%+.3f",
                              gg, unadj_bias,
                              ifelse(is.na(iv2sls2_bias), NA, iv2sls2_bias),
                              ifelse(is.na(pgc2gm_bias), NA, pgc2gm_bias)))
  }

  # Phase 2: prospective estimates at target
  lines <- c(lines, "", " Phase 2 -- Prospective estimates at target strength:")
  for (m in c("UNADJ", "IV2SLS", "IV2SLS2", "PGC2Gm")) {
    if (m %in% prospective$method) {
      r <- prospective[prospective$method == m, ]
      lines <- c(lines, sprintf(" %s: NDE=%.3f (bias=%+.3f), NIE=%.3f (bias=%+.3f)",
                                m, r$NDE_mean, r$NDE_bias, r$NIE_mean, r$NIE_bias))
    }
  }

  # Phase 3: exogeneity-robustness sweep over the joint rho x omega grid.
  # Report, for each leading estimator, the worst |NDE bias| across the whole
  # surface (and the omega cell where it occurs), plus a per-omega breakdown
  # of max bias and the earliest tipping rho.
  if (!is.null(rho_surface)) {
    has_omega <- all(c("omega_1", "omega_2") %in% names(rho_surface))
    hdr <- if (has_omega)
      " Phase 3 -- Robustness to instrument violation (rho_G1 x rho_G2 x omega):"
    else
      " Phase 3 -- Robustness to instrument violation (rho_G1 x rho_G2):"
    lines <- c(lines, "", hdr)
    for (m in c("UNADJ", "IV2SLS", "IV2SLS2", "PGC2", "PGC2Gm")) {
      if (m %in% rho_surface$method) {
        sub <- rho_surface[rho_surface$method == m, ]
        abs_b <- abs(sub$NDE_bias)
        max_b <- max(abs_b, na.rm = TRUE)
        # omega cell at which the worst bias occurs
        worst_txt <- ""
        if (has_omega) {
          iw <- which.max(abs_b)
          worst_txt <- sprintf(" (worst at omega_1=omega_2=%.1f)", sub$omega_1[iw])
        }
        lines <- c(lines, sprintf(" %s: max|NDE bias|=%.3f across grid%s",
                                  m, max_b, worst_txt))
        # per-omega breakdown: max bias and earliest tipping rho within each cell
        if (has_omega) {
          for (ov in sort(unique(sub$omega_1))) {
            so <- sub[sub$omega_1 == ov & sub$omega_2 == ov, ]
            if (!nrow(so)) next
            mb_o <- max(abs(so$NDE_bias), na.rm = TRUE)
            so0 <- so[so$rho_G1 == 0, ]
            tip <- so0$rho_G2[abs(so0$NDE_bias) > bias_threshold]
            tip_txt <- if (length(tip) > 0) sprintf("tips at rho_G2=%.1f", min(tip))
                       else "no tipping"
            lines <- c(lines, sprintf("    omega=%.1f: max|NDE bias|=%.3f (%s)",
                                      ov, mb_o, tip_txt))
          }
        } else {
          sub0 <- sub[sub$rho_G1 == 0, ]
          tip <- sub0$rho_G2[abs(sub0$NDE_bias) > bias_threshold]
          tip_txt <- if (length(tip) > 0) sprintf("tips at rho_G2=%.1f", min(tip))
                     else "no tipping"
          lines <- c(lines, sprintf("    (%s)", tip_txt))
        }
      }
    }
  }

  paste(lines, collapse = "\n")
}


#' Print method for iconic_prospect objects
#'
#' @param x An \code{iconic_prospect} object.
#' @param ... Unused.
#' @return Invisibly returns `x` (the `iconic_prospect` object); called for its side effect of printing a human-readable summary.
#' @export
#' @examples
#' if (check_torch_setup()) {
#'   data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
#'     M = rnorm(100))
#'   result <- iconic_prospect(data, n_iter = 2, gan_epochs = 5,
#'     gamma_G_grid = c(0.4, 0.8), run_rho_sweep = FALSE)
#'   print(result)
#' }
print.iconic_prospect <- function(x, ...) {
  cat("<iconic_prospect>\n")
  cat(x$summary, "\n")
  if (!is.null(x$texture_source))
    cat(" Texture:", x$texture_source, "\n")
  if (!is.null(x$inferred_confounding))
    cat(" Confounding: inferred from data\n")
  if (!is.null(x$recommendation_by_scenario)) {
    cat("\n Recommended estimator by collection scenario:\n")
    rbs <- x$recommendation_by_scenario
    for (i in seq_len(nrow(rbs))) {
      est <- rbs$estimator[i]
      est_txt <- if (is.na(est)) "none eligible" else est
      cat(sprintf("   if %-26s -> %s\n", rbs$scenario[i], est_txt))
    }
  } else if (!is.null(x$recommendation)) {
    cat("\n Recommended estimator:", x$recommendation$recommended, "\n")
  }
  invisible(x)
}

#' Summary method for iconic_prospect objects
#'
#' Prints the full prospective summary (same as \code{print()}).
#' @param object An \code{iconic_prospect} object.
#' @param ... Unused.
#' @return Invisibly returns \code{object}.
#' @export
#' @examples
#' if (check_torch_setup()) {
#'   data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
#'     M = rnorm(100))
#'   result <- iconic_prospect(data, n_iter = 2, gan_epochs = 5,
#'     gamma_G_grid = c(0.4, 0.8), run_rho_sweep = FALSE)
#'   summary(result)
#' }
summary.iconic_prospect <- function(object, ...) {
  print.iconic_prospect(object, ...)
  invisible(object)
}
