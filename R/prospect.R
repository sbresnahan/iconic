# ============================================================
# iconic_prospect: Prospective analysis for users without
# instruments or negative controls.
#
# When the user has Z, M, Y but no instruments (G, Gm) or negative
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
#' For users who have exposure (Z), mediator (M), and outcome (Y) but
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
#' @param data An \code{iconic_data} object (must have Z and Y;
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
#' @param separate_U Use separate confounders for Z->M and M->Y.
#' Default TRUE.
#' @param omega_1,omega_2 Assumed NC coverage. Default 0.7.
#' @param bias_threshold Tipping-point threshold. Default 0.10.
#' @param base_seed Base RNG seed. Default 500.
#' @param n_cores Number of parallel workers for simulation replicates.
#' Default 1 (sequential). Uses \code{parallel::mclapply} on Unix
#' and a PSOCK cluster on Windows.
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
#' \code{separate_U} \tab TRUE \tab Path-specific confounders \cr
#' \code{omega_1, omega_2} \tab 0.7 \tab NC coverage (simulation calibration) \cr
#' \code{bias_threshold} \tab 0.10 \tab Tipping-point threshold \cr
#' \code{allow_no_proxy} \tab TRUE \tab Proceed in prospective setting \cr
#' }
#'
#' @return An \code{iconic_prospect} S3 object: a named list with
#' \code{$strength_surface} (Phase 1: gamma_G x method estimates),
#' \code{$prospective} (Phase 2: full simulation at target strength),
#' \code{$summary}, \code{$recommendation},
#' \code{$texture_source}, and \code{$inferred_confounding} (when
#' \code{confounding = "inferred"}).
#' @export
#'
#' @examples
#' \dontrun{
#' data <- iconic_data(Z = rnorm(100), Y = matrix(rnorm(100*10), 10, 100),
#' M = rnorm(100))
#' result <- iconic_prospect(data, n_iter = 10, n_cores = 2)
#' print(result)
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
                            separate_U = TRUE,
                            omega_1 = 0.7, omega_2 = 0.7,
                            bias_threshold = 0.10,
                            base_seed = 500,
                            n_cores = 1,
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
  if (!has_iv && !has_nc) {
    message("No instruments or negative controls supplied: running the ",
            "prospective sweep (simulating what estimates you could expect ",
            "if you collected such data). Set allow_no_proxy = FALSE to ",
            "suppress this message.")
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
    if (inferred_conf$mo_confounding$available)
      mo_confounding <- inferred_conf$mo_confounding$estimate
    if (inferred_conf$omega_1$available)
      omega_1 <- inferred_conf$omega_1$estimate
    if (inferred_conf$omega_2$available)
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

  n_gamma <- length(gamma_G_grid)
  strength_rows <- lapply(seq_len(n_gamma), function(gi) {
    gg <- gamma_G_grid[gi]
    if (n_gamma > 1)
      message("Phase 1: gamma_G=", gg, " (", gi, "/", n_gamma, ")")
    worker <- function(i) {
      dat <- run_single_iteration(
        trained_gan = gan,
        n_synthetic_samples = n, n_features = n_features,
        mo_confounding = mo_confounding, phi = phi, gamma_G = gg,
        separate_U = separate_U, omega_1 = omega_1, omega_2 = omega_2,
        seed = base_seed + as.integer(gg * 10000) + i,
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
    combined <- do.call(rbind, .parallel_lapply(seq_len(n_iter), worker,
                                                n_cores = n_cores,
                                                progress = " Replicates"))
    s <- summarise_mediation_results(combined, true_NDE, true_NIE)
    s$gamma_G <- gg
    s
  })

  strength_surface <- do.call(rbind, strength_rows)

  # === Phase 2: Prospective simulation at target strength ===
  # Full simulation at the target instrument strength, showing what the
  # user could expect if they collected an instrument of that strength.

  prospect_worker <- function(i) {
    dat <- run_single_iteration(
      trained_gan = gan,
      n_synthetic_samples = n, n_features = n_features,
      mo_confounding = mo_confounding, phi = phi, gamma_G = target_gamma_G,
      separate_U = separate_U, omega_1 = omega_1, omega_2 = omega_2,
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
  prospect_combined <- do.call(rbind, .parallel_lapply(seq_len(n_iter),
                                                       prospect_worker,
                                                       n_cores = n_cores,
                                                       progress = "Phase 2 replicates"))
  prospective <- summarise_mediation_results(prospect_combined, true_NDE, true_NIE)

  # First-stage F at target strength (from one representative dataset)
  rep_dat <- run_single_iteration(
    trained_gan = gan,
    n_synthetic_samples = n, n_features = n_features,
    mo_confounding = mo_confounding, phi = phi, gamma_G = target_gamma_G,
    separate_U = separate_U, omega_1 = omega_1, omega_2 = omega_2,
    seed = base_seed + 77777,
    outcome_type = outcome_type, surv_h0 = surv_h0,
    surv_event_frac = surv_event_frac, surv_censor_rate = surv_censor_rate)
  if (outcome_type == "survival") {
    rep_idata <- iconic_data(
      Z = rep_dat$Z, outcome_type = "survival",
      surv_time = rep_dat$surv_time, surv_event = rep_dat$surv_event,
      M = rep_dat$M, W = t(rep_dat$W), W1 = t(rep_dat$W1), W2 = t(rep_dat$W2),
      G = rep_dat$G[, 1], Gm = rep_dat$Gm,
      covariates = rep_dat$synthetic_data)
  } else {
    rep_idata <- iconic_data(
      Z = rep_dat$Z, Y = t(rep_dat$Y), M = rep_dat$M,
      W = t(rep_dat$W), W1 = t(rep_dat$W1), W2 = t(rep_dat$W2),
      G = rep_dat$G[, 1], Gm = rep_dat$Gm,
      covariates = rep_dat$synthetic_data)
  }
  rep_diag <- iconic_diagnose(rep_idata)

  # === Summary ===
  summary_txt <- .build_prospect_summary(strength_surface, prospective,
                                         target_gamma_G, rep_diag, n,
                                         mo_confounding, phi)

  # Recommendation: which estimator would be best at the target strength?
  rec <- iconic_recommend(rep_idata, diagnosis = rep_diag)

  obj <- list(
    strength_surface = strength_surface,
    prospective = prospective,
    target_gamma_G = target_gamma_G,
    instrument_F = rep_diag$instrument_strength,
    recommendation = rec,
    n_iter = n_iter,
    n_samples = n,
    mo_confounding = mo_confounding,
    phi = phi,
    omega_1 = omega_1,
    omega_2 = omega_2,
    texture_source = texture_source,
    inferred_confounding = inferred_conf,
    summary = summary_txt
  )
  class(obj) <- c("iconic_prospect", "list")
  obj
}


#' Build prospect summary (internal)
#' @keywords internal
.build_prospect_summary <- function(strength_surface, prospective,
                                    target_gamma_G, rep_diag, n,
                                    mo_confounding, phi) {
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

  paste(lines, collapse = "\n")
}


#' Print method for iconic_prospect objects
#'
#' @param x An \code{iconic_prospect} object.
#' @param ... Unused.
#' @return Invisibly returns `x` (the `iconic_prospect` object); called for its side effect of printing a human-readable summary.
#' @export
print.iconic_prospect <- function(x, ...) {
  cat("<iconic_prospect>\n")
  cat(x$summary, "\n")
  if (!is.null(x$texture_source))
    cat(" Texture:", x$texture_source, "\n")
  if (!is.null(x$inferred_confounding))
    cat(" Confounding: inferred from data\n")
  if (!is.null(x$recommendation))
    cat("\n Recommended estimator if instruments collected:",
        x$recommendation$recommended, "\n")
  invisible(x)
}
