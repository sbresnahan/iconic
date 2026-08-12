# ============================================================
# Exported user-facing mediation simulation functions:
# run_mediation_sim() – single config, many iterations
# sweep_mediation_param() – grid sweep over one parameter
# run_null_mediation_sim() – Type I error for NDE and NIE
# sweep_mediation_null_by_conf() – Type I error vs confounding strength
#
# These mirror the total-effect simulation API (run_simulation, sweep_param,
# run_null_sim, sweep_null_by_conf) but estimate NDE and NIE instead of the
# total effect tau, and use the five (or six, when phi > 0) mediation
# estimators from mediation.R.
#
# the phi parameter (default 0) adds a mediator-specific genetic
# instrument (Gm). When phi > 0, the 2-stage MR estimator (IV2SLS2) is
# included in the results, and "phi" can be swept as a parameter.
#
# adds rho_G1, rho_G2, rho_pop, lambda_XM, lambda_MY, omega_1, omega_2
# parameters for the imperfect-independence DGP. When any is non-default,
# the two-stage proximal estimators (PGC2, PGC2Gm) are included in the
# results. rho_G1, rho_G2, rho_pop, omega_1, omega_2 can be swept.
# ============================================================

# Canonical method order
.mediation_method_order <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
                             "IV2SLS2", "PGC2", "PGC2Gm")


#' Run repeated mediation simulations for a single parameter configuration
#'
#' Generates \code{n_iter} synthetic datasets with mediator-outcome
#' confounding (via \code{mo_confounding}), runs all five (or six, when
#' \code{phi > 0}) mediation estimators, and summarises bias / RMSE / Type I
#' error for both NDE and NIE.
#'
#' When are non-default (rho_G1, rho_G2, rho_pop,
#' lambda_XM, lambda_MY, omega_1, omega_2), the two-stage proximal estimators (PGC2,
#' PGC2Gm) are also included.
#'
#' @param n_iter Number of simulation replicates. Default 100.
#' @param n_samples Observations per replicate. Default 500.
#' @param n_features Number of outcome and negative-control features. Default 20.
#' @param n_mediators Number of independent mediators. When > 1,
#' each mediator has its own genetic instrument Gm and
#' contributes additively to Y. Default 1 (single mediator,
#' backward compatible).
#' @param beta_X Direct effect of X on Y (true NDE). Default 0.10.
#' @param alpha_M Effect of X on mediator. Default 0.50.
#' @param beta_M Effect of mediator on Y (per-mediator true NIE = alpha_M * beta_M). Default 0.30.
#' @param conf_str Confounding strength delta. Default 0.80.
#' @param w_signal Proxy quality omega. Default 0.70.
#' @param mo_confounding Strength of U1 -> M (mediator-outcome confounding). Default 0.80.
#' @param phi Strength of the mediator instrument Gm -> M.
#' 0 = no mediator instrument (five estimators, backward
#' compatible). > 0 = generates Gm and includes the
#' 2-stage MR estimator (IV2SLS2). Default 0.
#' @param rho_G1 Correlation of G1 with conf_XM. Default 0.
#' @param rho_G2 Correlation of G2 with conf_MY. Default 0.
#' @param rho_pop Shared population structure. Default 0.
#' @param lambda_XM Optional per-path confounder loading vector (X->M path).
#' @param lambda_MY Optional per-path confounder loading vector (M->Y path).
#' @param omega_1 Coverage of conf_XM by W1. NULL = use w_signal.
#' @param omega_2 Coverage of conf_MY by W2. NULL = use w_signal.
#' @param feat_cor Within-module correlation for block-diagonal co-expression
#' modules in Y and W. 0 = independent features. Default 0.
#' @param base_seed Starting seed; replicate i uses base_seed + i. Default 100.
#' @param n_cores Number of parallel workers. Default 1.
#'
#' @return A list with \code{summary} (one row per method with NDE/NIE bias,
#' RMSE, Type I), \code{raw} (full per-feature results), \code{true_NDE},
#' \code{true_NIE}, and \code{params}.
#' @export
#'
#' @examples
#' res <- run_mediation_sim(n_iter = 3, n_samples = 100,
#'                          mo_confounding = 0.8)
#' res$summary
run_mediation_sim <- function(n_iter = 100,
                              n_samples = 500,
                              n_features = 20,
                              n_mediators = 1,
                              beta_X = 0.10,
                              alpha_M = 0.50,
                              beta_M = 0.30,
                              conf_str = 0.80,
                              w_signal = 0.70,
                              mo_confounding = 0.80,
                              phi = 0,
                              rho_G1 = 0,
                              rho_G2 = 0,
                              rho_pop = 0,
                              lambda_XM = NULL,
                              lambda_MY = NULL,
                              omega_1 = NULL,
                              omega_2 = NULL,
                              feat_cor = 0,
                              base_seed = 100,
                              n_cores = 1) {

  true_NDE <- beta_X
  true_NIE <- alpha_M * beta_M # per-mediator NIE

  worker <- function(i) {
    dat <- generate_toy_data(n = n_samples, n_features = n_features,
                             n_mediators = n_mediators,
                             beta_X = beta_X, alpha_M = alpha_M,
                             beta_M = beta_M, conf_str = conf_str,
                             w_signal = w_signal, mo_confounding = mo_confounding,
                             phi = phi, rho_G1 = rho_G1, rho_G2 = rho_G2,
                             rho_pop = rho_pop, lambda_XM = lambda_XM, lambda_MY = lambda_MY,
                             omega_1 = omega_1, omega_2 = omega_2,
                             feat_cor = feat_cor,
                             seed = base_seed + i)
    res <- run_mediation_methods(dat, n_features)
    res$iter <- i
    res
  }

  iter_results <- .parallel_lapply(seq_len(n_iter), worker, n_cores)
  combined <- do.call(rbind, iter_results)

  list(
    summary = summarise_mediation_results(combined, true_NDE, true_NIE),
    raw = combined,
    true_NDE = true_NDE,
    true_NIE = true_NIE,
    params = list(n_iter = n_iter, n_samples = n_samples,
                    n_features = n_features, n_mediators = n_mediators,
                    beta_X = beta_X,
                    alpha_M = alpha_M, beta_M = beta_M,
                    conf_str = conf_str, w_signal = w_signal,
                    mo_confounding = mo_confounding, phi = phi,
                    rho_G1 = rho_G1, rho_G2 = rho_G2,
                    rho_pop = rho_pop, lambda_XM = lambda_XM, lambda_MY = lambda_MY,
                    omega_1 = omega_1, omega_2 = omega_2,
                    feat_cor = feat_cor)
  )
}


#' Sweep a single mediation simulation parameter across a grid
#'
#' @param param Parameter to vary: one of "beta_X", "conf_str",
#' "w_signal", "alpha_M", "beta_M", "n_samples", "mo_confounding", "phi",
#' "rho_G1", "rho_G2", "rho_pop", "omega_1", "omega_2",
#' "feat_cor".
#' @param param_grid Numeric vector of values to sweep.
#' @param n_iter Replicates per grid point. Default 100.
#' @param n_samples Observations per replicate. Default 500.
#' @param n_features Features per replicate. Default 20.
#' @param n_mediators Number of independent mediators. Default 1.
#' @param beta_X Baseline direct effect. Default 0.10.
#' @param alpha_M Baseline mediator path. Default 0.50.
#' @param beta_M Baseline mediator effect. Default 0.30.
#' @param conf_str Baseline confounding strength. Default 0.80.
#' @param w_signal Baseline proxy quality. Default 0.70.
#' @param mo_confounding Baseline M-O confounding. Default 0.80.
#' @param phi Baseline mediator-instrument strength.
#' 0 = no mediator instrument. Default 0.
#' @param rho_G1 Baseline G1-conf_XM correlation. Default 0.
#' @param rho_G2 Baseline G2-conf_MY correlation. Default 0.
#' @param rho_pop Baseline population structure. Default 0.
#' @param lambda_XM Optional per-path confounder loading vector (X->M path).
#' @param lambda_MY Optional per-path confounder loading vector (M->Y path).
#' @param omega_1 Baseline W1 coverage. NULL = use w_signal.
#' @param omega_2 Baseline W2 coverage. NULL = use w_signal.
#' @param feat_cor Baseline within-module feature correlation. Default 0.
#' @param u_strength Numeric vector: per-confounder
#' strength scaling. \code{NULL} = uniform (legacy behavior). See
#' \code{\link{generate_toy_data}()}.
#' @param w_coverage_profile A list with \code{w1}/\code{w2} per-control
#' coverage vectors. \code{NULL} = uniform coverage.
#' See \code{\link{generate_toy_data}()}.
#' @param base_seed Seed offset. Default 0.
#' @param n_cores Parallel workers. Default 1.
#'
#' @return A list with \code{summary} (data frame) and \code{iter_bias}
#' (per-iteration NDE/NIE bias).
#' @export
#'
#' @examples
#' res <- sweep_mediation_param("conf_str", c(0.2, 0.8), n_iter = 3,
#'                              n_samples = 100)
#' res$summary
sweep_mediation_param <- function(param,
                                  param_grid,
                                  n_iter = 100,
                                  n_samples = 500,
                                  n_features = 20,
                                  n_mediators = 1,
                                  beta_X = 0.10,
                                  alpha_M = 0.50,
                                  beta_M = 0.30,
                                  conf_str = 0.80,
                                  w_signal = 0.70,
                                  mo_confounding = 0.80,
                                  phi = 0,
                                  rho_G1 = 0,
                                  rho_G2 = 0,
                                  rho_pop = 0,
                                  lambda_XM = NULL,
                                  lambda_MY = NULL,
                                  omega_1 = NULL,
                                  omega_2 = NULL,
                                  feat_cor = 0,
                                  u_strength = NULL,
                                  w_coverage_profile = NULL,
                                  base_seed = 0,
                                  n_cores = 1) {

  allowed <- c("beta_X", "conf_str", "w_signal", "alpha_M", "beta_M",
               "n_samples", "mo_confounding", "phi",
               "rho_G1", "rho_G2", "rho_pop", "omega_1", "omega_2",
               "feat_cor")
  param <- match.arg(param, allowed)

  # Build base args using 'n' (the generate_toy_data parameter name)
  base_args <- list(n = n_samples, n_features = n_features,
                    n_mediators = n_mediators,
                    beta_X = beta_X, alpha_M = alpha_M,
                    beta_M = beta_M, conf_str = conf_str, w_signal = w_signal,
                    mo_confounding = mo_confounding, phi = phi,
                    rho_G1 = rho_G1, rho_G2 = rho_G2, rho_pop = rho_pop,
                    lambda_XM = lambda_XM, lambda_MY = lambda_MY, omega_1 = omega_1, omega_2 = omega_2,
                    feat_cor = feat_cor,
                    u_strength = u_strength,
                    w_coverage_profile = w_coverage_profile)

  smry_list <- list()
  ibias_list <- list()

  for (gi in seq_along(param_grid)) {
    pval <- param_grid[gi]
    args <- base_args

    # Map the user-facing param name to the generate_toy_data arg name
    if (param == "n_samples") {
      args$n <- pval
    } else {
      args[[param]] <- pval
    }

    true_NDE <- args$beta_X
    true_NIE <- args$alpha_M * args$beta_M # per-mediator NIE

    worker <- function(i) {
      dat <- do.call(generate_toy_data,
                     c(args, list(seed = base_seed + gi * 1000L + i)))
      res <- run_mediation_methods(dat, args$n_features)
      res$iter <- i
      res
    }

    iter_results <- .parallel_lapply(seq_len(n_iter), worker, n_cores)
    combined <- do.call(rbind, iter_results)
    smry <- summarise_mediation_results(combined, true_NDE, true_NIE)
    smry$param_value <- pval
    smry$true_NDE <- true_NDE
    smry$true_NIE <- true_NIE
    smry_list[[gi]] <- smry

    # Per-iteration bias -- iterate over methods present in the data
    methods_present <- unique(combined$method)
    methods <- intersect(.mediation_method_order, methods_present)

    iters <- sort(unique(combined$iter))
    ibias <- do.call(rbind, lapply(iters, function(ii) {
      sub <- combined[combined$iter == ii, ]
      do.call(rbind, lapply(methods, function(m) {
        x <- sub[sub$method == m, ]
        data.frame(iter = ii, method = m,
                   NDE_bias = mean(x$NDE, na.rm = TRUE) - true_NDE,
                   NIE_bias = mean(x$NIE, na.rm = TRUE) - true_NIE,
                   stringsAsFactors = FALSE)
      }))
    }))
    ibias$pval <- pval
    ibias_list[[gi]] <- ibias
  }

  smry_all <- do.call(rbind, smry_list)
  smry_all$param <- param
  smry_all <- smry_all[, c("param", "param_value", "true_NDE", "true_NIE",
                           setdiff(names(smry_all),
                                   c("param", "param_value", "true_NDE", "true_NIE")))]
  list(
    summary = smry_all,
    iter_bias = do.call(rbind, ibias_list)
  )
}


#' Run null mediation simulations to estimate Type I error rates
#'
#' Sets \code{alpha_M = 0} and \code{beta_M = 0} (no true NIE) while
#' keeping \code{beta_X} active (NDE still present). Reports Type I
#' error for both NIE (should be 0 under the null) and NDE (should
#' reflect power for the direct effect).
#'
#' @param n_iter Number of replicates. Default 200.
#' @param n_samples Observations per replicate. Default 500.
#' @param n_features Features per replicate. Default 20.
#' @param beta_X Direct effect (NDE, still active under null NIE). Default 0.10.
#' @param conf_str Confounding strength delta. Default 0.80.
#' @param w_signal Proxy quality omega. Default 0.70.
#' @param mo_confounding Strength of U1 -> M. Default 0.80.
#' @param phi Strength of the mediator instrument Gm -> M.
#' 0 = no mediator instrument. Default 0.
#' @param rho_G1 Correlation of G1 with conf_XM. Default 0.
#' @param rho_G2 Correlation of G2 with conf_MY. Default 0.
#' @param rho_pop Shared population structure. Default 0.
#' @param lambda_XM Optional per-path confounder loading vector (X->M path).
#' @param lambda_MY Optional per-path confounder loading vector (M->Y path).
#' @param omega_1 Coverage of conf_XM by W1. NULL = use w_signal.
#' @param omega_2 Coverage of conf_MY by W2. NULL = use w_signal.
#' @param feat_cor Within-module correlation for block-diagonal co-expression
#' modules in Y and W. 0 = independent features. Default 0.
#' @param base_seed Seed offset. Default 300.
#' @param n_cores Parallel workers. Default 1.
#' @param alpha Significance threshold. Default 0.05.
#'
#' @return A list with \code{rates} (data frame: method, NIE_type1, NDE_type1)
#' and \code{raw} (full results).
#' @export
#' @examples
#' null <- run_null_mediation_sim(n_iter = 2, n_samples = 100,
#'   n_features = 5, mo_confounding = 0.8, phi = 0.8)
#' null$rates
run_null_mediation_sim <- function(n_iter = 200,
                                   n_samples = 500,
                                   n_features = 20,
                                   beta_X = 0.10,
                                   conf_str = 0.80,
                                   w_signal = 0.70,
                                   mo_confounding = 0.80,
                                   phi = 0,
                                   rho_G1 = 0,
                                   rho_G2 = 0,
                                   rho_pop = 0,
                                   lambda_XM = NULL,
                                   lambda_MY = NULL,
                                   omega_1 = NULL,
                                   omega_2 = NULL,
                                   feat_cor = 0,
                                   base_seed = 300,
                                   n_cores = 1,
                                   alpha = 0.05) {

  worker <- function(i) {
    dat <- generate_toy_data(n = n_samples, n_features = n_features,
                             beta_X = beta_X, alpha_M = 0, beta_M = 0,
                             conf_str = conf_str, w_signal = w_signal,
                             mo_confounding = mo_confounding,
                             phi = phi, rho_G1 = rho_G1, rho_G2 = rho_G2,
                             rho_pop = rho_pop, lambda_XM = lambda_XM, lambda_MY = lambda_MY,
                             omega_1 = omega_1, omega_2 = omega_2,
                             feat_cor = feat_cor,
                             seed = base_seed + i)
    res <- run_mediation_methods(dat, n_features)
    res$iter <- i
    res
  }

  iter_results <- .parallel_lapply(seq_len(n_iter), worker, n_cores)
  null_combined <- do.call(rbind, iter_results)

  # Iterate over methods present in the data
  methods_present <- unique(null_combined$method)
  methods <- intersect(.mediation_method_order, methods_present)

  rates <- do.call(rbind, lapply(methods, function(m) {
    sub <- null_combined[null_combined$method == m, ]
    data.frame(
      method = m,
      NIE_type1 = mean(sub$NIE_p < alpha, na.rm = TRUE),
      NDE_type1 = mean(sub$NDE_p < alpha, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  list(rates = rates, raw = null_combined)
}


#' Sweep mediation Type I error across confounding strength levels
#'
#' @param conf_grid Numeric vector of confounding strength values.
#' Default c(0.2, 0.4, 0.6, 0.8, 1.0).
#' @param n_iter Replicates per conf_str value. Default 100.
#' @param n_samples Observations per replicate. Default 500.
#' @param n_features Features per replicate. Default 20.
#' @param w_signal Proxy quality omega. Default 0.70.
#' @param mo_confounding Strength of U1 -> M. Default 0.80.
#' @param phi Strength of the mediator instrument Gm -> M.
#' 0 = no mediator instrument. Default 0.
#' @param rho_G1 Correlation of G1 with conf_XM. Default 0.
#' @param rho_G2 Correlation of G2 with conf_MY. Default 0.
#' @param rho_pop Shared population structure. Default 0.
#' @param lambda_XM Optional per-path confounder loading vector (X->M path).
#' @param lambda_MY Optional per-path confounder loading vector (M->Y path).
#' @param omega_1 Coverage of conf_XM by W1. NULL = use w_signal.
#' @param omega_2 Coverage of conf_MY by W2. NULL = use w_signal.
#' @param feat_cor Within-module feature correlation. Default 0.
#' @param base_seed Seed offset. Default 900.
#' @param n_cores Parallel workers. Default 1.
#' @param alpha Significance threshold. Default 0.05.
#'
#' @return A data frame with columns: conf_str, method, NIE_type1, NDE_type1.
#' @export
#'
#' @examples
#' t1e <- sweep_mediation_null_by_conf(c(0.2, 0.8), n_iter = 3,
#'                                     n_samples = 100)
#' head(t1e)
sweep_mediation_null_by_conf <- function(conf_grid = c(0.2, 0.4, 0.6, 0.8, 1.0),
                                         n_iter = 100,
                                         n_samples = 500,
                                         n_features = 20,
                                         w_signal = 0.70,
                                         mo_confounding = 0.80,
                                         phi = 0,
                                         rho_G1 = 0,
                                         rho_G2 = 0,
                                         rho_pop = 0,
                                         lambda_XM = NULL,
                                         lambda_MY = NULL,
                                         omega_1 = NULL,
                                         omega_2 = NULL,
                                         feat_cor = 0,
                                         base_seed = 900,
                                         n_cores = 1,
                                         alpha = 0.05) {

  # Determine methods from a pilot run (handles phi > 0 and cases)
  pilot <- generate_toy_data(n = n_samples, n_features = n_features,
                             beta_X = 0, alpha_M = 0, beta_M = 0,
                             conf_str = conf_grid[1], w_signal = w_signal,
                             mo_confounding = mo_confounding,
                             phi = phi, rho_G1 = rho_G1, rho_G2 = rho_G2,
                             rho_pop = rho_pop, lambda_XM = lambda_XM, lambda_MY = lambda_MY,
                             omega_1 = omega_1, omega_2 = omega_2,
                             feat_cor = feat_cor,
                             seed = base_seed)
  pilot_res <- run_mediation_methods(pilot, n_features)
  methods_present <- unique(pilot_res$method)
  methods <- intersect(.mediation_method_order, methods_present)

  do.call(rbind, lapply(seq_along(conf_grid), function(ci) {
    cs <- conf_grid[ci]

    worker <- function(i) {
      dat <- generate_toy_data(n = n_samples, n_features = n_features,
                               beta_X = 0, alpha_M = 0, beta_M = 0,
                               conf_str = cs, w_signal = w_signal,
                               mo_confounding = mo_confounding,
                               phi = phi, rho_G1 = rho_G1, rho_G2 = rho_G2,
                               rho_pop = rho_pop, lambda_XM = lambda_XM, lambda_MY = lambda_MY,
                               omega_1 = omega_1, omega_2 = omega_2,
                               feat_cor = feat_cor,
                               seed = base_seed + ci * 1000L + i)
      res <- run_mediation_methods(dat, n_features)
      res$iter <- i
      res
    }

    iter_results <- .parallel_lapply(seq_len(n_iter), worker, n_cores)
    combined <- do.call(rbind, iter_results)

    do.call(rbind, lapply(methods, function(m) {
      sub <- combined[combined$method == m, ]
      data.frame(conf_str = cs,
                 method = m,
                 NIE_type1 = mean(sub$NIE_p < alpha, na.rm = TRUE),
                 NDE_type1 = mean(sub$NDE_p < alpha, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
  }))
}
