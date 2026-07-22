# ============================================================
# Pleiotropy sensitivity analysis.
#
# gan_pleiotropy_sensitivity() benchmarks every estimator across a grid
# of horizontal-pleiotropy strengths (direct G -> Y paths that violate
# the exclusion restriction) and confounding strengths, running both an
# alternative arm (true effect > 0) and a null arm (true effect = 0) so
# that bias and Type I error can be read off the same sweep.
#
# This implements the pleiotropy sweep from critique #4 of the QED
# review, wrapping it into the package's sensitivity framework so it
# shares the same generator / parallel infrastructure as
# gan_sensitivity().
# ============================================================

#' Benchmark estimators across pleiotropy and confounding scenarios
#'
#' For each cell in the grid (pleiotropy strength x confounding strength),
#' generates `n_iter` synthetic datasets with [run_single_iteration()],
#' runs every estimator, and summarises bias / RMSE / power. Two arms are
#' run per cell:
#'
#' - **alternative** (`effect_size = tau`): the true total effect is `tau`,
#'   so `power` is the empirical power to detect it.
#' - **null** (`effect_size = 0`): the true total effect is 0, so `power`
#'   is the empirical Type I error rate.
#'
#' The `pleio` parameter adds a direct `G -> Y` path of the requested
#' strength, violating the exclusion restriction. IV/2SLS is consistent
#' only when `pleio = 0`; any `pleio > 0` introduces bias that does not
#' shrink with sample size.
#'
#' @param trained_gan   An `iconic_gan` (or `NULL` to use default texture).
#' @param pleio_grid    Horizontal-pleiotropy strengths (direct G -> Y
#'                      coefficients) to sweep. Default `c(0, 0.05, 0.10)`.
#' @param conf_grid     Confounding-strength values. Default `c(0.2, 0.5, 0.8)`.
#' @param tau           True total effect for the alternative arm. Default 0.25.
#' @param nc_model      Negative-control model (function or name). Default `"proxy"`.
#' @param n_iter        Replicates per cell per arm. Default 50.
#' @param n_samples     Samples per replicate. Default 500.
#' @param n_features    Features per replicate. Default 10.
#' @param coverage      Negative-control coverage. Default 0.7.
#' @param k             Number of latent confounders. Default 1.
#' @param base_seed     Base RNG seed. Default 900.
#' @param n_cores       Parallel workers across replicates. Default 1.
#'
#' @return A list with `summary` (one row per cell x arm x method, with
#'   `pleio`, `conf_strength`, `arm`, `true_total`, and the columns from
#'   `summarise_results()`) and `grid`.
#' @export
#'
#' @examples
#' \dontrun{
#' sens <- gan_pleiotropy_sensitivity(NULL,
#'         pleio_grid = c(0, 0.05, 0.10), conf_grid = c(0.5, 0.8),
#'         n_iter = 20, n_features = 10)
#' head(sens$summary)
#' }
gan_pleiotropy_sensitivity <- function(trained_gan  = NULL,
                                       pleio_grid   = c(0, 0.05, 0.10),
                                       conf_grid    = c(0.2, 0.5, 0.8),
                                       tau          = 0.25,
                                       nc_model     = "proxy",
                                       n_iter       = 50,
                                       n_samples    = 500,
                                       n_features   = 10,
                                       coverage     = 0.7,
                                       k            = 1,
                                       base_seed    = 900,
                                       n_cores      = 1) {

  grid <- expand.grid(pleio = pleio_grid, conf_strength = conf_grid,
                      KEEP.OUT.ATTRS = FALSE)
  arms <- c("alt", "null")

  smry <- lapply(seq_len(nrow(grid)), function(gi) {
    pl <- grid$pleio[gi]; cs <- grid$conf_strength[gi]

    arm_results <- lapply(arms, function(arm) {
      eff <- if (arm == "alt") tau else 0

      worker <- function(i) {
        dat <- run_single_iteration(
          trained_gan, n_synthetic_samples = n_samples,
          n_features = n_features, n_confounders = k,
          effect_size = eff, conf_strength = cs, coverage = coverage,
          nc_model = nc_model, pleio = pl,
          seed = base_seed + gi * 1000L + if (arm == "alt") 0L else 500000L + i)
        res <- run_methods(dat, n_features)
        res$iter <- i
        res
      }

      combined <- do.call(rbind, .parallel_lapply(seq_len(n_iter), worker, n_cores))
      s <- summarise_results(combined, eff)
      s$pleio         <- pl
      s$conf_strength <- cs
      s$arm           <- arm
      s$true_total    <- eff
      s
    })

    do.call(rbind, arm_results)
  })

  summary <- do.call(rbind, smry)
  front   <- c("pleio", "conf_strength", "arm", "true_total", "method")
  summary <- summary[, c(front, setdiff(names(summary), front))]
  list(summary = summary, grid = grid)
}
