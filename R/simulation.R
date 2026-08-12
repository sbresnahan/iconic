# ============================================================
# Exported user-facing simulation functions:
# run_simulation() – single parameter set, many iterations
# sweep_param() – grid sweep over one parameter
# run_null_sim() – Type I error experiment
# sweep_null_by_conf() – Type I error vs confounding strength
# compute_iter_bias() – per-seed bias vector for boxplots
#
# PGC now uses matrix bridge. Scalar variant (fit_pgc_scalar) is standalone-only.
# ============================================================


#' Run repeated simulations for a single parameter configuration
#'
#' @param n_iter Number of simulation replicates. Default 100.
#' @param n_samples Observations per replicate. Default 500.
#' @param n_features Number of outcome and negative-control features. Default 20.
#' @param beta_X Direct effect of X on Y. Default 0.10.
#' @param alpha_M Effect of X on mediator. Default 0.50.
#' @param beta_M Effect of mediator on Y. Default 0.30.
#' @param conf_str Confounding strength delta. Default 0.80.
#' @param w_signal Proxy quality omega. Default 0.70.
#' @param feat_cor Within-module correlation for block-diagonal co-expression
#' modules in Y and W. 0 = independent features. Default 0.
#' @param base_seed Starting seed; replicate i uses base_seed + i. Default 100.
#' @param n_cores Number of parallel workers. Default 1.
#'
#' @return A list with raw, summary, iter_bias, true_total, params.
#' @export
#'
#' @examples
#' res <- run_simulation(n_iter = 3, n_samples = 100, beta_X = 0.1,
#'                       conf_str = 0.8)
#' res$summary

run_simulation <- function(n_iter = 100,
                           n_samples = 500,
                           n_features = 20,
                           beta_X = 0.10,
                           alpha_M = 0.50,
                           beta_M = 0.30,
                           conf_str = 0.80,
                           w_signal= 0.70,
                           feat_cor = 0,
                           base_seed = 100,
                           n_cores= 1) {

  true_total <- beta_X + alpha_M * beta_M

  worker <- function(i) {
    dat <- generate_toy_data(n = n_samples, n_features = n_features,
                                  beta_X = beta_X, alpha_M = alpha_M,
                                  beta_M = beta_M, conf_str = conf_str,
                                  w_signal = w_signal, feat_cor = feat_cor,
                                  seed = base_seed + i)
    res <- run_methods(dat, n_features)
    res$iter <- i
    res
  }

  iter_results <- .parallel_lapply(seq_len(n_iter), worker, n_cores)
  combined <- do.call(rbind, iter_results)

  list(
    raw = combined,
    summary = summarise_results(combined, true_total),
    iter_bias = compute_iter_bias(combined, true_total),
    true_total = true_total,
    params = list(n_iter = n_iter, n_samples = n_samples,
                      n_features = n_features, beta_X = beta_X,
                      alpha_M = alpha_M, beta_M = beta_M,
                      conf_str = conf_str, w_signal = w_signal,
                      feat_cor = feat_cor)
  )
}


#' Sweep a single simulation parameter across a grid
#'
#' @param param Parameter to vary: one of "beta_X", "conf_str",
#' "w_signal", "alpha_M", "beta_M", "n_samples", "feat_cor".
#' @param param_grid Numeric vector of values to sweep.
#' @param n_iter Replicates per grid point. Default 100.
#' @param n_samples Observations per replicate. Default 500.
#' @param n_features Features per replicate. Default 20.
#' @param beta_X Baseline direct effect. Default 0.10.
#' @param alpha_M Baseline mediator path. Default 0.50.
#' @param beta_M Baseline mediator effect. Default 0.30.
#' @param conf_str Baseline confounding strength. Default 0.80.
#' @param w_signal Baseline proxy quality. Default 0.70.
#' @param feat_cor Baseline within-module feature correlation. Default 0.
#' @param base_seed Seed offset. Default 0.
#' @param n_cores Parallel workers. Default 1.
#'
#' @return A list with summary (data frame) and iter_bias (data frame).
#' @export
#'
#' @examples
#' res <- sweep_param("conf_str", c(0.2, 0.8), n_iter = 3, n_samples = 100)
#' res$summary
sweep_param <- function(param,
                        param_grid,
                        n_iter = 100,
                        n_samples = 500,
                        n_features = 20,
                        beta_X = 0.10,
                        alpha_M = 0.50,
                        beta_M = 0.30,
                        conf_str = 0.80,
                        w_signal = 0.70,
                        feat_cor = 0,
                        base_seed = 0,
                        n_cores = 1) {

  allowed <- c("beta_X", "conf_str", "w_signal", "alpha_M", "beta_M",
               "n_samples", "feat_cor")
  param <- match.arg(param, allowed)

  base_args <- list(n = n_samples, n_features = n_features,
                    beta_X = beta_X, alpha_M = alpha_M,
                    beta_M = beta_M, conf_str = conf_str, w_signal = w_signal,
                    feat_cor = feat_cor)

  smry_list <- list()
  ibias_list <- list()

  for (gi in seq_along(param_grid)) {
    pval <- param_grid[gi]
    args <- base_args
    if (param == "n_samples") { args$n <- pval } else { args[[param]] <- pval }
    true_total <- args$beta_X + args$alpha_M * args$beta_M

    worker <- function(i) {
      dat <- do.call(generate_toy_data,
                          c(args, list(seed = base_seed + gi * 1000L + i)))
      res <- run_methods(dat, args$n_features)
      res$iter <- i
      res
    }

    iter_results <- .parallel_lapply(seq_len(n_iter), worker, n_cores)
    combined <- do.call(rbind, iter_results)
    smry <- summarise_results(combined, true_total)
    smry$param_value <- pval
    smry$true_total <- true_total
    smry_list[[gi]] <- smry

    ibias <- compute_iter_bias(combined, true_total)
    ibias$pval <- pval
    ibias_list[[gi]] <- ibias
  }

  smry_all <- do.call(rbind, smry_list)
  smry_all$param <- param
  smry_all <- smry_all[, c("param", "param_value", "true_total",
                                 setdiff(names(smry_all),
                                         c("param", "param_value", "true_total")))]
  list(
    summary = smry_all,
    iter_bias = do.call(rbind, ibias_list)
  )
}


#' Run null simulations to estimate Type I error rates
#'
#' @param n_iter Number of replicates. Default 200.
#' @param n_samples Observations per replicate. Default 500.
#' @param n_features Features per replicate. Default 20.
#' @param conf_str Confounding strength delta. Default 0.80.
#' @param w_signal Proxy quality omega. Default 0.70.
#' @param feat_cor Within-module feature correlation. Default 0.
#' @param base_seed Seed offset. Default 300.
#' @param n_cores Parallel workers. Default 1.
#' @param alpha Significance threshold. Default 0.05.
#'
#' @return A list with rates (data frame) and raw (full results).
#' @export
#' @examples
#' null <- run_null_sim(n_iter = 2, n_samples = 100, n_features = 5)
#' null$rates
run_null_sim <- function(n_iter = 200,
                         n_samples = 500,
                         n_features = 20,
                         conf_str = 0.80,
                         w_signal = 0.70,
                         feat_cor = 0,
                         base_seed = 300,
                         n_cores = 1,
                         alpha = 0.05) {

  worker <- function(i) {
    dat <- generate_toy_data(n = n_samples, n_features = n_features,
                                  beta_X = 0, alpha_M = 0, beta_M = 0,
                                  conf_str = conf_str, w_signal = w_signal,
                                  feat_cor = feat_cor,
                                  seed = base_seed + i)
    res <- run_methods(dat, n_features)
    res$iter <- i
    res
  }

  iter_results <- .parallel_lapply(seq_len(n_iter), worker, n_cores)
  null_combined <- do.call(rbind, iter_results)

  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC")
  rates <- vapply(methods, function(m) {
    mean(null_combined$pvalue[null_combined$method == m] < alpha, na.rm = TRUE)
  }, numeric(1))

  list(
    rates = data.frame(
      method = methods,
      type1_error = rates,
      flag = ifelse(rates > 3 * alpha, "INFLATED", ifelse(rates < alpha / 2.5, "CONSERVATIVE", "OK")),
      stringsAsFactors = FALSE
    ),
    raw = null_combined
  )
}


#' Sweep Type I error rate across confounding strength levels
#'
#' @param conf_grid Numeric vector of confounding strength values. Default c(0.2, 0.4, 0.6, 0.8, 1.0).
#' @param n_iter Replicates per conf_str value. Default 100.
#' @param n_samples Observations per replicate. Default 500.
#' @param n_features Features per replicate. Default 20.
#' @param w_signal Proxy quality omega. Default 0.70.
#' @param feat_cor Within-module feature correlation. Default 0.
#' @param base_seed Seed offset. Default 900.
#' @param n_cores Parallel workers. Default 1.
#' @param alpha Significance threshold. Default 0.05.
#'
#' @return A data frame with columns: conf_str, method, type1_error.
#' @export
#'
#' @examples
#' t1e <- sweep_null_by_conf(c(0.2, 0.8), n_iter = 3, n_samples = 100)
#' plot_type1_vs_conf(t1e)
sweep_null_by_conf <- function(conf_grid = c(0.2, 0.4, 0.6, 0.8, 1.0),
                               n_iter = 100,
                               n_samples = 500,
                               n_features = 20,
                               w_signal = 0.70,
                               feat_cor = 0,
                               base_seed = 900,
                               n_cores = 1,
                               alpha = 0.05) {

  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC")

  do.call(rbind, lapply(seq_along(conf_grid), function(ci) {
    cs <- conf_grid[ci]

    worker <- function(i) {
      dat <- generate_toy_data(n = n_samples, n_features = n_features,
                                    beta_X = 0, alpha_M = 0, beta_M = 0,
                                    conf_str = cs, w_signal = w_signal,
                                    feat_cor = feat_cor,
                                    seed = base_seed + ci * 1000L + i)
      res <- run_methods(dat, n_features)
      res$iter <- i
      res
    }

    iter_results <- .parallel_lapply(seq_len(n_iter), worker, n_cores)
    combined <- do.call(rbind, iter_results)

    do.call(rbind, lapply(methods, function(m) {
      data.frame(conf_str = cs,
                 method = m,
                 type1_error = mean(combined$pvalue[combined$method == m] < alpha, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
  }))
}


#' Extract per-seed bias vector (internal)
#'
#' @param combined Data frame from run_methods().
#' @param true_total Scalar true total causal effect.
#' @return Data frame with columns: iter, method, bias.
#' @keywords internal
compute_iter_bias <- function(combined, true_total) {
  iters <- sort(unique(combined$iter))
  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC")
  do.call(rbind, lapply(iters, function(i) {
    sub <- combined[combined$iter == i, ]
    do.call(rbind, lapply(methods, function(m) {
      x <- sub$beta[sub$method == m]
      data.frame(iter = i, method = m,
                 bias = mean(x, na.rm = TRUE) - true_total,
                 stringsAsFactors = FALSE)
    }))
  }))
}


#' Parallel lapply dispatcher (internal)
#'
#' Cross-platform parallel lapply. Uses \code{parallel::mclapply}
#' (fork-based) on Unix-like systems and a PSOCK cluster
#' (\code{makeCluster}/\code{parLapply}) on Windows. Falls back to
#' sequential \code{lapply} when \code{n_cores <= 1} or the task list
#' has fewer than 2 elements.
#'
#' @param X List or vector to iterate over.
#' @param FUN Function applied to each element of X.
#' @param n_cores Number of parallel workers. Default 1 (sequential).
#' @param progress Optional character label for progress messages.
#' When non-NULL, prints a start message and then a
#' percentage-complete update at ~10% increments (plus
#' a final completion line with total elapsed time).
#' In parallel mode, tasks are split into ~10 chunks and
#' a progress message is printed after each chunk
#' completes. Default NULL (silent).
#' @keywords internal
#' @noRd
.parallel_lapply <- function(X, FUN, n_cores = 1, progress = NULL) {
  n <- length(X)
  if (n == 0) return(list())

  # ── Sequential path ──
  if (n_cores <= 1 || n < 2) {
    if (!is.null(progress)) {
      message(progress, ": ", n, " tasks (sequential)")
      start_time <- Sys.time()
      out <- vector("list", n)
      last_band <- 0L # 10% band already reported (skip the 0% band)
      for (i in seq_len(n)) {
        out[[i]] <- FUN(X[[i]])
        band <- as.integer(floor(10 * i / n)) # 0..10
        if (band > last_band) {
          elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
          pct <- round(100 * i / n)
          message(" ", progress, ": ", pct, "% (", i, "/", n, ") [",
                  round(elapsed, 1), "s]")
          last_band <- band
        }
      }
      return(out)
    }
    return(lapply(X, FUN))
  }

  # ── Parallel path (chunked for progress reporting) ──
  if (!is.null(progress))
    message(progress, ": ", n, " tasks on ", n_cores, " cores")

  # Split into ~10 chunks (one per 10% band). For small n, one chunk per
  # task. Breakpoints are inclusive starts; chunk k covers [bp[k], bp[k+1]-1].
  n_chunks <- min(10L, n)
  bp <- as.integer(seq.int(1L, n + 1L, length.out = n_chunks + 1L))
  bp[n_chunks + 1L] <- n + 1L # guard against rounding

  out <- vector("list", n)
  start_time <- Sys.time()

  if (.Platform$OS.type != "windows") {
    # Unix: fork-based mclapply, one call per chunk
    for (k in seq_len(n_chunks)) {
      idx <- bp[k]:(bp[k + 1L] - 1L)
      out[idx] <- parallel::mclapply(X[idx], FUN, mc.cores = n_cores)
      if (!is.null(progress)) {
        done <- bp[k + 1L] - 1L
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        pct <- round(100 * done / n)
        message(" ", progress, ": ", pct, "% (", done, "/", n, ") [",
                round(elapsed, 1), "s]")
      }
    }
  } else {
    # Windows: PSOCK cluster created ONCE, reused across all chunks
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Load the iconic namespace so internal functions resolve in workers
    parallel::clusterEvalQ(cl, {
      if (requireNamespace("iconic", quietly = TRUE))
        attachNamespace("iconic")
    })
    for (k in seq_len(n_chunks)) {
      idx <- bp[k]:(bp[k + 1L] - 1L)
      out[idx] <- parallel::parLapply(cl, X[idx], FUN)
      if (!is.null(progress)) {
        done <- bp[k + 1L] - 1L
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        pct <- round(100 * done / n)
        message(" ", progress, ": ", pct, "% (", done, "/", n, ") [",
                round(elapsed, 1), "s]")
      }
    }
  }

  out
}
