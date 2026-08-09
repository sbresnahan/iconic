# ============================================================
# Empirical negative-control (NC) validity diagnostics.
#
# These functions implement the data-driven screens specified in
# They operate on a fitted dataset
# (from run_single_iteration() or generate_toy_data()) and test the
# two identifying assumptions that negative controls must satisfy:
#
# (A1) W _|_ X | C -- controls are not affected by the exposure
# (A2) W _|_ G | C -- controls are independent of the instrument
#
# plus the proximal-inference completeness condition:
#
# (B) dim(W_valid) >= k -- enough valid controls to span the
# confounder subspace
#
# (A2') W _|_ Gm | C -- controls are independent of the mediator
# instrument (tested by nc_independence_check_gm)
#
# In the simulation, A1, A2, and A2' hold by construction (W = f(U), and
# U is independent of both X's causal path, G, and Gm). These screens are
# therefore most useful on REAL data, where they detect violations
# (e.g. a control that is actually downstream of the exposure, or a
# meQTL linking G or Gm to a methylation-based control).
# ============================================================


#' Screen negative controls for exposure dependence (W _|_ X | C)
#'
#' Regresses each negative-control feature on the exposure X (plus
#' observed covariates C) and flags controls whose association with X
#' survives Benjamini-Hochberg FDR control at the requested level.
#'
#' A control that is significantly associated with X after adjusting
#' for C violates the negative-control assumption (A1): it is either
# downstream of the exposure or shares a cause of X that is not in C.
#' Such controls should be dropped before running COCA / PGC.
#'
#' @section Assumption:
#' A1 states that negative-control outcomes are independent of the exposure
#' **conditional on observed covariates C and the unmeasured confounder U**
#' (U is unobserved). NC outcome proxies may cause or be caused by the
#' outcome or instruments; the screen tests the observable marginal
#' association with X given C.
#'
#' @section Criterion:
#' `criterion = "fdr"` (the legacy behavior) flags controls whose
#' association with X survives BH-FDR. `criterion = "magnitude"` flags
#' controls by partial-correlation *size* (|partial_r(W,X|C)| >
#' `magnitude_threshold`), which is less vulnerable to the "always
#' significant via U" problem: because W shares U with X by construction,
#' the FDR test can flag valid controls as violations purely from the
#' intended confounder-sharing signal. `criterion = "both"` (default)
#' requires a control to pass both branches to be deemed valid.
#'
#' Distinguishing "W downstream of X" (a true A1 violation) from "W shares
#' a cause with X" (the intended NC behavior) remains an open problem
#' without a clean test statistic. The `relative_effect` column
#' (|partial_r(W,X|C)| / |partial_r(W,U_proxy|C)|, when `u_proxy` is
#' supplied) is exposed as a *diagnostic*: a control downstream of X tends
#' to have a larger W-X association relative to its W-U association, but
#' the separation is not always clean. Analysts should inspect
#' `partial_r` and `relative_effect` together rather than relying on a
#' single gate.
#'
#' @param dat Dataset list from [run_single_iteration()] or
#' [generate_toy_data()], containing `W`, `X`, and
#' `synthetic_data` (covariates).
#' @param fdr_level Target false-discovery rate for BH correction.
#' Default 0.10.
#' @param alpha Per-test significance level used before FDR
#' adjustment (informational only). Default 0.05.
#' @param n_cores Number of parallel workers. Default 1 (sequential).
#' Uses \code{parallel::mclapply} on Unix and a PSOCK cluster on
#' Windows.
#' @param criterion Which criterion to use: `"fdr"` (legacy),
#' `"magnitude"`, or `"both"` (default).
#' @param magnitude_threshold |partial_r(W,X|C)| cutoff for the
#' magnitude branch. Default 0.10.
#' @param u_proxy Optional outcome-derived proxy for U (e.g. the
#' leading residualized-outcome principal component) used
#' to compute `relative_effect` = |partial_r(W,X|C)| /
#' |partial_r(W,U_proxy|C)|. When `NULL` (default),
#' `relative_effect` is `NA` and only the FDR/magnitude
#' branches are used.
#'
#' @return A data frame with one row per control feature:
#' `feature`, `p_value`, `fdr`, `partial_r`, `relative_effect`,
#' `significant` (composite, per `criterion`),
#' `verdict` (composite), `verdict_fdr`, `verdict_magnitude`.
#' @export
#' @importFrom stats lm as.formula p.adjust
#'
#' @examples
#' dat <- run_single_iteration(n_features = 10, seed = 1)
#' nc_validity_screen(dat)
nc_validity_screen <- function(dat, fdr_level = 0.10, alpha = 0.05,
                               n_cores = 1,
                               criterion = c("both", "fdr", "magnitude"),
                               magnitude_threshold = 0.10,
                               u_proxy = NULL) {
  criterion <- match.arg(criterion)
  W <- dat$W
  X <- dat$X
  cv <- dat$synthetic_data
  p <- ncol(W)
  n <- nrow(W)
  cnames <- if (!is.null(cv)) names(cv) else character(0)
  cs <- if (length(cnames)) paste0(" + ", paste(cnames, collapse = " + ")) else ""

  # Residualize X on covariates once (for the partial correlation).
  X_res <- if (!is.null(cv) && ncol(cv) > 0) {
    zf <- tryCatch(lm(X ~ ., data = cv), error = function(e) NULL)
    if (is.null(zf)) X else residuals(zf)
  } else X
  df_res <- n - 1 - if (!is.null(cv)) ncol(cv) else 0

  # Per-control: p-value (FDR branch) and partial correlation (magnitude branch).
  res_list <- .parallel_lapply(seq_len(p), function(f) {
    w <- W[, f]
    # FDR branch: p-value from w ~ X + C
    d <- if (!is.null(cv)) cbind(data.frame(w = w, X = X), cv) else
           data.frame(w = w, X = X)
    fit <- tryCatch(lm(as.formula(paste0("w ~ X", cs)), data = d),
                    error = function(e) NULL)
    pv <- if (is.null(fit)) NA_real_ else {
      sm <- summary(fit)$coefficients
      if (!"X" %in% rownames(sm)) NA_real_ else as.numeric(sm["X", 4])
    }
    # Magnitude branch: partial correlation of W with X given C.
    w_res <- if (!is.null(cv) && ncol(cv) > 0) {
      wf <- tryCatch(lm(w ~ ., data = cv), error = function(e) NULL)
      if (is.null(wf)) w else residuals(wf)
    } else w
    ok <- complete.cases(w_res, X_res)
    r <- if (sum(ok) < 5) NA_real_ else cor(w_res[ok], X_res[ok])
    # relative_effect vs u_proxy, if supplied.
    rel <- NA_real_
    if (!is.null(u_proxy)) {
      up_res <- if (!is.null(cv) && ncol(cv) > 0) {
        uf <- tryCatch(lm(u_proxy ~ ., data = cv), error = function(e) NULL)
        if (is.null(uf)) u_proxy else residuals(uf)
      } else u_proxy
      ok2 <- complete.cases(w_res, up_res)
      r_wu <- if (sum(ok2) < 5) NA_real_ else cor(w_res[ok2], up_res[ok2])
      if (!is.na(r) && !is.na(r_wu) && abs(r_wu) > 1e-8) rel <- abs(r) / abs(r_wu)
    }
    c(pv, r, rel)
  }, n_cores = n_cores, progress = "NC validity screen")
  res <- do.call(cbind, res_list)

  pvals <- res[1, ]
  partial_r <- res[2, ]
  rel_eff <- res[3, ]
  fdr <- p.adjust(pvals, method = "BH")
  sig_fdr <- !is.na(fdr) & fdr < fdr_level
  sig_mag <- !is.na(partial_r) & abs(partial_r) > magnitude_threshold

  # Composite per criterion.
  sig <- switch(criterion,
    fdr = sig_fdr,
    magnitude = sig_mag,
    both = sig_fdr & sig_mag)

  data.frame(
    feature = seq_len(p),
    p_value = pvals,
    fdr = fdr,
    partial_r = partial_r,
    relative_effect = rel_eff,
    significant = sig,
    verdict = ifelse(sig, "drop: associated with X", "valid"),
    verdict_fdr = ifelse(sig_fdr, "drop: associated with X (FDR)", "valid (FDR)"),
    verdict_magnitude= ifelse(sig_mag, "drop: associated with X (magnitude)", "valid (magnitude)"),
    stringsAsFactors = FALSE
  )
}


#' Test instrument-independence of negative controls (W _|_ G | C)
#'
#' For each control feature, computes the partial correlation with the
#' genetic instrument G after residualising both on observed covariates
#' C, and reports the p-value. Controls significantly associated with G
#' after FDR correction may carry meQTL / allele-specific effects that
#' violate the instrument-independence assumption (A2).
#'
#' @param dat Dataset list from [run_single_iteration()] or
#' [generate_toy_data()].
#' @param fdr_level Target FDR for BH correction. Default 0.10.
#' @param n_cores Number of parallel workers. Default 1 (sequential).
#' Uses \code{parallel::mclapply} on Unix and a PSOCK cluster on
#' Windows.
#'
#' @return A data frame with one row per control feature:
#' `feature`, `partial_r`, `p_value`, `fdr`, `significant`, `verdict`.
#' @export
#' @importFrom stats cor lm residuals complete.cases pt p.adjust
#'
#' @examples
#' dat <- run_single_iteration(n_features = 10, seed = 1)
#' nc_independence_check(dat)
nc_independence_check <- function(dat, fdr_level = 0.10, n_cores = 1) {
  W <- dat$W
  G <- dat$G[, 1]
  cv <- dat$synthetic_data
  p <- ncol(W)
  n <- nrow(W)

  # Residualise G on covariates once
  if (!is.null(cv) && ncol(cv) > 0) {
    g_fit <- tryCatch(lm(G ~ ., data = cv), error = function(e) NULL)
    G_res <- if (is.null(g_fit)) G else residuals(g_fit)
  } else {
    G_res <- G
  }
  df_res <- n - 1 - if (!is.null(cv)) ncol(cv) else 0

  res_list <- .parallel_lapply(seq_len(p), function(f) {
    w <- W[, f]
    if (!is.null(cv) && ncol(cv) > 0) {
      w_fit <- tryCatch(lm(w ~ ., data = cv), error = function(e) NULL)
      w_res <- if (is.null(w_fit)) w else residuals(w_fit)
    } else {
      w_res <- w
    }
    ok <- complete.cases(w_res, G_res)
    if (sum(ok) < 5) return(c(NA_real_, NA_real_))
    r <- cor(w_res[ok], G_res[ok])
    if (is.na(r) || df_res < 2) return(c(NA_real_, NA_real_))
    t_stat <- r * sqrt(df_res) / sqrt(1 - r^2)
    pv <- 2 * pt(-abs(t_stat), df = df_res)
    c(r, pv)
  }, n_cores = n_cores, progress = "NC independence (G)")

  res <- do.call(cbind, res_list)

  partial_r <- res[1, ]
  pvals <- res[2, ]
  fdr <- p.adjust(pvals, method = "BH")
  sig <- !is.na(fdr) & fdr < fdr_level

  data.frame(
    feature = seq_len(p),
    partial_r = partial_r,
    p_value = pvals,
    fdr = fdr,
    significant = sig,
    verdict = ifelse(sig, "drop: associated with G", "valid"),
    stringsAsFactors = FALSE
  )
}


#' Test mediator-instrument independence of negative controls (W _|_ Gm | C)
#'
#' For each control feature, computes the partial correlation with the
#' mediator-specific genetic instrument Gm after residualising both on
#' observed covariates C, and reports the p-value. Controls significantly
#' associated with Gm after FDR correction may carry eQTL / allele-specific
#' effects that violate the mediator-instrument-independence assumption
#' (A2'). This is the mediator-instrument analogue of
#' [nc_independence_check()]: just as the exposure instrument G must be
#' independent of the negative controls, so must the mediator instrument Gm.
#'
#' In the placental eQTL motivating example, this screen tests whether the
#' eQTL instrument (fetal genotype) is associated with the methylation-based
#' negative controls -- a violation would indicate shared genomic structure
#' between the eQTL SNPs and the control CpG sites.
#'
#' @param dat Dataset list from [run_single_iteration()] or
#' [generate_toy_data()], containing `Gm`, `W`, and
#' `synthetic_data`. If `Gm` is absent (i.e. no mediator
#' instrument was generated), the function returns `NULL`
#' with a message.
#' @param fdr_level Target FDR for BH correction. Default 0.10.
#' @param n_cores Number of parallel workers. Default 1 (sequential).
#' Uses \code{parallel::mclapply} on Unix and a PSOCK cluster on
#' Windows.
#'
#' @return A data frame with one row per control feature:
#' `feature`, `partial_r`, `p_value`, `fdr`, `significant`, `verdict`.
#' Returns `NULL` if `dat$Gm` is not present.
#' @export
#' @importFrom stats cor lm residuals complete.cases pt p.adjust
#'
#' @examples
#' dat <- run_single_iteration(n_features = 10, phi = 0.8, seed = 1)
#' nc_independence_check_gm(dat)
nc_independence_check_gm <- function(dat, fdr_level = 0.10, n_cores = 1) {
  if (is.null(dat$Gm)) {
    message("nc_independence_check_gm: dat$Gm is not present (no mediator ",
            "instrument was generated). Set phi > 0 in generate_toy_data() ",
            "or run_single_iteration() to generate Gm.")
    return(NULL)
  }

  W <- dat$W
  Gm <- dat$Gm
  cv <- dat$synthetic_data
  p <- ncol(W)
  n <- nrow(W)

  # Residualise Gm on covariates once
  if (!is.null(cv) && ncol(cv) > 0) {
    gm_fit <- tryCatch(lm(Gm ~ ., data = cv), error = function(e) NULL)
    Gm_res <- if (is.null(gm_fit)) Gm else residuals(gm_fit)
  } else {
    Gm_res <- Gm
  }
  df_res <- n - 1 - if (!is.null(cv)) ncol(cv) else 0

  res_list <- .parallel_lapply(seq_len(p), function(f) {
    w <- W[, f]
    if (!is.null(cv) && ncol(cv) > 0) {
      w_fit <- tryCatch(lm(w ~ ., data = cv), error = function(e) NULL)
      w_res <- if (is.null(w_fit)) w else residuals(w_fit)
    } else {
      w_res <- w
    }
    ok <- complete.cases(w_res, Gm_res)
    if (sum(ok) < 5) return(c(NA_real_, NA_real_))
    r <- cor(w_res[ok], Gm_res[ok])
    if (is.na(r) || df_res < 2) return(c(NA_real_, NA_real_))
    t_stat <- r * sqrt(df_res) / sqrt(1 - r^2)
    pv <- 2 * pt(-abs(t_stat), df = df_res)
    c(r, pv)
  }, n_cores = n_cores, progress = "NC independence (Gm)")

  res <- do.call(cbind, res_list)

  partial_r <- res[1, ]
  pvals <- res[2, ]
  fdr <- p.adjust(pvals, method = "BH")
  sig <- !is.na(fdr) & fdr < fdr_level

  data.frame(
    feature = seq_len(p),
    partial_r = partial_r,
    p_value = pvals,
    fdr = fdr,
    significant = sig,
    verdict = ifelse(sig, "drop: associated with Gm", "valid"),
    stringsAsFactors = FALSE
  )
}


#' Check negative-control completeness (dimensional + covariance-capture)
#'
#' Reports two components of the proximal-inference completeness condition:
#'
#' (1) **Dimensional** (the legacy count-based check): the number of valid
#' negative-control features vs the number of latent confounders k.
#' Bridge-function estimators (COCA, PGC) require at least as many
#' valid controls as confounders (Miao, Geng & Tchetgen Tchetgen,
#' 2018). When `dim(W_valid) < k`, no estimator built on those
#' controls can recover the causal effect, regardless of sample size.
#'
#' (2) **Covariance-capture**: whether the controls actually
#' *capture the confounder covariance* — i.e., whether adding W
#' reduces U's contribution to the outcome toward zero, operationalized
#' as the incremental R^2 of W for the outcome above covariates alone,
#' with a permutation null. This addresses the concern that
#' completeness is about covariance captured, not proxy number.
#'
#' The composite `completeness` verdict requires the dimensional component
#' to pass (satisfied/borderline) AND the capture component to be non-
#' negligible (strong/weak). When the dimensional component passes but
#' capture is negligible, the verdict is "weak-capture".
#'
#' @param dat Dataset list from [run_single_iteration()] or
#' [generate_toy_data()].
#' @param n_valid_controls Optional override: the number of valid
#' controls known from the study design. If
#' `NULL` (default), the function runs both
#' empirical screens and counts controls that
#' pass both.
#' @param fdr_level FDR level for the empirical screens (used
#' only when `n_valid_controls` is `NULL`).
#' @param n_cores Number of parallel workers for the empirical
#' screens. Default 1 (sequential).
#' @param outcome Outcome block for the capture test: `"Y"` (default)
#' or `"M"` (mediator). When `NULL`, the capture
#' component is skipped and only the dimensional
#' verdict is returned.
#' @param n_perm Number of permutations for the capture-test null.
#' Default 1000.
#' @param capture_thresholds Named list with `strong` and `weak` R^2 cutoffs
#' for the capture verdict. Default
#' `list(strong = 0.3, weak = 0.1)`.
#'
#' @return A list with:
#' `n_valid_controls` (count), `k` (number of confounders), `dim_W`,
#' `dimensional` ("satisfied", "borderline", "under-identified"),
#' `capture` (output of [nc_completeness_capture()], or NULL),
#' `completeness` (composite: "satisfied", "borderline",
#' "under-identified", or "weak-capture"),
#' `screen_X` (A1 screen results, if run), `screen_G` (A2 screen
#' results, if run).
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- run_single_iteration(n_features = 10, n_confounders = 1, seed = 1)
#' nc_completeness_check(dat)
#' }
nc_completeness_check <- function(dat, n_valid_controls = NULL,
                                  fdr_level = 0.10, n_cores = 1,
                                  outcome = "Y", n_perm = 1000,
                                  capture_thresholds = list(strong = 0.3, weak = 0.1)) {
  # Number of latent confounders k. dat$U may be a matrix (k columns),
  # a vector (k = 1), or NULL. Defend against all three.
  k <- if (!is.null(dat$U)) {
    if (is.matrix(dat$U)) ncol(dat$U) else 1L
  } else if (!is.null(dat$conf_XM) || !is.null(dat$conf_MY)) {
    # Path-specific loadings store conf_XM / conf_MY; each is a confounder composite.
    1L
  } else {
    1L
  }

  screen_X <- NULL
  screen_G <- NULL

  if (is.null(n_valid_controls)) {
    screen_X <- nc_validity_screen(dat, fdr_level = fdr_level, n_cores = n_cores)
    screen_G <- nc_independence_check(dat, fdr_level = fdr_level, n_cores = n_cores)
    valid_X <- !screen_X$significant
    valid_G <- !screen_G$significant
    valid_both <- valid_X & valid_G
    n_valid <- sum(valid_both, na.rm = TRUE)
  } else {
    n_valid <- n_valid_controls
  }

  dimensional <- if (n_valid > k) {
    "satisfied"
  } else if (n_valid == k) {
    "borderline"
  } else {
    "under-identified"
  }

  # Covariance-capture component. Skipped when outcome is NULL
  # (preserves default behavior for callers that only want the count check).
  capture <- NULL
  if (!is.null(outcome)) {
    capture <- tryCatch(
      nc_completeness_capture(dat, outcome = outcome, n_perm = n_perm,
                              n_cores = n_cores, thresholds = capture_thresholds),
      error = function(e) {
        message("nc_completeness_check: capture test skipped (", conditionMessage(e), ")")
        NULL
      }
    )
  }

  # Composite verdict.
  completeness <- dimensional
  if (!is.null(capture) && dimensional %in% c("satisfied", "borderline")) {
    cap <- capture$capture_verdict
    if (is.null(cap) || cap == "negligible") {
      completeness <- "weak-capture"
    }
    # else keep dimensional verdict (satisfied/borderline)
  }

  list(
    n_valid_controls = n_valid,
    k = k,
    dim_W = ncol(dat$W),
    dimensional = dimensional,
    capture = capture,
    completeness = completeness,
    screen_X = screen_X,
    screen_G = screen_G
  )
}


#' Covariance-capture completeness test
#'
#' Operationalizes the proximal completeness condition as whether the
#' negative-control panel W *captures the confounder covariance* — i.e.,
#' whether adding W reduces U's contribution to the outcome toward zero.
#'
#' The test computes the incremental R^2 of W for the outcome above
#' covariates C alone: `R^2(W | C) = R^2(Y ~ C + W) - R^2(Y ~ C)`, averaged
#' across outcome features. A permutation null is generated by permuting the
#' W-outcome association (holding the C-outcome association fixed) over
#' `n_perm` permutations, yielding a permutation p-value for "W captures
#' U-signal beyond chance."
#'
#' This addresses the concern that the count-based completeness check
#' (`dim(W_valid) >= k`) is necessary but not sufficient: completeness is
#' about the covariance the controls capture, not the number of proxies. A
#' panel can have many controls that each weakly capture U, or few controls
#' that together capture it well.
#'
#' @param dat Dataset list from [run_single_iteration()] or
#' [generate_toy_data()], containing `W`, an outcome block
#' (`Y` or `M`), and `synthetic_data` (covariates).
#' @param outcome Outcome block to test against: `"Y"` (default) or `"M"`.
#' @param n_perm Number of permutations for the null. Default 1000.
#' @param n_cores Number of parallel workers. Default 1.
#' @param thresholds Named list with `strong` and `weak` R^2 cutoffs for the
#' verdict. Default `list(strong = 0.3, weak = 0.1)`.
#'
#' @return A list with:
#' `capture_R2` (point estimate, averaged across outcome features),
#' `capture_pvalue` (permutation p-value),
#' `capture_verdict` ("strong", "weak", or "negligible"),
#' `null_distribution` (numeric vector of permuted R^2 values),
#' `n_features` (number of outcome features tested).
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- run_single_iteration(n_features = 10, n_confounders = 1, seed = 1)
#' nc_completeness_capture(dat)
#' }
nc_completeness_capture <- function(dat, outcome = "Y", n_perm = 1000,
                                    n_cores = 1,
                                    thresholds = list(strong = 0.3, weak = 0.1)) {
  W <- dat$W
  if (is.null(W)) stop("dat$W is not present.")
  cv <- dat$synthetic_data
  n <- nrow(W)
  pW <- ncol(W)

  Y_block <- switch(outcome,
    "Y" = dat$Y, "M" = dat$M,
    stop("outcome must be 'Y' or 'M'"))
  if (is.null(Y_block)) stop("dat$", outcome, " is not present.")
  # Normalize to (samples x features) to match W. The package stores both
  # Y and W as (samples x features) by convention, but defend against either
  # orientation by aligning the sample dimension to nrow(W).
  if (nrow(Y_block) != n) {
    if (ncol(Y_block) == n) {
      Y_block <- t(Y_block)
    } else {
      stop("dat$", outcome, " sample dimension (", nrow(Y_block), " or ",
           ncol(Y_block), ") does not match n (", n, ").")
    }
  }
  pY <- ncol(Y_block)

  cnames <- if (!is.null(cv)) names(cv) else character(0)
  cs <- if (length(cnames)) paste0(" + ", paste(cnames, collapse = " + ")) else ""
  w_str <- paste0(" + ", paste0("W", seq_len(pW), collapse = " + "))

  # R^2 of a fit: 1 - SS_res / SS_tot
  r2_of <- function(fit, y) {
    ss_res <- sum(residuals(fit)^2, na.rm = TRUE)
    ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
    if (ss_tot <= 0) return(NA_real_)
    1 - ss_res / ss_tot
  }

  # Per-feature incremental R^2 of W above C. Features are columns of Y_block.
  inc_r2_feature <- function(f) {
    y <- Y_block[, f]
    d_base <- if (!is.null(cv)) cbind(data.frame(y = y), cv) else data.frame(y = y)
    fit_base <- tryCatch(lm(as.formula(paste0("y ~ 1", cs)), data = d_base),
                         error = function(e) NULL)
    if (is.null(fit_base)) return(NA_real_)
    d_full <- d_base
    for (j in seq_len(pW)) d_full[[paste0("W", j)]] <- W[, j]
    fit_full <- tryCatch(lm(as.formula(paste0("y ~ 1", cs, w_str)), data = d_full),
                         error = function(e) NULL)
    if (is.null(fit_full)) return(NA_real_)
    r2_full <- r2_of(fit_full, y)
    r2_base <- r2_of(fit_base, y)
    r2_full - r2_base
  }

  obs_r2 <- vapply(seq_len(pY), inc_r2_feature, numeric(1))
  obs_r2_mean <- mean(obs_r2, na.rm = TRUE)

  # Permutation null: permute the W-outcome association (shuffle sample order
  # of W) holding the C-outcome association fixed. Parallelized across perms.
  perm_r2_one <- function(seed) {
    # Permute rows of W (breaks W-U/Y link, preserves W cross-feature corr).
    idx <- sample.int(n)
    Wp <- W[idx, , drop = FALSE]
    vals <- vapply(seq_len(pY), function(f) {
      y <- Y_block[, f]
      d_full <- if (!is.null(cv)) cbind(data.frame(y = y), cv) else data.frame(y = y)
      for (j in seq_len(pW)) d_full[[paste0("W", j)]] <- Wp[, j]
      fit_full <- tryCatch(lm(as.formula(paste0("y ~ 1", cs, w_str)), data = d_full),
                           error = function(e) NULL)
      if (is.null(fit_full)) return(NA_real_)
      d_base <- if (!is.null(cv)) cbind(data.frame(y = y), cv) else data.frame(y = y)
      fit_base <- tryCatch(lm(as.formula(paste0("y ~ 1", cs)), data = d_base),
                           error = function(e) NULL)
      if (is.null(fit_base)) return(NA_real_)
      r2_of(fit_full, y) - r2_of(fit_base, y)
    }, numeric(1))
    mean(vals, na.rm = TRUE)
  }

  null_dist <- .parallel_lapply(seq_len(n_perm), function(i) perm_r2_one(i),
                                n_cores = n_cores, progress = "NC capture null")
  null_dist <- unlist(null_dist)

  # Permutation p-value: fraction of null >= observed (one-sided).
  pval <- mean(null_dist >= obs_r2_mean, na.rm = TRUE)

  strong_cut <- if (!is.null(thresholds$strong)) thresholds$strong else 0.3
  weak_cut <- if (!is.null(thresholds$weak)) thresholds$weak else 0.1
  verdict <- if (obs_r2_mean > strong_cut && pval < 0.05) {
    "strong"
  } else if (obs_r2_mean >= weak_cut) {
    "weak"
  } else {
    "negligible"
  }

  list(
    capture_R2 = obs_r2_mean,
    capture_pvalue = pval,
    capture_verdict = verdict,
    null_distribution = null_dist,
    n_features = pY
  )
}


#' Negative-control support/range check
#'
#' Diagnostic for whether the negative-control panel captures the full
#' support of the confounder, or only part of it. Uses a confounder proxy
#' \code{U_tilde = resid(X ~ G + C)} (the instrument-purged exposure
#' residual, which carries the confounding signal) and asks how much of
#' \code{U_tilde} the NC panel explains, and whether each individual NC
#' covers a distinct share of that signal.
#'
#' A panel can pass the count-based completeness check
#' (\code{dim(W_valid) >= k}) and the covariance-capture test while still
#' covering only part of the confounder support — for example when every
#' control loads on the same single confounder direction. This diagnostic
#' reports the multivariate \code{R^2(U_tilde | W)} (how much of the
#' confounder proxy the panel explains) and a per-NC \code{support_ratio}
#' (the partial correlation of each control with \code{U_tilde} given the
#' other controls), which flags controls that add no unique coverage.
#'
#' This is a diagnostic, not a gate: \code{U_tilde} is itself an imperfect
#' proxy (it mixes the confounder with exposure noise), so the values are
#' interpretable only comparatively across controls and panels.
#'
#' @param dat Dataset list from [run_single_iteration()],
#' [generate_toy_data()], or the \code{.to_nc_dat()} bridge in
#' [iconic_diagnose()], containing \code{W}, \code{X}, and
#' \code{synthetic_data} (covariates). \code{G} is used when present.
#' @param fdr_level FDR level for flagging controls whose unique
#' contribution is indistinguishable from zero. Default 0.10.
#'
#' @return A list with:
#' \code{R2_utilde_given_W} (multivariate R^2 of the confounder proxy on
#' the full NC panel),
#' \code{support} (data frame with per-NC \code{support_ratio} = partial
#' correlation of the control with \code{U_tilde} given the other
#' controls, \code{p_value}, and \code{adds_coverage} flag),
#' \code{n_controls}, and \code{verdict} ("broad" when
#' \code{R2_utilde_given_W >= 0.5}, "partial" when >= 0.2, else "narrow").
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- run_single_iteration(n_features = 10, n_confounders = 1, seed = 1)
#' nc_support_check(dat)
#' }
nc_support_check <- function(dat, fdr_level = 0.10) {
  W <- dat$W
  if (is.null(W)) stop("dat$W is not present.")
  if (!is.matrix(W)) W <- as.matrix(W)
  n <- nrow(W)
  pW <- ncol(W)
  X <- dat$X
  if (is.null(X)) stop("dat$X is not present.")
  cv <- dat$synthetic_data
  G <- dat$G
  if (!is.null(G) && is.matrix(G)) G <- as.numeric(G[, 1])

  cnames <- if (!is.null(cv) && ncol(cv) > 0) names(cv) else character(0)
  cs <- if (length(cnames)) paste0(" + ", paste(cnames, collapse = " + ")) else ""

  # Confounder proxy: resid(X ~ G + C). Without an instrument, resid(X ~ C).
  d_x <- data.frame(X = as.numeric(X))
  if (!is.null(G)) d_x$G <- as.numeric(G)
  if (!is.null(cv)) d_x <- cbind(d_x, cv)
  rhs <- paste0(if (!is.null(G)) "G" else "1", cs)
  fit_x <- tryCatch(lm(as.formula(paste0("X ~ ", rhs)), data = d_x),
                    error = function(e) NULL)
  u_tilde <- if (!is.null(fit_x)) residuals(fit_x) else as.numeric(scale(X))

  # Multivariate R^2 of u_tilde on the full W panel (+ C)
  d_full <- data.frame(u = u_tilde)
  if (!is.null(cv)) d_full <- cbind(d_full, cv)
  for (j in seq_len(pW)) d_full[[paste0("W", j)]] <- W[, j]
  w_str <- paste0(" + ", paste0("W", seq_len(pW), collapse = " + "))
  fit_full <- tryCatch(lm(as.formula(paste0("u ~ 1", cs, w_str)), data = d_full),
                       error = function(e) NULL)
  fit_base <- tryCatch(lm(as.formula(paste0("u ~ 1", cs)), data = d_full),
                       error = function(e) NULL)
  r2_of <- function(fit, y) {
    if (is.null(fit)) return(NA_real_)
    ss_res <- sum(residuals(fit)^2, na.rm = TRUE)
    ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
    if (ss_tot <= 0) return(NA_real_)
    1 - ss_res / ss_tot
  }
  R2_full <- r2_of(fit_full, u_tilde)
  R2_base <- r2_of(fit_base, u_tilde)
  R2_utilde_W <- if (!is.na(R2_full) && !is.na(R2_base)) R2_full - R2_base else NA_real_

  # Per-NC unique contribution: partial cor of W_j with u_tilde given W_-j + C.
  # Computed as the drop in R^2 when control j is removed (unique R^2 share),
  # signed by the marginal correlation.
  support_rows <- lapply(seq_len(pW), function(j) {
    others <- setdiff(seq_len(pW), j)
    d_red <- data.frame(u = u_tilde)
    if (!is.null(cv)) d_red <- cbind(d_red, cv)
    for (jj in others) d_red[[paste0("W", jj)]] <- W[, jj]
    w_red <- if (length(others)) paste0(" + ", paste0("W", others, collapse = " + ")) else ""
    fit_red <- tryCatch(lm(as.formula(paste0("u ~ 1", cs, w_red)), data = d_red),
                        error = function(e) NULL)
    R2_red <- r2_of(fit_red, u_tilde)
    unique_r2 <- if (!is.na(R2_full) && !is.na(R2_red)) max(R2_full - R2_red, 0) else NA_real_
    # signed support ratio: sign by marginal cor, magnitude = sqrt(unique R2)
    marg_cor <- stats::cor(W[, j], u_tilde, use = "pairwise.complete.obs")
    sr <- sign(marg_cor) * sqrt(unique_r2)
    # p-value for the added variable: F-test comparing reduced vs full
    pval <- NA_real_
    if (!is.null(fit_red) && !is.null(fit_full)) {
      an <- tryCatch(stats::anova(fit_red, fit_full), error = function(e) NULL)
      if (!is.null(an) && nrow(an) >= 2) pval <- an$`Pr(>F)`[2]
    }
    data.frame(control = j, support_ratio = sr, unique_R2 = unique_r2,
               p_value = pval, stringsAsFactors = FALSE)
  })
  support <- do.call(rbind, support_rows)
  support$p_adj <- stats::p.adjust(support$p_value, method = "BH")
  support$adds_coverage <- !is.na(support$p_adj) & support$p_adj < fdr_level

  verdict <- if (is.na(R2_utilde_W)) "unknown" else if (R2_utilde_W >= 0.5) "broad" else if (R2_utilde_W >= 0.2) "partial" else "narrow"

  list(
    R2_utilde_given_W = R2_utilde_W,
    support = support,
    n_controls = pW,
    verdict = verdict
  )
}
