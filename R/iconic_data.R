# ============================================================
# iconic_data: Standardized data interface for the model
# selection workflow.
#
# iconic_data() – construct a validated S3 object from the
# user's real data (vectors or matrices).
# as_iconic_data() – convert a load_real_input_data() result.
# validate_iconic_data() – internal validation.
# print.iconic_data() – human-readable summary.
#
# This is the estimation-facing data constructor. It is distinct
# from load_real_input_data(), which is designed for GAN training
# (features x samples matrices -> tidy sample x variable frame).
# iconic_data() preserves the matrix structure so that estimation
# can loop over features and mediators.
# ============================================================

#' Construct a standardized data object for ICONIC model selection
#'
#' Creates an `iconic_data` S3 object from the user's real data,
#' standardizing vectors and matrices into a consistent format that
#' all downstream model selection functions consume.
#'
#' @param X Exposure: numeric vector (length n) or features x samples
#' matrix. If a matrix, column means are taken and scaled (one exposure
#' per sample).
#' @param Y Outcome: numeric vector (length n) or features x samples
#' matrix. When a matrix, estimation runs per-feature. When
#' \code{outcome_type = "binary"}, \code{Y} is the 0/1 outcome vector
#' (length n).
#' @param M Optional mediator: numeric vector (length n) or features x
#' samples matrix. When a matrix, estimation runs per-mediator x
#' per-outcome.
#' @param G Optional exposure instrument: numeric vector (length n) or
#' n x n_features matrix (as returned by [generate_toy_data()]). If a matrix,
#' the first column is extracted. E.g., a polygenic risk score.
#' @param Gm Optional mediator instrument: numeric vector (length n) or
#' matrix (n x n_mediators). When a matrix, each column is the instrument
#' for the corresponding mediator (e.g., per-isoform cis-eQTLs).
#' @param W Optional negative-control panel: features x samples matrix.
#' Single-panel NCs used for COCA, PGC.
#' @param W1 Optional path-specific NCs for the X->M path: features x
#' samples matrix. Captures conf_XM. When W1/W2 are absent but W is present,
#' W1 = W2 = W.
#' @param W2 Optional path-specific NCs for the M->Y path: features x
#' samples matrix. Captures conf_MY.
#' @param covariates Optional data frame of sample-level covariates (n rows).
#' Recognized columns `sex`, `GA`, `mother_ethnicity` are encoded; other
#' numeric columns are z-scored. Names colliding with estimator-reserved
#' tokens are renamed.
#' @param feature_names Optional character vector of outcome feature names.
#' @param mediator_names Optional character vector of mediator names.
#' @param trained_gan Optional \code{iconic_gan} from
#' \code{\link{train_gan_on_real_data}()}. When supplied, the GAN is
#' attached to the data object and reused by
#' \code{\link{iconic_sensitivity}()} and \code{\link{iconic_prospect}()}
#' instead of auto-training a new one. This avoids retraining when the
#' same data is used across multiple workflow steps.
#' @param scale Logical: center and scale all continuous inputs (X, Y, M,
#' G, Gm, W, W1, W2, and numeric covariates) to mean 0 / sd 1. Default
#' \code{TRUE}. Scaling parameters are recorded in \code{$scaling} for
#' back-transformation. Set \code{FALSE} to preserve the original scale.
#' @param recycle_lone_panel Logical: when exactly one of \code{W1} /
#' \code{W2} is supplied (no pooled \code{W}), use that lone panel as BOTH
#' path-specific bridges (\code{W1 = W2}), making the two-bridge estimators
#' PGC2 / PGC2Gm eligible. Default \code{FALSE}: a lone panel is retained
#' for IV2SLS2's path-specific augmentation and used to derive the pooled
#' \code{W} (so DIRECT / COCA / PGC and the NC validity screens run), but
#' \code{has_path_nc} stays \code{FALSE} and PGC2 / PGC2Gm remain
#' ineligible. Set to \code{TRUE} only when the single panel is assumed
#' complete for BOTH path confounder composites (the shared-panel special
#' case); IV2SLS2 is the more defensible primary estimator when coverage of
#' the other path's composite is in doubt. A warning is emitted when the
#' recycle is activated.
#' @param outcome_type Character: \code{"continuous"} (default,
#' backward-compatible), \code{"survival"}, or \code{"binary"}. When
#' \code{"survival"},
#' \code{Y} is not required; instead supply \code{surv_time} and
#' \code{surv_event}. Estimation uses Cox proportional-hazards
#' (\code{\link[survival]{coxph}}) or RMST pseudo-observation OLS
#' (see \code{effect_scale} in \code{\link{iconic_estimate}()}). When
#' \code{"binary"}, supply \code{Y} as the 0/1 outcome vector (length n);
#' estimation uses logistic regression (log-OR scale) or a linear
#' probability model (risk-difference scale).
#' @param surv_time Numeric follow-up time vector (length n). Required
#' when \code{outcome_type = "survival"}; ignored otherwise.
#' @param surv_event Numeric 0/1 event indicator (length n; 1 = event
#' observed, 0 = censored). Required when
#' \code{outcome_type = "survival"}; ignored otherwise.
#' @param Z Defunct. Renamed to \code{X}; passing a value errors
#'   with a message pointing to \code{X}. Retained in the signature only to
#'   catch and redirect old calls.
#'
#' @return An `iconic_data` S3 object: a named list with `$X`, `$Y`, `$M`,
#' `$G`, `$Gm`, `$W`, `$W1`, `$W2`, `$covariates`, `$n`, `$n_features`,
#' `$n_mediators`, `$has_instrument`, `$has_mediator_instrument`,
#' `$has_nc`, `$has_path_nc`, `$is_mediation`, `$feature_names`,
#' `$mediator_names`, `$trained_gan`, `$outcome_type`, and (when
#' survival) `$surv_time`, `$surv_event`, or (when binary) `$Y_bin`.
#' @export
#'
#' @examples
#' # Total-effect only (no mediation)
#' data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100))
#'
#' # Full mediation with instruments and NCs
#' data <- iconic_data(
#' X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
#' M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
#' W = matrix(rnorm(100 * 10), 10, 100)
#' )
#' print(data)
#'
#' # Survival outcome
#' data <- iconic_data(
#' X = rnorm(100), outcome_type = "survival",
#' surv_time = rexp(100), surv_event = rbinom(100, 1, 0.6),
#' M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
#' W = matrix(rnorm(100 * 10), 10, 100)
#' )
#' print(data)
#'
#' # Binary outcome (Y is the 0/1 outcome vector)
#' data <- iconic_data(
#' X = rnorm(100), Y = rbinom(100, 1, 0.4), outcome_type = "binary",
#' M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
#' W = matrix(rnorm(100 * 10), 10, 100)
#' )
#' print(data)
iconic_data <- function(X, Y = NULL, M = NULL, G = NULL, Gm = NULL, W = NULL,
                        W1 = NULL, W2 = NULL, covariates = NULL,
                        feature_names = NULL, mediator_names = NULL,
                        trained_gan = NULL,
                        outcome_type = c("continuous", "survival", "binary"),
                        surv_time = NULL, surv_event = NULL,
                        scale = TRUE,
                        recycle_lone_panel = FALSE,
                        Z = NULL) {
  outcome_type <- match.arg(outcome_type)

  ## Deprecated-argument trap: the exposure was renamed Z -> X.
  if (!is.null(Z))
    stop("argument `Z` was renamed to `X`; please use `X = ...`.",
         call. = FALSE)

  # Scaling helper: center and scale a numeric vector, recording the
  # parameters for back-transformation. Returns list(x, center, scale).
  .scale_vec <- function(v) {
    v <- as.numeric(v)
    ctr <- mean(v, na.rm = TRUE)
    s <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) s <- 1
    list(x = as.numeric((v - ctr) / s), center = ctr, scale = s)
  }
  # Scale a features x samples matrix column-wise (per feature).
  .scale_mat <- function(Mm) {
    ctr <- rowMeans(Mm, na.rm = TRUE)
    s <- apply(Mm, 1, stats::sd, na.rm = TRUE)
    s[!is.finite(s) | s == 0] <- 1
    list(x = (Mm - ctr) / s, center = ctr, scale = s)
  }
  scaling <- list(enabled = scale)

  ## --- X: exposure ---
  if (missing(X) || is.null(X)) stop("X (exposure) is required.")
  if (is.matrix(X)) {
    if (ncol(X) == 1) {
      X <- as.numeric(X)
    } else {
      X <- colMeans(X, na.rm = TRUE)
    }
  } else {
    X <- as.numeric(X)
  }
  if (scale) {
    sx <- .scale_vec(X)
    X <- sx$x
    scaling$X <- list(center = sx$center, scale = sx$scale)
  }
  n <- length(X)

  ## --- Y: outcome ---
  if (outcome_type == "survival") {
    # Survival outcome: Y is not required; surv_time + surv_event define
    # the outcome. n_features = 1 (single time-to-event outcome).
    if (is.null(surv_time) || is.null(surv_event))
      stop("When outcome_type = 'survival', surv_time and surv_event ",
           "are required.")
    surv_time <- as.numeric(surv_time)
    surv_event <- as.numeric(surv_event)
    if (length(surv_time) != n || length(surv_event) != n)
      stop("surv_time and surv_event must have length n (=", n, ").")
    if (any(is.na(surv_time)) || any(is.na(surv_event)))
      stop("surv_time and surv_event must not contain NA values.")
    if (any(surv_time <= 0))
      stop("surv_time must be strictly positive.")
    if (!all(surv_event %in% c(0, 1)))
      stop("surv_event must be 0/1 (1 = event observed, 0 = censored).")
    Y <- NULL
    n_features <- 1L
  } else if (outcome_type == "binary") {
    # Binary outcome: Y is the 0/1 outcome vector (length n). Stored as
    # the dedicated Y_bin field; n_features = 1 (single binary outcome).
    # Y_bin is NOT scaled (it is a 0/1 indicator).
    if (missing(Y) || is.null(Y))
      stop("When outcome_type = 'binary', Y must be the 0/1 outcome ",
           "vector (length n).")
    Y_bin <- as.numeric(Y)
    if (length(Y_bin) != n)
      stop("Y must have length n (=", n, ") when outcome_type = 'binary'.")
    if (any(is.na(Y_bin)))
      stop("Y must not contain NA values when outcome_type = 'binary'.")
    if (!all(Y_bin %in% c(0, 1)))
      stop("Y must be dichotomous (0/1) when outcome_type = 'binary'.")
    if (length(unique(Y_bin)) < 2L)
      stop("Y must contain both 0 and 1 when outcome_type = 'binary' ",
           "(the outcome has no variation).")
    Y <- NULL
    n_features <- 1L
  } else {
    if (missing(Y) || is.null(Y)) stop("Y (outcome) is required ",
         "(or set outcome_type = 'survival' with surv_time/surv_event, ",
         "or outcome_type = 'binary' with Y as the 0/1 outcome vector).")
    Y <- as.matrix(Y)
    if (nrow(Y) == n) {
      # User passed samples x features; transpose to features x samples
      Y <- t(Y)
    }
    if (ncol(Y) != n)
      stop("Y must have n samples. If a matrix, pass features x samples ",
           "(features in rows) or samples x features (samples in rows).")
    if (scale) {
      sy <- .scale_mat(Y)
      Y <- sy$x
      scaling$Y <- list(center = sy$center, scale = sy$scale)
    }
    n_features <- nrow(Y)
  }

  ## --- M: mediator ---
  is_mediation <- !is.null(M)
  n_mediators <- NULL
  if (is_mediation) {
    if (is.matrix(M)) {
      if (nrow(M) == n) M <- t(M) # samples x mediators -> mediators x samples
      if (ncol(M) != n)
        stop("M must have n samples.")
      n_mediators <- nrow(M)
    } else {
      M <- matrix(as.numeric(M), nrow = 1, ncol = n)
      n_mediators <- 1L
    }
    if (scale) {
      sm <- .scale_mat(M)
      M <- sm$x
      scaling$M <- list(center = sm$center, scale = sm$scale)
    }
  }

  ## --- G: exposure instrument ---
  has_instrument <- !is.null(G)
  if (has_instrument) {
    if (is.matrix(G)) {
      # G may be an n x n_features matrix (as returned by
      # generate_toy_data()); extract the first column.
      if (ncol(G) == 1) {
        G <- as.numeric(G)
      } else if (nrow(G) == n) {
        G <- as.numeric(G[, 1])
      } else {
        G <- as.numeric(G[, 1])
      }
    } else {
      G <- as.numeric(G)
    }
    if (length(G) != n) stop("G must have length n (or n rows if a matrix).")
    if (scale) {
      sg <- .scale_vec(G)
      G <- sg$x
      scaling$G <- list(center = sg$center, scale = sg$scale)
    }
  }

  ## --- Gm: mediator instrument ---
  has_mediator_instrument <- !is.null(Gm)
  if (has_mediator_instrument) {
    if (is.matrix(Gm)) {
      if (nrow(Gm) == n) Gm <- t(Gm)
      if (ncol(Gm) != n) stop("Gm must have n samples.")
      if (is_mediation && nrow(Gm) != n_mediators)
        warning("Gm has ", nrow(Gm), " rows but M has ", n_mediators,
                " mediators. Using row recycling.")
    } else {
      Gm <- matrix(as.numeric(Gm), nrow = 1, ncol = n)
    }
    if (scale) {
      sgm <- .scale_mat(Gm)
      Gm <- sgm$x
      scaling$Gm <- list(center = sgm$center, scale = sgm$scale)
    }
  }

  ## --- W: negative controls (single panel) ---
  has_nc <- !is.null(W)
  if (has_nc) {
    W <- as.matrix(W)
    if (nrow(W) == n) W <- t(W)
    if (ncol(W) != n) stop("W must have n samples.")
    if (nrow(W) != n_features && outcome_type == "continuous")
      warning("W has ", nrow(W), " features but Y has ", n_features,
              ". Using row recycling.")
    if (scale) {
      sw <- .scale_mat(W)
      W <- sw$x
      scaling$W <- list(center = sw$center, scale = sw$scale)
    }
  }

  ## --- W1, W2: path-specific NCs ---
  # has_path_nc (both panels present) gates the two-bridge estimators
  # PGC2 / PGC2Gm, which need W1 AND W2. A lone W1 or W2 panel is still
  # retained (below) so that IV2SLS2 can use it for path-specific
  # augmentation of the corresponding stage(s).
  recycled_lone_panel <- FALSE
  has_path_nc <- !is.null(W1) && !is.null(W2)
  if (!has_path_nc && has_nc) {
    # Backward-compatible: single-panel NCs used for both paths
    W1 <- W
    W2 <- W
    has_path_nc <- TRUE
  }
  if (has_path_nc) {
    W1 <- as.matrix(W1)
    W2 <- as.matrix(W2)
    if (nrow(W1) == n) W1 <- t(W1)
    if (nrow(W2) == n) W2 <- t(W2)
    if (ncol(W1) != n) stop("W1 must have n samples.")
    if (ncol(W2) != n) stop("W2 must have n samples.")
    # If W was not supplied but W1/W2 were, derive a combined W panel
    # so that single-panel estimators (COCA, PGC, IV2SLS) and the
    # instrument-strength check have a W matrix to work with.
    if (!has_nc) {
      W <- if (identical(W1, W2)) W1 else (W1 + W2) / 2
      has_nc <- TRUE
    }
    if (scale) {
      sw1 <- .scale_mat(W1)
      W1 <- sw1$x
      scaling$W1 <- list(center = sw1$center, scale = sw1$scale)
      sw2 <- .scale_mat(W2)
      W2 <- sw2$x
      scaling$W2 <- list(center = sw2$center, scale = sw2$scale)
      # re-derive combined W from scaled W1/W2 when it was derived above
      if (!is.null(scaling$W)) {
        # W was supplied and already scaled; leave it
      } else {
        W <- if (identical(W1, W2)) W1 else (W1 + W2) / 2
      }
    }
  } else if (!is.null(W1) || !is.null(W2)) {
    # Lone-panel case (exactly one of W1 / W2 supplied, no pooled W):
    # retain the supplied panel for IV2SLS2's path-specific augmentation.
    if (!is.null(W1)) {
      W1 <- as.matrix(W1)
      if (nrow(W1) == n) W1 <- t(W1)
      if (ncol(W1) != n) stop("W1 must have n samples.")
      if (scale) {
        sw1 <- .scale_mat(W1)
        W1 <- sw1$x
        scaling$W1 <- list(center = sw1$center, scale = sw1$scale)
      }
    }
    if (!is.null(W2)) {
      W2 <- as.matrix(W2)
      if (nrow(W2) == n) W2 <- t(W2)
      if (ncol(W2) != n) stop("W2 must have n samples.")
      if (scale) {
        sw2 <- .scale_mat(W2)
        W2 <- sw2$x
        scaling$W2 <- list(center = sw2$center, scale = sw2$scale)
      }
    }
    # Derive the pooled single panel from the lone path-specific panel so
    # the single-panel estimators (DIRECT, COCA, PGC) and the NC validity /
    # completeness screens have a W matrix to work with. This mirrors the
    # W1+W2 branch above, which derives W <- (W1 + W2) / 2.
    if (!has_nc) {
      W <- if (!is.null(W2)) W2 else W1
      has_nc <- TRUE
      scaling$W <- if (!is.null(W2)) scaling$W2 else scaling$W1
    }
    # Opt-in recycle: use the lone panel as BOTH path-specific bridges so
    # the two-bridge estimators (PGC2 / PGC2Gm) become eligible. This is the
    # "shared panel" special case (W1 = W2); it assumes the single panel is
    # complete for BOTH path confounder composites. Off by default because
    # IV2SLS2 is the more defensible primary estimator when coverage of the
    # other path's composite is in doubt.
    if (isTRUE(recycle_lone_panel) && !has_path_nc) {
      lone_label <- if (is.null(W1)) "W2" else "W1"
      if (is.null(W1)) {
        W1 <- W2
        scaling$W1 <- scaling$W2
      } else {
        W2 <- W1
        scaling$W2 <- scaling$W1
      }
      has_path_nc <- TRUE
      recycled_lone_panel <- TRUE
      warning("recycle_lone_panel = TRUE: using the lone ", lone_label,
              " panel as both the X->M (W1) and M->Y (W2) bridge. ",
              "PGC2/PGC2Gm then assume one panel is complete for BOTH path ",
              "confounder composites; IV2SLS2 is more defensible when ",
              "coverage of the X->M composite is in doubt.",
              call. = FALSE)
    }
  } else {
    W1 <- NULL
    W2 <- NULL
  }

  ## --- Covariates ---
  if (!is.null(covariates)) {
    covariates <- .encode_covariates(as.data.frame(covariates), n,
                                     paste0("S", seq_len(n)))
    if (scale && ncol(covariates) > 0) {
      num_cols <- vapply(covariates, is.numeric, logical(1))
      if (any(num_cols)) {
        cov_num <- as.matrix(covariates[, num_cols, drop = FALSE])
        sc <- .scale_mat(t(cov_num)) # scale per covariate (rows)
        covariates[, num_cols] <- as.data.frame(t(sc$x))
        scaling$covariates <- list(center = sc$center, scale = sc$scale,
                                   columns = names(covariates)[num_cols])
      }
    }
  } else {
    covariates <- data.frame(row.names = seq_len(n))[, 0, drop = FALSE]
  }

  ## --- Feature / mediator names ---
  if (is.null(feature_names)) {
    feature_names <- if (outcome_type == "survival") "survival"
                     else if (outcome_type == "binary") "binary"
                     else if (!is.null(rownames(Y))) rownames(Y)
                     else paste0("feature_", seq_len(n_features))
  }
  if (is_mediation && is.null(mediator_names)) {
    mediator_names <- if (!is.null(rownames(M))) rownames(M)
                      else paste0("mediator_", seq_len(n_mediators))
  }

  ## --- Validate ---
  obj <- list(
    X = X,
    Y = Y,
    M = if (is_mediation) M else NULL,
    G = if (has_instrument) G else NULL,
    Gm = if (has_mediator_instrument) Gm else NULL,
    W = if (has_nc) W else NULL,
    W1 = if (!is.null(W1)) W1 else NULL,
    W2 = if (!is.null(W2)) W2 else NULL,
    covariates = covariates,
    n = n,
    n_features = n_features,
    n_mediators = n_mediators,
    has_instrument = has_instrument,
    has_mediator_instrument = has_mediator_instrument,
    has_nc = has_nc,
    has_path_nc = has_path_nc,
    recycled_lone_panel = recycled_lone_panel,
    is_mediation = is_mediation,
    feature_names = feature_names,
    mediator_names = mediator_names,
    trained_gan = trained_gan,
    outcome_type = outcome_type,
    surv_time = if (outcome_type == "survival") surv_time else NULL,
    surv_event = if (outcome_type == "survival") surv_event else NULL,
    Y_bin = if (outcome_type == "binary") Y_bin else NULL,
    scaling = scaling
  )
  validate_iconic_data(obj)
  class(obj) <- c("iconic_data", "list")
  obj
}


#' Validate an iconic_data object (internal)
#'
#' Checks structural integrity: consistent n, minimum sample size,
#' required fields present.
#' @param obj An iconic_data list (pre-class assignment).
#' @keywords internal
#' @noRd
validate_iconic_data <- function(obj) {
  n <- obj$n
  if (n < 20)
    stop("At least 20 samples are required (got ", n, ").")
  if (obj$n_features < 1)
    stop("At least 1 outcome feature is required.")
  if (obj$is_mediation && !is.null(obj$n_mediators) && obj$n_mediators < 1)
    stop("At least 1 mediator is required when M is supplied.")
  if (obj$is_mediation && obj$has_mediator_instrument) {
    if (ncol(obj$Gm) != n)
      warning("Gm must have n samples (columns).")
    if (is.matrix(obj$Gm) && nrow(obj$Gm) != obj$n_mediators)
      warning("Gm has ", nrow(obj$Gm), " rows but M has ",
              obj$n_mediators, " mediators.")
  }
  invisible(TRUE)
}


#' Convert external data containers to iconic_data
#'
#' S3 generic bridging other data interfaces to the estimation
#' interface. Methods:
#'
#' * `default`: a list returned by [load_real_input_data()], an
#'   existing `iconic_data` object (returned as-is), or the exposure
#'   vector `X` with named arguments matching [iconic_data()] (Y, M,
#'   G, Gm, W, W1, W2, covariates, feature_names, mediator_names),
#'   which delegates to [iconic_data()].
#' * `SummarizedExperiment`: extracts the outcome panel from an assay
#'   and sample-level fields (exposure, instruments, negative
#'   controls, covariates) from `colData`. Requires the
#'   SummarizedExperiment package (listed under `Suggests`).
#'
#' @param input An object to convert: a list returned by
#' [load_real_input_data()], an `iconic_data` object, an exposure
#' vector (named-argument form), or a SummarizedExperiment.
#' @param ... Named arguments passed to [iconic_data()] when using
#' the named-argument form, or to the method.
#'
#' @return An `iconic_data` S3 object.
#' @export
#'
#' @examples
#' # From a load_real_input_data() result
#' input <- load_real_input_data(example = TRUE)
#' data <- as_iconic_data(input)
#' print(data)
#'
#' # From named components (delegates to iconic_data())
#' data <- as_iconic_data(rnorm(100), Y = matrix(rnorm(100*10), 10, 100),
#' G = rnorm(100), Gm = rnorm(100),
#' W = matrix(rnorm(100*10), 10, 100))
as_iconic_data <- function(input, ...) {
  UseMethod("as_iconic_data")
}


#' @rdname as_iconic_data
#' @export
as_iconic_data.default <- function(input, ...) {
  dots <- list(...)

  # Passthrough: if input is already an iconic_data object, return as-is
  if (inherits(input, "iconic_data")) {
    return(input)
  }

  # If input is a list with $original_matrices, use the original interface
  if (is.list(input) && !is.null(input$original_matrices)) {
    om <- input$original_matrices
    iconic_data(
      X = om$X,
      Y = om$Y,
      W = om$W,
      covariates = input$covariates,
      feature_names = input$feature_names
    )
  } else if (length(dots) > 0 || !is.null(input)) {
    # Named-argument form: delegate to iconic_data()
    # input is X, dots are the remaining named args
    do.call(iconic_data, c(list(X = input), dots))
  } else {
    stop("as_iconic_data() requires either a load_real_input_data() result ",
         "or named arguments (X=, Y=, ...).")
  }
}


#' @rdname as_iconic_data
#' @param assay Name or index of the assay of `input` holding the
#' outcome panel (features x samples). Default `1`. Set to `NULL` for
#' survival outcomes (when `surv_time`/`surv_event` are colData
#' columns and no continuous outcome panel is used).
#' @param mediator_assay Optional name or index of an assay holding
#' the mediator panel (mediators x samples).
#' @param exposure Character: name of the `colData` column holding the
#' exposure X.
#' @param instrument Optional character: name of the `colData` column
#' holding the exposure instrument G (e.g. a polygenic score).
#' @param mediator_instrument Optional character vector of `colData`
#' column names holding the mediator instrument(s) Gm (one column per
#' mediator).
#' @param negative_controls Optional character vector of `colData`
#' column names forming the negative-control panel W.
#' @param covariates Optional character vector of `colData` column
#' names to carry through as covariates.
#' @param surv_time,surv_event Optional character: names of `colData`
#' columns holding follow-up time and the 0/1 event indicator; set
#' together with `outcome_type = "survival"`.
#' @param bin_outcome Optional character: name of a `colData` column
#' holding the 0/1 outcome; set together with
#' `outcome_type = "binary"`.
#' @param outcome_type `\"continuous\"` (default),
#' `\"survival\"`, or `\"binary\"`.
#' @method as_iconic_data SummarizedExperiment
#' @export
#'
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'   se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(expr = matrix(rnorm(20 * 60), 20, 60,
#'                                 dimnames = list(paste0("gene", 1:20),
#'                                                 paste0("S", 1:60)))),
#'     colData = S4Vectors::DataFrame(
#'       bmi = rnorm(60), prs = rnorm(60),
#'       nc1 = rnorm(60), nc2 = rnorm(60), age = rnorm(60))
#'   )
#'   data <- as_iconic_data(se, assay = "expr", exposure = "bmi",
#'                          instrument = "prs",
#'                          negative_controls = c("nc1", "nc2"),
#'                          covariates = "age")
#'   print(data)
#' }
as_iconic_data.SummarizedExperiment <- function(input, assay = 1,
                                                mediator_assay = NULL,
                                                exposure,
                                                instrument = NULL,
                                                mediator_instrument = NULL,
                                                negative_controls = NULL,
                                                covariates = NULL,
                                                surv_time = NULL,
                                                surv_event = NULL,
                                                bin_outcome = NULL,
                                                outcome_type = c("continuous", "survival", "binary"),
                                                ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("The SummarizedExperiment method requires the ",
         "SummarizedExperiment package. Install with: ",
         "BiocManager::install('SummarizedExperiment')")
  outcome_type <- match.arg(outcome_type)
  cd <- SummarizedExperiment::colData(input)
  cd_names <- colnames(cd)

  get_col <- function(name, label) {
    if (!name %in% cd_names)
      stop(label, " column '", name, "' not found in colData. ",
           "Available: ", paste(cd_names, collapse = ", "))
    v <- cd[[name]]
    if (inherits(v, c("DataFrame", "DFrame")) && ncol(v) == 1) v <- v[[1]]
    as.numeric(v)
  }

  X <- get_col(exposure, "exposure")

  Y <- NULL
  if (outcome_type == "continuous") {
    if (is.null(assay))
      stop("assay = NULL is only valid with outcome_type = 'survival' ",
           "or 'binary'.")
    Y <- as.matrix(SummarizedExperiment::assay(input, assay))
  } else if (outcome_type == "binary") {
    if (is.null(bin_outcome))
      stop("When outcome_type = 'binary', supply bin_outcome = '<colData ",
           "column>' naming the 0/1 outcome column.")
    Y <- get_col(bin_outcome, "bin_outcome")
  }

  M <- NULL
  if (!is.null(mediator_assay))
    M <- as.matrix(SummarizedExperiment::assay(input, mediator_assay))

  G <- if (!is.null(instrument)) get_col(instrument, "instrument") else NULL

  Gm <- NULL
  if (!is.null(mediator_instrument)) {
    gm_cols <- vapply(mediator_instrument, get_col, numeric(length(X)),
                      label = "mediator_instrument")
    Gm <- t(gm_cols)   # mediators x samples
    rownames(Gm) <- mediator_instrument
  }

  W <- NULL
  if (!is.null(negative_controls)) {
    w_cols <- vapply(negative_controls, get_col, numeric(length(X)),
                     label = "negative_controls")
    W <- t(w_cols)     # NCs x samples
    rownames(W) <- negative_controls
  }

  cv <- NULL
  if (!is.null(covariates))
    cv <- as.data.frame(cd[, covariates, drop = FALSE])

  st <- if (!is.null(surv_time)) get_col(surv_time, "surv_time") else NULL
  se <- if (!is.null(surv_event)) get_col(surv_event, "surv_event") else NULL

  iconic_data(X = X, Y = Y, M = M, G = G, Gm = Gm, W = W,
              covariates = cv,
              feature_names = if (!is.null(Y)) rownames(Y),
              mediator_names = if (!is.null(M)) rownames(M),
              outcome_type = outcome_type,
              surv_time = st, surv_event = se, ...)
}


#' Print method for iconic_data objects
#'
#' @param x An `iconic_data` object.
#' @param ... Unused.
#' @return Invisibly returns `x` (the `iconic_data` object); called for its side effect of printing a human-readable summary.
#' @export
#' @examples
#' data <- iconic_data(X = rnorm(50), Y = matrix(rnorm(50 * 5), 5, 50),
#'   G = rnorm(50), W = matrix(rnorm(50 * 5), 5, 50))
#' print(data)
print.iconic_data <- function(x, ...) {
  if (x$outcome_type == "survival") {
    n_events <- sum(x$surv_event)
    cat("<iconic_data> ", x$n, " samples, survival outcome (",
        n_events, " events / ", x$n, ")", sep = "")
  } else if (x$outcome_type == "binary") {
    n_cases <- sum(x$Y_bin)
    cat("<iconic_data> ", x$n, " samples, binary outcome (",
        n_cases, " cases / ", x$n, ")", sep = "")
  } else {
    cat("<iconic_data> ", x$n, " samples, ", x$n_features, " outcome features",
        sep = "")
  }
  if (x$is_mediation)
    cat(", ", x$n_mediators, " mediator(s)", sep = "")
  cat("\n")

  present <- character(0)
  if (x$has_instrument) present <- c(present, "G (exposure instrument)")
  if (x$has_mediator_instrument) present <- c(present, "Gm (mediator instrument)")
  if (x$has_nc) present <- c(present, "W (negative controls)")
  if (isTRUE(x$recycled_lone_panel)) {
    present <- c(present, "W1/W2 (lone panel recycled to both paths)")
  } else if (x$has_path_nc &&
      !is.null(x$W1) && !identical(x$W1, x$W)) {
    present <- c(present, "W1/W2 (path-specific NCs)")
  } else if (!x$has_path_nc && (!is.null(x$W1) || !is.null(x$W2))) {
    present <- c(present, "W1/W2 (lone path-specific NC panel)")
  }

  if (length(present)) {
    cat(" Available:", paste(present, collapse = ", "), "\n")
  } else {
    cat(" No instruments or negative controls supplied.\n")
  }

  if (ncol(x$covariates) > 0)
    cat(" Covariates:", paste(names(x$covariates), collapse = ", "), "\n")

  if (!is.null(x$trained_gan))
    cat(" GAN: attached (", x$trained_gan$model_type, ")\n", sep = "")

  cat(" Mode:", if (x$is_mediation) "mediation" else "total effect", "\n")
  invisible(x)
}
