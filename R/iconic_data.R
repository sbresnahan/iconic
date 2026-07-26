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
#' @param Z Exposure: numeric vector (length n) or features x samples
#' matrix. If a matrix, column means are taken and scaled (one exposure
#' per sample).
#' @param Y Outcome: numeric vector (length n) or features x samples
#' matrix. When a matrix, estimation runs per-feature.
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
#' @param W1 Optional path-specific NCs for the Z->M path: features x
#' samples matrix. Captures U_XM. When W1/W2 are absent but W is present,
#' W1 = W2 = W.
#' @param W2 Optional path-specific NCs for the M->Y path: features x
#' samples matrix. Captures U_MY.
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
#' @param outcome_type Character: \code{"continuous"} (default,
#' backward-compatible) or \code{"survival"}. When \code{"survival"},
#' \code{Y} is not required; instead supply \code{surv_time} and
#' \code{surv_event}. Estimation uses Cox proportional-hazards
#' (\code{\link[survival]{coxph}}) or RMST pseudo-observation OLS
#' (see \code{effect_scale} in \code{\link{iconic_estimate}()}).
#' @param surv_time Numeric follow-up time vector (length n). Required
#' when \code{outcome_type = "survival"}; ignored otherwise.
#' @param surv_event Numeric 0/1 event indicator (length n; 1 = event
#' observed, 0 = censored). Required when
#' \code{outcome_type = "survival"}; ignored otherwise.
#'
#' @return An `iconic_data` S3 object: a named list with `$Z`, `$Y`, `$M`,
#' `$G`, `$Gm`, `$W`, `$W1`, `$W2`, `$covariates`, `$n`, `$n_features`,
#' `$n_mediators`, `$has_instrument`, `$has_mediator_instrument`,
#' `$has_nc`, `$has_path_nc`, `$is_mediation`, `$feature_names`,
#' `$mediator_names`, `$trained_gan`, `$outcome_type`, and (when
#' survival) `$surv_time`, `$surv_event`.
#' @export
#'
#' @examples
#' # Total-effect only (no mediation)
#' data <- iconic_data(Z = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100))
#'
#' # Full mediation with instruments and NCs
#' data <- iconic_data(
#' Z = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
#' M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
#' W = matrix(rnorm(100 * 10), 10, 100)
#' )
#' print(data)
#'
#' # Survival outcome
#' data <- iconic_data(
#' Z = rnorm(100), outcome_type = "survival",
#' surv_time = rexp(100), surv_event = rbinom(100, 1, 0.6),
#' M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
#' W = matrix(rnorm(100 * 10), 10, 100)
#' )
#' print(data)
iconic_data <- function(Z, Y = NULL, M = NULL, G = NULL, Gm = NULL, W = NULL,
                        W1 = NULL, W2 = NULL, covariates = NULL,
                        feature_names = NULL, mediator_names = NULL,
                        trained_gan = NULL,
                        outcome_type = c("continuous", "survival"),
                        surv_time = NULL, surv_event = NULL) {
  outcome_type <- match.arg(outcome_type)

  ## --- Z: exposure ---
  if (missing(Z) || is.null(Z)) stop("Z (exposure) is required.")
  if (is.matrix(Z)) {
    if (ncol(Z) == 1) {
      Z <- as.numeric(scale(as.numeric(Z)))
    } else {
      Z <- as.numeric(scale(colMeans(Z, na.rm = TRUE)))
    }
  } else {
    Z <- as.numeric(Z)
  }
  n <- length(Z)

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
  } else {
    if (missing(Y) || is.null(Y)) stop("Y (outcome) is required ",
         "(or set outcome_type = 'survival' with surv_time/surv_event).")
    Y <- as.matrix(Y)
    if (nrow(Y) == n) {
      # User passed samples x features; transpose to features x samples
      Y <- t(Y)
    }
    if (ncol(Y) != n)
      stop("Y must have n samples. If a matrix, pass features x samples ",
           "(features in rows) or samples x features (samples in rows).")
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
  }

  ## --- W1, W2: path-specific NCs ---
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
  } else {
    W1 <- NULL
    W2 <- NULL
  }

  ## --- Covariates ---
  if (!is.null(covariates)) {
    covariates <- .encode_covariates(as.data.frame(covariates), n,
                                     paste0("S", seq_len(n)))
  } else {
    covariates <- data.frame(row.names = seq_len(n))[, 0, drop = FALSE]
  }

  ## --- Feature / mediator names ---
  if (is.null(feature_names)) {
    feature_names <- if (outcome_type == "survival") "survival"
                     else if (!is.null(rownames(Y))) rownames(Y)
                     else paste0("feature_", seq_len(n_features))
  }
  if (is_mediation && is.null(mediator_names)) {
    mediator_names <- if (!is.null(rownames(M))) rownames(M)
                      else paste0("mediator_", seq_len(n_mediators))
  }

  ## --- Validate ---
  obj <- list(
    Z = Z,
    Y = Y,
    M = if (is_mediation) M else NULL,
    G = if (has_instrument) G else NULL,
    Gm = if (has_mediator_instrument) Gm else NULL,
    W = if (has_nc) W else NULL,
    W1 = if (has_path_nc) W1 else NULL,
    W2 = if (has_path_nc) W2 else NULL,
    covariates = covariates,
    n = n,
    n_features = n_features,
    n_mediators = n_mediators,
    has_instrument = has_instrument,
    has_mediator_instrument = has_mediator_instrument,
    has_nc = has_nc,
    has_path_nc = has_path_nc,
    is_mediation = is_mediation,
    feature_names = feature_names,
    mediator_names = mediator_names,
    trained_gan = trained_gan,
    outcome_type = outcome_type,
    surv_time = if (outcome_type == "survival") surv_time else NULL,
    surv_event = if (outcome_type == "survival") surv_event else NULL
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


#' Convert a load_real_input_data result (or named components) to iconic_data
#'
#' Bridges the GAN-training data interface to the estimation interface.
#' Can be called in two ways:
#'
#' 1. With a list returned by [load_real_input_data()] (the original
#' interface).
#' 2. With named arguments matching [iconic_data()] (Z, Y, M, G, Gm,
#' W, W1, W2, covariates, feature_names, mediator_names), which
#' simply delegates to [iconic_data()].
#'
#' @param input Either a list returned by [load_real_input_data()], or
#' the exposure vector Z (when using named-argument form).
#' @param ... Named arguments passed to [iconic_data()] when using
#' the named-argument form. Ignored when \code{input} is a list.
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
  dots <- list(...)

  # Passthrough: if input is already an iconic_data object, return as-is
  if (inherits(input, "iconic_data")) {
    return(input)
  }

  # If input is a list with $original_matrices, use the original interface
  if (is.list(input) && !is.null(input$original_matrices)) {
    om <- input$original_matrices
    iconic_data(
      Z = om$Z,
      Y = om$Y,
      W = om$W,
      covariates = input$covariates,
      feature_names = input$feature_names
    )
  } else if (length(dots) > 0 || !is.null(input)) {
    # Named-argument form: delegate to iconic_data()
    # input is Z, dots are the remaining named args
    do.call(iconic_data, c(list(Z = input), dots))
  } else {
    stop("as_iconic_data() requires either a load_real_input_data() result ",
         "or named arguments (Z=, Y=, ...).")
  }
}


#' Print method for iconic_data objects
#'
#' @param x An `iconic_data` object.
#' @param ... Unused.
#' @return Invisibly returns `x` (the `iconic_data` object); called for its side effect of printing a human-readable summary.
#' @export
print.iconic_data <- function(x, ...) {
  if (x$outcome_type == "survival") {
    n_events <- sum(x$surv_event)
    cat("<iconic_data> ", x$n, " samples, survival outcome (",
        n_events, " events / ", x$n, ")", sep = "")
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
  if (x$has_path_nc &&
      !is.null(x$W1) && !identical(x$W1, x$W))
    present <- c(present, "W1/W2 (path-specific NCs)")

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
