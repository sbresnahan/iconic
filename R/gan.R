# ============================================================
# TODO(v1.0): supplement — move generative pipeline detail to supplement
# (JYH #416). The comprehensive draft documents the GAN texture model
# inline; when condensed to a journal submission, this material moves to
# a supplementary methods section. No code change in v0.9.2.
# ============================================================
# Generative texture model trained on the user's real data.
#
# `train_gan_on_real_data()` learns the joint distribution of the tidy
# sample x variable frame produced by load_real_input_data() -- exposure
# level, outcome level and encoded covariates. If the torch package (and
# libtorch) is available it trains an adversarial network; otherwise it
# falls back to a multivariate-normal fit (MASS::mvrnorm). Either way the
# result exposes the same downstream contract via sample_texture(): draw n
# realistic base rows. The causal structure is layered on separately in
# run_single_iteration(), so the generator only supplies realism, never the
# ground-truth effect.
#
# Binary-column handling: columns whose training values are all in {0, 1}
# (e.g. encoded sex, one-hot ethnicity dummies) are flagged at training
# time. sample_texture() rounds them back to 0/1 after de-normalisation,
# and enforces mutual exclusivity within one-hot groups so the synthetic
# draws respect the categorical structure.
# ============================================================


#' Check whether a working torch installation is available
#'
#' Verifies the package is installed, libtorch is present, and a basic tensor
#' op succeeds. Used to gate the GAN path; a `FALSE` result triggers the
#' multivariate-normal fallback.
#'
#' @return `TRUE` if torch is usable, otherwise `FALSE`.
#' @export
check_torch_setup <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) return(FALSE)
  ok <- tryCatch({
    torch::torch_is_installed() &&
      is.numeric(as.numeric((torch::torch_ones(2) + torch::torch_ones(2))$sum()))
  }, error = function(e) FALSE)
  isTRUE(ok)
}


#' Detect binary (0/1) columns in a data frame (internal)
#'
#' Returns the names of columns whose unique values are a subset of \eqn{{0, 1}}.
#' Used to flag columns that must be rounded back to 0/1 after sampling.
#' @keywords internal
.detect_binary_cols <- function(X) {
  if (!ncol(X)) return(character(0))
  bin <- vapply(X, function(col) {
    vals <- unique(as.numeric(col[!is.na(col)]))
    length(vals) > 0 && all(vals %in% c(0, 1))
  }, logical(1))
  names(X)[bin]
}

#' Detect one-hot dummy groups from column names (internal)
#'
#' Columns sharing a common prefix followed by an underscore (e.g.
#' `mother_ethnicity_indian`, `mother_ethnicity_malay`) are grouped so that
#' mutual exclusivity can be enforced after rounding. Only groups with >= 2
#' members are returned.
#' @keywords internal
.detect_onehot_groups <- function(columns) {
  if (length(columns) < 2) return(list())
  parts <- strsplit(columns, "_")
  # A one-hot group shares everything except the last underscore-delimited
  # segment. Build prefix = paste(all but last segment, collapse = "_").
  prefixes <- vapply(parts, function(p) {
    if (length(p) <= 1) "" else paste(p[-length(p)], collapse = "_")
  }, character(1))
  groups <- split(columns, prefixes)
  groups <- groups[names(groups) != "" & lengths(groups) >= 2]
  unname(groups)
}

#' Round binary columns and enforce one-hot mutual exclusivity (internal)
#'
#' After de-normalisation, columns flagged as binary are rounded to the nearest
#' of \eqn{{0, 1}}. Within each one-hot group, the column with the highest
#' pre-rounding value wins (set to 1, others to 0), preserving the constraint
#' that exactly one level is active per row.
#' @keywords internal
.enforce_discrete <- function(df, binary_cols, onehot_groups) {
  if (!length(binary_cols)) return(df)

  # First, enforce one-hot mutual exclusivity using the raw (pre-round) values:
  # within each group, the column with the largest value is the "winner".
  for (grp in onehot_groups) {
    present <- intersect(grp, names(df))
    if (length(present) < 2) next
    mat <- as.matrix(df[, present, drop = FALSE])
    winner <- max.col(mat, ties.method = "first")
    for (j in seq_along(present))
      df[[present[j]]] <- as.numeric(winner == j)
  }

  # Round any remaining binary columns not part of a one-hot group.
  remaining <- setdiff(binary_cols, unlist(onehot_groups, use.names = FALSE))
  for (cl in remaining)
    if (cl %in% names(df))
      df[[cl]] <- as.numeric(round(pmin(pmax(df[[cl]], 0), 1)))

  df
}


#' Build the GAN generator network (internal to the torch path)
#'
#' MLP mapping a noise vector to a synthetic row: `noise_dim -> 256 -> 512 ->
#' 1024 -> output_dim`, ReLU + batch-norm + dropout, linear output.
#'
#' @param output_dim Number of output variables (columns of the training frame).
#' @param noise_dim  Latent noise dimension. Default 100.
#' @return A torch `nn_module` instance. Requires torch.
#' @keywords internal
create_generator <- function(output_dim, noise_dim = 100) {
  if (!check_torch_setup()) stop("torch is not available.")
  gen <- torch::nn_module(
    "iconic_generator",
    initialize = function(noise_dim, output_dim) {
      self$net <- torch::nn_sequential(
        torch::nn_linear(noise_dim, 256), torch::nn_relu(),
        torch::nn_batch_norm1d(256),      torch::nn_dropout(0.2),
        torch::nn_linear(256, 512),       torch::nn_relu(),
        torch::nn_batch_norm1d(512),      torch::nn_dropout(0.2),
        torch::nn_linear(512, 1024),      torch::nn_relu(),
        torch::nn_linear(1024, output_dim))
    },
    forward = function(z) self$net(z))
  gen(noise_dim = noise_dim, output_dim = output_dim)
}

#' Build the GAN discriminator network (internal to the torch path)
#'
#' MLP `input_dim -> 1024 -> 512 -> 256 -> 1` with leaky-ReLU and dropout,
#' returning a real-vs-fake logit (no sigmoid; paired with a BCE-with-logits loss).
#'
#' @param input_dim Number of input variables.
#' @return A torch `nn_module` instance. Requires torch.
#' @keywords internal
create_discriminator <- function(input_dim) {
  if (!check_torch_setup()) stop("torch is not available.")
  disc <- torch::nn_module(
    "iconic_discriminator",
    initialize = function(input_dim) {
      self$net <- torch::nn_sequential(
        torch::nn_linear(input_dim, 1024), torch::nn_leaky_relu(0.2), torch::nn_dropout(0.3),
        torch::nn_linear(1024, 512),       torch::nn_leaky_relu(0.2), torch::nn_dropout(0.3),
        torch::nn_linear(512, 256),        torch::nn_leaky_relu(0.2), torch::nn_dropout(0.3),
        torch::nn_linear(256, 1))
    },
    forward = function(x) self$net(x))
  disc(input_dim = input_dim)
}


#' Train a generative texture model on real data
#'
#' Fits an adversarial network (torch GAN) to the tidy training frame.
#' Requires the \pkg{torch} package; if torch is not available, an error is
#' raised with installation instructions.
#'
#' Binary columns (values in \eqn{{0, 1}}, e.g. encoded sex or one-hot ethnicity
#' dummies) are detected automatically and stored in the returned object.
#' [sample_texture()] rounds them back to 0/1 and enforces one-hot mutual
#' exclusivity, so synthetic draws respect the discrete structure.
#'
#' Feature-level residual correlation matrices (for the Y and W panels)
#' can be attached via the `feature_correlations` argument.  When present,
#' [run_single_iteration()] uses them to inject correlated noise into the
#' simulated outcome and negative-control panels.
#'
#' A feature-level copula texture model for the mediator (M) panel can be
#' attached via the `feature_texture` argument (an
#' `iconic_feature_texture` object from [train_feature_texture()]).  When
#' present, [run_single_iteration()] uses it to draw realistic mediator
#' feature vectors that preserve the marginal distributions and
#' cross-feature correlation structure of the user's mediator panel.
#'
#' @param real_data  Data frame of numeric variables (e.g. the
#'   `gan_training_data` element from [load_real_input_data()]). Complete cases
#'   are used.
#' @param feature_correlations Optional list with elements `Y`, `W`,
#'   each a p x p residual correlation matrix from
#'   [load_real_input_data()].  When supplied, the matrices are stored in the
#'   returned `iconic_gan` object and used by [run_single_iteration()] to
#'   generate correlated noise.  Default `NULL`.
#' @param feature_texture Optional `iconic_feature_texture` object from
#'   [train_feature_texture()].  When supplied, it is stored in the returned
#'   `iconic_gan` object and used by [run_single_iteration()] to inject
#'   realistic mediator texture.  Default `NULL`.
#' @param epochs     GAN training epochs. Default 300.
#' @param batch_size Mini-batch size. Default 32.
#' @param lr         Adam learning rate. Default 2e-4.
#' @param seed       Optional RNG seed (sets both R and torch seeds).
#' @param verbose    Print progress. Default TRUE.
#'
#' @return An `iconic_gan` object: a list with `model_type` (`"gan"`),
#'   `columns`, `norm` (per-column centre/scale), `binary_cols` (names of
#'   0/1 columns), `onehot_groups` (list of one-hot dummy groups),
#'   `feature_correlations` (list or NULL), `feature_texture` (object or
#'   NULL), training statistics, and the trained networks + loss history.
#' @export
#'
#' @examples
#' \donttest{
#' dat <- load_real_input_data(example = TRUE)
#' gan <- train_gan_on_real_data(dat$gan_training_data,
#'                               feature_correlations = dat$feature_correlations,
#'                               feature_texture = dat$feature_texture,
#'                               epochs = 50)
#' head(sample_texture(gan, 5))
#' }
train_gan_on_real_data <- function(real_data, feature_correlations = NULL,
                                   feature_texture = NULL,
                                   epochs = 300, batch_size = 32,
                                   lr = 2e-4, seed = NULL, verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  if (!check_torch_setup())
    stop("torch is required for the generative texture model. ",
         "Install with: install.packages('torch')")

  X <- as.data.frame(real_data)
  X <- X[stats::complete.cases(X), , drop = FALSE]
  num <- vapply(X, is.numeric, logical(1))
  X   <- X[, num, drop = FALSE]
  if (!nrow(X) || !ncol(X)) stop("real_data has no complete numeric rows/columns to train on.")

  columns      <- names(X)
  binary_cols  <- .detect_binary_cols(X)
  onehot_groups <- .detect_onehot_groups(columns)

  ctr     <- vapply(X, mean, numeric(1))
  scl     <- vapply(X, stats::sd, numeric(1)); scl[scl == 0 | is.na(scl)] <- 1
  Xn      <- as.matrix(sweep(sweep(X, 2, ctr, "-"), 2, scl, "/"))

  norm  <- list(center = ctr, scale = scl)
  stats <- list(mu = colMeans(as.matrix(X)), Sigma = stats::cov(as.matrix(X)),
                n_train = nrow(X))

  out <- .gan_train_loop(Xn, columns, norm, stats, epochs, batch_size, lr, seed, verbose,
                         binary_cols, onehot_groups, feature_correlations,
                         feature_texture)
  out
}


#' Torch GAN training loop (internal)
#' @keywords internal
.gan_train_loop <- function(Xn, columns, norm, stats, epochs, batch_size, lr,
                            seed, verbose, binary_cols, onehot_groups,
                            feature_correlations = NULL,
                            feature_texture = NULL) {
  if (!is.null(seed)) torch::torch_manual_seed(seed)
  noise_dim <- 100L
  d_in      <- ncol(Xn)
  n         <- nrow(Xn)
  bs        <- min(batch_size, n)

  G <- create_generator(d_in, noise_dim)
  D <- create_discriminator(d_in)
  g_opt <- torch::optim_adam(G$parameters, lr = lr, betas = c(0.5, 0.999))
  d_opt <- torch::optim_adam(D$parameters, lr = lr, betas = c(0.5, 0.999))
  bce   <- function(logit, target) torch::nnf_binary_cross_entropy_with_logits(logit, target)
  real_t <- torch::torch_tensor(Xn, dtype = torch::torch_float())

  d_hist <- numeric(epochs); g_hist <- numeric(epochs)
  for (ep in seq_len(epochs)) {
    idx   <- sample.int(n, bs, replace = n < bs)
    batch <- real_t[idx, , drop = FALSE]

    # Discriminator step (detach fake so the generator is not updated here).
    d_opt$zero_grad()
    z         <- torch::torch_randn(bs, noise_dim)
    fake      <- G(z)
    d_real    <- D(batch)
    d_fake    <- D(fake$detach())
    d_loss    <- bce(d_real, torch::torch_ones_like(d_real) * 0.9) +
                 bce(d_fake, torch::torch_zeros_like(d_fake))
    d_loss$backward(); d_opt$step()

    # Generator step (fresh noise/graph; non-saturating loss).
    g_opt$zero_grad()
    z2        <- torch::torch_randn(bs, noise_dim)
    gen_logit <- D(G(z2))
    g_loss    <- bce(gen_logit, torch::torch_ones_like(gen_logit))
    g_loss$backward(); g_opt$step()

    d_hist[ep] <- as.numeric(d_loss$item()); g_hist[ep] <- as.numeric(g_loss$item())
    if (verbose && ep %% max(1L, epochs %/% 5L) == 0)
      message(sprintf("epoch %d/%d  D=%.3f  G=%.3f", ep, epochs, d_hist[ep], g_hist[ep]))
  }
  G$eval()

  structure(list(model_type = "gan", columns = columns, norm = norm,
                 binary_cols = binary_cols, onehot_groups = onehot_groups,
                 feature_correlations = feature_correlations,
                 feature_texture = feature_texture,
                 generator = G, discriminator = D, noise_dim = noise_dim,
                 mu = stats$mu, Sigma = stats$Sigma, n_train = stats$n_train,
                 loss = list(d = d_hist, g = g_hist)),
            class = "iconic_gan")
}


#' Draw synthetic base rows from a trained texture model
#'
#' Binary columns (detected at training time) are rounded to \eqn{{0, 1}} after
#' de-normalisation, and one-hot dummy groups are made mutually exclusive
#' (the column with the highest pre-rounding value wins per row).
#'
#' @param trained_gan An `iconic_gan` object from [train_gan_on_real_data()].
#' @param n           Number of rows to draw.
#'
#' @return A data frame of `n` rows with the trained columns, on the original
#'   (de-normalised) scale. Binary columns contain only 0/1 values.
#' @export
sample_texture <- function(trained_gan, n) {
  stopifnot(inherits(trained_gan, "iconic_gan"))
  cols <- trained_gan$columns

  z   <- torch::torch_randn(n, trained_gan$noise_dim)
  out <- torch::with_no_grad(as.matrix(trained_gan$generator(z)$detach()))
  # De-normalise from the training z-scores back to the original scale.
  out <- sweep(sweep(out, 2, trained_gan$norm$scale, "*"),
               2, trained_gan$norm$center, "+")
  out <- as.data.frame(out)
  names(out) <- cols

  # Round binary columns back to 0/1 and enforce one-hot mutual exclusivity.
  bin  <- trained_gan$binary_cols %||% character(0)
  grps <- trained_gan$onehot_groups %||% list()
  out  <- .enforce_discrete(out, bin, grps)
  out
}


# Null-coalescing helper (internal): returns x if not NULL, otherwise y.
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Print method for iconic_gan objects
#'
#' @param x An `iconic_gan` object.
#' @param ... Unused.
#' @return Invisibly returns `x` (the `iconic_gan` object); called for its side effect of printing a human-readable summary.
#' @export
print.iconic_gan <- function(x, ...) {
  cat("<iconic_gan>\n")
  cat("  engine     : torch GAN\n")
  cat("  variables  :", paste(x$columns, collapse = ", "), "\n")
  cat("  trained on :", x$n_train, "samples\n")
  bin <- x$binary_cols %||% character(0)
  if (length(bin))
    cat("  binary cols:", paste(bin, collapse = ", "), "\n")
  fc <- x$feature_correlations %||% list()
  fc_present <- names(fc)[vapply(fc, function(m) !is.null(m), logical(1))]
  if (length(fc_present))
    cat("  feat. cors :", paste(fc_present, collapse = ", "), "\n")
  if (!is.null(x$feature_texture))
    cat("  feat. tex  : copula (", x$feature_texture$n_features, " features)\n", sep = "")
  if (!is.null(x$loss))
    cat(sprintf("  final loss : D=%.3f  G=%.3f\n",
                utils::tail(x$loss$d, 1), utils::tail(x$loss$g, 1)))
  invisible(x)
}
