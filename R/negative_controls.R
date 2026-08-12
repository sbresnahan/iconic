# ============================================================
# negative_controls.R: helpers for building the negative-control
# panel W from high-dimensional omics data.
#
# beta_to_m()           – methylation beta -> M-value transform
# residualize_matrix()  – regress covariates out of a feature matrix
# build_w_pcs()         – top PCs of a (residualized) matrix as W
# apply_fusion_weights()– FUSION/TWAS eQTL weights -> predicted
#                         expression (a germline-determined NC source)
#
# Distilled from the case-study pipelines (build_W_methylation.R,
# build_W_blood_predicted_CS2.R). All functions are pure R; irlba
# (Suggests) is used for fast truncated PCA when available.
# ============================================================


#' Convert methylation beta values to M-values
#'
#' Applies the logit transform \eqn{\log_2(\beta / (1 - \beta))} after
#' clipping betas away from 0 and 1. M-values are approximately
#' homoscedastic and are the recommended scale for linear modeling and
#' PCA of methylation data.
#'
#' @param beta Numeric matrix of beta values (probes x samples), in
#' \[0, 1\].
#' @param clip Numeric: betas are clipped to \[`clip`, `1 - clip`\]
#' before the transform. Default 1e-4.
#' @param drop_nonfinite Logical: drop probes (rows) with any
#' missing/non-finite values after the transform. Default `TRUE`.
#'
#' @return A numeric matrix of M-values, probes x samples.
#' @export
#'
#' @examples
#' b <- matrix(runif(200, 0.01, 0.99), nrow = 20,
#'             dimnames = list(paste0("cg", 1:20), paste0("S", 1:10)))
#' m <- beta_to_m(b)
#' range(m)
beta_to_m <- function(beta, clip = 1e-4, drop_nonfinite = TRUE) {
  beta <- as.matrix(beta)
  if (any(beta < 0 | beta > 1, na.rm = TRUE))
    stop("beta values must lie in [0, 1].")
  m <- log2(pmax(pmin(beta, 1 - clip), clip) /
              (1 - pmax(pmin(beta, 1 - clip), clip)))
  dimnames(m) <- dimnames(beta)
  if (isTRUE(drop_nonfinite)) {
    good <- rowSums(is.finite(m)) == ncol(m)
    if (any(!good))
      message("Dropped ", sum(!good), " probes with missing/non-finite ",
              "M-values.")
    m <- m[good, , drop = FALSE]
  }
  m
}


#' Residualize a feature matrix on covariates
#'
#' Regresses each feature (row) of `x` on the covariate design and
#' returns the residuals, using a precomputed projection for memory-
#' efficient chunked processing of large (e.g. 450k-probe) matrices.
#' This is the workhorse for removing cell composition, batch, and
#' technical factors before extracting negative-control PCs.
#'
#' @param x Numeric matrix, features x samples.
#' @param covariates data.frame or matrix with one row per sample
#' (i.e. `nrow(covariates) == ncol(x)`). An intercept is added
#' internally. Factors should already be expanded to dummies (see
#' [stats::model.matrix()]).
#' @param chunk_size Integer: number of features processed per block.
#' Default 5000.
#'
#' @return A numeric matrix of residuals, features x samples, centered
#' per feature.
#' @export
#'
#' @examples
#' x <- matrix(rnorm(100 * 40), nrow = 100,
#'             dimnames = list(paste0("f", 1:100), paste0("S", 1:40)))
#' batch <- factor(rep(c("A", "B"), each = 20))
#' cv <- model.matrix(~ batch)[, -1, drop = FALSE]
#' xr <- residualize_matrix(x, cv)
#' cor(as.numeric(xr[1, ]), as.numeric(cv))
residualize_matrix <- function(x, covariates, chunk_size = 5000) {
  x <- as.matrix(x)
  C <- as.matrix(as.data.frame(covariates))
  if (nrow(C) != ncol(x))
    stop("covariates must have one row per sample (ncol(x)).")
  if (any(!is.finite(C)))
    stop("covariates contain missing/non-finite values.")

  Xc <- scale(C, center = TRUE, scale = FALSE)
  XtX <- crossprod(Xc)
  XtX_inv <- tryCatch(solve(XtX), error = function(e) {
    message("Covariate crossproduct is singular; using pseudoinverse.")
    MASS::ginv(XtX)
  })
  proj <- XtX_inv %*% t(Xc)   # p x n

  n_feat <- nrow(x)
  out <- matrix(0, nrow = n_feat, ncol = ncol(x),
                dimnames = dimnames(x))
  for (start in seq(1, n_feat, by = chunk_size)) {
    end <- min(start + chunk_size - 1, n_feat)
    blk <- x[start:end, , drop = FALSE]
    blk <- blk - rowMeans(blk)
    B <- proj %*% t(blk)                 # p x k
    out[start:end, ] <- blk - tcrossprod(t(B), Xc)
  }
  out
}


#' Build a negative-control panel from principal components
#'
#' Computes the top `n_pcs` principal components of a features x
#' samples matrix (typically residualized with
#' [residualize_matrix()]) and returns the sample scores as a W panel
#' (PCs x samples) in the orientation expected by [iconic_data()].
#' Uses irlba for fast truncated PCA when available, falling back to
#' [stats::prcomp()].
#'
#' @param x Numeric matrix, features x samples (e.g. residualized
#' M-values or predicted expression).
#' @param n_pcs Integer: number of PCs to return. Default 20.
#' @param scale_features Logical: scale features to unit variance
#' before PCA. Default `FALSE` (appropriate for M-values); use `TRUE`
#' for predicted-expression panels where feature scales differ.
#' @param prefix Character: row-name prefix for the PCs. Default
#' `"PC"`.
#'
#' @return A list with:
#' \describe{
#'   \item{W}{numeric matrix, PCs x samples, suitable as the `W`
#'   argument of [iconic_data()].}
#'   \item{variance_explained}{numeric vector: percent of total
#'   variance explained by each returned PC.}
#' }
#' @export
#'
#' @examples
#' x <- matrix(rnorm(500 * 60), nrow = 500,
#'             dimnames = list(paste0("f", 1:500), paste0("S", 1:60)))
#' w <- build_w_pcs(x, n_pcs = 5)
#' dim(w$W)
#' w$variance_explained
build_w_pcs <- function(x, n_pcs = 20, scale_features = FALSE,
                        prefix = "PC") {
  x <- as.matrix(x)
  n <- ncol(x)
  if (n_pcs >= n)
    stop("n_pcs must be smaller than the number of samples.")
  xt <- t(x)   # samples x features
  if (isTRUE(scale_features)) xt <- scale(xt)

  if (requireNamespace("irlba", quietly = TRUE)) {
    pca <- irlba::irlba(xt, nu = n_pcs, nv = n_pcs, center = TRUE)
    pc_scores <- sweep(pca$u, 2, pca$d, `*`)
    pc_var <- pca$d^2 / (nrow(xt) - 1)
    total_var <- sum(apply(xt, 2, stats::var, na.rm = TRUE))
  } else {
    pca <- prcomp(xt, center = TRUE, scale. = FALSE, rank. = n_pcs)
    pc_scores <- pca$x
    pc_var <- pca$sdev^2
    total_var <- sum(apply(xt, 2, stats::var, na.rm = TRUE))
  }
  ve <- pc_var / total_var * 100

  W <- t(pc_scores)
  rownames(W) <- paste0(prefix, seq_len(n_pcs))
  colnames(W) <- colnames(x)
  list(W = W, variance_explained = ve)
}


#' Apply FUSION/TWAS eQTL weights to a dosage matrix
#'
#' Predicts genetically regulated expression by applying published
#' FUSION per-gene weight files (`.wgt.RDat`) to a genotype dosage
#' matrix. For each gene, the model with the best cross-validation
#' performance (`cv.performance["rsq", ]`) is used. The resulting
#' predicted-expression matrix is germline-determined and — after PCA
#' via [build_w_pcs()] — serves as a negative-control panel that
#' captures confounding structure without direct exposure mediation.
#'
#' @param genotypes Numeric dosage matrix, variants x samples, with
#' rownames matching the SNP IDs used in the weight files (usually
#' rsIDs).
#' @param pos A data.frame with (at least) columns `ID` (gene ID) and
#' `WGT` (weight-file name), as read from a FUSION `.pos` file; or a
#' path to such a file.
#' @param weights_dir Directory containing the `.wgt.RDat` files named
#' in `pos$WGT`. Ignored when `weights` is supplied.
#' @param weights Optional named list of pre-loaded weight objects,
#' each a list with elements `wgt.matrix`, `snps`, and
#' `cv.performance` (as stored in FUSION `.wgt.RDat` files). Names
#' must match `pos$ID`. Supplying this avoids repeated file I/O and
#' makes the function usable without the original files.
#' @param min_snps Integer: minimum number of matched SNPs required to
#' predict a gene. Default 2.
#' @param n_cores Integer: parallel workers via [parallel::mclapply()].
#' Default 1.
#'
#' @return A list with:
#' \describe{
#'   \item{predicted}{numeric matrix, genes x samples: predicted
#'   expression.}
#'   \item{info}{data.frame per predicted gene: `gene`, `n_snps`,
#'   `best_method`, `rsq` (cross-validation R2 of the selected
#'   model).}
#' }
#' @export
#'
#' @examples
#' # Simulate two FUSION-style weight objects in memory
#' dos <- matrix(rbinom(10 * 50, 2, 0.3), nrow = 10,
#'               dimnames = list(paste0("rs", 1:10), paste0("S", 1:50)))
#' mk_wgt <- function(snps, w) {
#'   list(wgt.matrix = matrix(w, ncol = 1, dimnames = list(NULL, "enet")),
#'        snps = data.frame(V2 = snps),
#'        cv.performance = matrix(0.2, nrow = 1,
#'                                dimnames = list("rsq", "enet")))
#' }
#' wlist <- list(GeneA = mk_wgt(paste0("rs", 1:5), rep(0.1, 5)),
#'               GeneB = mk_wgt(paste0("rs", 6:10), rep(-0.2, 5)))
#' pos <- data.frame(ID = c("GeneA", "GeneB"), WGT = c("a.wgt.RDat", "b.wgt.RDat"))
#' out <- apply_fusion_weights(dos, pos = pos, weights = wlist)
#' dim(out$predicted)
apply_fusion_weights <- function(genotypes, pos, weights_dir = NULL,
                                 weights = NULL, min_snps = 2,
                                 n_cores = 1) {
  G <- as.matrix(genotypes)
  if (is.null(rownames(G)))
    stop("genotypes must have SNP-ID rownames (usually rsIDs).")
  if (is.null(colnames(G)))
    colnames(G) <- paste0("sample", seq_len(ncol(G)))

  if (is.character(pos)) pos <- utils::read.delim(pos, stringsAsFactors = FALSE)
  pos <- as.data.frame(pos)
  if (!all(c("ID", "WGT") %in% colnames(pos)))
    stop("pos needs columns `ID` (gene) and `WGT` (weight file).")
  if (is.null(weights) && is.null(weights_dir))
    stop("Supply either `weights` (pre-loaded list) or `weights_dir`.")

  snp_lookup <- stats::setNames(seq_len(nrow(G)), rownames(G))

  predict_one <- function(gene_id, wgt_file) {
    wgt <- NULL
    if (!is.null(weights)) {
      wgt <- weights[[gene_id]]
    } else {
      path <- file.path(weights_dir, wgt_file)
      if (!file.exists(path)) return(NULL)
      wgt <- tryCatch({
        e <- new.env()
        load(path, envir = e)
        list(wgt.matrix = e$wgt.matrix, snps = e$snps,
             cv.performance = e$cv.performance)
      }, error = function(err) NULL)
    }
    if (is.null(wgt) || is.null(wgt$snps) || is.null(wgt$wgt.matrix))
      return(NULL)

    snp_ids <- wgt$snps$V2
    idx <- unname(snp_lookup[snp_ids])
    valid <- !is.na(idx)
    if (sum(valid) < min_snps) return(NULL)

    wgt_mat <- wgt$wgt.matrix[valid, , drop = FALSE]
    dosage_sub <- G[idx[valid], , drop = FALSE]

    available <- intersect(colnames(wgt_mat), colnames(wgt$cv.performance))
    if (length(available) == 0) return(NULL)
    rsq <- wgt$cv.performance["rsq", available]
    if (all(is.na(rsq))) return(NULL)
    best <- available[which.max(rsq)]

    pred <- as.numeric(crossprod(wgt_mat[, best], dosage_sub))
    list(gene = gene_id, pred = pred, n_snps = sum(valid),
         best_method = best, rsq = unname(rsq[best]))
  }

  worker <- function(i) predict_one(pos$ID[i], pos$WGT[i])
  res <- if (n_cores > 1) {
    parallel::mclapply(seq_len(nrow(pos)), worker, mc.cores = n_cores)
  } else {
    lapply(seq_len(nrow(pos)), worker)
  }
  res <- Filter(Negate(is.null), res)
  if (length(res) == 0)
    stop("No genes could be predicted (check SNP-ID matching and ",
         "weight files).")

  predicted <- do.call(rbind, lapply(res, `[[`, "pred"))
  rownames(predicted) <- vapply(res, `[[`, character(1), "gene")
  colnames(predicted) <- colnames(G)
  ## Drop genes with all-missing predictions
  predicted <- predicted[!apply(is.na(predicted), 1, all), , drop = FALSE]

  info <- data.frame(
    gene = vapply(res, `[[`, character(1), "gene"),
    n_snps = vapply(res, `[[`, integer(1), "n_snps"),
    best_method = vapply(res, `[[`, character(1), "best_method"),
    rsq = vapply(res, `[[`, numeric(1), "rsq"),
    stringsAsFactors = FALSE)
  info <- info[info$gene %in% rownames(predicted), , drop = FALSE]

  list(predicted = predicted, info = info)
}
