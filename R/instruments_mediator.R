# ============================================================
# instruments_mediator.R: helpers for building the mediator
# instrument G2 (cis-eQTL-based genetically predicted expression).
#
# call_cis_eqtls()             – per-gene cis-eQTL scan (single-SNP)
# build_mediator_instruments() – per-gene elastic-net composite
#                                instruments with out-of-fold QC
#
# Distilled from the case-study pipelines (call_GUSTO_eQTLs_gene.R,
# call_TCGA_eqtls_gene_CS2.R). Both functions residualize expression
# AND dosages on covariates (Frisch-Waugh-Lovell) so the fitted
# genetic effects are covariate-adjusted and the resulting
# instrument is orthogonal to the covariates ICONIC conditions on.
# glmnet is required only for build_mediator_instruments().
# ============================================================


#' Scan for cis-eQTLs of each gene
#'
#' For each gene, tests every cis SNP (within `cis_dist` of the gene
#' start, optionally restricted to a supplied candidate set) for
#' association with expression, after residualizing both expression
#' and dosages on `covariates` (Frisch-Waugh-Lovell). A lightweight,
#' pure-R alternative to MatrixEQTL for preparing mediator-instrument
#' candidate sets.
#'
#' @param expression Numeric matrix, genes x samples (e.g.
#' inverse-normal-transformed TPM).
#' @param genotypes Numeric dosage matrix, variants x samples (0-2).
#' @param gene_pos data.frame with columns `gene`, `chr`, and `tss`
#' (or `start`; the gene start is used as the cis-window anchor).
#' @param snp_pos data.frame with columns `snp`, `chr`, `pos`.
#' @param covariates Optional data.frame or matrix (n rows) of
#' covariates to partial out (e.g. ancestry PCs, technical factors).
#' @param cis_snps Optional named list: gene name -> character vector
#' of candidate SNP IDs (e.g. significant GTEx pairs). When supplied,
#' only these SNPs are tested for each gene (still intersected with
#' the cis window).
#' @param cis_dist Numeric: cis window half-width around the gene
#' start, in bp. Default 1e6.
#' @param min_maf Numeric: minimum minor allele frequency (computed
#' from mean dosage). Default 0.01.
#' @param fdr Numeric: BH FDR level for the per-gene `pass` flag in
#' the `best` table. Default 0.05.
#' @param n_cores Integer: parallel workers via [parallel::mclapply()].
#' Default 1.
#'
#' @return A list with:
#' \describe{
#'   \item{all}{data.frame of all tested gene-SNP pairs: `gene`,
#'   `snp`, `beta`, `se`, `t`, `p`.}
#'   \item{best}{data.frame with the top SNP per gene (`gene`, `snp`,
#'   `beta`, `p`, `q` (BH within gene), `pass`).}
#' }
#' @export
#'
#' @examples
#' n <- 80
#' dos <- matrix(rbinom(200 * n, 2, 0.3), nrow = 200,
#'               dimnames = list(paste0("1:", 1:200, ":A:G"), NULL))
#' expr <- matrix(rnorm(20 * n), nrow = 20,
#'                dimnames = list(paste0("Gene", 1:20), NULL))
#' expr[1, ] <- expr[1, ] + 0.5 * dos[5, ]   # one true cis-eQTL
#' gp <- data.frame(gene = paste0("Gene", 1:20), chr = "1",
#'                  tss = seq(1, 191, by = 10))
#' sp <- data.frame(snp = rownames(dos), chr = "1", pos = 1:200)
#' hits <- call_cis_eqtls(expr, dos, gp, sp)
#' head(hits$best)
call_cis_eqtls <- function(expression, genotypes, gene_pos, snp_pos,
                           covariates = NULL, cis_snps = NULL,
                           cis_dist = 1e6, min_maf = 0.01, fdr = 0.05,
                           n_cores = 1) {
  expression <- as.matrix(expression)
  genotypes  <- as.matrix(genotypes)
  if (ncol(expression) != ncol(genotypes))
    stop("expression and genotypes must have the same number of samples.")
  if (!all(c("gene", "chr") %in% colnames(gene_pos)))
    stop("gene_pos needs columns `gene` and `chr` (and `tss` or `start`).")
  if (!all(c("snp", "chr", "pos") %in% colnames(snp_pos)))
    stop("snp_pos needs columns `snp`, `chr`, `pos`.")
  if (!"tss" %in% colnames(gene_pos)) {
    if (!"start" %in% colnames(gene_pos))
      stop("gene_pos needs a `tss` or `start` column.")
    gene_pos$tss <- gene_pos$start
  }

  n <- ncol(expression)
  genes <- intersect(rownames(expression), gene_pos$gene)
  if (length(genes) == 0)
    stop("No overlap between expression rownames and gene_pos$gene.")
  snp_pos <- snp_pos[snp_pos$snp %in% rownames(genotypes), , drop = FALSE]

  ## FWL residualization on covariates
  qr_cvrt <- NULL
  if (!is.null(covariates)) {
    C <- as.matrix(as.data.frame(covariates))
    if (nrow(C) != n) stop("covariates must have n rows.")
    qr_cvrt <- qr(cbind(`(Intercept)` = 1, C))
    expr_r <- t(qr.resid(qr_cvrt, t(expression)))   # genes x samples
    geno_r <- qr.resid(qr_cvrt, t(genotypes))        # samples x snps
    geno_r <- t(geno_r)
  } else {
    expr_r <- expression - rowMeans(expression)
    geno_r <- genotypes - rowMeans(genotypes)
  }
  colnames(expr_r) <- colnames(expression)
  colnames(geno_r) <- colnames(genotypes)

  gp <- gene_pos[match(genes, gene_pos$gene), , drop = FALSE]

  scan_gene <- function(g) {
    gpos <- gp[gp$gene == g, , drop = FALSE]
    if (nrow(gpos) != 1) return(NULL)
    cis <- which(snp_pos$chr == gpos$chr &
                   abs(snp_pos$pos - gpos$tss) <= cis_dist)
    if (!is.null(cis_snps)) {
      cand <- cis_snps[[g]]
      if (is.null(cand)) return(NULL)
      cis <- cis[snp_pos$snp[cis] %in% cand]
    }
    if (length(cis) < 1) return(NULL)
    snps <- snp_pos$snp[cis]
    X <- geno_r[snps, , drop = FALSE]
    ## MAF must come from the RAW genotypes: geno_r is residualized
    ## (mean ~0), so allele frequencies computed on it are meaningless.
    Xraw <- genotypes[snps, , drop = FALSE]
    maf <- pmin(rowMeans(Xraw, na.rm = TRUE) / 2,
                1 - rowMeans(Xraw, na.rm = TRUE) / 2)
    keep <- is.finite(maf) & maf >= min_maf &
      apply(X, 1, stats::var, na.rm = TRUE) > 0
    if (sum(keep) < 1) return(NULL)
    X <- X[keep, , drop = FALSE]
    snps <- snps[keep]
    y <- expr_r[g, ]

    ## Correlation-based test on residualized data (equivalent to
    ## marginal OLS t-test after FWL).
    sx <- sqrt(rowSums((X - rowMeans(X))^2))
    sy <- sqrt(sum((y - mean(y))^2))
    if (sy == 0) return(NULL)
    r <- rowSums(sweep(X - rowMeans(X), 2, y - mean(y), `*`)) / (sx * sy)
    dfree <- n - 2 - if (!is.null(qr_cvrt)) qr_cvrt$rank - 1L else 0L
    tstat <- r * sqrt(dfree / pmax(1e-12, 1 - r^2))
    p <- 2 * stats::pt(-abs(tstat), df = dfree)
    beta <- r * (sy / sx)
    se <- abs(beta / tstat)
    data.frame(gene = g, snp = snps, beta = beta, se = se,
               t = tstat, p = p, stringsAsFactors = FALSE)
  }

  res_list <- if (n_cores > 1) {
    parallel::mclapply(genes, scan_gene, mc.cores = n_cores)
  } else {
    lapply(genes, scan_gene)
  }
  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) == 0)
    stop("No testable gene-SNP pairs (check cis window, MAF, and IDs).")
  all_hits <- do.call(rbind, res_list)
  rownames(all_hits) <- NULL

  ## Best SNP per gene with within-gene BH q-value
  best <- do.call(rbind, lapply(split(all_hits, all_hits$gene), function(d) {
    d <- d[order(d$p), , drop = FALSE]
    d$q <- stats::p.adjust(d$p, method = "BH")
    d[1, , drop = FALSE]
  }))
  best$pass <- best$q <= fdr
  rownames(best) <- NULL

  list(all = all_hits, best = best)
}


#' Build mediator instruments via per-gene elastic net
#'
#' Trains one cross-validated elastic-net model per gene predicting
#' expression from cis SNPs, after residualizing both expression and
#' dosages on `covariates` (Frisch-Waugh-Lovell). Genes are kept when
#' the model is non-degenerate (at least one non-zero SNP weight), the
#' signed out-of-fold R2 exceeds `cv_r2_min`, and the one-sided
#' out-of-fold correlation p-value is below `cv_p_max`. The in-sample
#' predicted expression of passing genes is the composite mediator
#' instrument `Gm` consumed by [iconic_data()].
#'
#' Requires the glmnet package (listed under `Suggests`).
#'
#' @param expression Numeric matrix, genes x samples (e.g.
#' inverse-normal-transformed TPM).
#' @param genotypes Numeric dosage matrix, variants x samples (0-2).
#' @param gene_pos data.frame with columns `gene`, `chr`, and `tss`
#' (or `start`).
#' @param snp_pos data.frame with columns `snp`, `chr`, `pos`.
#' @param covariates Optional data.frame or matrix (n rows) of
#' covariates to partial out (ancestry PCs, technical factors, ...).
#' @param cis_snps Optional named list: gene name -> character vector
#' of candidate SNP IDs (e.g. significant GTEx cis-eQTL pairs, as from
#' [call_cis_eqtls()] or an external resource). When supplied, the
#' elastic net for each gene sees only these SNPs (intersected with
#' the cis window).
#' @param alpha Numeric: elastic-net mixing (0 = ridge, 1 = lasso).
#' Default 0.3.
#' @param nfolds Integer: cross-validation folds. Default 5.
#' @param cis_dist Numeric: cis window half-width around the gene
#' start, in bp. Default 1e6.
#' @param cv_r2_min Numeric: minimum signed out-of-fold R2. Default
#' 0.005.
#' @param cv_p_max Numeric: maximum one-sided out-of-fold correlation
#' p-value (nominal screen only; the load-bearing gates are
#' `cv_r2_min` and the non-degeneracy requirement). Default 0.10.
#' @param n_cores Integer: parallel workers via [parallel::mclapply()].
#' Default 1.
#' @param seed Optional integer for reproducible CV fold assignments.
#' @param verbose Logical: print progress. Default `TRUE`.
#'
#' @return A list with:
#' \describe{
#'   \item{Gm}{numeric matrix, genes x samples: in-sample predicted
#'   expression for genes passing the quality filter. Pass to
#'   [iconic_data()] as `Gm` (it is transposed internally).}
#'   \item{qc}{data.frame of per-gene quality metrics: `gene`,
#'   `n_cis_snps`, `n_nonzero`, `cv_r2`, `cor_oof`, `cv_p`,
#'   `lambda_min`, `pass`.}
#'   \item{weights}{named list (all trained genes) of data.frames with
#'   `snp` and `weight` (non-zero coefficients at lambda.min).}
#' }
#' @export
#'
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(42)
#'   n <- 150; p <- 80
#'   dos <- matrix(rbinom(p * n, 2, 0.3), nrow = p,
#'                 dimnames = list(paste0("1:", 1:p, ":A:G"),
#'                                 paste0("S", 1:n)))
#'   expr <- matrix(rnorm(4 * n), nrow = 4,
#'                  dimnames = list(paste0("Gene", 1:4), paste0("S", 1:n)))
#'   expr[1, ] <- expr[1, ] + 0.9 * scale(dos[5, ]) + 0.9 * scale(dos[15, ])
#'   gp <- data.frame(gene = paste0("Gene", 1:4), chr = "1",
#'                    tss = c(5, 25, 45, 65))
#'   sp <- data.frame(snp = rownames(dos), chr = "1", pos = 1:p)
#'   fit <- build_mediator_instruments(expr, dos, gp, sp, seed = 1)
#'   fit$qc
#' }
build_mediator_instruments <- function(expression, genotypes, gene_pos,
                                       snp_pos, covariates = NULL,
                                       cis_snps = NULL, alpha = 0.3,
                                       nfolds = 5, cis_dist = 1e6,
                                       cv_r2_min = 0.005, cv_p_max = 0.10,
                                       n_cores = 1, seed = NULL,
                                       verbose = TRUE) {
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("build_mediator_instruments() requires the glmnet package. ",
         "Install with: install.packages('glmnet')")
  if (!is.null(seed)) withr::local_seed(seed)

  expression <- as.matrix(expression)
  genotypes  <- as.matrix(genotypes)
  if (ncol(expression) != ncol(genotypes))
    stop("expression and genotypes must have the same number of samples.")
  if (!all(c("gene", "chr") %in% colnames(gene_pos)))
    stop("gene_pos needs columns `gene` and `chr` (and `tss` or `start`).")
  if (!all(c("snp", "chr", "pos") %in% colnames(snp_pos)))
    stop("snp_pos needs columns `snp`, `chr`, `pos`.")
  if (!"tss" %in% colnames(gene_pos)) {
    if (!"start" %in% colnames(gene_pos))
      stop("gene_pos needs a `tss` or `start` column.")
    gene_pos$tss <- gene_pos$start
  }

  n <- ncol(expression)
  genes <- intersect(rownames(expression), gene_pos$gene)
  if (!is.null(cis_snps))
    genes <- intersect(genes, names(cis_snps))
  if (length(genes) == 0)
    stop("No genes to test (check expression rownames, gene_pos, and ",
         "cis_snps names).")
  snp_pos <- snp_pos[snp_pos$snp %in% rownames(genotypes), , drop = FALSE]

  ## FWL residualization: expression AND dosages on covariates
  qr_cvrt <- NULL
  if (!is.null(covariates)) {
    C <- as.matrix(as.data.frame(covariates))
    if (nrow(C) != n) stop("covariates must have n rows.")
    qr_cvrt <- qr(cbind(`(Intercept)` = 1, C))
    expr_resid <- qr.resid(qr_cvrt, t(expression))   # samples x genes
    colnames(expr_resid) <- rownames(expression)
  } else {
    expr_resid <- t(expression - rowMeans(expression))
  }

  train_gene <- function(g) {
    gpos <- gene_pos[gene_pos$gene == g, , drop = FALSE]
    if (nrow(gpos) != 1) return(NULL)
    cis <- which(snp_pos$chr == gpos$chr &
                   abs(snp_pos$pos - gpos$tss) <= cis_dist)
    if (!is.null(cis_snps)) {
      cand <- cis_snps[[g]]
      if (is.null(cand)) return(NULL)
      cis <- cis[snp_pos$snp[cis] %in% cand]
    }
    if (length(cis) < 2) return(NULL)
    snps <- snp_pos$snp[cis]
    X <- t(genotypes[snps, , drop = FALSE])          # samples x snps
    y_resid <- expr_resid[, g]

    if (!is.null(qr_cvrt)) X <- qr.resid(qr_cvrt, X)
    snp_var <- apply(X, 2, stats::var, na.rm = TRUE)
    keep <- is.finite(snp_var) & snp_var > 0
    if (sum(keep) < 2) return(NULL)
    X <- X[, keep, drop = FALSE]

    cv_fit <- glmnet::cv.glmnet(X, y_resid, alpha = alpha, nfolds = nfolds,
                                type.measure = "mse", standardize = TRUE,
                                keep = TRUE)
    lam_idx  <- which.min(abs(cv_fit$lambda - cv_fit$lambda.min))
    pred_oof <- as.numeric(cv_fit$fit.preval[, lam_idx])

    coefs    <- glmnet::coef.glmnet(cv_fit, s = "lambda.min")
    nz_names <- setdiff(rownames(coefs)[as.numeric(coefs) != 0], "(Intercept)")
    n_nonzero <- length(nz_names)

    ## Signed out-of-fold R2 (negative for null/anti-predictive models)
    tss     <- sum((y_resid - mean(y_resid))^2)
    mse_oof <- sum((y_resid - pred_oof)^2)
    cv_r2   <- if (tss > 0) 1 - mse_oof / tss else NA_real_
    cor_oof <- if (stats::var(pred_oof) > 0)
      stats::cor(pred_oof, y_resid, use = "complete.obs") else NA_real_
    cv_p <- if (is.finite(cor_oof))
      stats::cor.test(pred_oof, y_resid, alternative = "greater")$p.value
    else NA_real_

    pred_insample <- as.numeric(stats::predict(cv_fit, newx = X,
                                               s = "lambda.min",
                                               type = "response"))
    names(pred_insample) <- colnames(expression)

    list(gene = g, pred_expr = pred_insample, cv_r2 = cv_r2,
         cor_oof = cor_oof, cv_p = cv_p, lambda_min = cv_fit$lambda.min,
         n_cis_snps = ncol(X), n_nonzero = n_nonzero,
         weights = data.frame(snp = nz_names,
                              weight = as.numeric(coefs[nz_names, 1]),
                              stringsAsFactors = FALSE))
  }

  if (isTRUE(verbose))
    message("Training elastic-net instruments for ", length(genes), " genes...")
  res_list <- if (n_cores > 1) {
    parallel::mclapply(genes, function(g)
      tryCatch(train_gene(g), error = function(e) NULL), mc.cores = n_cores)
  } else {
    lapply(genes, function(g)
      tryCatch(train_gene(g), error = function(e) NULL))
  }
  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) == 0)
    stop("No genes could be trained (need >= 2 non-constant cis SNPs).")
  names(res_list) <- vapply(res_list, function(r) r$gene, character(1))

  qc <- do.call(rbind, lapply(res_list, function(r) data.frame(
    gene = r$gene, n_cis_snps = r$n_cis_snps, n_nonzero = r$n_nonzero,
    cv_r2 = r$cv_r2, cor_oof = r$cor_oof, cv_p = r$cv_p,
    lambda_min = r$lambda_min, stringsAsFactors = FALSE)))
  rownames(qc) <- NULL
  qc$pass <- qc$n_nonzero >= 1L &
    is.finite(qc$cv_r2) & qc$cv_r2 > cv_r2_min &
    is.finite(qc$cv_p)  & qc$cv_p  < cv_p_max

  if (isTRUE(verbose)) {
    message("Quality filter (OOF R2 > ", cv_r2_min, ", one-sided p < ",
            cv_p_max, ", >= 1 non-zero SNP):")
    message("  Intercept-only fits (excluded): ", sum(qc$n_nonzero == 0L))
    message("  Passing: ", sum(qc$pass), " / ", nrow(qc), " genes")
  }

  pass_genes <- qc$gene[qc$pass]
  Gm <- NULL
  if (length(pass_genes) > 0) {
    Gm <- do.call(rbind, lapply(res_list[pass_genes], `[[`, "pred_expr"))
    rownames(Gm) <- pass_genes
    if (!is.null(colnames(expression)))
      Gm <- Gm[, colnames(expression), drop = FALSE]
    row_sd <- apply(Gm, 1, stats::sd)
    if (any(row_sd < 1e-12))
      warning(sum(row_sd < 1e-12),
              " passing genes have ~zero-variance instruments; ",
              "inspect qc before use.")
  } else {
    warning("No genes passed the quality filter; Gm is NULL. ",
            "Inspect the qc table and consider relaxing cv_r2_min.")
  }

  list(Gm = Gm, qc = qc,
       weights = lapply(res_list, `[[`, "weights"))
}
