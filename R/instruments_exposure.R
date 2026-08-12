# ============================================================
# instruments_exposure.R: helpers for building the exposure
# instrument G1 (polygenic score) from GWAS summary statistics.
#
# qc_gwas_sumstats()          – clean + QC GWAS summary statistics
# build_prs_ldpred2()         – LDpred2-auto Bayesian PRS (bigsnpr)
# score_pgs_panel()           – score a published weight panel
#                               (e.g. PGS Catalog) on a dosage matrix
# check_instrument_strength() – first-stage partial F / R2 for G1
#
# These distill the case-study pipelines (build_G1_prs_ldpred2.R,
# build_G1_prs_CS2.R) into reusable functions. bigsnpr is required
# only for build_prs_ldpred2(); the rest are pure R.
# ============================================================


#' Standardize and QC GWAS summary statistics
#'
#' Renames common GWAS summary-statistic column conventions to the
#' bigsnpr convention (`chr`, `pos`, `a0`, `a1`, `beta`, `beta_se`,
#' `n_eff`, `freq`) and applies the standard LDpred2 QC filters:
#' missing/invalid values, extreme effect sizes, ambiguous-strand
#' (A/T and C/G) SNPs, and — when effect allele frequencies are
#' available — the standard-deviation ratio check from the LDpred2
#' tutorial.
#'
#' @param sumstats A data.frame of GWAS summary statistics. Recognized
#' column aliases: `chr`/`chromosome`/`CHR`; `pos`/`base_pair_location`/
#' `BP`; `a1`/`effect_allele`/`ALT`; `a0`/`other_allele`/`REF`;
#' `beta`/`BETA`/`effect`/`logOR`; `beta_se`/`standard_error`/`se`/`SE`;
#' `freq`/`effect_allele_frequency`/`EAF`/`FRQ`; `n_eff`/`n`/`N`;
#' `p`/`p_value`/`P`; `rsid`/`rs_id`/`SNP`/`ID`.
#' @param n_eff Optional scalar: effective sample size, used when the
#' input has no sample-size column.
#' @param drop_ambiguous Logical: drop A/T and C/G SNPs whose strand
#' cannot be resolved. Default `TRUE`.
#' @param beta_iqr_mult Numeric: drop SNPs with
#' \eqn{|\beta - \mathrm{median}(\beta)| >} `beta_iqr_mult`
#' \eqn{\times \mathrm{IQR}(\beta)}. Default 10. Set to `Inf` to skip.
#' @param sd_ratio_range Numeric length-2: acceptable range for the
#' ratio of the summary-statistic-implied genotype SD to the
#' frequency-implied SD, \eqn{\sqrt{2 f (1-f)}}. SNPs outside the range
#' are dropped. Default `c(0.5, 2)` (LDpred2 tutorial). Requires `freq`;
#' skipped with a message when frequencies are absent.
#' @param sd_y Optional scalar: phenotypic SD used in the SD-ratio check.
#' When `NULL` (default), estimated as the 1st percentile of
#' \eqn{\sqrt{0.5 (n_{eff} \cdot se^2 + \beta^2)}} across SNPs, the
#' continuous-trait approximation used in the LDpred2 tutorial.
#'
#' @return A list with:
#' \describe{
#'   \item{sumstats}{data.frame of QC-passing variants with standardized
#'   columns `chr`, `pos`, `a0`, `a1`, `beta`, `beta_se`, and (when
#'   available) `n_eff`, `freq`, `p`, `rsid`.}
#'   \item{qc}{data.frame with counts of input variants and variants
#'   dropped by each filter.}
#'   \item{sd_ratio}{numeric vector of SD ratios for retained variants,
#'   or `NULL` when the check was skipped.}
#' }
#' @export
#'
#' @examples
#' ss <- data.frame(
#'   chromosome = 1, base_pair_location = 1000 + 0:9 * 1000L,
#'   effect_allele = rep(c("A", "C"), 5), other_allele = rep(c("G", "T"), 5),
#'   beta = rnorm(10, 0, 0.05), standard_error = 0.01,
#'   effect_allele_frequency = runif(10, 0.1, 0.5)
#' )
#' out <- qc_gwas_sumstats(ss, n_eff = 50000)
#' str(out$sumstats)
#' out$qc
qc_gwas_sumstats <- function(sumstats, n_eff = NULL,
                             drop_ambiguous = TRUE,
                             beta_iqr_mult = 10,
                             sd_ratio_range = c(0.5, 2),
                             sd_y = NULL) {
  sumstats <- as.data.frame(sumstats)
  n_in <- nrow(sumstats)

  ## --- Column standardization -------------------------------------------
  aliases <- list(
    chr   = c("chr", "chromosome", "CHR", "chrom", "CHROM"),
    pos   = c("pos", "base_pair_location", "BP", "bp", "POS", "position"),
    a1    = c("a1", "effect_allele", "A1", "ALT", "alt", "ea"),
    a0    = c("a0", "other_allele", "A0", "REF", "ref", "nea", "a2"),
    beta  = c("beta", "BETA", "effect", "b", "logOR"),
    beta_se = c("beta_se", "standard_error", "se", "SE", "stderr"),
    freq  = c("freq", "effect_allele_frequency", "EAF", "eaf", "FRQ",
              "frq", "maf_a1"),
    n_eff = c("n_eff", "n", "N", "neff", "samplesize"),
    p     = c("p", "p_value", "P", "pval", "p_value_nominal"),
    rsid  = c("rsid", "rs_id", "SNP", "ID", "variant_id", "marker")
  )
  out_list <- list()
  for (col in names(aliases)) {
    hit <- intersect(aliases[[col]], colnames(sumstats))
    if (length(hit) > 0) {
      out_list[[col]] <- sumstats[[hit[1]]]
    }
  }
  out <- as.data.frame(out_list, stringsAsFactors = FALSE)
  required <- c("chr", "pos", "a0", "a1", "beta", "beta_se")
  missing_cols <- setdiff(required, colnames(out))
  if (length(missing_cols) > 0)
    stop("Could not find required columns: ", paste(missing_cols, collapse = ", "),
         ". Supply columns named like ", paste(vapply(missing_cols, function(m) aliases[[m]][1], character(1)), collapse = ", "), ".")

  if (is.null(out$n_eff)) {
    if (is.null(n_eff))
      stop("No sample-size column found; supply `n_eff` (scalar).")
    out$n_eff <- n_eff
  }

  out$chr <- as.character(out$chr)
  out$pos <- as.integer(out$pos)
  out$a0  <- toupper(as.character(out$a0))
  out$a1  <- toupper(as.character(out$a1))
  out$beta    <- as.numeric(out$beta)
  out$beta_se <- as.numeric(out$beta_se)
  out$n_eff   <- as.numeric(out$n_eff)
  if (!is.null(out$freq)) out$freq <- as.numeric(out$freq)

  ## --- Filter 1: missing / invalid --------------------------------------
  keep <- stats::complete.cases(out[, required]) &
    is.finite(out$beta) & is.finite(out$beta_se) & out$beta_se > 0
  n_missing <- sum(!keep)
  out <- out[keep, , drop = FALSE]

  ## --- Filter 2: extreme betas ------------------------------------------
  n_extreme <- 0L
  if (is.finite(beta_iqr_mult) && nrow(out) > 10) {
    beta_iqr <- stats::IQR(out$beta, na.rm = TRUE)
    beta_med <- stats::median(out$beta, na.rm = TRUE)
    if (beta_iqr > 0) {
      keep <- abs(out$beta - beta_med) <= beta_iqr_mult * beta_iqr
      n_extreme <- sum(!keep)
      out <- out[keep, , drop = FALSE]
    }
  }

  ## --- Filter 3: ambiguous strand (A/T, C/G) -----------------------------
  n_ambig <- 0L
  if (isTRUE(drop_ambiguous)) {
    ambig <- (out$a0 == "A" & out$a1 == "T") |
             (out$a0 == "T" & out$a1 == "A") |
             (out$a0 == "C" & out$a1 == "G") |
             (out$a0 == "G" & out$a1 == "C")
    n_ambig <- sum(ambig)
    out <- out[!ambig, , drop = FALSE]
  }

  ## --- Filter 4: SD-ratio check (needs freq) -----------------------------
  sd_ratio <- NULL
  n_sdratio <- 0L
  if (!is.null(out$freq) && !is.null(sd_ratio_range)) {
    if (is.null(sd_y))
      sd_y <- as.numeric(stats::quantile(
        sqrt(0.5 * (out$n_eff * out$beta_se^2 + out$beta^2)), 0.01,
        na.rm = TRUE))
    sd_ss   <- sd_y / sqrt(out$n_eff * out$beta_se^2 + out$beta^2)
    sd_geno <- sqrt(2 * out$freq * (1 - out$freq))
    sd_ratio <- sd_ss / sd_geno
    bad <- sd_ratio < sd_ratio_range[1] | sd_ratio > sd_ratio_range[2] |
      is.na(sd_ratio)
    n_sdratio <- sum(bad)
    if (n_sdratio > 0.1 * nrow(out))
      message("More than 10% of variants flagged by the SD-ratio check; ",
              "consider shrink_corr = 0.4 in build_prs_ldpred2().")
    out <- out[!bad, , drop = FALSE]
    sd_ratio <- sd_ratio[!bad]
  } else if (is.null(out$freq)) {
    message("No effect-allele-frequency column found; skipping the ",
            "SD-ratio QC step.")
  }

  qc <- data.frame(
    n_input = n_in,
    dropped_missing = n_missing,
    dropped_extreme_beta = n_extreme,
    dropped_ambiguous = n_ambig,
    dropped_sd_ratio = n_sdratio,
    n_output = nrow(out)
  )
  list(sumstats = out, qc = qc, sd_ratio = sd_ratio)
}


#' Build an exposure polygenic score with LDpred2-auto
#'
#' Runs the LDpred2-auto Bayesian shrinkage pipeline (Privé et al.) via
#' the bigsnpr package: match summary statistics to an LD reference,
#' estimate SNP heritability by LD score regression, run a grid of
#' auto chains, filter chains by convergence, and average the
#' per-chain effects. Optionally scores a target genotype set.
#'
#' This function requires the bigsnpr package (listed under `Suggests`)
#' and a PLINK/bed-format LD reference panel converted with
#' [bigsnpr::snp_readBed()].
#'
#' @param sumstats A data.frame of GWAS summary statistics, ideally
#' pre-processed by [qc_gwas_sumstats()]. Must contain (after
#' standardization) `chr`, `pos`, `a0`, `a1`, `beta`, `beta_se`,
#' `n_eff`; `freq` is used for the SD-ratio QC when present.
#' @param ld_ref Either a path to a bigsnpr `.rds` backing file (from
#' [bigsnpr::snp_readBed()]) or an attached `bigSNP` object. Should
#' match the ancestry of the GWAS.
#' @param genotypes Optional: path to a bigsnpr `.rds` or an attached
#' `bigSNP` object for the target cohort to be scored. When `NULL`
#' (default), only the LDpred2 effect sizes are returned.
#' @param hapmap3 Optional character vector of variant IDs (rsIDs or
#' `chr:pos:a0:a1` keys) to restrict the analysis to (e.g. HapMap3).
#' Restriction improves chain stability; `NULL` (default) uses all
#' matched variants.
#' @param max_variants Optional integer: if the matched variant count
#' exceeds this, systematically thin to approximately `max_variants`
#' variants (evenly spaced along the genome) before computing LD.
#' Default `NULL` (no thinning). Use e.g. `250000` to bound memory.
#' @param n_chains Integer: number of LDpred2-auto chains. Default 30.
#' @param p_range Numeric length-2: range of the initial per-chain
#' proportion-of-causal-variants grid, log-spaced. Default
#' `c(1e-4, 0.2)`.
#' @param burn_in Integer: burn-in iterations per chain. Default 500.
#' @param num_iter Integer: sampling iterations per chain (after
#' burn-in). Default 500.
#' @param shrink_corr Numeric: LD shrinkage parameter. Default 0.95;
#' reduce to ~0.4 when the SD-ratio QC flags >10% of variants.
#' @param allow_jump_sign Logical: allow sign jumps across chains.
#' Default `FALSE` (recommended for auto mode).
#' @param use_mle Logical: use MLE for the hyperparameter updates.
#' Default `FALSE` (appropriate for moderate GWAS sample sizes).
#' @param chain_filter_quantile Numeric: chains are kept when the range
#' of their posterior `corr_est` exceeds
#' `chain_filter_quantile`-quantile of all chain ranges times 0.95.
#' Default 0.95.
#' @param impute_mean Logical: mean-impute missing target genotypes
#' before scoring. Default `TRUE`.
#' @param ncores Integer: threads for bigsnpr. Default 1.
#' @param verbose Logical: print progress. Default `TRUE`.
#'
#' @return A list with:
#' \describe{
#'   \item{beta}{data.frame of LDpred2 effect sizes (`chr`, `pos`,
#'   `a0`, `a1`, `beta`) in LD-reference orientation.}
#'   \item{score}{numeric vector of per-sample polygenic scores for
#'   `genotypes` (scaled to unit variance), or `NULL` when `genotypes`
#'   was not supplied.}
#'   \item{h2}{LD-score-regression SNP heritability estimate used to
#'   initialize the chains.}
#'   \item{chains_kept}{integer indices of chains passing the
#'   convergence filter.}
#'   \item{qc}{list with variant counts at each stage.}
#' }
#' @export
#'
#' @examples
#' if (requireNamespace("bigsnpr", quietly = TRUE)) {
#'   # Tiny fake genotype panel standing in for an LD reference
#'   # (snp_fake genotypes are all-missing mocks; fill with random dosages)
#'   fake <- bigsnpr::snp_fake(100, 500)
#'   fake$genotypes[] <- rbinom(100 * 500, 2, 0.3)
#'   fake$map$chromosome <- 1L
#'   fake$map$physical.pos <- sort(sample(1:1e6, 500))
#'   ss <- data.frame(
#'     chr = fake$map$chromosome, pos = fake$map$physical.pos,
#'     a0 = fake$map$allele2, a1 = fake$map$allele1,
#'     beta = rnorm(500, 0, 0.05), beta_se = 0.02, n_eff = 20000
#'   )
#'   prs <- build_prs_ldpred2(ss, ld_ref = fake, n_chains = 3,
#'                            burn_in = 50, num_iter = 50, verbose = FALSE)
#'   str(prs$beta)
#' }
build_prs_ldpred2 <- function(sumstats, ld_ref, genotypes = NULL,
                              hapmap3 = NULL, max_variants = NULL,
                              n_chains = 30, p_range = c(1e-4, 0.2),
                              burn_in = 500, num_iter = 500,
                              shrink_corr = 0.95, allow_jump_sign = FALSE,
                              use_mle = FALSE, chain_filter_quantile = 0.95,
                              impute_mean = TRUE, ncores = 1,
                              verbose = TRUE) {
  if (!requireNamespace("bigsnpr", quietly = TRUE))
    stop("build_prs_ldpred2() requires the bigsnpr package. ",
         "Install with: install.packages('bigsnpr')")

  ## --- Attach LD reference ----------------------------------------------
  attach_ref <- function(x, label) {
    if (is.character(x)) {
      rds <- if (grepl("\\.rds$", x, ignore.case = TRUE)) x else paste0(x, ".rds")
      if (!file.exists(rds))
        stop(label, " backing file not found: ", rds,
             ". Convert PLINK files first with bigsnpr::snp_readBed().")
      return(bigsnpr::snp_attach(rds))
    }
    if (inherits(x, "bigSNP")) return(x)
    stop(label, " must be a path to a bigsnpr .rds file or an attached bigSNP object.")
  }
  obj.ref <- attach_ref(ld_ref, "ld_ref")

  ## --- Standardize sumstats ---------------------------------------------
  if (!all(c("chr", "pos", "a0", "a1", "beta", "beta_se") %in% colnames(sumstats)))
    sumstats <- qc_gwas_sumstats(sumstats)$sumstats
  if (is.null(sumstats$n_eff))
    stop("sumstats needs an `n_eff` column (effective sample size).")

  normalize_chr <- function(x) {
    x <- as.character(x)
    ifelse(grepl("^chr", x), x, paste0("chr", x))
  }
  sumstats$chr <- normalize_chr(sumstats$chr)

  ## --- Match to LD reference ---------------------------------------------
  map.ref <- stats::setNames(obj.ref$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
  map.ref$chr <- normalize_chr(map.ref$chr)
  if (isTRUE(verbose)) message("Matching ", nrow(sumstats), " sumstats to ",
                               nrow(map.ref), " LD-reference variants...")
  info <- bigsnpr::snp_match(sumstats, map.ref, return_flip_and_rev = TRUE)

  ## Optional HapMap3-style restriction
  if (!is.null(hapmap3)) {
    key <- paste0(info$chr, ":", info$pos, ":", info$a0, ":", info$a1)
    in_rsid <- info$rsid %in% hapmap3
    in_key  <- key %in% hapmap3
    keep <- if (sum(in_key) > sum(in_rsid)) in_key else in_rsid
    if (isTRUE(verbose))
      message("Variant-list restriction: ", nrow(info), " -> ", sum(keep))
    info <- info[keep, , drop = FALSE]
    if (nrow(info) == 0)
      stop("Variant restriction removed all variants. Check ID format ",
           "(expected rsIDs or chr:pos:a0:a1 keys).")
  }

  df_beta <- info[, intersect(c("chr", "pos", "a0", "a1", "beta", "beta_se",
                                "n_eff", "freq", "_NUM_ID_"), colnames(info))]
  n_matched <- nrow(df_beta)

  ## --- SD-ratio QC (when freq available) ---------------------------------
  n_sd_flagged <- 0L
  if ("freq" %in% colnames(df_beta)) {
    sd_y_est <- as.numeric(stats::quantile(
      sqrt(0.5 * (df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2)),
      0.01, na.rm = TRUE))
    sd_ss   <- sd_y_est / sqrt(df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2)
    sd_geno <- sqrt(2 * df_beta$freq * (1 - df_beta$freq))
    ratio   <- sd_ss / sd_geno
    bad <- ratio < 0.5 | ratio > 2 | is.na(ratio)
    n_sd_flagged <- sum(bad)
    if (isTRUE(verbose))
      message("SD-ratio QC: flagged ", n_sd_flagged, " / ", nrow(df_beta),
              " variants outside [0.5, 2]")
    if (n_sd_flagged > 0.1 * nrow(df_beta))
      message(">10% flagged; consider shrink_corr = 0.4.")
    df_beta <- df_beta[!bad, , drop = FALSE]
  }

  ## --- Optional systematic thinning ---------------------------------------
  df_beta <- df_beta[order(df_beta$chr, df_beta$pos), ]
  if (!is.null(max_variants) && nrow(df_beta) > max_variants) {
    step <- ceiling(nrow(df_beta) / max_variants)
    df_beta <- df_beta[seq(1, nrow(df_beta), by = step), ]
    if (isTRUE(verbose))
      message("Thinned to ", nrow(df_beta), " variants (every ", step, ").")
  }

  ## Drop monomorphic LD-reference variants (zero-variance columns make
  ## the LD matrix undefined)
  maf_ref <- bigsnpr::snp_MAF(obj.ref$genotypes, ind.col = df_beta$`_NUM_ID_`,
                              ncores = ncores)
  keep_poly <- is.finite(maf_ref) & maf_ref > 0
  if (any(!keep_poly)) {
    if (isTRUE(verbose))
      message("Dropping ", sum(!keep_poly), " monomorphic variants.")
    df_beta <- df_beta[keep_poly, , drop = FALSE]
  }
  if (nrow(df_beta) < 10)
    stop("Fewer than 10 polymorphic matched variants remain; cannot ",
         "estimate LD. Check the LD reference and sumstats build.")

  ## --- LD matrix -----------------------------------------------------------
  if (isTRUE(verbose)) message("Computing LD matrix (", nrow(df_beta),
                               " variants)...")
  corr <- bigsnpr::snp_cor(obj.ref$genotypes, ind.col = df_beta$`_NUM_ID_`,
                           ncores = ncores)

  ## --- LDSC for h2 init -----------------------------------------------------
  if (isTRUE(verbose)) message("LD score regression for h2 init...")
  ldsc <- tryCatch(
    bigsnpr::snp_ldsc2(corr, df_beta, blocks = 200, ncores = ncores),
    error = function(e) NULL)
  if (is.null(ldsc) || !is.finite(ldsc[["h2"]])) {
    ## Fallback moment estimator: E[chi2] ~ 1 + n_eff * h2 / M under no LD
    chi2 <- (df_beta$beta / df_beta$beta_se)^2
    h2_init <- max(1e-4, nrow(df_beta) * (mean(chi2, na.rm = TRUE) - 1) /
                     mean(df_beta$n_eff, na.rm = TRUE))
    message("LD score regression failed on these data; falling back to a ",
            "moment estimate of h2 (", signif(h2_init, 3), ").")
  } else {
    h2_init <- ldsc[["h2"]]
    if (isTRUE(verbose))
      message(sprintf("  LDSC h2 = %.4f (SE = %.4f)", h2_init, ldsc[["h2_se"]]))
  }

  ## --- LDpred2-auto chains ---------------------------------------------------
  if (isTRUE(verbose)) message("Running LDpred2-auto (", n_chains, " chains)...")
  corr_sfbm <- bigsnpr::as_SFBM(corr)
  multi_auto <- bigsnpr::snp_ldpred2_auto(
    corr_sfbm, df_beta,
    h2_init         = h2_init,
    vec_p_init      = bigsnpr::seq_log(p_range[1], p_range[2], n_chains),
    burn_in         = burn_in,
    num_iter        = num_iter,
    allow_jump_sign = allow_jump_sign,
    shrink_corr     = shrink_corr,
    use_MLE         = use_mle,
    ncores          = ncores
  )

  ## --- Chain convergence filter ----------------------------------------------
  range_per_chain <- vapply(multi_auto,
                            function(auto) diff(range(auto$corr_est)),
                            numeric(1))
  thr <- chain_filter_quantile *
    stats::quantile(range_per_chain, 0.95, na.rm = TRUE)
  keep <- which(range_per_chain > thr)
  if (length(keep) == 0) {
    warning("No LDpred2-auto chains passed the convergence filter; ",
            "returning the average across all chains. Inspect chain ",
            "diagnostics before trusting these effects (consider ",
            "HapMap3 restriction, larger burn_in/num_iter, or an ",
            "ancestry-matched LD reference).")
    keep <- seq_along(multi_auto)
  } else if (length(keep) < 3) {
    warning("Fewer than 3 LDpred2-auto chains passed the convergence ",
            "filter. Consider restricting to HapMap3 variants, increasing ",
            "burn_in/num_iter, or checking LD-reference ancestry match.")
  }
  beta_auto <- rowMeans(vapply(multi_auto[keep],
                               function(auto) auto$beta_est,
                               numeric(nrow(df_beta))))
  if (isTRUE(verbose))
    message("Kept ", length(keep), " / ", length(multi_auto), " chains.")

  beta_df <- data.frame(chr = df_beta$chr, pos = df_beta$pos,
                        a0 = df_beta$a0, a1 = df_beta$a1, beta = beta_auto,
                        stringsAsFactors = FALSE)

  ## --- Optional scoring of a target cohort ------------------------------------
  score <- NULL
  if (!is.null(genotypes)) {
    obj.tgt <- attach_ref(genotypes, "genotypes")
    map.tgt <- stats::setNames(obj.tgt$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
    map.tgt$chr <- normalize_chr(map.tgt$chr)
    matched <- bigsnpr::snp_match(
      data.frame(chr = df_beta$chr, pos = df_beta$pos, a0 = df_beta$a0,
                 a1 = df_beta$a1, beta = beta_auto),
      map.tgt)
    beta_full <- rep(0, ncol(obj.tgt$genotypes))
    beta_full[matched$`_NUM_ID_`] <- matched$beta
    if (isTRUE(verbose))
      message("Scoring ", nrow(obj.tgt$fam), " target samples on ",
              nrow(matched), " matched variants...")
    if (isTRUE(impute_mean)) {
      ## Mean-impute in chunks to bound memory
      G <- obj.tgt$genotypes
      chunk <- 10000L
      for (start in seq(1, ncol(G), by = chunk)) {
        end <- min(start + chunk - 1, ncol(G))
        block <- G[, start:end, drop = FALSE]
        if (anyNA(block)) {
          cm <- colMeans(block, na.rm = TRUE)
          for (j in seq_len(ncol(block))) {
            idx <- is.na(block[, j])
            if (any(idx)) block[idx, j] <- cm[j]
          }
          G[, start:end] <- block
        }
      }
    }
    raw <- as.numeric(bigstatsr::big_prodVec(obj.tgt$genotypes, beta_full,
                                           ncores = ncores))
    score <- as.numeric(scale(raw))
    names(score) <- obj.tgt$fam$sample.ID
  }

  list(beta = beta_df, score = score, h2 = h2_init, chains_kept = keep,
       qc = list(n_matched = n_matched, n_sd_flagged = n_sd_flagged,
                 n_final = nrow(df_beta)))
}


#' Score a published polygenic-score panel on a dosage matrix
#'
#' Applies a published set of per-variant weights (e.g. a PGS Catalog
#' scoring file, which is already LD-aware) directly to a genotype
#' dosage matrix. Matching is by `chr:pos:ref:alt` key (with automatic
#' allele-flip handling: when the genotype stores the opposite
#' orientation, the weight is negated) or by rsID when present.
#'
#' @param weights A data.frame of variant weights. Recognized columns:
#' `chr`/`chr_name`, `pos`/`chr_position`, `a1`/`effect_allele`,
#' `a0`/`other_allele`, `weight`/`effect_weight`, and optionally
#' `rsid`. PGS Catalog scoring files can be read directly with
#' [utils::read.delim()] (skip the header comment lines).
#' @param genotypes A numeric dosage matrix (variants x samples, dosages
#' in 0-2) whose rownames are either `chr:pos:ref:alt` keys or rsIDs
#' matching `weights`. A [bigsnpr::snp_attach()] `bigSNP` object is also
#' accepted (requires bigsnpr).
#' @param impute_mean Logical: mean-impute missing dosages per variant.
#' Default `TRUE`.
#' @param scale_score Logical: scale the resulting score to mean 0 /
#' sd 1. Default `TRUE`.
#'
#' @return A list with:
#' \describe{
#'   \item{score}{named numeric vector of per-sample scores.}
#'   \item{matched}{data.frame of matched variants with the weights
#'   used (after any flip correction).}
#'   \item{n_input}{number of weight-panel variants.}
#' }
#' @export
#'
#' @examples
#' # Toy dosage matrix: 4 variants x 50 samples
#' dos <- matrix(rbinom(4 * 50, 2, 0.3), nrow = 4,
#'               dimnames = list(c("1:100:A:G", "1:200:C:T",
#'                                 "1:300:G:A", "2:400:T:C"),
#'                               paste0("S", 1:50)))
#' wts <- data.frame(chr = c(1, 1, 2), pos = c(100, 300, 400),
#'                   a1 = c("A", "G", "T"), a0 = c("G", "A", "C"),
#'                   weight = c(0.1, -0.2, 0.05))
#' out <- score_pgs_panel(wts, dos)
#' plot(out$score)
score_pgs_panel <- function(weights, genotypes, impute_mean = TRUE,
                            scale_score = TRUE) {
  weights <- as.data.frame(weights)

  ## --- Standardize weight columns ------------------------------------------
  aliases <- list(
    chr    = c("chr", "chr_name", "chromosome", "CHR"),
    pos    = c("pos", "chr_position", "base_pair_location", "BP", "POS"),
    a1     = c("a1", "effect_allele", "ALT", "alt"),
    a0     = c("a0", "other_allele", "REF", "ref"),
    weight = c("weight", "effect_weight", "beta", "BETA"),
    rsid   = c("rsid", "rs_id", "SNP", "ID")
  )
  w_list <- list()
  for (col in names(aliases)) {
    hit <- intersect(aliases[[col]], colnames(weights))
    if (length(hit) > 0) w_list[[col]] <- weights[[hit[1]]]
  }
  w <- as.data.frame(w_list, stringsAsFactors = FALSE)
  if (is.null(w$weight))
    stop("weights needs an effect-weight column (weight/effect_weight/beta).")
  has_pos <- !is.null(w$chr) && !is.null(w$pos) &&
    !is.null(w$a1) && !is.null(w$a0)
  if (!has_pos && is.null(w$rsid))
    stop("weights needs either chr/pos/a0/a1 columns or an rsid column.")
  if (has_pos) {
    w$chr <- as.character(w$chr)
    w$pos <- as.integer(w$pos)
    w$a0  <- toupper(as.character(w$a0))
    w$a1  <- toupper(as.character(w$a1))
  }
  w$weight <- as.numeric(w$weight)
  w <- w[is.finite(w$weight), , drop = FALSE]
  n_input <- nrow(w)

  ## --- bigSNP input: extract map + delegate to bigsnpr scoring --------------
  if (inherits(genotypes, "bigSNP")) {
    if (!requireNamespace("bigsnpr", quietly = TRUE))
      stop("Scoring a bigSNP object requires the bigsnpr package.")
    map <- stats::setNames(genotypes$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
    ss <- data.frame(chr = w$chr, pos = w$pos, a0 = w$a0, a1 = w$a1,
                     beta = w$weight)
    matched <- bigsnpr::snp_match(ss, map)
    beta_full <- rep(0, ncol(genotypes$genotypes))
    beta_full[matched$`_NUM_ID_`] <- matched$beta
    raw <- as.numeric(bigstatsr::big_prodVec(genotypes$genotypes, beta_full))
    sc <- if (isTRUE(scale_score)) as.numeric(scale(raw)) else raw
    names(sc) <- genotypes$fam$sample.ID
    return(list(score = sc, matched = matched, n_input = n_input))
  }

  ## --- Matrix input -----------------------------------------------------------
  G <- as.matrix(genotypes)
  if (is.null(rownames(G)))
    stop("genotypes must have variant-ID rownames (chr:pos:ref:alt or rsID).")
  if (is.null(colnames(G)))
    colnames(G) <- paste0("sample", seq_len(ncol(G)))

  ## Build keys
  key_of <- function(chr, pos, ref, alt) paste0(chr, ":", pos, ":", ref, ":", alt)
  g_ids <- rownames(G)
  g_has_colons <- grepl(":", g_ids, fixed = TRUE)

  matched_idx <- integer(0)
  matched_flip <- logical(0)

  if (has_pos && any(g_has_colons)) {
    w_key      <- key_of(w$chr, w$pos, w$a0, w$a1)          # ref:alt = a0:a1
    w_key_flip <- key_of(w$chr, w$pos, w$a1, w$a0)
    idx_same <- match(g_ids, w_key)
    idx_flip <- match(g_ids, w_key_flip)
    ## Per weight row, find genotype row (same orientation preferred)
    gi_same <- match(w_key, g_ids)
    gi_flip <- match(w_key_flip, g_ids)
    use_flip <- is.na(gi_same) & !is.na(gi_flip)
    gi <- ifelse(is.na(gi_same), gi_flip, gi_same)
    keep <- !is.na(gi)
    matched_idx <- gi[keep]
    matched_flip <- use_flip[keep]
    w_m <- w[keep, , drop = FALSE]
  } else if (!is.null(w$rsid)) {
    gi <- match(w$rsid, g_ids)
    keep <- !is.na(gi)
    matched_idx <- gi[keep]
    matched_flip <- rep(FALSE, sum(keep))
    w_m <- w[keep, , drop = FALSE]
  } else {
    stop("Cannot match: weights lack chr/pos/a0/a1 and genotype rownames ",
         "are not chr:pos:ref:alt keys, and no rsid column is available.")
  }

  if (nrow(w_m) == 0)
    stop("No weight-panel variants matched genotype rownames. Check ID ",
         "formats (expected chr:pos:ref:alt or rsIDs).")
  if (nrow(w_m) < 0.5 * n_input)
    message("Only ", nrow(w_m), " / ", n_input,
            " weight-panel variants matched the genotype matrix.")

  eff <- w_m$weight
  eff[matched_flip] <- -eff[matched_flip]

  Gsub <- G[matched_idx, , drop = FALSE]
  if (isTRUE(impute_mean) && anyNA(Gsub)) {
    cm <- rowMeans(Gsub, na.rm = TRUE)
    for (j in seq_len(nrow(Gsub))) {
      idx <- is.na(Gsub[j, ])
      if (any(idx)) Gsub[j, idx] <- cm[j]
    }
  } else if (anyNA(Gsub)) {
    stop("genotypes contain missing dosages; set impute_mean = TRUE or ",
         "impute beforehand.")
  }

  raw <- as.numeric(crossprod(eff, Gsub))
  sc <- if (isTRUE(scale_score)) as.numeric(scale(raw)) else raw
  names(sc) <- colnames(G)

  w_m$weight_used <- eff
  w_m$flipped <- matched_flip
  list(score = sc, matched = w_m, n_input = n_input)
}


#' Check first-stage instrument strength (partial F)
#'
#' Regresses the exposure `X` on the instrument `G` (plus optional
#' covariates) and reports the partial F statistic and partial R2 for
#' the instrument — the weak-instrument diagnostic. The partial F
#' (not the overall model F) is the relevant quantity when covariates
#' such as ancestry PCs explain much of the exposure.
#'
#' @param G Numeric vector (length n): the exposure instrument
#' (e.g. a polygenic score from [build_prs_ldpred2()] or
#' [score_pgs_panel()]).
#' @param X Numeric vector (length n): the exposure.
#' @param covariates Optional data.frame or matrix (n rows) of
#' covariates to partial out (e.g. ancestry PCs).
#' @param min_f Numeric: weak-instrument threshold on the partial F.
#' Default 10 (Staiger-Stock rule of thumb).
#'
#' @return A list with `F` (partial F statistic), `df1`, `df2`,
#' `pvalue`, `partial_r2`, `n`, and `weak` (`TRUE` when F < `min_f`).
#' @export
#'
#' @examples
#' n <- 300
#' G <- rnorm(n)
#' X <- 0.3 * G + rnorm(n)
#' pcs <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, paste0("PC", 1:3)))
#' check_instrument_strength(G, X, covariates = pcs)
check_instrument_strength <- function(G, X, covariates = NULL, min_f = 10) {
  G <- as.numeric(G)
  X <- as.numeric(X)
  n <- length(X)
  if (length(G) != n) stop("G and X must have the same length.")

  if (!is.null(covariates)) {
    C <- as.matrix(as.data.frame(covariates))
    if (nrow(C) != n) stop("covariates must have n rows.")
    keep <- stats::complete.cases(data.frame(G = G, X = X, C))
    if (sum(keep) < n) {
      message("Dropping ", n - sum(keep), " samples with missing values.")
      G <- G[keep]; X <- X[keep]; C <- C[keep, , drop = FALSE]
      n <- sum(keep)
    }
    fit_full <- stats::lm(X ~ G + C)
    fit_null <- stats::lm(X ~ C)
  } else {
    keep <- stats::complete.cases(data.frame(G = G, X = X))
    G <- G[keep]; X <- X[keep]; n <- sum(keep)
    fit_full <- stats::lm(X ~ G)
    fit_null <- stats::lm(X ~ 1)
  }

  af <- stats::anova(fit_null, fit_full)
  f_stat <- af[["F"]][2]
  pval   <- af[["Pr(>F)"]][2]
  rss0 <- sum(stats::resid(fit_null)^2)
  rss1 <- sum(stats::resid(fit_full)^2)
  partial_r2 <- (rss0 - rss1) / rss0

  if (is.finite(f_stat) && f_stat < min_f)
    warning("First-stage partial F = ", round(f_stat, 2), " < ", min_f,
            ": weak instrument. Estimates relying on G may be biased.")

  list(F = f_stat, df1 = af[["Df"]][2], df2 = fit_full$df.residual,
       pvalue = pval, partial_r2 = partial_r2, n = n,
       weak = is.finite(f_stat) && f_stat < min_f)
}
