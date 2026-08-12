# Tests for R/instruments_exposure.R:
# qc_gwas_sumstats, build_prs_ldpred2, score_pgs_panel, check_instrument_strength

test_that("qc_gwas_sumstats standardises column aliases", {
  ss <- data.frame(
    chromosome = 1, base_pair_location = 1000 + 0:9 * 1000L,
    effect_allele = rep(c("A", "C"), 5), other_allele = rep(c("G", "T"), 5),
    beta = rnorm(10, 0, 0.05), standard_error = 0.01,
    effect_allele_frequency = runif(10, 0.1, 0.5)
  )
  out <- qc_gwas_sumstats(ss, n_eff = 50000)
  expect_true(all(c("chr", "pos", "a0", "a1", "beta", "beta_se") %in%
                    names(out$sumstats)))
  expect_equal(nrow(out$sumstats), 10)
  expect_true(is.data.frame(out$qc))
})

test_that("qc_gwas_sumstats drops ambiguous and extreme-beta variants", {
  # >10 variants so the IQR-based extreme-beta filter is active
  ss <- data.frame(
    chr = 1, pos = 1000 + 0:14 * 1000L,
    a1 = c("A", "C", "G", rep("G", 12)),
    a0 = c("T", "G", "C", rep("A", 12)),   # rows 1-3 ambiguous (AT/CG/GC)
    beta = c(rep(c(0.01, 0.02, 0.03, 0.04), length.out = 14), 50),  # last outlier
    beta_se = 0.01, n_eff = 50000
  )
  out <- qc_gwas_sumstats(ss)
  kept <- out$sumstats
  # ambiguous A/T, C/G, G/C variants removed
  expect_false(any(kept$a1 == "A" & kept$a0 == "T"))
  expect_false(any(kept$a1 == "C" & kept$a0 == "G"))
  # extreme beta outlier removed
  expect_true(all(abs(kept$beta) < 1))
  expect_lt(nrow(kept), 15)
})

test_that("qc_gwas_sumstats handles an empty result gracefully", {
  ss <- data.frame(chr = 1, pos = 1:3, a1 = "A", a0 = "T",
                   beta = 0.1, beta_se = 0.01, n_eff = 1000)
  out <- qc_gwas_sumstats(ss)  # all A/T ambiguous -> 0 rows
  expect_equal(nrow(out$sumstats), 0)
  expect_true(all(c("chr", "pos", "a0", "a1") %in% names(out$sumstats)))
})

test_that("score_pgs_panel scores match a manual weighted sum", {
  set.seed(1)
  dos <- matrix(rbinom(4 * 50, 2, 0.3), nrow = 4,
                dimnames = list(c("1:100:A:G", "1:200:C:T",
                                  "1:300:G:A", "2:400:T:C"),
                                paste0("S", 1:50)))
  # Genotype key is chr:pos:ref:alt with dosage counting the alt allele.
  # A weight is used as-is when a0=ref/a1=alt, negated when swapped.
  wts <- data.frame(chr = c(1, 1, 2), pos = c(100, 300, 400),
                    a1 = c("G", "G", "C"), a0 = c("A", "A", "T"),
                    weight = c(0.1, -0.2, 0.05))
  out <- score_pgs_panel(wts, dos, scale_score = FALSE)
  expect_length(out$score, 50)
  expect_equal(nrow(out$matched), 3)
  # manual check on the raw (unscaled) score:
  #  1:100:A:G -> a0=A(ref),a1=G(alt) same       -> +0.1 * dos
  #  1:300:G:A -> a0=A(alt),a1=G(ref) flip       -> -(-0.2)=+0.2 * dos
  #  2:400:T:C -> a0=T(ref),a1=C(alt) same       -> +0.05 * dos
  manual <- 0.1 * dos["1:100:A:G", ] + 0.2 * dos["1:300:G:A", ] +
    0.05 * dos["2:400:T:C", ]
  expect_equal(as.numeric(out$score), as.numeric(manual), tolerance = 1e-8)
})

test_that("score_pgs_panel negates weights on strand flip", {
  set.seed(2)
  dos <- matrix(rbinom(2 * 40, 2, 0.3), nrow = 2,
                dimnames = list(c("1:100:A:G", "1:200:C:T"), NULL))
  # effect allele (a1) = the ref allele of the genotype key -> flip -> negate
  wts <- data.frame(chr = 1, pos = 100, a1 = "A", a0 = "G", weight = 0.5)
  out <- score_pgs_panel(wts, dos, scale_score = FALSE)
  expect_equal(nrow(out$matched), 1)
  expect_equal(as.numeric(out$score), as.numeric(-0.5 * dos["1:100:A:G", ]),
               tolerance = 1e-8)
})

test_that("score_pgs_panel falls back to rsID matching", {
  set.seed(3)
  dos <- matrix(rbinom(2 * 30, 2, 0.3), nrow = 2,
                dimnames = list(c("rs1", "rs2"), NULL))
  wts <- data.frame(rsid = c("rs1", "rs2"), weight = c(1, 2))
  out <- score_pgs_panel(wts, dos, scale_score = FALSE)
  expect_equal(as.numeric(out$score),
               as.numeric(1 * dos["rs1", ] + 2 * dos["rs2", ]),
               tolerance = 1e-8)
})

test_that("check_instrument_strength flags weak instruments", {
  set.seed(4)
  n <- 300
  G <- rnorm(n)
  X_strong <- 0.6 * G + rnorm(n)
  s <- check_instrument_strength(G, X_strong)
  expect_true(s$F > 10)
  expect_false(s$weak)
  expect_gt(s$partial_r2, 0.05)

  X_weak <- 0.01 * G + rnorm(n)
  expect_warning(w <- check_instrument_strength(G, X_weak),
                 "weak instrument")
  expect_true(w$weak)
})

test_that("check_instrument_strength adjusts for covariates", {
  set.seed(5)
  n <- 200
  pcs <- matrix(rnorm(n * 2), n, 2, dimnames = list(NULL, c("PC1", "PC2")))
  G <- rnorm(n)
  X <- 0.5 * G + 0.8 * pcs[, 1] + rnorm(n)
  s <- check_instrument_strength(G, X, covariates = pcs)
  expect_equal(s$n, n)
  expect_true(is.finite(s$F))
})

test_that("build_prs_ldpred2 runs on a fake panel and scores a target", {
  skip_if_not_installed("bigsnpr")
  skip_if_not_installed("bigstatsr")
  set.seed(6)
  fake <- bigsnpr::snp_fake(100, 500)
  fake$genotypes[] <- rbinom(100 * 500, 2, 0.3)
  fake$map$chromosome <- 1L
  fake$map$physical.pos <- sort(sample(1:1e6, 500))
  ss <- data.frame(
    chr = fake$map$chromosome, pos = fake$map$physical.pos,
    a0 = fake$map$allele2, a1 = fake$map$allele1,
    beta = rnorm(500, 0, 0.05), beta_se = 0.02, n_eff = 20000
  )
  prs <- suppressWarnings(
    build_prs_ldpred2(ss, ld_ref = fake, genotypes = fake,
                      n_chains = 3, burn_in = 50, num_iter = 50,
                      verbose = FALSE)
  )
  expect_true(is.data.frame(prs$beta))
  expect_gt(nrow(prs$beta), 0)
  expect_true(!is.null(prs$score))
  expect_length(prs$score, 100)
  expect_true(is.finite(prs$h2))
})
