# Tests for R/instruments_mediator.R:
# call_cis_eqtls, build_mediator_instruments

test_that("call_cis_eqtls recovers a planted cis-eQTL", {
  set.seed(11)
  n <- 80
  dos <- matrix(rbinom(200 * n, 2, 0.3), nrow = 200,
                dimnames = list(paste0("1:", 1:200, ":A:G"), NULL))
  expr <- matrix(rnorm(20 * n), nrow = 20,
                 dimnames = list(paste0("Gene", 1:20), NULL))
  expr[1, ] <- expr[1, ] + 0.8 * scale(dos[5, ])   # strong cis-eQTL
  gp <- data.frame(gene = paste0("Gene", 1:20), chr = "1",
                   tss = seq(1, 191, by = 10))
  sp <- data.frame(snp = rownames(dos), chr = "1", pos = 1:200)
  hits <- call_cis_eqtls(expr, dos, gp, sp)
  expect_true(all(c("all", "best") %in% names(hits)))
  best1 <- hits$best[hits$best$gene == "Gene1", ]
  expect_equal(best1$snp, "1:5:A:G")
  expect_lt(best1$p, 1e-3)
})

test_that("call_cis_eqtls respects the cis window", {
  set.seed(12)
  n <- 60
  dos <- matrix(rbinom(50 * n, 2, 0.3), nrow = 50,
                dimnames = list(paste0("1:", 1:50, ":A:G"), NULL))
  expr <- matrix(rnorm(2 * n), nrow = 2,
                 dimnames = list(c("Gene1", "Gene2"), NULL))
  gp <- data.frame(gene = c("Gene1", "Gene2"), chr = "1", tss = c(10, 40))
  sp <- data.frame(snp = rownames(dos), chr = "1", pos = 1:50)
  # window of 2 bp around tss=10 covers SNPs at pos 8..12 only
  hits <- call_cis_eqtls(expr, dos, gp, sp, cis_dist = 2)
  g1 <- hits$all[hits$all$gene == "Gene1", ]
  expect_true(all(g1$snp %in% paste0("1:", 8:12, ":A:G")))
})

test_that("call_cis_eqtls errors on mismatched inputs", {
  dos <- matrix(rbinom(20 * 30, 2, 0.3), nrow = 20,
                dimnames = list(paste0("1:", 1:20, ":A:G"), NULL))
  expr <- matrix(rnorm(5 * 25), nrow = 5,
                 dimnames = list(paste0("Gene", 1:5), NULL))
  gp <- data.frame(gene = paste0("Gene", 1:5), chr = "1", tss = 1:5)
  sp <- data.frame(snp = rownames(dos), chr = "1", pos = 1:20)
  expect_error(call_cis_eqtls(expr, dos, gp, sp),
               "same number of samples")
})

test_that("build_mediator_instruments passes genes with real cis signal", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 150; p <- 80
  dos <- matrix(rbinom(p * n, 2, 0.3), nrow = p,
                dimnames = list(paste0("1:", 1:p, ":A:G"),
                                paste0("S", 1:n)))
  expr <- matrix(rnorm(4 * n), nrow = 4,
                 dimnames = list(paste0("Gene", 1:4), paste0("S", 1:n)))
  expr[1, ] <- expr[1, ] + 0.9 * scale(dos[5, ]) + 0.9 * scale(dos[15, ])
  gp <- data.frame(gene = paste0("Gene", 1:4), chr = "1",
                   tss = c(5, 25, 45, 65))
  sp <- data.frame(snp = rownames(dos), chr = "1", pos = 1:p)
  fit <- build_mediator_instruments(expr, dos, gp, sp, seed = 1)
  expect_true(all(c("Gm", "qc", "weights") %in% names(fit)))
  qc <- fit$qc
  expect_true(qc$pass[qc$gene == "Gene1"])
  expect_false(any(qc$pass[qc$gene != "Gene1"]))
  # Gm is genes x samples, restricted to passing genes
  expect_equal(ncol(fit$Gm), n)
  expect_equal(rownames(fit$Gm), "Gene1")
})

test_that("build_mediator_instruments is reproducible under the same seed", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 150; p <- 80
  dos <- matrix(rbinom(p * n, 2, 0.3), nrow = p,
                dimnames = list(paste0("1:", 1:p, ":A:G"),
                                paste0("S", 1:n)))
  expr <- matrix(rnorm(4 * n), nrow = 4,
                 dimnames = list(paste0("Gene", 1:4), paste0("S", 1:n)))
  expr[1, ] <- expr[1, ] + 0.9 * scale(dos[5, ]) + 0.9 * scale(dos[15, ])
  gp <- data.frame(gene = paste0("Gene", 1:4), chr = "1",
                   tss = c(5, 25, 45, 65))
  sp <- data.frame(snp = rownames(dos), chr = "1", pos = 1:p)
  f1 <- build_mediator_instruments(expr, dos, gp, sp, seed = 7)
  f2 <- build_mediator_instruments(expr, dos, gp, sp, seed = 7)
  expect_identical(f1$Gm, f2$Gm)
  expect_identical(f1$qc$cv_r2, f2$qc$cv_r2)
})

test_that("build_mediator_instruments warns when no gene passes", {
  skip_if_not_installed("glmnet")
  set.seed(13)
  n <- 100; p <- 40
  dos <- matrix(rbinom(p * n, 2, 0.3), nrow = p,
                dimnames = list(paste0("1:", 1:p, ":A:G"), NULL))
  expr <- matrix(rnorm(3 * n), nrow = 3,
                 dimnames = list(paste0("Gene", 1:3), NULL))
  gp <- data.frame(gene = paste0("Gene", 1:3), chr = "1", tss = c(5, 20, 35))
  sp <- data.frame(snp = rownames(dos), chr = "1", pos = 1:p)
  expect_warning(
    fit <- build_mediator_instruments(expr, dos, gp, sp, seed = 1),
    "No genes passed"
  )
  expect_null(fit$Gm)
})
