# Tests for R/negative_controls.R:
# beta_to_m, residualize_matrix, build_w_pcs, apply_fusion_weights

test_that("beta_to_m applies the logit transform with clipping", {
  b <- matrix(c(0.5, 0.25, 0.0, 1.0), nrow = 2, byrow = TRUE,
              dimnames = list(c("cg1", "cg2"), c("S1", "S2")))
  m <- beta_to_m(b)
  expect_equal(m[1, 1], log2(0.5 / 0.5), tolerance = 1e-8)
  expect_equal(m[1, 2], log2(0.25 / 0.75), tolerance = 1e-8)
  # 0 and 1 are clipped to [clip, 1 - clip] before the transform
  expect_true(all(is.finite(m[2, ])))
  expect_equal(m[2, 1], log2(1e-4 / (1 - 1e-4)), tolerance = 1e-6)
})

test_that("beta_to_m drops non-finite rows when requested", {
  b <- matrix(runif(30), nrow = 3,
              dimnames = list(c("cg1", "cg2", "cg3"), NULL))
  b[2, 3] <- NA
  m <- beta_to_m(b)
  expect_equal(nrow(m), 2)
  expect_false("cg2" %in% rownames(m))
  m_keep <- beta_to_m(b, drop_nonfinite = FALSE)
  expect_equal(nrow(m_keep), 3)
})

test_that("residualize_matrix removes covariate effects", {
  set.seed(21)
  n <- 40
  batch <- factor(rep(c("A", "B"), each = n / 2))
  cv <- model.matrix(~ batch)[, -1, drop = FALSE]
  x <- matrix(rnorm(50 * n), nrow = 50,
              dimnames = list(paste0("f", 1:50), paste0("S", 1:n)))
  x <- x + matrix(rep(2 * as.numeric(batch == "B"), 50), nrow = 50,
                  byrow = TRUE)
  xr <- residualize_matrix(x, cv)
  expect_equal(dim(xr), dim(x))
  # after residualization the batch covariate explains ~nothing
  cors <- apply(xr, 1, function(r) cor(r, as.numeric(batch == "B")))
  expect_true(all(abs(cors) < 0.05))
})

test_that("residualize_matrix chunked path matches the one-shot path", {
  set.seed(22)
  n <- 25
  cv <- data.frame(batch = rnorm(n))
  x <- matrix(rnorm(30 * n), nrow = 30)
  xr_one <- residualize_matrix(x, cv, chunk_size = 1000)
  xr_chunked <- residualize_matrix(x, cv, chunk_size = 7)
  expect_equal(xr_chunked, xr_one, tolerance = 1e-10)
})

test_that("build_w_pcs returns PCs in samples orientation", {
  set.seed(23)
  x <- matrix(rnorm(500 * 60), nrow = 500,
              dimnames = list(paste0("f", 1:500), paste0("S", 1:60)))
  w <- build_w_pcs(x, n_pcs = 5)
  expect_equal(dim(w$W), c(5, 60))
  expect_length(w$variance_explained, 5)
  expect_true(all(w$variance_explained >= 0))
  expect_lte(sum(w$variance_explained), 100 + 1e-6)  # percentages
})

test_that("build_w_pcs honours n_pcs and prefix", {
  set.seed(24)
  x <- matrix(rnorm(200 * 30), nrow = 200)
  w <- build_w_pcs(x, n_pcs = 4, prefix = "NC")
  expect_equal(dim(w$W), c(4, 30))
  expect_true(all(grepl("^NC", rownames(w$W))))
  expect_error(build_w_pcs(x, n_pcs = 30), "n_pcs must be smaller")
})

test_that("apply_fusion_weights predicts from in-memory weights", {
  set.seed(25)
  dos <- matrix(rbinom(10 * 50, 2, 0.3), nrow = 10,
                dimnames = list(paste0("rs", 1:10), paste0("S", 1:50)))
  mk_wgt <- function(snps, w) {
    list(wgt.matrix = matrix(w, ncol = 1, dimnames = list(NULL, "enet")),
         snps = data.frame(V2 = snps),
         cv.performance = matrix(0.2, nrow = 1,
                                 dimnames = list("rsq", "enet")))
  }
  wlist <- list(GeneA = mk_wgt(paste0("rs", 1:5), rep(0.1, 5)),
                GeneB = mk_wgt(paste0("rs", 6:10), rep(-0.2, 5)))
  pos <- data.frame(ID = c("GeneA", "GeneB"), WGT = c("a.wgt.RDat", "b.wgt.RDat"))
  out <- apply_fusion_weights(dos, pos = pos, weights = wlist)
  expect_equal(dim(out$predicted), c(2, 50))
  expect_equal(rownames(out$predicted), c("GeneA", "GeneB"))
  # manual check: prediction is the weighted sum of dosages
  expect_equal(as.numeric(out$predicted["GeneA", ]),
               as.numeric(0.1 * colSums(dos[1:5, ])), tolerance = 1e-8)
})

test_that("apply_fusion_weights picks the best model by cv rsq", {
  set.seed(26)
  dos <- matrix(rbinom(4 * 30, 2, 0.3), nrow = 4,
                dimnames = list(paste0("rs", 1:4), NULL))
  wgt <- list(
    wgt.matrix = matrix(c(1, 1, 1, 1,   2, 2, 2, 2), nrow = 4,
                        dimnames = list(NULL, c("lasso", "enet"))),
    snps = data.frame(V2 = paste0("rs", 1:4)),
    cv.performance = rbind(rsq = c(lasso = 0.1, enet = 0.4))
  )
  pos <- data.frame(ID = "GeneA", WGT = "a.wgt.RDat")
  out <- apply_fusion_weights(dos, pos = pos, weights = list(GeneA = wgt))
  # enet has the higher cv rsq -> weights of 2 used
  expect_equal(as.numeric(out$predicted["GeneA", ]),
               as.numeric(2 * colSums(dos)), tolerance = 1e-8)
  expect_equal(out$info$best_method[out$info$gene == "GeneA"], "enet")
})
