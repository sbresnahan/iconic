## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.5, comment = "#>"
)
set.seed(1)
library(iconic)

## ----install-bioc, eval=FALSE-------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("iconic")

## ----install-github, eval=FALSE-----------------------------------------------
# BiocManager::install("sbresnahan/iconic")

## ----qc-sumstats--------------------------------------------------------------
ss <- data.frame(
  chromosome = 1, base_pair_location = 1000 + 0:9 * 1000L,
  effect_allele = rep(c("A", "C"), 5), other_allele = rep(c("G", "T"), 5),
  beta = rnorm(10, 0, 0.05), standard_error = 0.01,
  effect_allele_frequency = runif(10, 0.1, 0.5)
)
out <- qc_gwas_sumstats(ss, n_eff = 50000)
out$qc

## ----ldpred2, eval=requireNamespace("bigsnpr", quietly = TRUE)----------------
# Tiny fake genotype panel standing in for an LD reference
# (snp_fake genotypes are all-missing mocks; fill with random dosages)
fake <- bigsnpr::snp_fake(100, 500)
fake$genotypes[] <- rbinom(100 * 500, 2, 0.3)
fake$map$chromosome <- 1L
fake$map$physical.pos <- sort(sample(1:1e6, 500))
ss <- data.frame(
  chr = fake$map$chromosome, pos = fake$map$physical.pos,
  a0 = fake$map$allele2, a1 = fake$map$allele1,
  beta = rnorm(500, 0, 0.05), beta_se = 0.02, n_eff = 20000
)
prs <- build_prs_ldpred2(ss, ld_ref = fake, n_chains = 3,
                         burn_in = 50, num_iter = 50, verbose = FALSE)
str(prs$beta)

## ----score-panel--------------------------------------------------------------
dos <- matrix(rbinom(4 * 50, 2, 0.3), nrow = 4,
              dimnames = list(c("1:100:A:G", "1:200:C:T",
                                "1:300:G:A", "2:400:T:C"),
                              paste0("S", 1:50)))
wts <- data.frame(chr = c(1, 1, 2), pos = c(100, 300, 400),
                  a1 = c("A", "G", "T"), a0 = c("G", "A", "C"),
                  weight = c(0.1, -0.2, 0.05))
out <- score_pgs_panel(wts, dos)
head(out$score)

## ----instrument-strength------------------------------------------------------
n <- 300
G <- rnorm(n)
X <- 0.3 * G + rnorm(n)
pcs <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, paste0("PC", 1:3)))
check_instrument_strength(G, X, covariates = pcs)

## ----cis-eqtls----------------------------------------------------------------
n <- 80
dos <- matrix(rbinom(200 * n, 2, 0.3), nrow = 200,
              dimnames = list(paste0("1:", 1:200, ":A:G"), NULL))
expr <- matrix(rnorm(20 * n), nrow = 20,
               dimnames = list(paste0("Gene", 1:20), NULL))
expr[1, ] <- expr[1, ] + 0.5 * dos[5, ]   # one true cis-eQTL
gp <- data.frame(gene = paste0("Gene", 1:20), chr = "1",
                 tss = seq(1, 191, by = 10))
sp <- data.frame(snp = rownames(dos), chr = "1", pos = 1:200)
hits <- call_cis_eqtls(expr, dos, gp, sp)
head(hits$best)

## ----mediator-instruments, eval=requireNamespace("glmnet", quietly = TRUE)----
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
fit$qc

## ----beta-to-m----------------------------------------------------------------
b <- matrix(runif(200, 0.01, 0.99), nrow = 20,
            dimnames = list(paste0("cg", 1:20), paste0("S", 1:10)))
m <- beta_to_m(b)
range(m)

## ----residualize--------------------------------------------------------------
x <- matrix(rnorm(100 * 40), nrow = 100,
            dimnames = list(paste0("f", 1:100), paste0("S", 1:40)))
batch <- factor(rep(c("A", "B"), each = 20))
cv <- model.matrix(~ batch)[, -1, drop = FALSE]
xr <- residualize_matrix(x, cv)
cor(as.numeric(xr[1, ]), as.numeric(cv))

## ----build-w------------------------------------------------------------------
x <- matrix(rnorm(500 * 60), nrow = 500,
            dimnames = list(paste0("f", 1:500), paste0("S", 1:60)))
w <- build_w_pcs(x, n_pcs = 5)
dim(w$W)
w$variance_explained

## ----fusion-weights-----------------------------------------------------------
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
dim(out$predicted)

## ----summarized-experiment, eval=requireNamespace("SummarizedExperiment", quietly = TRUE)----
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(expr = matrix(rnorm(20 * 60), 20, 60,
                              dimnames = list(paste0("gene", 1:20),
                                              paste0("S", 1:60)))),
  colData = S4Vectors::DataFrame(
    bmi = rnorm(60), prs = rnorm(60),
    nc1 = rnorm(60), nc2 = rnorm(60), age = rnorm(60))
)
data <- as_iconic_data(se, assay = "expr", exposure = "bmi",
                       instrument = "prs",
                       negative_controls = c("nc1", "nc2"),
                       covariates = "age")
print(data)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

