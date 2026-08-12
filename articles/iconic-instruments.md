# Genetic Instruments and Negative Controls from External Data

## Introduction

`iconic` identifies natural direct and indirect effects (NDE/NIE) in
observational omics data by combining *genetic instruments* (for the
exposure and, optionally, the mediator) with *negative controls* that
proxy the unmeasured confounders. In practice these inputs come from
external resources: GWAS summary statistics for exposure instruments,
cis-eQTL scans or transcriptome-wide association study (TWAS) weights
for mediator instruments, and high-dimensional omics panels (for example
DNA methylation) for negative-control construction.

This vignette demonstrates the helper functions that bridge those
external resources into an `iconic_data` object:

- [`qc_gwas_sumstats()`](https://seantbresnahan.com/iconic/reference/qc_gwas_sumstats.md),
  [`build_prs_ldpred2()`](https://seantbresnahan.com/iconic/reference/build_prs_ldpred2.md),
  [`score_pgs_panel()`](https://seantbresnahan.com/iconic/reference/score_pgs_panel.md),
  and
  [`check_instrument_strength()`](https://seantbresnahan.com/iconic/reference/check_instrument_strength.md)
  for exposure instruments;
- [`call_cis_eqtls()`](https://seantbresnahan.com/iconic/reference/call_cis_eqtls.md)
  and
  [`build_mediator_instruments()`](https://seantbresnahan.com/iconic/reference/build_mediator_instruments.md)
  for mediator instruments;
- [`beta_to_m()`](https://seantbresnahan.com/iconic/reference/beta_to_m.md),
  [`residualize_matrix()`](https://seantbresnahan.com/iconic/reference/residualize_matrix.md),
  [`build_w_pcs()`](https://seantbresnahan.com/iconic/reference/build_w_pcs.md),
  and
  [`apply_fusion_weights()`](https://seantbresnahan.com/iconic/reference/apply_fusion_weights.md)
  for negative-control panel construction;
- [`as_iconic_data()`](https://seantbresnahan.com/iconic/reference/as_iconic_data.md)
  for importing a `SummarizedExperiment` directly.

All examples use small simulated inputs so the vignette builds quickly;
on real data the same calls apply unchanged.

## Installation

`iconic` is being submitted to Bioconductor. Once accepted, install the
release version with:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("iconic")
```

The development version is available from GitHub:

``` r

BiocManager::install("sbresnahan/iconic")
```

## Exposure instruments from GWAS summary statistics

### Quality-control the summary statistics

[`qc_gwas_sumstats()`](https://seantbresnahan.com/iconic/reference/qc_gwas_sumstats.md)
standardises the column names of a GWAS summary statistics table
(recognising common aliases such as
`chromosome`/`base_pair_location`/`effect_allele`), and applies the
usual filters: missing values, extreme effect sizes, strand-ambiguous
A/T and C/G variants, and the standard-deviation-ratio check that
compares the implied phenotypic variance of each variant against the
bulk distribution.

``` r

ss <- data.frame(
  chromosome = 1, base_pair_location = 1000 + 0:9 * 1000L,
  effect_allele = rep(c("A", "C"), 5), other_allele = rep(c("G", "T"), 5),
  beta = rnorm(10, 0, 0.05), standard_error = 0.01,
  effect_allele_frequency = runif(10, 0.1, 0.5)
)
out <- qc_gwas_sumstats(ss, n_eff = 50000)
out$qc
```

    #>   n_input dropped_missing dropped_extreme_beta dropped_ambiguous
    #> 1      10               0                    0                 0
    #>   dropped_sd_ratio n_output
    #> 1                0       10

### Build a polygenic score with LDpred2

When the `bigsnpr` package is available,
[`build_prs_ldpred2()`](https://seantbresnahan.com/iconic/reference/build_prs_ldpred2.md)
matches the summary statistics to an LD reference panel, runs
LDpred2-auto (with a moment-estimator fallback for the heritability),
filters out non-convergent chains, and can optionally score a target
genotype panel.

``` r

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
```

    #> 'data.frame':    500 obs. of  5 variables:
    #>  $ chr : chr  "chr1" "chr1" "chr1" "chr1" ...
    #>  $ pos : int  641 2214 3291 3328 4408 5202 16397 18023 18235 19677 ...
    #>  $ a0  : chr  "C" "T" "C" "C" ...
    #>  $ a1  : chr  "T" "C" "T" "T" ...
    #>  $ beta: num  -0.01351 -0.01108 -0.02115 0.00702 0.00673 ...

On real data the chain-convergence filter and the sd-ratio QC do the
heavy lifting; on this tiny fake panel the warning about non-convergent
chains is expected.

### Score a panel of pre-computed weights

If you already have effect weights (for example from the PGS Catalog or
a published score),
[`score_pgs_panel()`](https://seantbresnahan.com/iconic/reference/score_pgs_panel.md)
applies them to a dosage matrix directly, matching variants by
`chr:pos:ref:alt` key with automatic strand-flip handling.

``` r

dos <- matrix(rbinom(4 * 50, 2, 0.3), nrow = 4,
              dimnames = list(c("1:100:A:G", "1:200:C:T",
                                "1:300:G:A", "2:400:T:C"),
                              paste0("S", 1:50)))
wts <- data.frame(chr = c(1, 1, 2), pos = c(100, 300, 400),
                  a1 = c("A", "G", "T"), a0 = c("G", "A", "C"),
                  weight = c(0.1, -0.2, 0.05))
out <- score_pgs_panel(wts, dos)
head(out$score)
```

    #>          S1          S2          S3          S4          S5          S6 
    #>  0.66042273 -2.06296996 -0.02042545  0.31999864 -0.70127362 -1.38212179

### Verify instrument strength

Whatever the source of the instrument, check the first-stage partial
F-statistic before relying on it. The conventional rule of thumb is F
\>= 10.

``` r

n <- 300
G <- rnorm(n)
X <- 0.3 * G + rnorm(n)
pcs <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, paste0("PC", 1:3)))
check_instrument_strength(G, X, covariates = pcs)
```

    #> $F
    #> [1] 24.52549
    #> 
    #> $df1
    #> [1] 1
    #> 
    #> $df2
    #> [1] 295
    #> 
    #> $pvalue
    #> [1] 1.236496e-06
    #> 
    #> $partial_r2
    #> [1] 0.07675598
    #> 
    #> $n
    #> [1] 300
    #> 
    #> $weak
    #> [1] FALSE

## Mediator instruments

### Scan for cis-eQTLs

[`call_cis_eqtls()`](https://seantbresnahan.com/iconic/reference/call_cis_eqtls.md)
residualises expression and genotypes on covariates (FWL) and runs a
per-gene cis scan with Benjamini-Hochberg correction within gene.

``` r

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
```

    #>     gene       snp       beta        se         t            p          q  pass
    #> 1  Gene1   1:5:A:G  0.7381040 0.1822037  4.050982 0.0001195414 0.02390828  TRUE
    #> 2 Gene10  1:16:A:G -0.3992653 0.1416802 -2.818075 0.0061192857 0.65270286 FALSE
    #> 3 Gene11   1:7:A:G  0.4951092 0.1624197  3.048331 0.0031411347 0.62822693 FALSE
    #> 4 Gene12  1:18:A:G -0.6368031 0.2365183 -2.692406 0.0086788359 0.48041135 FALSE
    #> 5 Gene13 1:175:A:G -0.6659835 0.2183781 -3.049682 0.0031285656 0.57001683 FALSE
    #> 6 Gene14 1:115:A:G  0.5094369 0.1871389  2.722240 0.0079954784 0.61274755 FALSE

### Build a genetically predicted mediator (GReX-style)

[`build_mediator_instruments()`](https://seantbresnahan.com/iconic/reference/build_mediator_instruments.md)
fits an elastic-net cis-prediction model per gene (via `glmnet`), keeps
genes whose cross-validated out-of-fold prediction passes quality gates,
and returns the predicted mediator matrix `Gm` in the genes x samples
orientation expected by
[`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md).

``` r

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
```

    #>    gene n_cis_snps n_nonzero        cv_r2     cor_oof         cv_p lambda_min
    #> 1 Gene1         80        14  0.603994566  0.79234657 7.091717e-34  0.3699326
    #> 2 Gene2         80         0 -0.016596372 -0.14033085 9.566256e-01  0.7635586
    #> 3 Gene3         80        10 -0.004780502  0.07266243 1.884416e-01  0.4009905
    #> 4 Gene4         80         0 -0.014790698 -0.16869659 9.804737e-01  0.5855206
    #>    pass
    #> 1  TRUE
    #> 2 FALSE
    #> 3 FALSE
    #> 4 FALSE

Only Gene1, which carries the two planted cis-eQTLs, passes the
cross-validation gates.

## Negative controls from high-dimensional panels

### Methylation beta to M-values

Methylation beta values are heteroscedastic at the extremes;
[`beta_to_m()`](https://seantbresnahan.com/iconic/reference/beta_to_m.md)
applies the logit transform with clipping.

``` r

b <- matrix(runif(200, 0.01, 0.99), nrow = 20,
            dimnames = list(paste0("cg", 1:20), paste0("S", 1:10)))
m <- beta_to_m(b)
range(m)
```

    #> [1] -5.243843  6.478614

### Residualise on covariates

Before extracting negative-control factors, remove technical covariates
(for example batch) from the panel.

``` r

x <- matrix(rnorm(100 * 40), nrow = 100,
            dimnames = list(paste0("f", 1:100), paste0("S", 1:40)))
batch <- factor(rep(c("A", "B"), each = 20))
cv <- model.matrix(~ batch)[, -1, drop = FALSE]
xr <- residualize_matrix(x, cv)
cor(as.numeric(xr[1, ]), as.numeric(cv))
```

    #> [1] -1.581456e-17

### Extract negative-control factors with PCA

[`build_w_pcs()`](https://seantbresnahan.com/iconic/reference/build_w_pcs.md)
computes the top principal components of the (residualised) panel — the
`W` matrix used by the proximal and COCA estimators.

``` r

x <- matrix(rnorm(500 * 60), nrow = 500,
            dimnames = list(paste0("f", 1:500), paste0("S", 1:60)))
w <- build_w_pcs(x, n_pcs = 5)
dim(w$W)
```

    #> [1]  5 60

``` r

w$variance_explained
```

    #> [1] 3.120368 2.967639 2.823363 2.772791 2.757102

### Apply FUSION-style TWAS weights

Pre-computed transcriptomic weights (FUSION `.wgt.RDat` files, or an
in-memory list with the same structure) can be applied to a dosage panel
to build a genetically predicted expression panel, which can itself
serve as a negative-control or mediator-instrument panel.

``` r

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
```

    #> [1]  2 50

## Importing a SummarizedExperiment

For Bioconductor-native workflows,
[`as_iconic_data()`](https://seantbresnahan.com/iconic/reference/as_iconic_data.md)
dispatches on `SummarizedExperiment`: the primary assay becomes the
outcome panel `Y`, an optional second assay becomes the mediator panel
`M`, and `colData` columns supply the exposure, instruments, negative
controls, covariates, and survival endpoints.

``` r

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
```

    #> <iconic_data> 60 samples, 20 outcome features
    #>  Available: G (exposure instrument), W (negative controls), W1/W2 (path-specific NCs) 
    #>  Covariates: age 
    #>  Mode: total effect

From here the standard workflow applies:
[`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md),
[`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md),
and the sensitivity machinery described in
[`vignette("iconic-walkthrough")`](https://seantbresnahan.com/iconic/articles/iconic-walkthrough.md).

## Session information

``` r

sessionInfo()
```

    #> R version 4.6.1 (2026-06-24)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.4 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    #> 
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] doRNG_1.8.6.3    rngtools_1.5.2   foreach_1.5.2    iconic_0.99.0   
    #> [5] BiocStyle_2.40.0
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] SummarizedExperiment_1.42.0 gtable_0.3.6               
    #>  [3] shape_1.4.6.1               xfun_0.60                  
    #>  [5] bslib_0.12.0                ggplot2_4.0.3              
    #>  [7] Biobase_2.72.0              lattice_0.22-9             
    #>  [9] bigassertr_0.2.0            ps_1.9.3                   
    #> [11] vctrs_0.7.3                 tools_4.6.1                
    #> [13] generics_0.1.4              stats4_4.6.1               
    #> [15] parallel_4.6.1              tibble_3.3.1               
    #> [17] pkgconfig_2.0.3             Matrix_1.7-5               
    #> [19] data.table_1.18.4           RColorBrewer_1.1-3         
    #> [21] bigstatsr_1.6.2             S7_0.2.2                   
    #> [23] desc_1.4.3                  S4Vectors_0.50.1           
    #> [25] lifecycle_1.0.5             compiler_4.6.1             
    #> [27] farver_2.1.2                textshaping_1.0.5          
    #> [29] bigparallelr_0.3.2          Seqinfo_1.2.0              
    #> [31] codetools_0.2-20            htmltools_0.5.9            
    #> [33] sass_0.4.10                 yaml_2.3.12                
    #> [35] glmnet_5.0                  pillar_1.11.1              
    #> [37] pkgdown_2.2.1               jquerylib_0.1.4            
    #> [39] cachem_1.1.0                DelayedArray_0.38.2        
    #> [41] iterators_1.0.14            abind_1.4-8                
    #> [43] tidyselect_1.2.1            digest_0.6.39              
    #> [45] dplyr_1.2.1                 bookdown_0.47              
    #> [47] splines_4.6.1               cowplot_1.2.0              
    #> [49] fastmap_1.2.0               grid_4.6.1                 
    #> [51] cli_3.6.6                   SparseArray_1.12.2         
    #> [53] magrittr_2.0.5              S4Arrays_1.12.0            
    #> [55] survival_3.8-6              withr_3.0.3                
    #> [57] scales_1.4.0                XVector_0.52.0             
    #> [59] rmarkdown_2.31              bigsparser_0.7.3           
    #> [61] matrixStats_1.5.0           rmio_0.4.0                 
    #> [63] bit_4.6.0                   otel_0.2.0                 
    #> [65] ragg_1.5.2                  evaluate_1.0.5             
    #> [67] ff_4.5.3                    knitr_1.51                 
    #> [69] GenomicRanges_1.64.0        IRanges_2.46.0             
    #> [71] doParallel_1.0.17           irlba_2.3.7                
    #> [73] rlang_1.3.0                 Rcpp_1.1.2                 
    #> [75] glue_1.8.1                  BiocManager_1.30.27        
    #> [77] BiocGenerics_0.58.1         jsonlite_2.0.0             
    #> [79] R6_2.6.1                    bigsnpr_1.12.21            
    #> [81] systemfonts_1.3.2           MatrixGenerics_1.24.0      
    #> [83] fs_2.1.0                    flock_0.7
