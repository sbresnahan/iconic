# Build mediator instruments via per-gene elastic net

Trains one cross-validated elastic-net model per gene predicting
expression from cis SNPs, after residualizing both expression and
dosages on `covariates` (Frisch-Waugh-Lovell). Genes are kept when the
model is non-degenerate (at least one non-zero SNP weight), the signed
out-of-fold R2 exceeds `cv_r2_min`, and the one-sided out-of-fold
correlation p-value is below `cv_p_max`. The in-sample predicted
expression of passing genes is the composite mediator instrument `Gm`
consumed by
[`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md).

## Usage

``` r
build_mediator_instruments(
  expression,
  genotypes,
  gene_pos,
  snp_pos,
  covariates = NULL,
  cis_snps = NULL,
  alpha = 0.3,
  nfolds = 5,
  cis_dist = 1e+06,
  cv_r2_min = 0.005,
  cv_p_max = 0.1,
  n_cores = 1,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- expression:

  Numeric matrix, genes x samples (e.g. inverse-normal-transformed TPM).

- genotypes:

  Numeric dosage matrix, variants x samples (0-2).

- gene_pos:

  data.frame with columns `gene`, `chr`, and `tss` (or `start`).

- snp_pos:

  data.frame with columns `snp`, `chr`, `pos`.

- covariates:

  Optional data.frame or matrix (n rows) of covariates to partial out
  (ancestry PCs, technical factors, ...).

- cis_snps:

  Optional named list: gene name -\> character vector of candidate SNP
  IDs (e.g. significant GTEx cis-eQTL pairs, as from
  [`call_cis_eqtls()`](https://seantbresnahan.com/iconic/reference/call_cis_eqtls.md)
  or an external resource). When supplied, the elastic net for each gene
  sees only these SNPs (intersected with the cis window).

- alpha:

  Numeric: elastic-net mixing (0 = ridge, 1 = lasso). Default 0.3.

- nfolds:

  Integer: cross-validation folds. Default 5.

- cis_dist:

  Numeric: cis window half-width around the gene start, in bp. Default
  1e6.

- cv_r2_min:

  Numeric: minimum signed out-of-fold R2. Default 0.005.

- cv_p_max:

  Numeric: maximum one-sided out-of-fold correlation p-value (nominal
  screen only; the load-bearing gates are `cv_r2_min` and the
  non-degeneracy requirement). Default 0.10.

- n_cores:

  Integer: parallel workers via
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html).
  Default 1.

- seed:

  Optional integer for reproducible CV fold assignments.

- verbose:

  Logical: print progress. Default `TRUE`.

## Value

A list with:

- Gm:

  numeric matrix, genes x samples: in-sample predicted expression for
  genes passing the quality filter. Pass to
  [`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md)
  as `Gm` (it is transposed internally).

- qc:

  data.frame of per-gene quality metrics: `gene`, `n_cis_snps`,
  `n_nonzero`, `cv_r2`, `cor_oof`, `cv_p`, `lambda_min`, `pass`.

- weights:

  named list (all trained genes) of data.frames with `snp` and `weight`
  (non-zero coefficients at lambda.min).

## Details

Requires the glmnet package (listed under `Suggests`).

## Examples

``` r
if (requireNamespace("glmnet", quietly = TRUE)) {
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
}
#> Training elastic-net instruments for 4 genes...
#> Quality filter (OOF R2 > 0.005, one-sided p < 0.1, >= 1 non-zero SNP):
#>   Intercept-only fits (excluded): 2
#>   Passing: 1 / 4 genes
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
```
