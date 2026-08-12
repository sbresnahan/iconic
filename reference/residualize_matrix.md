# Residualize a feature matrix on covariates

Regresses each feature (row) of `x` on the covariate design and returns
the residuals, using a precomputed projection for memory- efficient
chunked processing of large (e.g. 450k-probe) matrices. This is the
workhorse for removing cell composition, batch, and technical factors
before extracting negative-control PCs.

## Usage

``` r
residualize_matrix(x, covariates, chunk_size = 5000)
```

## Arguments

- x:

  Numeric matrix, features x samples.

- covariates:

  data.frame or matrix with one row per sample (i.e.
  `nrow(covariates) == ncol(x)`). An intercept is added internally.
  Factors should already be expanded to dummies (see
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)).

- chunk_size:

  Integer: number of features processed per block. Default 5000.

## Value

A numeric matrix of residuals, features x samples, centered per feature.

## Examples

``` r
x <- matrix(rnorm(100 * 40), nrow = 100,
            dimnames = list(paste0("f", 1:100), paste0("S", 1:40)))
batch <- factor(rep(c("A", "B"), each = 20))
cv <- model.matrix(~ batch)[, -1, drop = FALSE]
xr <- residualize_matrix(x, cv)
cor(as.numeric(xr[1, ]), as.numeric(cv))
#> [1] 2.89022e-17
```
