# PGC estimator: Proxy G-Component Correction (matrix bridge)

A three-step bridge-function estimator:

1.  Residualises X on G to isolate the U-driven component X_resid.

2.  Regresses X_resid on the FULL W matrix to construct W_hat, a proxy
    for unmeasured confounding. This step requires `ncol(W) >= k` (the
    proximal completeness condition): if W has fewer valid columns than
    confounders, the bridge cannot span the confounder subspace and the
    estimator is under-identified.

3.  Fits Y ~ X + W_hat to absorb confounding bias.

## Usage

``` r
fit_pgc(y, X, g, W, covars = NULL)
```

## Arguments

- y:

  Numeric outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- g:

  Numeric instrument vector (length n).

- W:

  Numeric negative-control matrix (n x q) or vector (length n). If a
  matrix, the bridge uses all q columns. Pass only validity-screened
  columns for the completeness condition to be meaningful.

- covars:

  Optional data frame of additional covariates (n rows).

## Value

Named list: `beta`, `se`, `pvalue`.

## Details

Unlike the scalar version
([`fit_pgc_scalar`](https://seantbresnahan.com/iconic/reference/fit_pgc_scalar.md)),
which collapses W to `rowMeans(W)` and is algebraically equivalent to
IV/2SLS, the matrix bridge preserves the dimensional structure of W and
is the estimator for which the completeness condition is binding.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, seed = 1)
fit_pgc(dat$Y[, 1], dat$X, dat$G[, 1], dat$W)
#> $beta
#> [1] 0.3440353
#> 
#> $se
#> [1] 0.02441775
#> 
#> $pvalue
#> [1] 1.219016e-31
#> 
```
