# IV2SLS estimator: Two-Stage Least Squares with genetic instrument

Uses the genetic instrument G to instrument for the exposure X,
controlling for the negative-control W and any additional covariates.
Requires AER.

## Usage

``` r
fit_iv2sls(y, X, g, w, covars = NULL, min_f = 10)
```

## Arguments

- y:

  Numeric outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- g:

  Numeric instrument vector (length n).

- w:

  Numeric negative-control vector (length n).

- covars:

  Optional data frame of additional covariates (n rows).

- min_f:

  Minimum acceptable partial F-statistic for the excluded instrument.
  Default 10.

## Value

Named list: `beta`, `se`, `pvalue`.

## Details

A weak-instrument check is applied using the partial F-statistic for the
excluded instrument G (testing G conditional on W and covariates in the
first stage), following Stock & Yogo (2005). If the partial F is below
`min_f`, the function returns `NA`.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, seed = 1)
fit_iv2sls(dat$Y[, 1], dat$X, dat$G[, 1], dat$W[, 1])
#> $beta
#> [1] 0.2796918
#> 
#> $se
#> [1] 0.03449663
#> 
#> $pvalue
#> [1] 5.381187e-14
#> 
```
