# DIRECT estimator: OLS with instrument and negative-control as covariates

Regresses Y on X plus the genetic instrument G, the negative-control W,
and any additional covariates. This is a "naive" adjustment that uses
whatever observables are available but does NOT correct for unmeasured
confounding via a ratio or IV approach.

## Usage

``` r
fit_direct(y, X, g, w, covars = NULL)
```

## Arguments

- y:

  Numeric outcome vector (length n).

- X:

  Numeric exposure vector (length n), assumed pre-scaled.

- g:

  Numeric instrument vector (length n).

- w:

  Numeric negative-control vector (length n) or matrix (n x q). When a
  matrix, all q columns are included as separate covariates.

- covars:

  Optional data frame of additional covariates (n rows).

## Value

Named list: `beta`, `se`, `pvalue`.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, seed = 1)
fit_direct(dat$Y[, 1], dat$X, dat$G[, 1], dat$W[, 1])
#> $beta
#> [1] 0.4543358
#> 
#> $se
#> [1] 0.02923705
#> 
#> $pvalue
#> [1] 5.121819e-36
#> 
```
