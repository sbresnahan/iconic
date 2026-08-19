# DIRECT binary estimator: logistic / LPM with instrument and NC covariates

Regresses the binary outcome on X plus the genetic instrument G, the
negative-control panel W, and covariates. Naive adjustment; does not
correct for unmeasured confounding via a ratio or IV approach.

## Usage

``` r
fit_direct_bin(
  y,
  X,
  g,
  w,
  covars = NULL,
  effect_scale = c("logor", "riskdiff")
)
```

## Arguments

- y:

  Numeric 0/1 outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- g:

  Numeric instrument vector (length n).

- w:

  Numeric NC vector (length n) or matrix (n x q).

- covars:

  Optional data frame of covariates (n rows).

- effect_scale:

  Character: `"logor"` or `"riskdiff"`.

## Value

Named list: `beta`, `se`, `pvalue`.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "binary", seed = 1)
fit_direct_bin(dat$y_bin, dat$X, dat$G[, 1], dat$W[, 1])
#> $beta
#> [1] 0.5398956
#> 
#> $se
#> [1] 0.2533408
#> 
#> $pvalue
#> [1] 0.03308055
#> 
```
