# IV2SLS binary estimator: two-stage predictor substitution with logistic / LPM

Two-stage predictor substitution (2SPS): the first stage regresses X on
the instrument G (plus W and covariates) via OLS, producing fitted
\\\hat X\\; the second stage regresses the binary outcome on \\\hat X\\
(plus W and covariates) via logistic regression (log-OR) or a linear
probability model (risk difference).

## Usage

``` r
fit_iv2sls_bin(
  y,
  X,
  g,
  w,
  covars = NULL,
  min_f = 10,
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

- min_f:

  Minimum partial F for the excluded instrument. Default 10.

- effect_scale:

  Character: `"logor"` or `"riskdiff"`.

## Value

Named list: `beta`, `se`, `pvalue`.

## Details

A weak-instrument check (partial F for the excluded instrument G, Stock
& Yogo 2005) is applied to the OLS first stage. If the partial F is
below `min_f`, the function returns NA.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "binary", seed = 1)
fit_iv2sls_bin(dat$y_bin, dat$X, dat$G[, 1], dat$W[, 1])
#> $beta
#> [1] -0.01277931
#> 
#> $se
#> [1] 0.2630275
#> 
#> $pvalue
#> [1] 0.9612497
#> 
```
