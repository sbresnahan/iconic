# PGC binary estimator: proxy G-component correction with logistic / LPM

Three-step bridge-function estimator with a binary outcome stage:

1.  Residualise X on G -\> X_resid (OLS).

2.  Bridge X_resid on the FULL W matrix -\> W_hat (OLS).

3.  Regress the binary outcome on X + W_hat via logistic regression
    (log-OR) or a linear probability model (risk difference).

## Usage

``` r
fit_pgc_bin(y, X, g, W, covars = NULL, effect_scale = c("logor", "riskdiff"))
```

## Arguments

- y:

  Numeric 0/1 outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- g:

  Numeric instrument vector (length n).

- W:

  Numeric NC matrix (n x q) or vector (length n).

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
fit_pgc_bin(dat$y_bin, dat$X, dat$G[, 1], dat$W)
#> $beta
#> [1] 0.2603184
#> 
#> $se
#> [1] 0.1982229
#> 
#> $pvalue
#> [1] 0.189095
#> 
```
