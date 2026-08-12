# PGC survival estimator: proxy G-component correction with Cox / RMST

Three-step bridge-function estimator with a survival outcome stage:

1.  Residualise X on G -\> X_resid (OLS).

2.  Bridge X_resid on the FULL W matrix -\> W_hat (OLS).

3.  Regress the survival outcome on X + W_hat via Cox (log-HR) or RMST
    pseudo-observation OLS.

## Usage

``` r
fit_pgc_surv(
  time,
  event,
  X,
  g,
  W,
  covars = NULL,
  effect_scale = c("loghr", "rmst"),
  tau = NULL
)
```

## Arguments

- time:

  Numeric follow-up time vector (length n).

- event:

  Numeric 0/1 event indicator (length n).

- X:

  Numeric exposure vector (length n).

- g:

  Numeric instrument vector (length n).

- W:

  Numeric NC matrix (n x q) or vector (length n).

- covars:

  Optional data frame of covariates (n rows).

- effect_scale:

  Character: `"loghr"` or `"rmst"`.

- tau:

  RMST horizon (rmst only).

## Value

Named list: `beta`, `se`, `pvalue`.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
fit_pgc_surv(dat$surv_time, dat$surv_event, dat$X, dat$G[, 1], dat$W)
#> $beta
#> [1] 0.3859335
#> 
#> $se
#> [1] 1.470987
#> 
#> $pvalue
#> [1] 0.002221541
#> 
```
