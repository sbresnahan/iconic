# IV2SLS survival estimator: two-stage predictor substitution with Cox / RMST

Two-stage predictor substitution (2SPS): the first stage regresses X on
the instrument G (plus W and covariates) via OLS, producing fitted
\\\hat X\\; the second stage regresses the survival outcome on \\\hat
X\\ (plus W and covariates) via Cox (log-HR) or RMST pseudo-observation
OLS.

## Usage

``` r
fit_iv2sls_surv(
  time,
  event,
  X,
  g,
  w,
  covars = NULL,
  min_f = 10,
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

- w:

  Numeric NC vector (length n) or matrix (n x q).

- covars:

  Optional data frame of covariates (n rows).

- min_f:

  Minimum partial F for the excluded instrument. Default 10.

- effect_scale:

  Character: `"loghr"` or `"rmst"`.

- tau:

  RMST horizon (rmst only).

## Value

Named list: `beta`, `se`, `pvalue`.

## Details

A weak-instrument check (partial F for the excluded instrument G, Stock
& Yogo 2005) is applied to the OLS first stage. If the partial F is
below `min_f`, the function returns NA.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
fit_iv2sls_surv(dat$surv_time, dat$surv_event, dat$X, dat$G[, 1], dat$W[, 1])
#> $beta
#> [1] 0.3740847
#> 
#> $se
#> [1] 1.45366
#> 
#> $pvalue
#> [1] 0.03551488
#> 
```
