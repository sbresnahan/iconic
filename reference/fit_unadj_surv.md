# UNADJ survival estimator: unadjusted Cox / RMST regression

Regresses the survival outcome on the exposure X (plus covariates) with
no instrument or negative-control adjustment. Bias reference.

## Usage

``` r
fit_unadj_surv(
  time,
  event,
  X,
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

- covars:

  Optional data frame of covariates (n rows).

- effect_scale:

  Character: `"loghr"` (Cox) or `"rmst"` (pseudo-observation OLS).
  Default `"loghr"`.

- tau:

  RMST restriction time horizon. Default `NULL` (90th percentile of
  follow-up). Used only when `effect_scale = "rmst"`.

## Value

Named list: `beta`, `se`, `pvalue`.

## Details

When `effect_scale = "loghr"` (default), fits a Cox proportional-hazards
model and returns the log-hazard ratio for X. When
`effect_scale = "rmst"`, regresses leave-one-out RMST
pseudo-observations (Graw et al. 2009) on X via OLS, returning an effect
on the restricted-mean-survival-time (time) scale — a collapsible,
linear alternative to the non-collapsible hazard ratio.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
fit_unadj_surv(dat$surv_time, dat$surv_event, dat$X)
#> $beta
#> [1] 0.5547407
#> 
#> $se
#> [1] 1.741489
#> 
#> $pvalue
#> [1] 5.581735e-08
#> 
fit_unadj_surv(dat$surv_time, dat$surv_event, dat$X, effect_scale = "rmst")
#> $beta
#> [1] -1.64363
#> 
#> $se
#> [1] 0.3177295
#> 
#> $pvalue
#> [1] 5.622797e-07
#> 
```
