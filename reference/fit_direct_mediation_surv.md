# DIRECT survival mediation estimator: Cox / RMST with G and W covariates

Adjusts for the instrument G and negative-control W in both the mediator
(OLS) and outcome (Cox / RMST) stages. Naive adjustment.

## Usage

``` r
fit_direct_mediation_surv(
  time,
  event,
  X,
  M,
  g,
  w,
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

- M:

  Numeric mediator vector (length n).

- g:

  Numeric instrument vector (length n).

- w:

  Numeric NC vector (length n) or matrix (n x q).

- covars:

  Optional data frame of covariates (n rows).

- effect_scale:

  Character: `"loghr"` or `"rmst"`.

- tau:

  RMST horizon (rmst only).

## Value

Named list (same fields as
[`fit_unadj_mediation_surv`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation_surv.md)).

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "survival",
mo_confounding = 0.8, seed = 1)
fit_direct_mediation_surv(dat$surv_time, dat$surv_event, dat$X, dat$M,
dat$G[, 1], dat$W[, 1])
#> $NDE
#> [1] -0.4158999
#> 
#> $NDE_se
#> [1] 0.6597463
#> 
#> $NDE_p
#> [1] 0.3497413
#> 
#> $NIE
#> [1] 0.9305827
#> 
#> $NIE_se
#> [1] 2.619077
#> 
#> $NIE_p
#> [1] 0.7223579
#> 
#> $alpha_M
#> [1] 0.7229852
#> 
#> $alpha_se
#> [1] 0.02021802
#> 
#> $beta_M
#> [1] 1.287139
#> 
#> $beta_M_se
#> [1] 3.622409
#> 
```
