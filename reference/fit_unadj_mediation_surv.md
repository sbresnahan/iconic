# UNADJ survival mediation estimator: naive Baron-Kenny with Cox / RMST

Stage 1 (OLS): `M ~ X` -\> alpha_M. Stage 2 (Cox / RMST): survival
outcome ~ X + M -\> NDE (coef on X), beta_M (coef on M).
`NIE = alpha_M * beta_M`. No confounding adjustment; bias reference.

## Usage

``` r
fit_unadj_mediation_surv(
  time,
  event,
  X,
  M,
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

- covars:

  Optional data frame of covariates (n rows).

- effect_scale:

  Character: `"loghr"` or `"rmst"`.

- tau:

  RMST horizon (rmst only).

## Value

Named list: `NDE`, `NDE_se`, `NDE_p`, `NIE`, `NIE_se`, `NIE_p`,
`alpha_M`, `alpha_se`, `beta_M`, `beta_M_se`.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "survival",
mo_confounding = 0.8, seed = 1)
fit_unadj_mediation_surv(dat$surv_time, dat$surv_event, dat$X, dat$M)
#> $NDE
#> [1] -0.3825676
#> 
#> $NDE_se
#> [1] 0.6821078
#> 
#> $NDE_p
#> [1] 0.1793237
#> 
#> $NIE
#> [1] 1.002574
#> 
#> $NIE_se
#> [1] 2.836666
#> 
#> $NIE_p
#> [1] 0.7237632
#> 
#> $alpha_M
#> [1] 0.7652063
#> 
#> $alpha_se
#> [1] 0.01901884
#> 
#> $beta_M
#> [1] 1.310201
#> 
#> $beta_M_se
#> [1] 3.706918
#> 
```
