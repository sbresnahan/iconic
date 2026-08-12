# DIRECT survival estimator: Cox / RMST with instrument and NC covariates

Regresses the survival outcome on X plus the genetic instrument G, the
negative-control panel W, and covariates. Naive adjustment; does not
correct for unmeasured confounding via a ratio or IV approach.

## Usage

``` r
fit_direct_surv(
  time,
  event,
  X,
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

Named list: `beta`, `se`, `pvalue`.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
fit_direct_surv(dat$surv_time, dat$surv_event, dat$X, dat$G[, 1], dat$W[, 1])
#> $beta
#> [1] 0.4221789
#> 
#> $se
#> [1] 1.525281
#> 
#> $pvalue
#> [1] 0.005884582
#> 
```
