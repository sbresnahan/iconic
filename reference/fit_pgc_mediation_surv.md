# PGC survival mediation estimator: single-panel bridge with Cox / RMST

Bridge-function-adjusted natural direct and indirect effects with a
survival outcome stage, using a single negative-control panel W for both
the X-\>M and M-\>Y confounding paths. This is the survival analogue of
[`fit_pgc_mediation`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation.md)
(continuous outcome):

1.  OLS: residualise X on G -\> X_resid.

2.  OLS: bridge X_resid on the FULL W matrix -\> W_hat (proxy for U).

3.  OLS: `M ~ X + W_hat + covars` -\> alpha_M.

4.  Cox / RMST: `Surv(t,e) ~ X + M + W_hat + covars` -\> NDE (coef on
    X), beta_M (coef on M).

`NIE = alpha_M * beta_M`. Unlike
[`fit_pgc_mediation2_surv`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation2_surv.md),
which uses path-specific W1/W2 bridges, this estimator uses a single
combined W panel and is appropriate when separate conf_XM / conf_MY
confounders are not assumed.

## Usage

``` r
fit_pgc_mediation_surv(
  time,
  event,
  X,
  M,
  g,
  W,
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

- M:

  Numeric mediator vector (length n).

- g:

  Numeric instrument for X (length n).

- W:

  Numeric NC matrix (n x q) or vector (length n).

- covars:

  Optional data frame of covariates (n rows).

- min_f:

  Minimum partial F for G. Default 10.

- effect_scale:

  Character: `"loghr"` or `"rmst"`.

- tau:

  RMST horizon (rmst only).

## Value

Named list (same fields as
[`fit_unadj_mediation_surv`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation_surv.md)).
Returns all-NA if the first-stage partial F for G is below `min_f`.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
fit_pgc_mediation_surv(dat$surv_time, dat$surv_event, dat$X, dat$M,
dat$G[, 1], dat$W)
#> $NDE
#> [1] 0.2197512
#> 
#> $NDE_se
#> [1] 1.245767
#> 
#> $NDE_p
#> [1] 0.7995041
#> 
#> $NIE
#> [1] 0.1678713
#> 
#> $NIE_se
#> [1] 0.697834
#> 
#> $NIE_p
#> [1] 0.8098957
#> 
#> $alpha_M
#> [1] 0.4982158
#> 
#> $alpha_se
#> [1] 0.005077079
#> 
#> $beta_M
#> [1] 0.336945
#> 
#> $beta_M_se
#> [1] 1.400662
#> 
```
