# IV2SLS survival mediation estimator: single-instrument 2SPS with Cox / RMST

Two-stage predictor substitution with a single instrument (G for X, no
mediator instrument):

1.  OLS: `X ~ g + W + covars` -\> X_hat (purge U1 from X).

2.  OLS: `M ~ X_hat + covars` -\> alpha_M.

3.  Cox / RMST: `Surv(t,e) ~ X_hat + M + W + covars` -\> NDE (coef on
    X_hat), beta_M (coef on M).

`NIE = alpha_M * beta_M`. Weak-instrument gate (partial F for G) applies
to the OLS first stage. Unlike
[`fit_iv2sls_mediation2_surv`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation2_surv.md),
this estimator does not instrument the mediator and is NOT
point-identified under M-O confounding — it is included for parity with
the continuous
[`fit_iv2sls_mediation`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation.md).

## Usage

``` r
fit_iv2sls_mediation_surv(
  time,
  event,
  X,
  M,
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

- M:

  Numeric mediator vector (length n).

- g:

  Numeric instrument for X (length n).

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

Named list (same fields as
[`fit_unadj_mediation_surv`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation_surv.md)).

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
fit_iv2sls_mediation_surv(dat$surv_time, dat$surv_event, dat$X, dat$M,
dat$G[, 1], dat$W[, 1])
#> $NDE
#> [1] -0.03876021
#> 
#> $NDE_se
#> [1] 0.9619814
#> 
#> $NDE_p
#> [1] 0.8662954
#> 
#> $NIE
#> [1] 0.4245667
#> 
#> $NIE_se
#> [1] 1.169592
#> 
#> $NIE_p
#> [1] 0.7166018
#> 
#> $alpha_M
#> [1] 0.5015545
#> 
#> $alpha_se
#> [1] 0.02736802
#> 
#> $beta_M
#> [1] 0.8465018
#> 
#> $beta_M_se
#> [1] 2.331477
#> 
```
