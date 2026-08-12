# PGC2 / PGC2Gm survival mediation estimator: path-specific bridges with Cox / RMST

Two-stage proximal mediation with path-specific negative controls and a
survival outcome stage. Stages 1-2 (X bridge on W1, M bridge on W2,
X_hat / M_hat construction) are identical to
[`fit_pgc_mediation2`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation2.md)
and remain OLS. Only stage 3 switches to Cox (log-HR) or RMST
pseudo-observation OLS:
`Surv(t,e) ~ X_hat + M_hat + W_hat_X + W_hat_M + covars`.

## Usage

``` r
fit_pgc_mediation2_surv(
  time,
  event,
  X,
  M,
  g,
  W1,
  W2,
  gm = NULL,
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

- W1:

  Numeric NC matrix (n x q) or vector for the X-\>M path.

- W2:

  Numeric NC matrix (n x q) or vector for the M-\>Y path.

- gm:

  Optional numeric mediator instrument (length n).

- covars:

  Optional data frame of covariates (n rows).

- min_f:

  Minimum partial F for G1. Default 10.

- effect_scale:

  Character: `"loghr"` or `"rmst"`.

- tau:

  RMST horizon (rmst only).

## Value

Named list (same fields as
[`fit_unadj_mediation_surv`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation_surv.md)).
Returns all-NA if the first-stage partial F for G1 is below `min_f`.

## Details

When `gm = NULL` (PGC2), stage 2 uses pure NC identification. When `gm`
is supplied (PGC2Gm), the mediator instrument helps isolate conf_MY
before bridging W2.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 500, outcome_type = "survival",
mo_confounding = 0.8, rho_G2 = 0.3,
lambda_XM = c(1, 0), lambda_MY = c(0, 1), seed = 1)
fit_pgc_mediation2_surv(dat$surv_time, dat$surv_event, dat$X, dat$M,
dat$G[, 1], dat$W1, dat$W2, gm = dat$Gm)
#> $NDE
#> [1] -0.6529615
#> 
#> $NDE_se
#> [1] 0.520502
#> 
#> $NDE_p
#> [1] 2.288635e-05
#> 
#> $NIE
#> [1] 0.8102306
#> 
#> $NIE_se
#> [1] 2.526972
#> 
#> $NIE_p
#> [1] 0.7484885
#> 
#> $alpha_M
#> [1] 0.5002583
#> 
#> $alpha_se
#> [1] 0.0116492
#> 
#> $beta_M
#> [1] 1.619625
#> 
#> $beta_M_se
#> [1] 5.051194
#> 
```
