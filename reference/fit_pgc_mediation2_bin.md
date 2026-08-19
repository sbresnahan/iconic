# PGC2 / PGC2Gm binary mediation estimator: path-specific bridges with logistic / LPM

Two-stage proximal mediation with path-specific negative controls and a
binary outcome stage. Stages 1-2 (X bridge on W1, M bridge on W2, X_hat
/ M_hat construction) are identical to
[`fit_pgc_mediation2`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation2.md)
and remain OLS. Only stage 3 switches to logistic regression (log-OR) or
a linear probability model (risk difference):
`y ~ X_hat + M_hat + W_hat_X + W_hat_M + covars`.

## Usage

``` r
fit_pgc_mediation2_bin(
  y,
  X,
  M,
  g,
  W1,
  W2,
  gm = NULL,
  covars = NULL,
  min_f = 10,
  effect_scale = c("logor", "riskdiff")
)
```

## Arguments

- y:

  Numeric 0/1 outcome vector (length n).

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

  Character: `"logor"` or `"riskdiff"`.

## Value

Named list (same fields as
[`fit_unadj_mediation_bin`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation_bin.md)).
Returns all-NA if the first-stage partial F for G1 is below `min_f`.

## Details

When `gm = NULL` (PGC2), stage 2 uses pure NC identification. When `gm`
is supplied (PGC2Gm), the mediator instrument helps isolate conf_MY
before bridging W2.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 500, outcome_type = "binary",
mo_confounding = 0.8, rho_G2 = 0.3,
lambda_XM = c(1, 0), lambda_MY = c(0, 1), seed = 1)
fit_pgc_mediation2_bin(dat$y_bin, dat$X, dat$M,
dat$G[, 1], dat$W1, dat$W2, gm = dat$Gm)
#> $NDE
#> [1] -0.5286415
#> 
#> $NDE_se
#> [1] 0.2240008
#> 
#> $NDE_p
#> [1] 0.018275
#> 
#> $NIE
#> [1] 0.733721
#> 
#> $NIE_se
#> [1] 0.1223357
#> 
#> $NIE_p
#> [1] 2.002509e-09
#> 
#> $alpha_M
#> [1] 0.5002583
#> 
#> $alpha_se
#> [1] 0.0116492
#> 
#> $beta_M
#> [1] 1.466684
#> 
#> $beta_M_se
#> [1] 0.2421483
#> 
```
