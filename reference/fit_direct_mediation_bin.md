# DIRECT binary mediation estimator: logistic / LPM with G and W covariates

Adjusts for the instrument G and negative-control W in both the mediator
(OLS) and outcome (logistic / LPM) stages. Naive adjustment.

## Usage

``` r
fit_direct_mediation_bin(
  y,
  X,
  M,
  g,
  w,
  covars = NULL,
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

  Numeric instrument vector (length n).

- w:

  Numeric NC vector (length n) or matrix (n x q).

- covars:

  Optional data frame of covariates (n rows).

- effect_scale:

  Character: `"logor"` or `"riskdiff"`.

## Value

Named list (same fields as
[`fit_unadj_mediation_bin`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation_bin.md)).

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "binary",
mo_confounding = 0.8, seed = 1)
fit_direct_mediation_bin(dat$y_bin, dat$X, dat$M,
dat$G[, 1], dat$W[, 1])
#> $NDE
#> [1] -1.516507
#> 
#> $NDE_se
#> [1] 0.7137921
#> 
#> $NDE_p
#> [1] 0.03362188
#> 
#> $NIE
#> [1] 2.016358
#> 
#> $NIE_se
#> [1] 0.6771035
#> 
#> $NIE_p
#> [1] 0.002902143
#> 
#> $alpha_M
#> [1] 0.7229852
#> 
#> $alpha_se
#> [1] 0.02021802
#> 
#> $beta_M
#> [1] 2.788934
#> 
#> $beta_M_se
#> [1] 0.9332854
#> 
```
