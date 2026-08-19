# IV2SLS binary mediation estimator: single-instrument 2SPS with logistic / LPM

Two-stage predictor substitution with a single instrument (G for X, no
mediator instrument):

1.  OLS: `X ~ g + W + covars` -\> X_hat (purge U1 from X).

2.  OLS: `M ~ X_hat + covars` -\> alpha_M.

3.  Logistic / LPM: `y ~ X_hat + M + W + covars` -\> NDE (coef on
    X_hat), beta_M (coef on M).

`NIE = alpha_M * beta_M`. Weak-instrument gate (partial F for G) applies
to the OLS first stage. Unlike
[`fit_iv2sls_mediation2_bin`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation2_bin.md),
this estimator does not instrument the mediator and is NOT
point-identified under M-O confounding — it is included for parity with
the continuous
[`fit_iv2sls_mediation`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation.md).

## Usage

``` r
fit_iv2sls_mediation_bin(
  y,
  X,
  M,
  g,
  w,
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

- w:

  Numeric NC vector (length n) or matrix (n x q).

- covars:

  Optional data frame of covariates (n rows).

- min_f:

  Minimum partial F for the excluded instrument. Default 10.

- effect_scale:

  Character: `"logor"` or `"riskdiff"`.

## Value

Named list (same fields as
[`fit_unadj_mediation_bin`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation_bin.md)).
Returns all-NA if the first-stage partial F for G is below `min_f`.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 500, outcome_type = "binary",
mo_confounding = 0.8, seed = 1)
fit_iv2sls_mediation_bin(dat$y_bin, dat$X, dat$M, dat$G[, 1], dat$W[, 1])
#> $NDE
#> [1] -0.2332658
#> 
#> $NDE_se
#> [1] 0.2096974
#> 
#> $NDE_p
#> [1] 0.2659696
#> 
#> $NIE
#> [1] 0.7960189
#> 
#> $NIE_se
#> [1] 0.1735209
#> 
#> $NIE_p
#> [1] 4.48686e-06
#> 
#> $alpha_M
#> [1] 0.820922
#> 
#> $alpha_se
#> [1] 0.02989651
#> 
#> $beta_M
#> [1] 0.9696644
#> 
#> $beta_M_se
#> [1] 0.2084024
#> 
```
