# PGC binary mediation estimator: single-panel bridge with logistic / LPM

Bridge-function-adjusted natural direct and indirect effects with a
binary outcome stage, using a single negative-control panel W for both
the X-\>M and M-\>Y confounding paths. This is the binary analogue of
[`fit_pgc_mediation`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation.md)
(continuous outcome):

1.  OLS: residualise X on G -\> X_resid.

2.  OLS: bridge X_resid on the FULL W matrix -\> W_hat (proxy for U).

3.  OLS: `M ~ X + W_hat + covars` -\> alpha_M.

4.  Logistic / LPM: `y ~ X + M + W_hat + covars` -\> NDE (coef on X),
    beta_M (coef on M).

`NIE = alpha_M * beta_M`. Unlike
[`fit_pgc_mediation2_bin`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation2_bin.md),
which uses path-specific W1/W2 bridges, this estimator uses a single
combined W panel and is appropriate when separate conf_XM / conf_MY
confounders are not assumed.

## Usage

``` r
fit_pgc_mediation_bin(
  y,
  X,
  M,
  g,
  W,
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

- W:

  Numeric NC matrix (n x q) or vector (length n).

- covars:

  Optional data frame of covariates (n rows).

- min_f:

  Minimum partial F for G. Default 10.

- effect_scale:

  Character: `"logor"` or `"riskdiff"`.

## Value

Named list (same fields as
[`fit_unadj_mediation_bin`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation_bin.md)).

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 500, outcome_type = "binary",
mo_confounding = 0.8, seed = 1)
fit_pgc_mediation_bin(dat$y_bin, dat$X, dat$M, dat$G[, 1], dat$W)
#> $NDE
#> [1] -0.5448331
#> 
#> $NDE_se
#> [1] 0.3615117
#> 
#> $NDE_p
#> [1] 0.1317858
#> 
#> $NIE
#> [1] 0.985109
#> 
#> $NIE_se
#> [1] 0.3428178
#> 
#> $NIE_p
#> [1] 0.004058685
#> 
#> $alpha_M
#> [1] 0.5912409
#> 
#> $alpha_se
#> [1] 0.01002849
#> 
#> $beta_M
#> [1] 1.666172
#> 
#> $beta_M_se
#> [1] 0.5791386
#> 
```
