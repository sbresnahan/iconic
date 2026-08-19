# UNADJ binary mediation estimator: naive Baron-Kenny with logistic / LPM

Stage 1 (OLS): `M ~ X` -\> alpha_M. Stage 2 (logistic / LPM): binary
outcome ~ X + M -\> NDE (coef on X), beta_M (coef on M).
`NIE = alpha_M * beta_M`. No confounding adjustment; bias reference.

## Usage

``` r
fit_unadj_mediation_bin(
  y,
  X,
  M,
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

- covars:

  Optional data frame of covariates (n rows).

- effect_scale:

  Character: `"logor"` or `"riskdiff"`.

## Value

Named list: `NDE`, `NDE_se`, `NDE_p`, `NIE`, `NIE_se`, `NIE_p`,
`alpha_M`, `alpha_se`, `beta_M`, `beta_M_se`.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "binary",
mo_confounding = 0.8, seed = 1)
fit_unadj_mediation_bin(dat$y_bin, dat$X, dat$M)
#> $NDE
#> [1] -1.325125
#> 
#> $NDE_se
#> [1] 0.4800089
#> 
#> $NDE_p
#> [1] 0.005769064
#> 
#> $NIE
#> [1] 1.853043
#> 
#> $NIE_se
#> [1] 0.4723713
#> 
#> $NIE_p
#> [1] 8.750686e-05
#> 
#> $alpha_M
#> [1] 0.7652063
#> 
#> $alpha_se
#> [1] 0.01901884
#> 
#> $beta_M
#> [1] 2.421625
#> 
#> $beta_M_se
#> [1] 0.6143711
#> 
```
