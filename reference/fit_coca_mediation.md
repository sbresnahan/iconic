# COCA mediation estimator: negative-control calibration of both stages

Uses the negative-control outcome W to calibrate both the mediator and
outcome regressions, extending
[`fit_coca`](https://seantbresnahan.com/iconic/reference/fit_coca.md) to
the mediation setting.

## Usage

``` r
fit_coca_mediation(y, X, M, w, covars = NULL, ratio_cap = 10)
```

## Arguments

- y:

  Numeric primary outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- M:

  Numeric mediator vector (length n).

- w:

  Numeric negative-control outcome vector (length n). Recommended: pass
  `rowMeans(W_matrix)` for stability.

- covars:

  Optional data frame of additional covariates (n rows).

- ratio_cap:

  Maximum absolute value of any ratio estimate before flagging as
  unstable. Default 10.

## Value

Named list: `NDE`, `NDE_se`, `NDE_p`, `NIE`, `NIE_se`, `NIE_p`.

## Details

Stage 1: `W ~ M + X` -\> calibrated alpha_M = -beta_X / beta_M. Stage 2:
`W ~ Y + X + M` -\> calibrated NDE = -beta_X / beta_Y, calibrated beta_M
= -beta_M / beta_Y. NIE = alpha_M \* beta_M (both calibrated).

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, mo_confounding = 0.8, seed = 1)
fit_coca_mediation(dat$Y[, 1], dat$X, dat$M, rowMeans(dat$W))
#> $NDE
#>         X 
#> -5.482667 
#> 
#> $NDE_se
#> [1] 3.891432
#> 
#> $NDE_p
#>         X 
#> 0.1588626 
#> 
#> $NIE
#>        X 
#> 5.734214 
#> 
#> $NIE_se
#> [1] 3.882897
#> 
#> $NIE_p
#>         X 
#> 0.1397326 
#> 
#> $alpha_M
#>         X 
#> 0.4962959 
#> 
#> $alpha_se
#> [1] 0.01958304
#> 
#> $beta_M
#>        M 
#> 11.55402 
#> 
#> $beta_M_se
#> [1] 7.81046
#> 
```
