# UNADJ mediation estimator: naive Baron-Kenny style

Stage 1: `M ~ X` (estimate alpha_M). Stage 2: `Y ~ X + M` (estimate NDE
= beta_X, beta_M). NIE = alpha_M \* beta_M.

## Usage

``` r
fit_unadj_mediation(y, X, M, covars = NULL)
```

## Arguments

- y:

  Numeric outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- M:

  Numeric mediator vector (length n).

- covars:

  Optional data frame of additional covariates (n rows).

## Value

Named list: `NDE`, `NDE_se`, `NDE_p`, `NIE`, `NIE_se`, `NIE_p`.

## Details

Does not adjust for unmeasured confounding. Provided as a bias reference
floor, analogous to UNADJ in the total-effect setting.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, mo_confounding = 0.8, seed = 1)
fit_unadj_mediation(dat$Y[, 1], dat$X, dat$M)
#> $NDE
#> [1] -0.3296719
#> 
#> $NDE_se
#> [1] 0.04116539
#> 
#> $NDE_p
#> [1] 9.912469e-14
#> 
#> $NIE
#> [1] 0.896143
#> 
#> $NIE_se
#> [1] 0.04478842
#> 
#> $NIE_p
#> [1] 4.656993e-89
#> 
#> $alpha_M
#> [1] 0.7652063
#> 
#> $alpha_se
#> [1] 0.01901884
#> 
#> $beta_M
#> [1] 1.171113
#> 
#> $beta_M_se
#> [1] 0.05078045
#> 
```
