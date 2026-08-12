# DIRECT mediation estimator: OLS with instrument and NC as covariates

Adjusts for the genetic instrument G and negative-control W in both the
mediator and outcome regressions. Like
[`fit_direct`](https://seantbresnahan.com/iconic/reference/fit_direct.md),
this is a naive adjustment that does not correct for unmeasured
confounding via a ratio or IV approach.

## Usage

``` r
fit_direct_mediation(y, X, M, g, w, covars = NULL)
```

## Arguments

- y:

  Numeric outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- M:

  Numeric mediator vector (length n).

- g:

  Numeric instrument vector (length n).

- w:

  Numeric negative-control vector (length n).

- covars:

  Optional data frame of additional covariates (n rows).

## Value

Named list: `NDE`, `NDE_se`, `NDE_p`, `NIE`, `NIE_se`, `NIE_p`.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, mo_confounding = 0.8, seed = 1)
fit_direct_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$W[, 1])
#> $NDE
#> [1] -0.3486436
#> 
#> $NDE_se
#> [1] 0.06224249
#> 
#> $NDE_p
#> [1] 7.165418e-08
#> 
#> $NIE
#> [1] 0.8687889
#> 
#> $NIE_se
#> [1] 0.06284506
#> 
#> $NIE_p
#> [1] 1.818792e-43
#> 
#> $alpha_M
#> [1] 0.7229852
#> 
#> $alpha_se
#> [1] 0.02021802
#> 
#> $beta_M
#> [1] 1.201669
#> 
#> $beta_M_se
#> [1] 0.08016613
#> 
```
