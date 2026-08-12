# PGC mediation estimator: bridge-function-adjusted natural effects (scalar bridge)

The original ICONIC PGC mediation implementation, which summarises the
negative-control panel as a scalar (`rowMeans(W)`) before bridging. Like
[`fit_pgc_scalar`](https://seantbresnahan.com/iconic/reference/fit_pgc_scalar.md),
this version is algebraically equivalent to the IV2SLS mediation
estimator when the instrument is valid. Use
[`fit_pgc_mediation`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation.md)
(matrix bridge) when the completeness condition is of interest.

## Usage

``` r
fit_pgc_scalar_mediation(y, X, M, g, w, covars = NULL)
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

  Numeric negative-control vector (length n). Pass `rowMeans(W_matrix)`
  for stability.

- covars:

  Optional data frame of additional covariates (n rows).

## Value

Named list: `NDE`, `NDE_se`, `NDE_p`, `NIE`, `NIE_se`, `NIE_p`.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, mo_confounding = 0.8, seed = 1)
fit_pgc_scalar_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], rowMeans(dat$W))
#> $NDE
#> [1] -0.3304709
#> 
#> $NDE_se
#> [1] 0.04147903
#> 
#> $NDE_p
#> [1] 1.300326e-13
#> 
#> $NIE
#> [1] 0.6425387
#> 
#> $NIE_se
#> [1] 0.04652626
#> 
#> $NIE_p
#> [1] 2.211053e-43
#> 
#> $alpha_M
#> [1] 0.5454283
#> 
#> $alpha_se
#> [1] 0.02684247
#> 
#> $beta_M
#> [1] 1.178044
#> 
#> $beta_M_se
#> [1] 0.06257225
#> 
```
