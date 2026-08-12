# IV2SLS mediation estimator: instrumented exposure in both stages

Uses the genetic instrument G to purge U1 from X, then estimates the
mediator and outcome regressions with the cleaned exposure.

## Usage

``` r
fit_iv2sls_mediation(y, X, M, g, w, covars = NULL, min_f = 10)
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

- min_f:

  Minimum acceptable partial F-statistic for the excluded instrument.
  Default 10.

## Value

Named list: `NDE`, `NDE_se`, `NDE_p`, `NIE`, `NIE_se`, `NIE_p`.

## Details

Strategy:

1.  `X ~ G + W` -\> X_hat (purge U1 from X).

2.  `M ~ X_hat` -\> alpha_M (clean effect of X on M).

3.  `Y ~ X_hat + M + W` -\> NDE = beta_X, beta_M (OLS).

The IV cleans X of U1 confounding, but M remains endogenous via U1 -\>
M. With a single instrument, natural effects are not fully identified
(Rudolph et al., 2024); NDE and NIE are approximations whose bias from
M-O confounding is the key finding the simulation demonstrates.

## References

Rudolph, K. E., et al. (2024). Natural direct and indirect effects with
an instrumental variable. *Biometrics*.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 300, mo_confounding = 0.8, seed = 1)
fit_iv2sls_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$W[, 1])
#> $NDE
#> [1] -0.1252387
#> 
#> $NDE_se
#> [1] 0.02585698
#> 
#> $NDE_p
#> [1] 2.059719e-06
#> 
#> $NIE
#> [1] 0.6381454
#> 
#> $NIE_se
#> [1] 0.03637006
#> 
#> $NIE_p
#> [1] 6.393409e-69
#> 
#> $alpha_M
#> [1] 0.7874435
#> 
#> $alpha_se
#> [1] 0.03647012
#> 
#> $beta_M
#> [1] 0.8104016
#> 
#> $beta_M_se
#> [1] 0.0269171
#> 
```
