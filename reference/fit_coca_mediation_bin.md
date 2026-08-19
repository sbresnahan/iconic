# COCA binary mediation estimator: NOT supported (returns NA)

COCA mediation regresses the negative control W on the outcome Y in both
stages (`W ~ M + X`, `W ~ y + X + M`) and recovers the effect as a ratio
of regression coefficients. That identification argument assumes a
linear structural outcome model; with a binary (nonlinear, logistic)
outcome the ratio recovers neither the causal log-odds ratio nor the
risk difference. COCA mediation is therefore unsupported for binary
outcomes and always returns NA.

## Usage

``` r
fit_coca_mediation_bin(y, X, M, w, covars = NULL, ...)
```

## Arguments

- y:

  Numeric 0/1 outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- M:

  Numeric mediator vector (length n).

- w:

  Numeric NC vector (length n).

- covars:

  Optional data frame of covariates (n rows).

- ...:

  Ignored (signature compatibility).

## Value

All-NA list with an informative `"reason"` attribute.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "binary", seed = 1)
fit_coca_mediation_bin(dat$y_bin, dat$X, dat$M, dat$W[, 1])
#> $NDE
#> [1] NA
#> 
#> $NDE_se
#> [1] NA
#> 
#> $NDE_p
#> [1] NA
#> 
#> $NIE
#> [1] NA
#> 
#> $NIE_se
#> [1] NA
#> 
#> $NIE_p
#> [1] NA
#> 
#> $alpha_M
#> [1] NA
#> 
#> $alpha_se
#> [1] NA
#> 
#> $beta_M
#> [1] NA
#> 
#> $beta_M_se
#> [1] NA
#> 
#> attr(,"reason")
#> [1] "COCA mediation regresses W on Y (W ~ y + X + M) and recovers the effect as a ratio, an identification argument that assumes a linear structural outcome model. With a binary (nonlinear, logistic) outcome the ratio recovers neither the log-OR nor the risk difference. COCA mediation is unsupported for binary outcomes."
# all NA (COCA mediation unsupported for binary outcomes)
```
