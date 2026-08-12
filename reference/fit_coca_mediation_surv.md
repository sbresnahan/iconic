# COCA survival mediation estimator: NOT supported (returns NA)

COCA mediation regresses the negative control W on the outcome Y in both
stages (`W ~ M + X`, `W ~ y + X + M`), placing the outcome on the
right-hand side. This is structurally impossible when the outcome is a
[`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object.
COCA mediation is therefore unsupported for survival outcomes and always
returns NA.

## Usage

``` r
fit_coca_mediation_surv(time, event, X, M, w, covars = NULL, ...)
```

## Arguments

- time:

  Numeric follow-up time vector (length n).

- event:

  Numeric 0/1 event indicator (length n).

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
dat <- generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
fit_coca_mediation_surv(dat$surv_time, dat$surv_event, dat$X, dat$M, dat$W[, 1])
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
#> [1] "COCA mediation regresses W on Y (W ~ y + X + M), placing the outcome on the RHS, which is impossible with a Surv object. COCA mediation is unsupported for survival outcomes."
# all NA (COCA mediation unsupported for survival outcomes)
```
