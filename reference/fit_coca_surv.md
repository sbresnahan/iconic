# COCA survival estimator: NOT supported (returns NA)

COCA (Correlated Outcome Control Approach) regresses the negative
control W on the outcome Y (`W ~ y + X`) and recovers the causal effect
as a ratio \\-\hat\beta_X / \hat\beta_Y\\. This places the outcome on
the right-hand side of the regression, which is structurally impossible
when the outcome is a
[`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html) object
(time-to-event). COCA is therefore unsupported for survival outcomes and
always returns `list(beta=NA, se=NA, pvalue=NA)`.

## Usage

``` r
fit_coca_surv(time, event, X, w, covars = NULL, ...)
```

## Arguments

- time:

  Numeric follow-up time vector (length n).

- event:

  Numeric 0/1 event indicator (length n).

- X:

  Numeric exposure vector (length n).

- w:

  Numeric NC vector (length n).

- covars:

  Optional data frame of covariates (n rows).

- ...:

  Ignored (accepted for signature compatibility).

## Value

`list(beta = NA, se = NA, pvalue = NA)` with an informative `"reason"`
attribute.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "survival", seed = 1)
fit_coca_surv(dat$surv_time, dat$surv_event, dat$X, dat$W[, 1])
#> $beta
#> [1] NA
#> 
#> $se
#> [1] NA
#> 
#> $pvalue
#> [1] NA
#> 
#> attr(,"reason")
#> [1] "COCA regresses W on Y (W ~ y + X), placing the outcome on the RHS, which is impossible with a Surv object. COCA is unsupported for survival outcomes."
# $beta [1] NA (COCA unsupported for survival outcomes)
```
