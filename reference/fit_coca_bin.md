# COCA binary estimator: NOT supported (returns NA)

COCA (Correlated Outcome Control Approach) regresses the negative
control W on the outcome Y (`W ~ y + X`) and recovers the causal effect
as a ratio \\-\hat\beta_X / \hat\beta_Y\\. That identification argument
assumes a linear structural outcome model: the W-Y regression
coefficient must be proportional to the confounder-Outcome association
on the same scale as the causal effect. With a binary outcome the
structural model is nonlinear (logistic), so the linear COCA ratio
recovers neither the causal log-odds ratio nor the risk difference. COCA
is therefore unsupported for binary outcomes and always returns
`list(beta=NA, se=NA, pvalue=NA)`.

## Usage

``` r
fit_coca_bin(y, X, w, covars = NULL, ...)
```

## Arguments

- y:

  Numeric 0/1 outcome vector (length n).

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
dat <- generate_toy_data(n = 200, outcome_type = "binary", seed = 1)
fit_coca_bin(dat$y_bin, dat$X, dat$W[, 1])
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
#> [1] "COCA regresses W on Y (W ~ y + X) and recovers the effect as a ratio, an identification argument that assumes a linear structural outcome model. With a binary (nonlinear, logistic) outcome the ratio recovers neither the log-OR nor the risk difference. COCA is unsupported for binary outcomes."
# $beta [1] NA (COCA unsupported for binary outcomes)
```
