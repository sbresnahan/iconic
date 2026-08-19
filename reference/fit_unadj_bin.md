# UNADJ binary estimator: unadjusted logistic / linear-probability regression

Regresses the binary outcome on the exposure X (plus covariates) with no
instrument or negative-control adjustment. Bias reference.

## Usage

``` r
fit_unadj_bin(y, X, covars = NULL, effect_scale = c("logor", "riskdiff"))
```

## Arguments

- y:

  Numeric 0/1 outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- covars:

  Optional data frame of covariates (n rows).

- effect_scale:

  Character: `"logor"` (logistic) or `"riskdiff"` (linear probability
  model). Default `"logor"`.

## Value

Named list: `beta`, `se`, `pvalue`.

## Details

When `effect_scale = "logor"` (default), fits a logistic regression and
returns the conditional log-odds ratio for X. When
`effect_scale = "riskdiff"`, fits a linear probability model (OLS on the
0/1 outcome), returning a risk difference — a collapsible, linear
alternative to the non-collapsible odds ratio.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, outcome_type = "binary", seed = 1)
fit_unadj_bin(dat$y_bin, dat$X)
#> $beta
#> [1] 0.4066461
#> 
#> $se
#> [1] 0.1527228
#> 
#> $pvalue
#> [1] 0.007752997
#> 
fit_unadj_bin(dat$y_bin, dat$X, effect_scale = "riskdiff")
#> $beta
#> [1] 0.09560984
#> 
#> $se
#> [1] 0.03461115
#> 
#> $pvalue
#> [1] 0.006278234
#> 
```
