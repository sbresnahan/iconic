# Convert a linear predictor to a binary outcome

Given an n-vector (or n x p matrix) of linear predictors `eta` and a
target prevalence, returns a list with the 0/1 outcome vector `y_bin`.
When `eta` is a matrix, the first column is used (binary outcomes are
scalar — `n_features = 1`).

## Usage

``` r
.linpred_to_binary(eta, prev = 0.5)
```

## Arguments

- eta:

  numeric vector or matrix of linear predictors.

- prev:

  target prevalence (marginal event probability); the intercept is
  solved so that approximately this fraction of subjects are cases.
  Default 0.5.

## Value

list(y_bin, p, intercept, prev)
