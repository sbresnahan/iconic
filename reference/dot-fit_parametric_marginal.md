# Fit a parametric distribution to a data vector (internal)

Attempts to fit normal, log-normal, gamma, and beta distributions (where
applicable) by maximum likelihood, and returns the best fit by AIC along
with the KS test p-value.

## Usage

``` r
.fit_parametric_marginal(x)
```

## Arguments

- x:

  Numeric vector (the feature's observed values).

## Value

A list with `family`, `params`, `aic`, `ks_p`, or NULL if no parametric
fit succeeds.
