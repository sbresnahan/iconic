# Inverse CDF for a fitted marginal (internal)

Returns the quantile function for a parametric or empirical marginal.

## Usage

``` r
.marginal_quantile(marginal, p)
```

## Arguments

- marginal:

  A marginal specification list from .fit_parametric_marginal or an
  empirical CDF object.

- p:

  Vector of probabilities in (0, 1).

## Value

Numeric vector of quantiles.
