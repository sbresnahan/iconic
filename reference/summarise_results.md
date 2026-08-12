# Summarise simulation results across features (internal)

Summarise simulation results across features (internal)

## Usage

``` r
summarise_results(combined, true_total)
```

## Arguments

- combined:

  Data frame from
  [`run_methods()`](https://seantbresnahan.com/iconic/reference/run_methods.md).

- true_total:

  Scalar true total causal effect.

## Value

Data frame with one row per method: mean, median, sd, bias, abs_bias,
rmse, power, n.
