# Summarise mediation simulation results across features (internal)

Iterates over the methods present in `combined` (rather than a hardcoded
list), so it handles both the 5-method (no Gm) and 6-method (with Gm)
cases.

## Usage

``` r
summarise_mediation_results(combined, true_NDE, true_NIE)
```

## Arguments

- combined:

  Data frame from
  [`run_mediation_methods()`](https://seantbresnahan.com/iconic/reference/run_mediation_methods.md).

- true_NDE:

  Scalar true natural direct effect.

- true_NIE:

  Scalar true natural indirect effect.

## Value

Data frame with one row per method: NDE/NIE mean, bias, sd, rmse, Type I
error rates, and counts.
