# Extract per-seed bias vector (internal)

Extract per-seed bias vector (internal)

## Usage

``` r
compute_iter_bias(combined, true_total)
```

## Arguments

- combined:

  Data frame from run_methods().

- true_total:

  Scalar true total causal effect.

## Value

Data frame with columns: iter, method, bias.
