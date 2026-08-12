# Generate correlated noise for negative-control models (internal)

When `noise_cor` is a valid `p x p` correlation matrix, draws `n`
samples from `MVN(0, noise_cor)` and scales by `sd0`. Returns an
`n x (2*p)` matrix: columns `1:p` are the "uncaptured" noise component
and columns `(p+1):(2*p)` are the idiosyncratic noise component (both
with the same correlation structure). When `noise_cor` is NULL, uses
independent `rnorm` (backward compatible).

## Usage

``` r
.generate_nc_noise(n, p, noise_cor, sd0)
```

## Arguments

- n:

  Number of samples.

- p:

  Number of features.

- noise_cor:

  A `p x p` correlation matrix, or NULL.

- sd0:

  Idiosyncratic noise SD.

## Value

An `n x (2*p)` numeric matrix.
