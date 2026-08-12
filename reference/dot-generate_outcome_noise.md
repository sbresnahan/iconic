# Generate outcome noise, optionally correlated across features (internal)

When `noise_cor` is a valid `p x p` correlation matrix, draws `n`
samples from `MVN(0, 0.2^2 * noise_cor)`. Otherwise uses independent
`rnorm(n, 0, 0.2)` (backward compatible).

## Usage

``` r
.generate_outcome_noise(n, p, noise_cor)
```

## Arguments

- n:

  Number of samples.

- p:

  Number of features.

- noise_cor:

  A `p x p` correlation matrix, or NULL.

## Value

An `n x p` numeric matrix of outcome noise.
