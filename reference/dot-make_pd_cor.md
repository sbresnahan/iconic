# Near-PD correction for a correlation matrix (internal)

Clips eigenvalues to a small positive floor and renormalizes to unit
diagonal. Same approach as .residual_correlation() in load_data.R.

## Usage

``` r
.make_pd_cor(cor_mat)
```

## Arguments

- cor_mat:

  A symmetric matrix (intended to be a correlation matrix).

## Value

A valid correlation matrix (symmetric, positive-definite, unit
diagonal).
