# Convert methylation beta values to M-values

Applies the logit transform \\\log_2(\beta / (1 - \beta))\\ after
clipping betas away from 0 and 1. M-values are approximately
homoscedastic and are the recommended scale for linear modeling and PCA
of methylation data.

## Usage

``` r
beta_to_m(beta, clip = 1e-04, drop_nonfinite = TRUE)
```

## Arguments

- beta:

  Numeric matrix of beta values (probes x samples), in \[0, 1\].

- clip:

  Numeric: betas are clipped to \[`clip`, `1 - clip`\] before the
  transform. Default 1e-4.

- drop_nonfinite:

  Logical: drop probes (rows) with any missing/non-finite values after
  the transform. Default `TRUE`.

## Value

A numeric matrix of M-values, probes x samples.

## Examples

``` r
b <- matrix(runif(200, 0.01, 0.99), nrow = 20,
            dimnames = list(paste0("cg", 1:20), paste0("S", 1:10)))
m <- beta_to_m(b)
range(m)
#> [1] -6.061127  6.545274
```
