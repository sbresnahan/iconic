# Compute the residual correlation matrix of a features x samples panel (internal)

Residualizes each feature (row) on the sample-level exposure and
covariates, then computes the p x p correlation matrix of the residuals.
Returns NULL if the matrix has fewer than 2 features or the correlation
is degenerate.

## Usage

``` r
.residual_correlation(
  mat,
  exposure,
  covars,
  residualize_on = c("XC", "XCW"),
  W_mat = NULL
)
```

## Arguments

- mat:

  A `features x samples` numeric matrix.

- exposure:

  A length-n_samples numeric vector (per-sample exposure).

- covars:

  A data frame of encoded covariates (n_samples rows).

- residualize_on:

  Character: `"XC"` (default, legacy) residualizes on exposure +
  covariates only; `"XCW"` additionally residualizes on the
  negative-control bridge proxy Ŵ, obtained by regressing `W_mat` on
  exposure + covariates. This partials out the U-signature that W
  captures, reducing the double-counting of U's cross-feature signature
  when the residual matrix is later injected as correlated noise
  alongside synthetic U loadings. Tradeoff: can over-partial if W is a
  weak proxy — recommended only when the completeness-capture test (§1)
  reports "strong".

- W_mat:

  Optional `features x samples` NC matrix used to construct the bridge
  proxy Ŵ when `residualize_on = "XCW"`. Ignored otherwise.

## Value

A p x p correlation matrix, or NULL.
