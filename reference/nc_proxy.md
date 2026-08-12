# Direct-proxy negative-control model

Each control is a coverage-weighted mixture of the captured confounders
plus noise. With one confounder and `coverage = omega` this reduces to
the classic `w_signal` proxy in
[`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md).
Setting `captured` to a strict subset of the confounders models controls
that miss part of the confounder subspace (invalidating the
negative-control assumption).

## Usage

``` r
nc_proxy(U, covariates, params)
```

## Arguments

- U:

  `n x k` confounder matrix.

- covariates:

  Covariate data frame (unused; kept for the NC contract).

- params:

  List with `n_features`, and optionally `coverage` (scalar in `[0, 1]`,
  default 0.7), `captured` (integer confounder indices the controls see,
  default all), `noise_sd` (default 0.3), `MMCon` (loading multiplier,
  default 1), `mode` (default `"shared"`), and `noise_cor` (a `p x p`
  correlation matrix for correlated noise, or `NULL` for independent
  noise).

## Value

`n x n_features` matrix of negative controls.

## Details

When `noise_cor` is supplied (a `p x p` correlation matrix), the noise
component is drawn from a multivariate normal with that correlation
structure, so the negative controls retain realistic cross-feature
correlations conditional on the confounder. When `noise_cor` is `NULL`
(default), the noise is independent across features.

## Modes

- `"shared"`:

  All columns carry the same `rowMeans(U[, captured])` signal. This is
  the original behaviour (backward-compatible) and is numerically
  stable, but the W matrix has only one effective dimension regardless
  of `n_features`, so the proximal completeness condition is never
  binding for the matrix-bridge PGC.

- `"distinct"`:

  Column `f` captures confounder
  `captured[((f - 1) \%\% length(captured)) + 1]`. Different columns
  therefore carry signals from different confounders, giving the W
  matrix genuine dimensional structure. This is the mode to use when
  benchmarking the completeness cliff: the matrix-bridge PGC is
  identified only when `n_features >= k`.

## Examples

``` r
U <- matrix(rnorm(100), 100, 1)
W <- nc_proxy(U, covariates = NULL,
  params = list(n_features = 10, coverage = 0.7))
dim(W)
#> [1] 100  10
```
