# Build a negative-control panel from principal components

Computes the top `n_pcs` principal components of a features x samples
matrix (typically residualized with
[`residualize_matrix()`](https://seantbresnahan.com/iconic/reference/residualize_matrix.md))
and returns the sample scores as a W panel (PCs x samples) in the
orientation expected by
[`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md).
Uses irlba for fast truncated PCA when available, falling back to
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html).

## Usage

``` r
build_w_pcs(x, n_pcs = 20, scale_features = FALSE, prefix = "PC")
```

## Arguments

- x:

  Numeric matrix, features x samples (e.g. residualized M-values or
  predicted expression).

- n_pcs:

  Integer: number of PCs to return. Default 20.

- scale_features:

  Logical: scale features to unit variance before PCA. Default `FALSE`
  (appropriate for M-values); use `TRUE` for predicted-expression panels
  where feature scales differ.

- prefix:

  Character: row-name prefix for the PCs. Default `"PC"`.

## Value

A list with:

- W:

  numeric matrix, PCs x samples, suitable as the `W` argument of
  [`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md).

- variance_explained:

  numeric vector: percent of total variance explained by each returned
  PC.

## Examples

``` r
x <- matrix(rnorm(500 * 60), nrow = 500,
            dimnames = list(paste0("f", 1:500), paste0("S", 1:60)))
w <- build_w_pcs(x, n_pcs = 5)
dim(w$W)
#> [1]  5 60
w$variance_explained
#> [1] 3.027690 2.944088 2.717830 2.669895 2.621035
```
