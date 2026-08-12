# PGC estimator: Proxy G-Component Correction (scalar bridge)

The original ICONIC PGC implementation, which summarises the
negative-control panel as a scalar (`rowMeans(W)`) before bridging. This
version is numerically stable and works in small samples, but the scalar
bridge produces \\\hat W\\ proportional to the G-residualised exposure
by construction, making the estimator algebraically equivalent to
IV/2SLS when the instrument is valid. The proximal completeness
condition (`dim(W_valid) >= k`) is NOT binding for this estimator.

## Usage

``` r
fit_pgc_scalar(y, X, g, w, covars = NULL)
```

## Arguments

- y:

  Numeric outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- g:

  Numeric instrument vector (length n).

- w:

  Numeric negative-control vector (length n). Pass `rowMeans(W_matrix)`
  for stability.

- covars:

  Optional data frame of additional covariates (n rows).

## Value

Named list: `beta`, `se`, `pvalue`.

## Details

Use [`fit_pgc`](https://seantbresnahan.com/iconic/reference/fit_pgc.md)
(matrix bridge) when the completeness condition is of interest. Use this
function as a stable fallback when `n` is small relative to the number
of valid controls, or as an IV-equivalent benchmark.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, seed = 1)
fit_pgc_scalar(dat$Y[, 1], dat$X, dat$G[, 1], rowMeans(dat$W))
#> $beta
#> [1] 0.2989159
#> 
#> $se
#> [1] 0.03344448
#> 
#> $pvalue
#> [1] 2.869071e-16
#> 
```
