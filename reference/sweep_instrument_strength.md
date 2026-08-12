# Sweep instrument strength

Runs a custom DGP that varies the G-\>X coefficient (pi_GX) and records
IV2SLS performance binned by first-stage partial F. Uses
[`fit_iv2sls()`](https://seantbresnahan.com/iconic/reference/fit_iv2sls.md)
with `min_f = 0` (no weak-IV guard) so bias is visible across all
instrument strengths.

## Usage

``` r
sweep_instrument_strength(
  pi_GX_grid = c(0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
  n_iter = 50,
  n = 200,
  k = 1,
  conf_str = 0.8,
  tau = 0.25,
  coverage = 0.7,
  n_cores = 1
)
```

## Arguments

- pi_GX_grid:

  Numeric vector of G-\>X coefficients to sweep.

- n_iter:

  Number of replications per grid point.

- n:

  Sample size.

- k:

  Number of confounders.

- conf_str:

  Confounding strength.

- tau:

  True total effect.

- coverage:

  Negative-control coverage.

- n_cores:

  Number of parallel workers for the replicate loops. Default 1
  (sequential). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix and a PSOCK cluster on Windows.

## Value

A data frame with columns `pi_GX`, `arm`, `iter`, `partial_F`, `beta`,
`se`, `pvalue`, `rejected`.

## Examples

``` r
res <- sweep_instrument_strength(pi_GX_grid = c(0.05, 0.2), n_iter = 2)
#> Running instrument-strength sweep...
#>  pi_GX=0.05: 4 tasks (sequential)
#>   pi_GX=0.05: 25% (1/4) [0s]
#>   pi_GX=0.05: 50% (2/4) [0s]
#>   pi_GX=0.05: 75% (3/4) [0s]
#>   pi_GX=0.05: 100% (4/4) [0s]
#>  pi_GX = 0.05 done (mean F = 1.9)
#>  pi_GX=0.2: 4 tasks (sequential)
#>   pi_GX=0.2: 25% (1/4) [0s]
#>   pi_GX=0.2: 50% (2/4) [0s]
#>   pi_GX=0.2: 75% (3/4) [0s]
#>   pi_GX=0.2: 100% (4/4) [0s]
#>  pi_GX = 0.20 done (mean F = 37.4)
head(res)
#>   pi_GX  arm iter  partial_F       beta        se    pvalue rejected
#> 1  0.05  alt    1  4.4765442 0.35939447 0.3463920 0.3007575    FALSE
#> 2  0.05 null    1  2.3100322 0.60615488 0.4923835 0.2197672    FALSE
#> 3  0.05  alt    2  0.2283631 1.46302545 2.5320641 0.5640593    FALSE
#> 4  0.05 null    2  0.4112861 1.71448760 2.6166955 0.5130963    FALSE
#> 5  0.20  alt    1 16.5591634 0.32253791 0.1957043 0.1009286    FALSE
#> 6  0.20 null    1 48.5586431 0.01714582 0.1021277 0.8668451    FALSE
```
