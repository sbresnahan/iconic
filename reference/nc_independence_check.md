# Test instrument-independence of negative controls (W *\|* G \| C)

For each control feature, computes the partial correlation with the
genetic instrument G after residualising both on observed covariates C,
and reports the p-value. Controls significantly associated with G after
FDR correction may carry meQTL / allele-specific effects that violate
the instrument-independence assumption (A2).

## Usage

``` r
nc_independence_check(dat, fdr_level = 0.1, n_cores = 1)
```

## Arguments

- dat:

  Dataset list from
  [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
  or
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md).

- fdr_level:

  Target FDR for BH correction. Default 0.10.

- n_cores:

  Number of parallel workers. Default 1 (sequential). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix and a PSOCK cluster on Windows.

## Value

A data frame with one row per control feature: `feature`, `partial_r`,
`p_value`, `fdr`, `significant`, `verdict`.

## Examples

``` r
dat <- run_single_iteration(n_features = 10, seed = 1)
nc_independence_check(dat)
#> NC independence (G): 10 tasks (sequential)
#>  NC independence (G): 10% (1/10) [0s]
#>  NC independence (G): 20% (2/10) [0s]
#>  NC independence (G): 30% (3/10) [0s]
#>  NC independence (G): 40% (4/10) [0s]
#>  NC independence (G): 50% (5/10) [0s]
#>  NC independence (G): 60% (6/10) [0s]
#>  NC independence (G): 70% (7/10) [0s]
#>  NC independence (G): 80% (8/10) [0s]
#>  NC independence (G): 90% (9/10) [0s]
#>  NC independence (G): 100% (10/10) [0s]
#>    feature    partial_r   p_value       fdr significant verdict
#> 1        1 -0.012333113 0.7834542 0.9349553       FALSE   valid
#> 2        2 -0.024678277 0.5823385 0.9349553       FALSE   valid
#> 3        3 -0.027609123 0.5383489 0.9349553       FALSE   valid
#> 4        4 -0.017619234 0.6945948 0.9349553       FALSE   valid
#> 5        5  0.003662625 0.9349553 0.9349553       FALSE   valid
#> 6        6 -0.013731240 0.7596215 0.9349553       FALSE   valid
#> 7        7 -0.009543787 0.8315904 0.9349553       FALSE   valid
#> 8        8 -0.007105246 0.8742022 0.9349553       FALSE   valid
#> 9        9 -0.021990376 0.6240949 0.9349553       FALSE   valid
#> 10      10 -0.014680248 0.7435717 0.9349553       FALSE   valid
```
