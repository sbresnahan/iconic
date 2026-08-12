# Test mediator-instrument independence of negative controls (W *\|* Gm \| C)

For each control feature, computes the partial correlation with the
mediator-specific genetic instrument Gm after residualising both on
observed covariates C, and reports the p-value. Controls significantly
associated with Gm after FDR correction may carry eQTL / allele-specific
effects that violate the mediator-instrument-independence assumption
(A2'). This is the mediator-instrument analogue of
[`nc_independence_check()`](https://seantbresnahan.com/iconic/reference/nc_independence_check.md):
just as the exposure instrument G must be independent of the negative
controls, so must the mediator instrument Gm.

## Usage

``` r
nc_independence_check_gm(dat, fdr_level = 0.1, n_cores = 1)
```

## Arguments

- dat:

  Dataset list from
  [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
  or
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md),
  containing `Gm`, `W`, and `synthetic_data`. If `Gm` is absent (i.e. no
  mediator instrument was generated), the function returns `NULL` with a
  message.

- fdr_level:

  Target FDR for BH correction. Default 0.10.

- n_cores:

  Number of parallel workers. Default 1 (sequential). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix and a PSOCK cluster on Windows.

## Value

A data frame with one row per control feature: `feature`, `partial_r`,
`p_value`, `fdr`, `significant`, `verdict`. Returns `NULL` if `dat$Gm`
is not present.

## Details

For a mediator cis-eQTL instrument, this screen tests whether the eQTL
instrument is associated with the negative-control panel – a violation
would indicate shared genomic structure between the instrument SNPs and
the control features.

## Examples

``` r
dat <- run_single_iteration(n_features = 10, phi = 0.8, seed = 1)
nc_independence_check_gm(dat)
#> NC independence (Gm): 10 tasks (sequential)
#>  NC independence (Gm): 10% (1/10) [0s]
#>  NC independence (Gm): 20% (2/10) [0s]
#>  NC independence (Gm): 30% (3/10) [0s]
#>  NC independence (Gm): 40% (4/10) [0s]
#>  NC independence (Gm): 50% (5/10) [0s]
#>  NC independence (Gm): 60% (6/10) [0s]
#>  NC independence (Gm): 70% (7/10) [0s]
#>  NC independence (Gm): 80% (8/10) [0s]
#>  NC independence (Gm): 90% (9/10) [0s]
#>  NC independence (Gm): 100% (10/10) [0s]
#>    feature    partial_r   p_value       fdr significant verdict
#> 1        1 -0.005942095 0.8946642 0.9568965       FALSE   valid
#> 2        2 -0.002425623 0.9568965 0.9568965       FALSE   valid
#> 3        3  0.014979713 0.7385298 0.9568965       FALSE   valid
#> 4        4 -0.003474061 0.9382972 0.9568965       FALSE   valid
#> 5        5 -0.019308759 0.6669895 0.9568965       FALSE   valid
#> 6        6  0.003262888 0.9420410 0.9568965       FALSE   valid
#> 7        7  0.011658359 0.7950306 0.9568965       FALSE   valid
#> 8        8 -0.009218868 0.8372426 0.9568965       FALSE   valid
#> 9        9  0.003499733 0.9378421 0.9568965       FALSE   valid
#> 10      10  0.002933582 0.9478818 0.9568965       FALSE   valid
```
