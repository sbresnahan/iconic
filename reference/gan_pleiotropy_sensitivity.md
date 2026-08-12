# Benchmark estimators across pleiotropy and confounding scenarios

For each cell in the grid (pleiotropy strength x confounding strength),
generates `n_iter` synthetic datasets with
[`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md),
runs every estimator, and summarises bias / RMSE / power. Two arms are
run per cell:

## Usage

``` r
gan_pleiotropy_sensitivity(
  trained_gan = NULL,
  pleio_grid = c(0, 0.05, 0.1),
  conf_grid = c(0.2, 0.5, 0.8),
  tau = 0.25,
  nc_model = "proxy",
  n_iter = 50,
  n_samples = 500,
  n_features = 10,
  coverage = 0.7,
  k = 1,
  base_seed = 900,
  n_cores = 1
)
```

## Arguments

- trained_gan:

  An `iconic_gan` (or `NULL` to use default texture).

- pleio_grid:

  Horizontal-pleiotropy strengths (direct G -\> Y coefficients) to
  sweep. Default `c(0, 0.05, 0.10)`.

- conf_grid:

  Confounding-strength values. Default `c(0.2, 0.5, 0.8)`.

- tau:

  True total effect for the alternative arm. Default 0.25.

- nc_model:

  Negative-control model (function or name). Default `"proxy"`.

- n_iter:

  Replicates per cell per arm. Default 50.

- n_samples:

  Samples per replicate. Default 500.

- n_features:

  Features per replicate. Default 10.

- coverage:

  Negative-control coverage. Default 0.7.

- k:

  Number of latent confounders. Default 1.

- base_seed:

  Base RNG seed. Default 900.

- n_cores:

  Parallel workers across replicates. Default 1.

## Value

A list with `summary` (one row per cell x arm x method, with `pleio`,
`conf_strength`, `arm`, `true_total`, and the columns from
[`summarise_results()`](https://seantbresnahan.com/iconic/reference/summarise_results.md))
and `grid`.

## Details

- **alternative** (`effect_size = tau`): the true total effect is `tau`,
  so `power` is the empirical power to detect it.

- **null** (`effect_size = 0`): the true total effect is 0, so `power`
  is the empirical Type I error rate.

The `pleio` parameter adds a direct `G -> Y` path of the requested
strength, violating the exclusion restriction. IV/2SLS is consistent
only when `pleio = 0`; any `pleio > 0` introduces bias that does not
shrink with sample size.

## Examples

``` r
sens <- gan_pleiotropy_sensitivity(NULL,
  pleio_grid = c(0, 0.10), conf_grid = 0.8,
  n_iter = 2, n_samples = 100, n_features = 5)
head(sens$summary)
#>   pleio conf_strength  arm true_total method       mean     median          sd
#> 1     0           0.8  alt       0.25  UNADJ  0.5600946  0.5828345 0.072883744
#> 2     0           0.8  alt       0.25 DIRECT  0.3634369  0.3843584 0.041171431
#> 3     0           0.8  alt       0.25   COCA -0.2458452 -0.1827060 0.146524016
#> 4     0           0.8  alt       0.25 IV2SLS  0.1550511  0.1586644 0.009041299
#> 5     0           0.8  alt       0.25    PGC  0.2625661  0.2631845 0.017282523
#> 6     0           0.8 null       0.00  UNADJ  0.2969582  0.2921351 0.047463252
#>          bias   abs_bias      rmse power  n
#> 1  0.31009459 0.31009459 0.3177098   1.0 10
#> 2  0.11343695 0.11343695 0.1199730   1.0 10
#> 3 -0.49584516 0.49584516 0.5149610   0.4 10
#> 4 -0.09494887 0.09494887 0.0953355   0.2 10
#> 5  0.01256611 0.01256611 0.0206573   1.0 10
#> 6  0.29695818 0.29695818 0.3003525   1.0 10
```
