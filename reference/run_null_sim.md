# Run null simulations to estimate Type I error rates

Run null simulations to estimate Type I error rates

## Usage

``` r
run_null_sim(
  n_iter = 200,
  n_samples = 500,
  n_features = 20,
  conf_str = 0.8,
  w_signal = 0.7,
  feat_cor = 0,
  base_seed = 300,
  n_cores = 1,
  alpha = 0.05
)
```

## Arguments

- n_iter:

  Number of replicates. Default 200.

- n_samples:

  Observations per replicate. Default 500.

- n_features:

  Features per replicate. Default 20.

- conf_str:

  Confounding strength delta. Default 0.80.

- w_signal:

  Proxy quality omega. Default 0.70.

- feat_cor:

  Within-module feature correlation. Default 0.

- base_seed:

  Seed offset. Default 300.

- n_cores:

  Parallel workers. Default 1.

- alpha:

  Significance threshold. Default 0.05.

## Value

A list with rates (data frame) and raw (full results).

## Examples

``` r
null <- run_null_sim(n_iter = 2, n_samples = 100, n_features = 5)
null$rates
#>        method type1_error     flag
#> UNADJ   UNADJ         1.0 INFLATED
#> DIRECT DIRECT         1.0 INFLATED
#> COCA     COCA         0.2 INFLATED
#> IV2SLS IV2SLS         0.1       OK
#> PGC       PGC         0.9 INFLATED
```
