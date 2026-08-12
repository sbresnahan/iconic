# Run repeated simulations for a single parameter configuration

Run repeated simulations for a single parameter configuration

## Usage

``` r
run_simulation(
  n_iter = 100,
  n_samples = 500,
  n_features = 20,
  beta_X = 0.1,
  alpha_M = 0.5,
  beta_M = 0.3,
  conf_str = 0.8,
  w_signal = 0.7,
  feat_cor = 0,
  base_seed = 100,
  n_cores = 1
)
```

## Arguments

- n_iter:

  Number of simulation replicates. Default 100.

- n_samples:

  Observations per replicate. Default 500.

- n_features:

  Number of outcome and negative-control features. Default 20.

- beta_X:

  Direct effect of X on Y. Default 0.10.

- alpha_M:

  Effect of X on mediator. Default 0.50.

- beta_M:

  Effect of mediator on Y. Default 0.30.

- conf_str:

  Confounding strength delta. Default 0.80.

- w_signal:

  Proxy quality omega. Default 0.70.

- feat_cor:

  Within-module correlation for block-diagonal co-expression modules in
  Y and W. 0 = independent features. Default 0.

- base_seed:

  Starting seed; replicate i uses base_seed + i. Default 100.

- n_cores:

  Number of parallel workers. Default 1.

## Value

A list with raw, summary, iter_bias, true_total, params.

## Examples

``` r
res <- run_simulation(n_iter = 3, n_samples = 100, beta_X = 0.1,
                      conf_str = 0.8)
res$summary
#>   method       mean     median         sd        bias   abs_bias       rmse
#> 1  UNADJ 0.57509719 0.56731766 0.07836763  0.32509719 0.32509719 0.33425635
#> 2 DIRECT 0.43484712 0.43395739 0.08633872  0.18484712 0.18484712 0.20371203
#> 3   COCA 0.07020157 0.06864553 0.07711674 -0.17979843 0.17979843 0.19538513
#> 4 IV2SLS 0.21230851 0.21629677 0.05063672 -0.03769149 0.03769149 0.06278528
#> 5    PGC 0.32168376 0.31735060 0.04050204  0.07168376 0.07168376 0.08216835
#>       power  n
#> 1 1.0000000 60
#> 2 1.0000000 60
#> 3 0.2500000 60
#> 4 0.9166667 60
#> 5 1.0000000 60
```
