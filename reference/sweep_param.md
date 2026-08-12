# Sweep a single simulation parameter across a grid

Sweep a single simulation parameter across a grid

## Usage

``` r
sweep_param(
  param,
  param_grid,
  n_iter = 100,
  n_samples = 500,
  n_features = 20,
  beta_X = 0.1,
  alpha_M = 0.5,
  beta_M = 0.3,
  conf_str = 0.8,
  w_signal = 0.7,
  feat_cor = 0,
  base_seed = 0,
  n_cores = 1
)
```

## Arguments

- param:

  Parameter to vary: one of "beta_X", "conf_str", "w_signal", "alpha_M",
  "beta_M", "n_samples", "feat_cor".

- param_grid:

  Numeric vector of values to sweep.

- n_iter:

  Replicates per grid point. Default 100.

- n_samples:

  Observations per replicate. Default 500.

- n_features:

  Features per replicate. Default 20.

- beta_X:

  Baseline direct effect. Default 0.10.

- alpha_M:

  Baseline mediator path. Default 0.50.

- beta_M:

  Baseline mediator effect. Default 0.30.

- conf_str:

  Baseline confounding strength. Default 0.80.

- w_signal:

  Baseline proxy quality. Default 0.70.

- feat_cor:

  Baseline within-module feature correlation. Default 0.

- base_seed:

  Seed offset. Default 0.

- n_cores:

  Parallel workers. Default 1.

## Value

A list with summary (data frame) and iter_bias (data frame).

## Examples

``` r
res <- sweep_param("conf_str", c(0.2, 0.8), n_iter = 3, n_samples = 100)
res$summary
#>       param param_value true_total method       mean    median         sd
#> 1  conf_str         0.2       0.25  UNADJ 0.27414069 0.2738119 0.01981494
#> 2  conf_str         0.2       0.25 DIRECT 0.25208653 0.2517659 0.03398081
#> 3  conf_str         0.2       0.25   COCA 0.14966349 0.1541537 0.05256443
#> 4  conf_str         0.2       0.25 IV2SLS 0.24758739 0.2483571 0.03091597
#> 5  conf_str         0.2       0.25    PGC 0.24791252 0.2458816 0.02057153
#> 6  conf_str         0.8       0.25  UNADJ 0.60478040 0.6114569 0.08029133
#> 7  conf_str         0.8       0.25 DIRECT 0.46760502 0.4782414 0.08292496
#> 8  conf_str         0.8       0.25   COCA 0.09518619 0.1098339 0.08089260
#> 9  conf_str         0.8       0.25 IV2SLS 0.26672002 0.2715802 0.06235480
#> 10 conf_str         0.8       0.25    PGC 0.36326156 0.3723539 0.04141499
#>            bias    abs_bias       rmse     power  n
#> 1   0.024140695 0.024140695 0.03112653 1.0000000 60
#> 2   0.002086525 0.002086525 0.03376098 1.0000000 60
#> 3  -0.100336515 0.100336515 0.11306806 0.7500000 60
#> 4  -0.002412607 0.002412607 0.03075204 1.0000000 60
#> 5  -0.002087484 0.002087484 0.02050591 1.0000000 60
#> 6   0.354780404 0.354780404 0.36360471 1.0000000 60
#> 7   0.217605024 0.217605024 0.23262392 1.0000000 60
#> 8  -0.154813813 0.154813813 0.17436132 0.4166667 60
#> 9   0.016720025 0.016720025 0.06405371 0.9333333 60
#> 10  0.113261556 0.113261556 0.12047736 1.0000000 60
```
