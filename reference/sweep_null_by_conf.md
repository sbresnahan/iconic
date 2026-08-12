# Sweep Type I error rate across confounding strength levels

Sweep Type I error rate across confounding strength levels

## Usage

``` r
sweep_null_by_conf(
  conf_grid = c(0.2, 0.4, 0.6, 0.8, 1),
  n_iter = 100,
  n_samples = 500,
  n_features = 20,
  w_signal = 0.7,
  feat_cor = 0,
  base_seed = 900,
  n_cores = 1,
  alpha = 0.05
)
```

## Arguments

- conf_grid:

  Numeric vector of confounding strength values. Default c(0.2, 0.4,
  0.6, 0.8, 1.0).

- n_iter:

  Replicates per conf_str value. Default 100.

- n_samples:

  Observations per replicate. Default 500.

- n_features:

  Features per replicate. Default 20.

- w_signal:

  Proxy quality omega. Default 0.70.

- feat_cor:

  Within-module feature correlation. Default 0.

- base_seed:

  Seed offset. Default 900.

- n_cores:

  Parallel workers. Default 1.

- alpha:

  Significance threshold. Default 0.05.

## Value

A data frame with columns: conf_str, method, type1_error.

## Examples

``` r
t1e <- sweep_null_by_conf(c(0.2, 0.8), n_iter = 3, n_samples = 100)
plot_type1_vs_conf(t1e)
```
