# Baseline bias distribution (single-setting hero plot)

Baseline bias distribution (single-setting hero plot)

## Usage

``` r
plot_bias_distribution(sim_result, methods = iconic_method_order, title = NULL)
```

## Arguments

- sim_result:

  Object returned by run_simulation().

- methods:

  Methods to include. Default: all five.

- title:

  Plot title. If NULL a default is constructed.

## Value

A `ggplot` object showing the distribution of bias estimates by method.

## Examples

``` r
sim <- run_simulation(n_iter = 2, n_samples = 100, n_features = 5)
plot_bias_distribution(sim)
```
