# Boxplot of estimate distributions from run_simulation()

Boxplot of estimate distributions from run_simulation()

## Usage

``` r
plot_estimate_distribution(
  sim_result,
  methods = iconic_method_order,
  title = NULL
)
```

## Arguments

- sim_result:

  Object returned by run_simulation().

- methods:

  Methods to include. Default: all five.

- title:

  Plot title.

## Value

A `ggplot` object showing the distribution of estimates by method.

## Examples

``` r
sim <- run_simulation(n_iter = 2, n_samples = 100, n_features = 5)
plot_estimate_distribution(sim)
```
