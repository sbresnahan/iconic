# Plot absolute bias vs a swept parameter

Plot absolute bias vs a swept parameter

## Usage

``` r
plot_bias(
  sweep_df,
  param_label = "Parameter Value",
  methods = iconic_method_order,
  title = "Absolute Bias vs Parameter"
)
```

## Arguments

- sweep_df:

  Data frame from sweep_param()\$summary.

- param_label:

  X-axis label. Default "Parameter Value".

- methods:

  Methods to include. Default: all five.

- title:

  Plot title.

## Value

Invisibly returns `NULL`; called for the side effect of plotting
absolute bias vs a parameter.

## Examples

``` r
sweep <- sweep_param("conf_str", c(0.5, 0.8), n_iter = 2,
  n_samples = 100, n_features = 5)
plot_bias(sweep$summary)
```
