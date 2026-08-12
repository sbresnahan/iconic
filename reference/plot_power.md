# Plot detection rate (power) vs a swept parameter

Plot detection rate (power) vs a swept parameter

## Usage

``` r
plot_power(
  sweep_df,
  param_label = "Parameter Value",
  methods = iconic_method_order,
  title = "Detection Rate vs Parameter",
  legend_pos = "bottomleft"
)
```

## Arguments

- sweep_df:

  Data frame from sweep_param()\$summary.

- param_label:

  X-axis label.

- methods:

  Methods to include.

- title:

  Plot title.

- legend_pos:

  Legend position. Default "bottomleft".

## Value

Invisibly returns `NULL`; called for the side effect of plotting
detection rate vs a parameter.

## Examples

``` r
sweep <- sweep_param("conf_str", c(0.5, 0.8), n_iter = 2,
  n_samples = 100, n_features = 5)
plot_power(sweep$summary)
```
