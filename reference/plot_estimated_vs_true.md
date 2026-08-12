# Plot estimated effect vs true effect

Plot estimated effect vs true effect

## Usage

``` r
plot_estimated_vs_true(
  sweep_df,
  methods = iconic_method_order,
  title = "Estimated vs True Effect"
)
```

## Arguments

- sweep_df:

  Data frame from sweep_param()\$summary.

- methods:

  Methods to plot. Default: all five.

- title:

  Plot title.

## Value

Invisibly returns `NULL`; called for the side effect of plotting
estimated vs true effects.

## Examples

``` r
sweep <- sweep_param("beta_X", c(0.1, 0.3), n_iter = 2,
  n_samples = 100, n_features = 5)
plot_estimated_vs_true(sweep$summary)
```
