# Heatmap of a sensitivity metric across the scenario grid

Draws a coverage x confounding-strength heatmap of a chosen metric for
one estimator, at a fixed number of confounders. Useful for reading off
where a method's bias/RMSE blows up.

## Usage

``` r
plot_sensitivity_heatmap(
  sens,
  metric = "rmse",
  method = "IV2SLS",
  k = NULL,
  title = NULL
)
```

## Arguments

- sens:

  Object returned by
  [`gan_sensitivity()`](https://seantbresnahan.com/iconic/reference/gan_sensitivity.md).

- metric:

  Column of `sens$summary` to map (e.g. "rmse", "abs_bias", "power").
  Default "rmse".

- method:

  Estimator to display. Default "IV2SLS".

- k:

  Number of confounders to slice on. Default: first in the grid.

- title:

  Optional plot title.

## Value

Invisibly returns `NULL`; called for the side effect of drawing a
sensitivity heatmap.

## Examples

``` r
sens <- gan_sensitivity(NULL, conf_grid = c(0.5, 0.8),
  coverage_grid = c(0.3, 0.7), n_iter = 2, n_samples = 100,
  n_features = 5)
plot_sensitivity_heatmap(sens)
```
