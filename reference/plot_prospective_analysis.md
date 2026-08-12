# Prospective analysis figure

Produces a 2-panel figure showing expected gains from adding genetic
instruments. Panel A: NDE bias vs instrument strength. Panel B:
prospective NDE/NIE estimates at target instrument strength.

## Usage

``` r
plot_prospective_analysis(prospect, file = NULL, width = 12, height = 8)
```

## Arguments

- prospect:

  Result of
  [`iconic_prospect()`](https://seantbresnahan.com/iconic/reference/iconic_prospect.md).

- file:

  Optional file path to save the figure.

- width:

  Figure width in inches.

- height:

  Figure height in inches.

## Value

A `patchwork` ggplot object.

## Examples

``` r
if (check_torch_setup()) {
  data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
    M = rnorm(100))
  pros <- iconic_prospect(data, n_iter = 2, gan_epochs = 5,
    gamma_G_grid = c(0.4, 0.8), run_rho_sweep = FALSE)
  plot_prospective_analysis(pros)
}
```
