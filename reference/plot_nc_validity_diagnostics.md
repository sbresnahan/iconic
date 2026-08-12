# NC validity diagnostics figure

Produces a 4-panel figure sweeping the four empirical NC validity
diagnostics. Panels: (A) A1 W perp X\|C, (B) A2 W perp G\|C, (C) A3
dim(W_valid) \>= k completeness grid, (D) A2' W perp Gm\|C.

## Usage

``` r
plot_nc_validity_diagnostics(panels, file = NULL, width = 8, height = 6)
```

## Arguments

- panels:

  List returned by
  [`sweep_nc_validity()`](https://seantbresnahan.com/iconic/reference/sweep_nc_validity.md).

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
# Toy panels standing in for sweep_nc_validity() output
panels <- list(
  panel_a = data.frame(contamination = c(0, 0.1),
    violated_mean = c(0.1, 0.3), violated_sd = 0.05,
    confounding_mean = c(0.1, 0.12), confounding_sd = 0.05),
  panel_b = data.frame(meqtl = c(0, 0.1),
    violated_mean = c(0.1, 0.3), violated_sd = 0.05,
    clean_mean = c(0.1, 0.12), clean_sd = 0.05),
  panel_c = expand.grid(n_valid = 1:2, k = 1:2),
  panel_d = data.frame(eqtl = c(0, 0.1),
    violated_mean = c(0.1, 0.3), violated_sd = 0.05,
    clean_mean = c(0.1, 0.12), clean_sd = 0.05))
panels$panel_c$pgc_bias <- c(0.05, 0.20, 0.08, 0.25)
panels$panel_c$completeness <- c("satisfied", "under-identified",
  "satisfied", "under-identified")
plot_nc_validity_diagnostics(panels)
```
