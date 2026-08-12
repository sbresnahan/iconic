# Degradation surface figure

Produces a 3-panel figure showing estimator bias under simultaneous
instrument exogeneity violations. Panels: (A) IV2SLS2 NDE bias surface,
(B) PGC2Gm NDE bias surface, (C) crossover map showing which estimator
has lower \|NDE bias\| at each grid cell.

## Usage

``` r
plot_degradation_surface(
  results,
  rho_G1_grid = c(0, 0.1, 0.2, 0.3, 0.5),
  rho_G2_grid = c(0, 0.1, 0.2, 0.3, 0.5),
  omega_facet = TRUE,
  file = NULL,
  width = 15,
  height = 5.5
)
```

## Arguments

- results:

  Data frame with columns `rho_G1`, `rho_G2`, `method`, `NDE_bias`,
  `NIE_bias`, and derived `NDE_abs`, `NIE_abs`. Typically built by
  sweeping
  [`sweep_mediation_param()`](https://seantbresnahan.com/iconic/reference/sweep_mediation_param.md)
  across a rho_G1 x rho_G2 grid.

- rho_G1_grid:

  Numeric vector of exposure-instrument violation values.

- rho_G2_grid:

  Numeric vector of mediator-instrument violation values.

- omega_facet:

  Logical: when the sensitivity surface swept `omega_1` / `omega_2`
  (more than one distinct value), facet the rho_G1 x rho_G2 heatmaps by
  the omega cell. Default `TRUE`. Ignored when omega was not swept.

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
results <- expand.grid(rho_G1 = c(0, 0.2), rho_G2 = c(0, 0.2),
  method = c("IV2SLS2", "PGC2Gm"))
results$NDE_bias <- rnorm(nrow(results), 0, 0.05)
results$NIE_bias <- rnorm(nrow(results), 0, 0.05)
plot_degradation_surface(results, rho_G1_grid = c(0, 0.2),
  rho_G2_grid = c(0, 0.2))
```
