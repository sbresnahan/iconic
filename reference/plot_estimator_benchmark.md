# Estimator benchmark figure

Produces a 4-panel publication figure benchmarking the eight ICONIC
estimators under unmeasured confounding and mediator-outcome
confounding. Panels: (A) total-effect bias vs confounding strength, (B)
NDE/NIE bias vs confounding strength, (C) NDE/NIE bias vs sample size,
(D) NIE Type I error vs confounding strength.

## Usage

``` r
plot_estimator_benchmark(
  panel_a,
  panel_b,
  panel_c,
  panel_d,
  conf_grid = c(0.2, 0.4, 0.6, 0.8, 1),
  n_grid = c(100, 200, 500, 1000),
  file = NULL,
  width = 8,
  height = 6
)
```

## Arguments

- panel_a:

  Result of `sweep_param("conf_str", ...)` for the total-effect panel.

- panel_b:

  Result of `sweep_mediation_param("conf_str", ...)` for the mediation
  confounding-strength panel.

- panel_c:

  Result of `sweep_mediation_param("n_samples", ...)` for the mediation
  sample-size panel.

- panel_d:

  Result of `sweep_mediation_null_by_conf(...)` for the Type I error
  panel.

- conf_grid:

  Numeric vector of confounding-strength values.

- n_grid:

  Numeric vector of sample sizes.

- file:

  Optional file path to save the figure (PDF or PNG).

- width:

  Figure width in inches.

- height:

  Figure height in inches.

## Value

A `patchwork` ggplot object.

## Examples

``` r
# Toy inputs standing in for sweep_param()/sweep_mediation_param() output
methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
  "IV2SLS2", "PGC2", "PGC2Gm")
it <- expand.grid(method = methods, pval = c(0.5, 0.8), rep = 1:3)
it$bias <- rnorm(nrow(it), 0, 0.05)
it$NDE_bias <- rnorm(nrow(it), 0, 0.05)
it$NIE_bias <- rnorm(nrow(it), 0, 0.05)
panel <- list(iter_bias = it)
panel_d <- expand.grid(method = methods, conf_str = c(0.5, 0.8))
panel_d$NIE_type1 <- runif(nrow(panel_d), 0, 0.1)
plot_estimator_benchmark(panel, panel, panel, panel_d,
  conf_grid = c(0.5, 0.8), n_grid = c(0.5, 0.8))
```
