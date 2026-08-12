# Sweep negative-control validity diagnostics

Runs simulation sweeps for the four empirical NC validity diagnostics
(A1, A2, A2', A3) and returns data frames ready for plotting with
[`plot_nc_validity_diagnostics()`](https://seantbresnahan.com/iconic/reference/plot_nc_validity_diagnostics.md).

## Usage

``` r
sweep_nc_validity(
  n_samples = 500,
  n_iter = 50,
  phi_val = 0.8,
  contam_grid = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
  meqtl_grid = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
  eqtl_grid = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
  k_grid = c(1, 2, 3),
  n_valid_grid = c(1, 2, 3),
  n_cores = 1
)
```

## Arguments

- n_samples:

  Number of synthetic samples per iteration.

- n_iter:

  Number of replications per sweep point.

- phi_val:

  Mediator-instrument strength (for A2' panel).

- contam_grid:

  X-\>W contamination strength grid (A1).

- meqtl_grid:

  G-\>W (meQTL) strength grid (A2).

- eqtl_grid:

  Gm-\>W (eQTL) strength grid (A2').

- k_grid:

  Number of confounders grid (A3).

- n_valid_grid:

  Number of valid controls grid (A3).

- n_cores:

  Number of parallel workers for the replicate loops. Default 1
  (sequential). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix and a PSOCK cluster on Windows.

## Value

A list with elements `panel_a`, `panel_b`, `panel_c`, `panel_d` (data
frames).

## Examples

``` r
# \donttest{
panels <- sweep_nc_validity(n_samples = 100, n_iter = 2, k_grid = 1,
  n_valid_grid = 1, contam_grid = c(0, 0.1), meqtl_grid = c(0, 0.1),
  eqtl_grid = c(0, 0.1))
#> Running NC validity sweeps (n_iter = 2, n = 100) ...
#>  Panel A: nc_validity_screen (A1)
#>  Panel A replicates: 2 tasks (sequential)
#> NC validity screen: 10 tasks (sequential)
#>  NC validity screen: 10% (1/10) [0s]
#>  NC validity screen: 20% (2/10) [0s]
#>  NC validity screen: 30% (3/10) [0s]
#>  NC validity screen: 40% (4/10) [0s]
#>  NC validity screen: 50% (5/10) [0s]
#>  NC validity screen: 60% (6/10) [0s]
#>  NC validity screen: 70% (7/10) [0s]
#>  NC validity screen: 80% (8/10) [0s]
#>  NC validity screen: 90% (9/10) [0s]
#>  NC validity screen: 100% (10/10) [0s]
#>   Panel A replicates: 50% (1/2) [0s]
#> NC validity screen: 10 tasks (sequential)
#>  NC validity screen: 10% (1/10) [0s]
#>  NC validity screen: 20% (2/10) [0s]
#>  NC validity screen: 30% (3/10) [0s]
#>  NC validity screen: 40% (4/10) [0s]
#>  NC validity screen: 50% (5/10) [0s]
#>  NC validity screen: 60% (6/10) [0s]
#>  NC validity screen: 70% (7/10) [0s]
#>  NC validity screen: 80% (8/10) [0s]
#>  NC validity screen: 90% (9/10) [0s]
#>  NC validity screen: 100% (10/10) [0s]
#>   Panel A replicates: 100% (2/2) [0s]
#>  Panel A replicates: 2 tasks (sequential)
#> NC validity screen: 10 tasks (sequential)
#>  NC validity screen: 10% (1/10) [0s]
#>  NC validity screen: 20% (2/10) [0s]
#>  NC validity screen: 30% (3/10) [0s]
#>  NC validity screen: 40% (4/10) [0s]
#>  NC validity screen: 50% (5/10) [0s]
#>  NC validity screen: 60% (6/10) [0s]
#>  NC validity screen: 70% (7/10) [0s]
#>  NC validity screen: 80% (8/10) [0s]
#>  NC validity screen: 90% (9/10) [0s]
#>  NC validity screen: 100% (10/10) [0s]
#>   Panel A replicates: 50% (1/2) [0s]
#> NC validity screen: 10 tasks (sequential)
#>  NC validity screen: 10% (1/10) [0s]
#>  NC validity screen: 20% (2/10) [0s]
#>  NC validity screen: 30% (3/10) [0s]
#>  NC validity screen: 40% (4/10) [0s]
#>  NC validity screen: 50% (5/10) [0s]
#>  NC validity screen: 60% (6/10) [0s]
#>  NC validity screen: 70% (7/10) [0s]
#>  NC validity screen: 80% (8/10) [0s]
#>  NC validity screen: 90% (9/10) [0s]
#>  NC validity screen: 100% (10/10) [0s]
#>   Panel A replicates: 100% (2/2) [0s]
#>  Panel B: nc_independence_check (A2)
#>  Panel B replicates: 2 tasks (sequential)
#> NC independence (G): 10 tasks (sequential)
#>  NC independence (G): 10% (1/10) [0s]
#>  NC independence (G): 20% (2/10) [0s]
#>  NC independence (G): 30% (3/10) [0s]
#>  NC independence (G): 40% (4/10) [0s]
#>  NC independence (G): 50% (5/10) [0s]
#>  NC independence (G): 60% (6/10) [0s]
#>  NC independence (G): 70% (7/10) [0s]
#>  NC independence (G): 80% (8/10) [0s]
#>  NC independence (G): 90% (9/10) [0s]
#>  NC independence (G): 100% (10/10) [0s]
#>   Panel B replicates: 50% (1/2) [0s]
#> NC independence (G): 10 tasks (sequential)
#>  NC independence (G): 10% (1/10) [0s]
#>  NC independence (G): 20% (2/10) [0s]
#>  NC independence (G): 30% (3/10) [0s]
#>  NC independence (G): 40% (4/10) [0s]
#>  NC independence (G): 50% (5/10) [0s]
#>  NC independence (G): 60% (6/10) [0s]
#>  NC independence (G): 70% (7/10) [0s]
#>  NC independence (G): 80% (8/10) [0s]
#>  NC independence (G): 90% (9/10) [0s]
#>  NC independence (G): 100% (10/10) [0s]
#>   Panel B replicates: 100% (2/2) [0s]
#>  Panel B replicates: 2 tasks (sequential)
#> NC independence (G): 10 tasks (sequential)
#>  NC independence (G): 10% (1/10) [0s]
#>  NC independence (G): 20% (2/10) [0s]
#>  NC independence (G): 30% (3/10) [0s]
#>  NC independence (G): 40% (4/10) [0s]
#>  NC independence (G): 50% (5/10) [0s]
#>  NC independence (G): 60% (6/10) [0s]
#>  NC independence (G): 70% (7/10) [0s]
#>  NC independence (G): 80% (8/10) [0s]
#>  NC independence (G): 90% (9/10) [0s]
#>  NC independence (G): 100% (10/10) [0s]
#>   Panel B replicates: 50% (1/2) [0s]
#> NC independence (G): 10 tasks (sequential)
#>  NC independence (G): 10% (1/10) [0s]
#>  NC independence (G): 20% (2/10) [0s]
#>  NC independence (G): 30% (3/10) [0s]
#>  NC independence (G): 40% (4/10) [0s]
#>  NC independence (G): 50% (5/10) [0s]
#>  NC independence (G): 60% (6/10) [0s]
#>  NC independence (G): 70% (7/10) [0s]
#>  NC independence (G): 80% (8/10) [0s]
#>  NC independence (G): 90% (9/10) [0s]
#>  NC independence (G): 100% (10/10) [0s]
#>   Panel B replicates: 100% (2/2) [0s]
#>  Panel C: nc_completeness_check (A3)
#>  Panel C replicates: 2 tasks (sequential)
#>   Panel C replicates: 50% (1/2) [0s]
#>   Panel C replicates: 100% (2/2) [0s]
#> NC capture null: 1000 tasks (sequential)
#>  NC capture null: 10% (100/1000) [2.2s]
#>  NC capture null: 20% (200/1000) [4.4s]
#>  NC capture null: 30% (300/1000) [6.6s]
#>  NC capture null: 40% (400/1000) [8.8s]
#>  NC capture null: 50% (500/1000) [11s]
#>  NC capture null: 60% (600/1000) [13.2s]
#>  NC capture null: 70% (700/1000) [15.4s]
#>  NC capture null: 80% (800/1000) [17.6s]
#>  NC capture null: 90% (900/1000) [19.8s]
#>  NC capture null: 100% (1000/1000) [22s]
#>  Panel D: nc_independence_check_gm (A2')
#>  Panel D replicates: 2 tasks (sequential)
#> NC independence (Gm): 10 tasks (sequential)
#>  NC independence (Gm): 10% (1/10) [0s]
#>  NC independence (Gm): 20% (2/10) [0s]
#>  NC independence (Gm): 30% (3/10) [0s]
#>  NC independence (Gm): 40% (4/10) [0s]
#>  NC independence (Gm): 50% (5/10) [0s]
#>  NC independence (Gm): 60% (6/10) [0s]
#>  NC independence (Gm): 70% (7/10) [0s]
#>  NC independence (Gm): 80% (8/10) [0s]
#>  NC independence (Gm): 90% (9/10) [0s]
#>  NC independence (Gm): 100% (10/10) [0s]
#>   Panel D replicates: 50% (1/2) [0s]
#> NC independence (Gm): 10 tasks (sequential)
#>  NC independence (Gm): 10% (1/10) [0s]
#>  NC independence (Gm): 20% (2/10) [0s]
#>  NC independence (Gm): 30% (3/10) [0s]
#>  NC independence (Gm): 40% (4/10) [0s]
#>  NC independence (Gm): 50% (5/10) [0s]
#>  NC independence (Gm): 60% (6/10) [0s]
#>  NC independence (Gm): 70% (7/10) [0s]
#>  NC independence (Gm): 80% (8/10) [0s]
#>  NC independence (Gm): 90% (9/10) [0s]
#>  NC independence (Gm): 100% (10/10) [0s]
#>   Panel D replicates: 100% (2/2) [0s]
#>  Panel D replicates: 2 tasks (sequential)
#> NC independence (Gm): 10 tasks (sequential)
#>  NC independence (Gm): 10% (1/10) [0s]
#>  NC independence (Gm): 20% (2/10) [0s]
#>  NC independence (Gm): 30% (3/10) [0s]
#>  NC independence (Gm): 40% (4/10) [0s]
#>  NC independence (Gm): 50% (5/10) [0s]
#>  NC independence (Gm): 60% (6/10) [0s]
#>  NC independence (Gm): 70% (7/10) [0s]
#>  NC independence (Gm): 80% (8/10) [0s]
#>  NC independence (Gm): 90% (9/10) [0s]
#>  NC independence (Gm): 100% (10/10) [0s]
#>   Panel D replicates: 50% (1/2) [0s]
#> NC independence (Gm): 10 tasks (sequential)
#>  NC independence (Gm): 10% (1/10) [0s]
#>  NC independence (Gm): 20% (2/10) [0s]
#>  NC independence (Gm): 30% (3/10) [0s]
#>  NC independence (Gm): 40% (4/10) [0s]
#>  NC independence (Gm): 50% (5/10) [0s]
#>  NC independence (Gm): 60% (6/10) [0s]
#>  NC independence (Gm): 70% (7/10) [0s]
#>  NC independence (Gm): 80% (8/10) [0s]
#>  NC independence (Gm): 90% (9/10) [0s]
#>  NC independence (Gm): 100% (10/10) [0s]
#>   Panel D replicates: 100% (2/2) [0s]
#>  Done in 22.2 s
panels$panel_a
#>   contamination violated_mean violated_sd confounding_mean confounding_sd
#> 1           0.0           0.0   0.0000000                1              0
#> 2           0.1           0.3   0.1414214                1              0
# }
```
