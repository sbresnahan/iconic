# Sweep negative-control validity diagnostics

Runs simulation sweeps for the empirical NC validity diagnostics and
returns data frames ready for plotting with
[`plot_nc_validity_diagnostics()`](https://seantbresnahan.com/iconic/reference/plot_nc_validity_diagnostics.md).
Panels: (A) A1 W perp X \| C, (B) A2 W perp G \| C, (C) A2' W perp Gm \|
C, (D) A3 covariance-capture versus true coverage omega, (E) A3 support
R2(Utilde \| W) versus omega.

## Usage

``` r
sweep_nc_validity(
  n_samples = 500,
  n_iter = 50,
  phi_val = 0.8,
  contam_grid = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
  meqtl_grid = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
  eqtl_grid = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
  omega_grid = c(0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1),
  n_perm = 200,
  cs_confounders = 2,
  include_a3_grid = FALSE,
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

  Mediator-instrument strength (for the A2' panel and the
  capture/support sweeps).

- contam_grid:

  X-\>W contamination strength grid (A1).

- meqtl_grid:

  G-\>W (meQTL) strength grid (A2).

- eqtl_grid:

  Gm-\>W (eQTL) strength grid (A2').

- omega_grid:

  True negative-control coverage grid for the capture (D) and
  support (E) sweeps. Includes 0 by default so the permutation-null
  calibration is visible.

- n_perm:

  Permutations per replicate for the capture test. Default 200, matching
  [`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md).

- cs_confounders:

  Number of latent confounders for the capture/support sweeps (distinct
  loadings). Default 2.

- include_a3_grid:

  Logical; also run the legacy A3 dimensional grid (PGC bias over
  n_valid x k). Default FALSE (no longer plotted).

- k_grid, n_valid_grid:

  Grids for the legacy A3 dimensional sweep.

- n_cores:

  Number of parallel workers for the replicate loops. Default 1
  (sequential). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix and a PSOCK cluster on Windows.

## Value

A list with elements `panel_a`, `panel_b`, `panel_d`,
`panel_capture_support` (per-replicate raw), and optionally `panel_c`
when `include_a3_grid = TRUE`.

## Examples

``` r
# \donttest{
panels <- sweep_nc_validity(n_samples = 100, n_iter = 2,
  contam_grid = c(0, 0.1), meqtl_grid = c(0, 0.1),
  eqtl_grid = c(0, 0.1), omega_grid = c(0, 0.5), n_perm = 20)
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
#>  Panel C: nc_independence_check_gm (A2')
#>  Panel C replicates: 2 tasks (sequential)
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
#>   Panel C replicates: 50% (1/2) [0s]
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
#>   Panel C replicates: 100% (2/2) [0s]
#>  Panel C replicates: 2 tasks (sequential)
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
#>   Panel C replicates: 50% (1/2) [0s]
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
#>   Panel C replicates: 100% (2/2) [0s]
#>  Panels D-E: nc_completeness_capture + nc_support_check vs omega
#>   omega = 0: 2 tasks (sequential)
#> NC capture null: 20 tasks (sequential)
#>  NC capture null: 10% (2/20) [0.1s]
#>  NC capture null: 20% (4/20) [0.1s]
#>  NC capture null: 30% (6/20) [0.1s]
#>  NC capture null: 40% (8/20) [0.2s]
#>  NC capture null: 50% (10/20) [0.2s]
#>  NC capture null: 60% (12/20) [0.3s]
#>  NC capture null: 70% (14/20) [0.3s]
#>  NC capture null: 80% (16/20) [0.4s]
#>  NC capture null: 90% (18/20) [0.4s]
#>  NC capture null: 100% (20/20) [0.5s]
#>    omega = 0: 50% (1/2) [0.5s]
#> NC capture null: 20 tasks (sequential)
#>  NC capture null: 10% (2/20) [0s]
#>  NC capture null: 20% (4/20) [0.1s]
#>  NC capture null: 30% (6/20) [0.1s]
#>  NC capture null: 40% (8/20) [0.2s]
#>  NC capture null: 50% (10/20) [0.2s]
#>  NC capture null: 60% (12/20) [0.3s]
#>  NC capture null: 70% (14/20) [0.3s]
#>  NC capture null: 80% (16/20) [0.4s]
#>  NC capture null: 90% (18/20) [0.4s]
#>  NC capture null: 100% (20/20) [0.5s]
#>    omega = 0: 100% (2/2) [1.1s]
#>   omega = 0.5: 2 tasks (sequential)
#> NC capture null: 20 tasks (sequential)
#>  NC capture null: 10% (2/20) [0s]
#>  NC capture null: 20% (4/20) [0.1s]
#>  NC capture null: 30% (6/20) [0.1s]
#>  NC capture null: 40% (8/20) [0.2s]
#>  NC capture null: 50% (10/20) [0.2s]
#>  NC capture null: 60% (12/20) [0.3s]
#>  NC capture null: 70% (14/20) [0.3s]
#>  NC capture null: 80% (16/20) [0.4s]
#>  NC capture null: 90% (18/20) [0.4s]
#>  NC capture null: 100% (20/20) [0.5s]
#>    omega = 0.5: 50% (1/2) [0.5s]
#> NC capture null: 20 tasks (sequential)
#>  NC capture null: 10% (2/20) [0s]
#>  NC capture null: 20% (4/20) [0.1s]
#>  NC capture null: 30% (6/20) [0.1s]
#>  NC capture null: 40% (8/20) [0.2s]
#>  NC capture null: 50% (10/20) [0.2s]
#>  NC capture null: 60% (12/20) [0.3s]
#>  NC capture null: 70% (14/20) [0.3s]
#>  NC capture null: 80% (16/20) [0.4s]
#>  NC capture null: 90% (18/20) [0.4s]
#>  NC capture null: 100% (20/20) [0.5s]
#>    omega = 0.5: 100% (2/2) [1s]
#>  Done in 2.3 s
panels$panel_a
#>   contamination violated_mean violated_sd confounding_mean confounding_sd
#> 1           0.0           0.0   0.0000000                1              0
#> 2           0.1           0.3   0.1414214                1              0
panels$panel_capture_support
#>   omega rep capture_R2 capture_p    null_R2 support_R2 frac_adds_coverage
#> 1   0.0   1 0.07755442      0.75 0.10999436 0.11589077                0.0
#> 2   0.0   2 0.04619164      0.90 0.09367601 0.05081462                0.0
#> 3   0.5   1 0.49594649      0.00 0.08636304 0.65417453                0.4
#> 4   0.5   2 0.48886836      0.00 0.10329703 0.63588698                0.2
# }
```
