# Diagnose data and determine estimator eligibility

Runs a battery of diagnostic checks on the user's data and returns an
eligibility report for all eight ICONIC estimators. The checks include:

## Usage

``` r
iconic_diagnose(
  data,
  fdr_level = 0.1,
  min_f = 10,
  k = NULL,
  g_threshold = NULL,
  gm_threshold = NULL,
  n_cores = 1,
  allow_no_proxy = TRUE
)
```

## Arguments

- data:

  An `iconic_data` object from
  [`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md).

- fdr_level:

  FDR level for NC screens. Default 0.10.

- min_f:

  Minimum partial F for instrument strength. Default 10. Used as the
  scalar threshold when `g_threshold`/`gm_threshold` are NULL (legacy
  behavior). When thresholds are supplied, `min_f` is unused for the
  panel decision (the threshold's `R` governs).

- k:

  Number of latent confounders assumed for the completeness check.
  Default `NULL`: infer from the data via Horn parallel analysis on the
  residualized-outcome correlation matrix
  ([`infer_confounding()`](https://seantbresnahan.com/iconic/reference/infer_confounding.md));
  falls back to 1 when inference is not possible (fewer than 5 outcome
  features). Supply an integer to override inference.

- g_threshold:

  Optional list `list(E = 0.5, R = 10)` controlling G-dependent method
  eligibility via the panel distribution. A method is eligible if at
  least fraction `E` of instruments have F_G \>= `R`. Default NULL:
  legacy scalar behavior (median F_G vs `min_f`).

- gm_threshold:

  Optional list `list(E = 0.5, R = 10)` controlling Gm-dependent method
  eligibility (IV2SLS2, PGC2Gm) via the panel distribution. A method is
  eligible if at least fraction `E` of mediators have F_Gm \>= `R`.
  Default NULL: legacy scalar behavior (median F_Gm vs `min_f`).

- n_cores:

  Number of parallel workers for NC validity screens and panel
  instrument-strength computation. Default 1 (sequential). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix and a PSOCK cluster on Windows.

- allow_no_proxy:

  Logical: when `TRUE` (default), proceed with a message if no
  instruments and no NCs are supplied (only UNADJ is eligible). When
  `FALSE`, error instead.

## Value

An `iconic_diagnosis` S3 object: a named list with
`$instrument_strength`, `$nc_validity`, `$nc_independence`,
`$nc_independence_gm`, `$completeness`, `$eligibility`, and `$summary`.
`$instrument_strength$F_G` and `$instrument_strength$F_Gm` are numeric
vectors containing the full panel distributions (one entry per
instrument / mediator).

## Details

- **Instrument strength**: first-stage partial F-statistics for G
  (exposure instrument) and Gm (mediator instrument), with the
  Stock-Yogo weak-instrument threshold (F \>= 10).

- **NC validity (A1)**: screens each negative-control feature for
  association with the exposure X (the empirically testable projection
  of "W independent of X given C, U";
  [`nc_validity_screen()`](https://seantbresnahan.com/iconic/reference/nc_validity_screen.md)).

- **NC independence (A2)**: screens each NC for association with the
  instrument G (the empirically testable projection of "W independent of
  G given C, U";
  [`nc_independence_check()`](https://seantbresnahan.com/iconic/reference/nc_independence_check.md)).

- **NC independence (A2')**: screens each NC for association with the
  mediator instrument Gm (the empirically testable projection of "W
  independent of G_m given C, U";
  [`nc_independence_check_gm()`](https://seantbresnahan.com/iconic/reference/nc_independence_check_gm.md)).

- **Path completeness**: checks whether the valid NC panel has enough
  features to span the confounder subspace and captures the confounder
  covariance
  ([`nc_completeness_check()`](https://seantbresnahan.com/iconic/reference/nc_completeness_check.md)).

## Defaults

|  |  |  |
|----|----|----|
| **Parameter** | **Default** | **Source** |
| `fdr_level` | 0.10 | BH-FDR level for NC screens |
| `min_f` | 10 | Stock-Yogo weak-instrument threshold |
| `k` | 1 | Single-confounder assumption (typical mediation) |
| `g_threshold` | NULL | Legacy scalar (median F_G vs min_f) |
| `gm_threshold` | NULL | Legacy scalar (median F_Gm vs min_f) |
| `n_cores` | 1 | Sequential |
| `allow_no_proxy` | TRUE | Proceed with UNADJ-only when no IV/NC |

## Examples

``` r
set.seed(1)
data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
G = rnorm(100), W = matrix(rnorm(100 * 10), 10, 100))
diag <- iconic_diagnose(data)
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
#> Estimating features: 10 tasks (sequential)
#>  Estimating features: 10% (1/10) [0s]
#>  Estimating features: 20% (2/10) [0s]
#>  Estimating features: 30% (3/10) [0s]
#>  Estimating features: 40% (4/10) [0s]
#>  Estimating features: 50% (5/10) [0s]
#>  Estimating features: 60% (6/10) [0.1s]
#>  Estimating features: 70% (7/10) [0.1s]
#>  Estimating features: 80% (8/10) [0.1s]
#>  Estimating features: 90% (9/10) [0.1s]
#>  Estimating features: 100% (10/10) [0.1s]
#> omega_1 NC coverage: 10 tasks (sequential)
#>  omega_1 NC coverage: 10% (1/10) [0s]
#>  omega_1 NC coverage: 20% (2/10) [0s]
#>  omega_1 NC coverage: 30% (3/10) [0s]
#>  omega_1 NC coverage: 40% (4/10) [0s]
#>  omega_1 NC coverage: 50% (5/10) [0s]
#>  omega_1 NC coverage: 60% (6/10) [0s]
#>  omega_1 NC coverage: 70% (7/10) [0s]
#>  omega_1 NC coverage: 80% (8/10) [0.1s]
#>  omega_1 NC coverage: 90% (9/10) [0.1s]
#>  omega_1 NC coverage: 100% (10/10) [0.1s]
#> omega_2 NC coverage: 10 tasks (sequential)
#>  omega_2 NC coverage: 10% (1/10) [0s]
#>  omega_2 NC coverage: 20% (2/10) [0s]
#>  omega_2 NC coverage: 30% (3/10) [0s]
#>  omega_2 NC coverage: 40% (4/10) [0s]
#>  omega_2 NC coverage: 50% (5/10) [0s]
#>  omega_2 NC coverage: 60% (6/10) [0s]
#>  omega_2 NC coverage: 70% (7/10) [0s]
#>  omega_2 NC coverage: 80% (8/10) [0.1s]
#>  omega_2 NC coverage: 90% (9/10) [0.1s]
#>  omega_2 NC coverage: 100% (10/10) [0.1s]
#> k permutation analysis: 100 tasks (sequential)
#>  k permutation analysis: 10% (10/100) [0s]
#>  k permutation analysis: 20% (20/100) [0s]
#>  k permutation analysis: 30% (30/100) [0s]
#>  k permutation analysis: 40% (40/100) [0s]
#>  k permutation analysis: 50% (50/100) [0s]
#>  k permutation analysis: 60% (60/100) [0s]
#>  k permutation analysis: 70% (70/100) [0s]
#>  k permutation analysis: 80% (80/100) [0s]
#>  k permutation analysis: 90% (90/100) [0s]
#>  k permutation analysis: 100% (100/100) [0s]
#> NC capture null: 200 tasks (sequential)
#>  NC capture null: 10% (20/200) [0.3s]
#>  NC capture null: 20% (40/200) [0.7s]
#>  NC capture null: 30% (60/200) [1.1s]
#>  NC capture null: 40% (80/200) [1.4s]
#>  NC capture null: 50% (100/200) [1.8s]
#>  NC capture null: 60% (120/200) [2.1s]
#>  NC capture null: 70% (140/200) [2.5s]
#>  NC capture null: 80% (160/200) [2.8s]
#>  NC capture null: 90% (180/200) [3.2s]
#>  NC capture null: 100% (200/200) [3.6s]
#> iconic_diagnose complete. Call summary() or print() on the result for the full diagnosis.
print(diag)
#> <iconic_diagnosis>
#> Diagnostic summary:
#>  G (exposure instrument): partial F = 0.7 (WEAK)
#>  NC validity (A1): 10/10 controls valid (0 flagged)
#>  NC independence (A2): 10/10 controls valid (0 flagged)
#>  Completeness: 10 valid NCs vs k=3 (inferred) -> satisfied
#>   Capture: incremental R^2 = 0.103 (p = 0.455) -> weak
#>   Support: R^2(U~|W) = 0.097 -> narrow coverage (0/10 controls add unique coverage)
#> 
#>  Eligible estimators: 3/8
#>   UNADJ, DIRECT, COCA 
```
