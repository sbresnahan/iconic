# Recommend the best causal estimator for the user's data

Ranks all eligible estimators directly on per-estimand robustness (NDE
and NIE separately) from a sensitivity analysis. Returns the top-ranked
estimator with a transparent rationale. Every estimator is ranked on its
measured robustness to assumption violations, not on an a priori
identification-strength class.

## Usage

``` r
iconic_recommend(
  data,
  diagnosis = NULL,
  estimate = NULL,
  sensitivity = NULL,
  criterion = c("combined", "minimax_bias", "ci_coverage"),
  completeness_penalty = c(satisfied = 1, borderline = 0.7, `weak-capture` = 0.5,
    `under-identified` = 0),
  min_f = 10,
  g_threshold = NULL,
  gm_threshold = NULL,
  auto_sensitivity = TRUE,
  rho_G1_grid = c(0, 0.1, 0.2, 0.3, 0.5),
  rho_G2_grid = c(0, 0.1, 0.2, 0.3, 0.5),
  omega_1 = c(0.3, 0.7, 1),
  omega_2 = c(0.3, 0.7, 1),
  n_iter_sens = 30,
  gan_epochs = 100,
  n_cores = 1,
  verbose = FALSE
)
```

## Arguments

- data:

  An `iconic_data` object.

- diagnosis:

  Optional `iconic_diagnosis` from
  [`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md).
  If `NULL`, auto-eligibility is computed from the data.

- estimate:

  Optional estimate data frame from
  [`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md).
  Used to report point estimates alongside the recommendation.

- sensitivity:

  Optional `iconic_sensitivity` from
  [`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md).
  Used for robustness ranking. When `NULL` and `auto_sensitivity = TRUE`
  (the default), the sensitivity suite is run automatically so the
  recommendation is robustness-based out of the box.

- criterion:

  Character: `"combined"` (default), `"minimax_bias"`, or
  `"ci_coverage"`. Controls how robustness is ranked. `"combined"`
  blends bias and CI coverage; `"minimax_bias"` ranks by worst-case bias
  across the violation grid; `"ci_coverage"` ranks by CI coverage of the
  true effect. When `sensitivity` is supplied, also populates
  `$per_scenario` (best estimator at the origin vs at violation cells).

- completeness_penalty:

  Named numeric vector with values in \\\[0, 1\]\\ mapping the graded
  completeness verdict to a confidence multiplier applied to
  bridge-dependent estimators (DIRECT, COCA, PGC, PGC2, PGC2Gm). Default
  `c(satisfied = 1.0, borderline = 0.7, "weak-capture" = 0.5, "under-identified" = 0)`.
  Instrument-only estimators (IV2SLS, IV2SLS2) are not discounted.
  Override to sensitivity-check the discount.

- min_f:

  Instrument-strength gate (first-stage F threshold) used when
  `diagnosis` is `NULL` and
  [`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md)
  is run internally. Default 10. Also used for the requirement labels
  when the supplied `diagnosis` predates the stored `min_f` field.

- g_threshold, gm_threshold:

  Optional instrument-selection thresholds passed to
  [`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md)
  when `diagnosis` is `NULL`. Default `NULL` (use the `iconic_diagnose`
  defaults).

- auto_sensitivity:

  Logical: when `TRUE` (default) and `sensitivity = NULL`, run
  [`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md)
  internally to obtain the robustness surface. Requires the torch
  backend; when torch is unavailable the function falls back to
  eligibility-only ranking with a message. Set `FALSE` to skip the
  auto-run.

- rho_G1_grid, rho_G2_grid:

  Instrument-exogeneity violation grid used for the auto-run sensitivity
  suite. Default `c(0, 0.1, 0.2, 0.3, 0.5)`.

- omega_1, omega_2:

  Negative-control coverage grid for the auto-run. Default
  `c(0.3, 0.7, 1.0)`, swept on the diagonal (`omega_1 == omega_2`).

- n_iter_sens:

  Replicates per grid cell for the auto-run. Default 30.

- gan_epochs:

  GAN training epochs for the auto-run. Default 100.

- n_cores:

  Cores for the auto-run sensitivity sweep. Default 1.

- verbose:

  Logical: print progress messages. Default `FALSE` (quiet). Also
  silences the auto-run sensitivity suite.

## Value

An `iconic_recommendation` S3 object: a named list with `$ranking` (data
frame: estimator, eligible, rank, per-estimand robustness scores,
`composite`, `confidence_mult`, `final_score`, rationale),
`$recommended` (top eligible, non-naive estimator by the data-driven
composite score), `$recommended_NDE` and `$recommended_NIE` (top
estimator per estimand), `$per_scenario` (when `sensitivity` is
supplied), `$completeness_penalty`, and `$summary`.

## Details

When `sensitivity` is supplied (from
[`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md)),
eligible estimators are ranked by robustness: the estimator whose
estimates degrade least across the assumption-violation surface ranks
higher. NDE and NIE robustness are computed separately, so a mediation
estimator is ranked on the estimand of interest rather than a pooled
maximum.

## Examples

``` r
data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100*10), 10, 100),
G = rnorm(100), W = matrix(rnorm(100*10), 10, 100))
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
#>  omega_1 NC coverage: 70% (7/10) [0.1s]
#>  omega_1 NC coverage: 80% (8/10) [0.1s]
#>  omega_1 NC coverage: 90% (9/10) [0.1s]
#>  omega_1 NC coverage: 100% (10/10) [0.1s]
#> omega_2 NC coverage: 10 tasks (sequential)
#>  omega_2 NC coverage: 10% (1/10) [0s]
#>  omega_2 NC coverage: 20% (2/10) [0s]
#>  omega_2 NC coverage: 30% (3/10) [0s]
#>  omega_2 NC coverage: 40% (4/10) [0s]
#>  omega_2 NC coverage: 50% (5/10) [0.1s]
#>  omega_2 NC coverage: 60% (6/10) [0.1s]
#>  omega_2 NC coverage: 70% (7/10) [0.1s]
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
#>  NC capture null: 10% (20/200) [0.4s]
#>  NC capture null: 20% (40/200) [0.7s]
#>  NC capture null: 30% (60/200) [1.1s]
#>  NC capture null: 40% (80/200) [1.5s]
#>  NC capture null: 50% (100/200) [1.8s]
#>  NC capture null: 60% (120/200) [2.2s]
#>  NC capture null: 70% (140/200) [2.6s]
#>  NC capture null: 80% (160/200) [2.9s]
#>  NC capture null: 90% (180/200) [3.3s]
#>  NC capture null: 100% (200/200) [3.7s]
#> iconic_diagnose complete. Call summary() or print() on the result for the full diagnosis.
est <- iconic_estimate(data, diagnosis = diag)
#> Estimating features: 10 tasks (sequential)
#>  Estimating features: 10% (1/10) [0s]
#>  Estimating features: 20% (2/10) [0s]
#>  Estimating features: 30% (3/10) [0s]
#>  Estimating features: 40% (4/10) [0s]
#>  Estimating features: 50% (5/10) [0s]
#>  Estimating features: 60% (6/10) [0s]
#>  Estimating features: 70% (7/10) [0s]
#>  Estimating features: 80% (8/10) [0s]
#>  Estimating features: 90% (9/10) [0s]
#>  Estimating features: 100% (10/10) [0s]
rec <- iconic_recommend(data, diagnosis = diag, estimate = est,
                        auto_sensitivity = FALSE)
print(rec)
#> <iconic_recommendation>
#>  Recommended (composite): COCA
#>  Estimated total effect (mean): -0.4582
#> 
#>  Eligible alternatives: UNADJ, DIRECT 
#> 
#>  Full ranking:
#>  1. COCA -- requires: valid NCs (A1), completeness (A2 not required)
#>  2. UNADJ -- requires: no assumptions (naive OLS)
#>  3. DIRECT -- requires: G + W as covariates (no causal identification)
#> 
#>  Ineligible:
#>  IV2SLS -- ineligible -- requires G + F_G>=10 (F_G=0.9)
#>  PGC -- ineligible -- requires G + W + F_G>=10 + completeness (completeness: satisfied)
#>  IV2SLS2 -- ineligible -- requires mediation data (supply M)
#>  PGC2 -- ineligible -- requires mediation data (supply M)
#>  PGC2Gm -- ineligible -- requires mediation data (supply M)
```
