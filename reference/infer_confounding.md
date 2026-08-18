# Infer confounding parameters from the user's data

Estimates the held-fixed confounding parameters that
[`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md)
and
[`iconic_prospect()`](https://seantbresnahan.com/iconic/reference/iconic_prospect.md)
use when `confounding = "inferred"`.

## Usage

``` r
infer_confounding(
  data,
  diagnosis = NULL,
  estimate = NULL,
  n_cores = 1,
  max_infer_tasks = 50
)
```

## Arguments

- data:

  An `iconic_data` object.

- diagnosis:

  Optional `iconic_diagnosis` from
  [`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md).
  Used to check instrument strength before using estimator gaps.

- estimate:

  Optional estimate data frame from
  [`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md).
  If `NULL`, estimates are computed internally.

- n_cores:

  Number of parallel workers for omega inference and k permutation.
  Default 1 (sequential). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix and a PSOCK cluster on Windows.

- max_infer_tasks:

  Cap on the number of mediators and of outcome features used when
  `estimate = NULL` and estimates must be computed internally. The
  confounding-strength gaps are averages across the mediator x feature
  grid, so a random subset of at most `max_infer_tasks` mediators and
  `max_infer_tasks` features gives an unbiased Monte Carlo estimate at
  much lower cost. Default 50. Panels smaller than the cap are used in
  full.

## Value

An `iconic_confounding` S3 object: a named list with `$conf_strength`,
`$mo_confounding`, `$omega_1`, `$omega_2`, `$k`, `$unavailable`
(character vector of parameters that could not be inferred), and
`$warnings` (character vector of accumulated warnings). Each parameter
slot is a list with `estimate`, `method`, `assumption`, `available`
(logical), and `warning` (character or NULL).

## Inference methods

- **conf_strength** (delta): the gap between the unadjusted OLS estimate
  and the IV2SLS estimate, averaged across features. Requires a valid
  exposure instrument G (F \>= 10). Assumes the exclusion restriction
  holds – using IV2SLS validity to calibrate a benchmark that tests
  IV2SLS validity is circular, so the estimate should be interpreted as
  a best-case calibration.

- **mo_confounding** (delta_mo): the gap between the IV2SLS NIE and the
  IV2SLS2 NIE, averaged across features. Requires valid G + Gm (both F
  \>= 10). Same circularity caveat.

- **omega** (omega_1, omega_2): the square root of the R-squared from
  regressing each negative-control feature on the outcome residualized
  on X + C, averaged across features. Conflates NC coverage with
  confounder strength – reported as a composite, not pure coverage.

- **k**: the number of latent confounders, estimated via parallel
  analysis (Horn, 1965) on the correlation matrix of outcomes
  residualized on X + C. Requires at least 5 outcome features. Returns a
  point estimate and a bootstrap confidence interval.

## Defaults (when inference is unavailable)

When a parameter cannot be inferred from the data, the following
defaults are used (from the simulation calibration):

|  |  |  |
|----|----|----|
| **Parameter** | **Default** | **Source** |
| `conf_strength` (delta) | 0.8 | Simulation calibration |
| `mo_confounding` (delta_mo) | 0.8 | Simulation calibration |
| `omega_1, omega_2` | 0.7 | NC coverage (simulation calibration) |
| `k` | 1 | Single-confounder assumption (typical mediation) |
| `phi` | 0.8 | Strong mediator instrument assumption |
| `lambda_XM`, `lambda_MY` | shared | Per-path confounder loadings |

These defaults are reported with a warning so the user knows the value
was not inferred from their data.

## Examples

``` r
data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100*10), 10, 100),
M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
W = matrix(rnorm(100*10), 10, 100))
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
#> Estimating mediation effects: 10 tasks (sequential)
#>  Estimating mediation effects: 10% (1/10) [0s]
#>  Estimating mediation effects: 20% (2/10) [0s]
#>  Estimating mediation effects: 30% (3/10) [0.1s]
#>  Estimating mediation effects: 40% (4/10) [0.1s]
#>  Estimating mediation effects: 50% (5/10) [0.1s]
#>  Estimating mediation effects: 60% (6/10) [0.1s]
#>  Estimating mediation effects: 70% (7/10) [0.2s]
#>  Estimating mediation effects: 80% (8/10) [0.2s]
#>  Estimating mediation effects: 90% (9/10) [0.2s]
#>  Estimating mediation effects: 100% (10/10) [0.3s]
#> omega_1 NC coverage: 10 tasks (sequential)
#>  omega_1 NC coverage: 10% (1/10) [0s]
#>  omega_1 NC coverage: 20% (2/10) [0s]
#>  omega_1 NC coverage: 30% (3/10) [0s]
#>  omega_1 NC coverage: 40% (4/10) [0s]
#>  omega_1 NC coverage: 50% (5/10) [0s]
#>  omega_1 NC coverage: 60% (6/10) [0.1s]
#>  omega_1 NC coverage: 70% (7/10) [0.1s]
#>  omega_1 NC coverage: 80% (8/10) [0.1s]
#>  omega_1 NC coverage: 90% (9/10) [0.1s]
#>  omega_1 NC coverage: 100% (10/10) [0.1s]
#> omega_2 NC coverage: 10 tasks (sequential)
#>  omega_2 NC coverage: 10% (1/10) [0s]
#>  omega_2 NC coverage: 20% (2/10) [0s]
#>  omega_2 NC coverage: 30% (3/10) [0s]
#>  omega_2 NC coverage: 40% (4/10) [0s]
#>  omega_2 NC coverage: 50% (5/10) [0s]
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
#>  NC capture null: 20% (40/200) [0.8s]
#>  NC capture null: 30% (60/200) [1.2s]
#>  NC capture null: 40% (80/200) [1.6s]
#>  NC capture null: 50% (100/200) [2s]
#>  NC capture null: 60% (120/200) [2.4s]
#>  NC capture null: 70% (140/200) [2.8s]
#>  NC capture null: 80% (160/200) [3.2s]
#>  NC capture null: 90% (180/200) [3.6s]
#>  NC capture null: 100% (200/200) [4s]
#> iconic_diagnose complete. Call summary() or print() on the result for the full diagnosis.
est <- iconic_estimate(data, diagnosis = diag)
#> Estimating mediation effects: 10 tasks (sequential)
#>  Estimating mediation effects: 10% (1/10) [0s]
#>  Estimating mediation effects: 20% (2/10) [0s]
#>  Estimating mediation effects: 30% (3/10) [0s]
#>  Estimating mediation effects: 40% (4/10) [0s]
#>  Estimating mediation effects: 50% (5/10) [0.1s]
#>  Estimating mediation effects: 60% (6/10) [0.1s]
#>  Estimating mediation effects: 70% (7/10) [0.1s]
#>  Estimating mediation effects: 80% (8/10) [0.1s]
#>  Estimating mediation effects: 90% (9/10) [0.1s]
#>  Estimating mediation effects: 100% (10/10) [0.1s]
conf <- infer_confounding(data, diagnosis = diag, estimate = est,
                          max_infer_tasks = 5)
#> omega_1 NC coverage: 10 tasks (sequential)
#>  omega_1 NC coverage: 10% (1/10) [0s]
#>  omega_1 NC coverage: 20% (2/10) [0s]
#>  omega_1 NC coverage: 30% (3/10) [0s]
#>  omega_1 NC coverage: 40% (4/10) [0s]
#>  omega_1 NC coverage: 50% (5/10) [0s]
#>  omega_1 NC coverage: 60% (6/10) [0.1s]
#>  omega_1 NC coverage: 70% (7/10) [0.1s]
#>  omega_1 NC coverage: 80% (8/10) [0.1s]
#>  omega_1 NC coverage: 90% (9/10) [0.1s]
#>  omega_1 NC coverage: 100% (10/10) [0.1s]
#> omega_2 NC coverage: 10 tasks (sequential)
#>  omega_2 NC coverage: 10% (1/10) [0s]
#>  omega_2 NC coverage: 20% (2/10) [0s]
#>  omega_2 NC coverage: 30% (3/10) [0s]
#>  omega_2 NC coverage: 40% (4/10) [0s]
#>  omega_2 NC coverage: 50% (5/10) [0s]
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
print(conf)
#> <iconic_confounding>
#>  conf_strength    default (0.8) -- weak instrument
#>  mo_confounding   default (0.8) -- weak mediator instrument
#>  omega_1          0.192 (sqrt(R^2) of W on Y residualized on X+C)
#>  warning: composite: coverage x confounder strength, not pure coverage
#>  omega_2          0.192 (sqrt(R^2) of W on Y residualized on X+C)
#>  warning: composite: coverage x confounder strength, not pure coverage
#>  k                5 [CI: 4, 6] (parallel analysis (Horn, 1965))
#> 
#>  Unavailable: rho_G1, rho_G2, conf_strength, mo_confounding 
#> 
#>  Warnings:
#>   conf_strength: weak instrument (F_G=2.2 < 10), inference unreliable, using default 0.8. 
#>   mo_confounding: weak mediator instrument (F_Gm < 10), inference unreliable, using default 0.8. 
```
