# Fit all eligible estimators on real data

Applies all eligible causal estimators to the user's real data (supplied
as an
[`iconic_data`](https://seantbresnahan.com/iconic/reference/iconic_data.md)
object) and returns tidy per-feature (and per-mediator) point estimates,
standard errors, and p-values.

## Usage

``` r
iconic_estimate(
  data,
  methods = NULL,
  diagnosis = NULL,
  alpha = 0.05,
  n_cores = 1,
  min_f = NULL,
  run_all = FALSE,
  se_method = c("delta", "bootstrap", "composite"),
  n_boot = 500,
  effect_scale = c("loghr", "rmst"),
  tau = NULL
)
```

## Arguments

- data:

  An `iconic_data` object from
  [`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md).

- methods:

  Optional character vector of estimator names to run (e.g.
  `c("IV2SLS", "PGC")`). Default `NULL`: run all eligible.

- diagnosis:

  Optional `iconic_diagnosis` object from
  [`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md).
  When supplied, only eligible estimators are run (unless `methods` or
  `run_all` overrides).

- alpha:

  Significance threshold for significance flags. Default 0.05.

- n_cores:

  Number of parallel workers for per-feature (and per-mediator)
  estimation. Default 1 (sequential). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix and a PSOCK cluster on Windows.

- min_f:

  Minimum partial F for the per-transcript weak-instrument gate inside
  IV2SLS, IV2SLS2, and PGC2Gm. Default `NULL`: inherits from
  `diagnosis$min_f` when a diagnosis is supplied, otherwise 10. Pass an
  explicit value to override.

- run_all:

  Logical. Default `FALSE`. When `TRUE`, overrides diagnosis-based
  eligibility and runs every method whose required data exists (via
  [`.auto_eligible_methods()`](https://seantbresnahan.com/iconic/reference/dot-auto_eligible_methods.md)).
  The per-transcript `min_f` gate still applies. This is the "force run"
  escape hatch for exploratory analysis.

- se_method:

  Character: `"delta"` (default), `"bootstrap"`, or `"composite"`. When
  `"bootstrap"`, mediation NDE_se/NIE_se are replaced by the SD of
  `n_boot` nonparametric bootstrap resamples. When `"composite"`,
  mediation NIE_p is replaced by the Huang (2019) JT-comp composite null
  p-value, which accounts for the three-case structure of H0:
  alpha\*beta=0 and provides higher power than the Sobel/Wald test when
  signals are sparse. Only applies in mediation mode. NDE_p and NIE_se
  are unchanged.

- n_boot:

  Integer: number of bootstrap resamples when `se_method = "bootstrap"`.
  Default 500.

- effect_scale:

  Character: `"loghr"` (default) or `"rmst"`. Only used when
  `data$outcome_type = "survival"`. `"loghr"` fits Cox
  proportional-hazards models and reports log-hazard ratios. `"rmst"`
  regresses leave-one-out RMST pseudo-observations (Graw et al. 2009)
  via OLS, reporting effects on the restricted-mean-survival-time (time)
  scale — a collapsible alternative where the NDE/NIE product
  decomposition is exact. Ignored (with a message) when
  `outcome_type = "continuous"`.

- tau:

  Numeric: RMST restriction time horizon. Default `NULL` (90th
  percentile of follow-up). Used only when `effect_scale = "rmst"`.

## Value

A data frame. For total-effect mode: `feature`, `method`, `beta`, `se`,
`pvalue`, `significant`. For mediation mode: `feature`, `mediator`,
`method`, `NDE`, `NDE_se`, `NDE_p`, `NIE`, `NIE_se`, `NIE_p`,
`NDE_significant`, `NIE_significant`. When `outcome_type = "survival"`,
estimates are on the log-HR scale (`effect_scale = "loghr"`) or the
RMST/time scale (`effect_scale = "rmst"`).

## Details

If a `diagnosis` from
[`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md)
is supplied, only eligible estimators are run by default. If `methods`
is supplied, only those methods are run (user override). If neither is
supplied, all estimators whose required inputs are present are run.

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
#>  Estimating features: 50% (5/10) [0.1s]
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
head(est)
#>           feature method       beta        se    pvalue significant
#> UNADJ   feature_1  UNADJ 0.04862181 0.1008958 0.6309506       FALSE
#> DIRECT  feature_1 DIRECT 0.05973886 0.1105759 0.5904037       FALSE
#> COCA    feature_1   COCA         NA        NA        NA          NA
#> UNADJ1  feature_2  UNADJ 0.10539909 0.1004526 0.2966474       FALSE
#> DIRECT1 feature_2 DIRECT 0.13696567 0.1078149 0.2073367       FALSE
#> COCA1   feature_2   COCA 0.87041158 0.8546164 0.3084489       FALSE

# Survival outcome
sdat <- iconic_data(X = rnorm(100), outcome_type = "survival",
surv_time = rexp(100), surv_event = rbinom(100, 1, 0.6),
G = rnorm(100), W = matrix(rnorm(100 * 10), 10, 100))
est <- iconic_estimate(sdat, effect_scale = "loghr")
est_rmst <- iconic_estimate(sdat, effect_scale = "rmst")
```
